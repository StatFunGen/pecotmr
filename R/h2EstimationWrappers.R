# =============================================================================
# Heritability Estimation Wrappers
# -----------------------------------------------------------------------------
# Consolidated entry points for univariate heritability estimation. Three
# methods are exposed via the `estimateH2` S4 generic:
#
#   * gLDSC (Generalized LD Score Regression, Xiong et al. 2024)
#   * HDL/sHDL (High-Definition Likelihood, Ning et al. 2020 + Zhao 2023)
#   * LDER (LD Eigenvalue Regression, Song et al. 2022)
#
# Plus shared utilities (block ops, WLS, jackknife SE, enrichment,
# genome-build checks, tauStar standardization, meta-analysis) and the
# converter that bridges H2Estimate results into the sLDSC postprocessing
# pipeline.
#
# The H2Estimate result class itself lives in R/h2Estimate.R.
# =============================================================================

#' @title Shared Utilities for Heritability Estimation
#' @description Internal helper functions for block operations, regression,
#'   jackknife SE, enrichment computation, and meta-analysis.
#' @name pecotmr-h2-utils
#' @keywords internal
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
NULL

# =============================================================================
# Block-level operations
# =============================================================================

#' @title Get SNP Indices Per Block
#' @description For each LD block, find the SNP indices from a reference that
#'   fall within the block boundaries.
#' @param snpInfo A data.frame with columns CHR, BP.
#' @param ldBlocks The LD block intervals: a \code{GRanges}, or the path
#'   of a BED file, read as BED with a 0-based start.
#' @return A list of integer vectors, one per block.
#' @keywords internal
snpsPerBlock <- function(snpInfo, ldBlocks) {
    blocksGr <- .asLdBlockRanges(ldBlocks)
    snpGr <- GRanges(
        seqnames = snpInfo$CHR,
        ranges = IRanges(start = snpInfo$BP, width = 1L)
    )
    hits <- findOverlaps(snpGr, blocksGr)
    split(queryHits(hits), subjectHits(hits))
}

# =============================================================================
# Baseline enrichment computation
# =============================================================================

# Enrichment SE: jackknife over per-block enrichment (preferred) or the
# delta-method fallback tauSe * M / |h2| when no blocks are available.
# @noRd
.baselineEnrichmentSe <- function(annotationJackknifeCoefs, tauSe, M_a, M, h2) {
    if (is.null(annotationJackknifeCoefs)) {
        return(tauSe * M / abs(h2))
    }
    nBlocks <- nrow(annotationJackknifeCoefs)
    rawH2 <- as.vector(annotationJackknifeCoefs %*% M_a)
    # A zero block h2 would divide the enrichment to Inf; NA drops the block
    # from the jackknife mean instead.
    h2Blocks <- replace(rawH2, rawH2 == 0, NA)
    enrichmentBlocks <- sweep(
        annotationJackknifeCoefs,
        1,
        h2Blocks,
        FUN = "/"
    ) *
        M
    enrichmentMean <- colMeans(enrichmentBlocks, na.rm = TRUE)
    enrichmentVar <- (nBlocks - 1) /
        nBlocks *
        colSums(sweep(enrichmentBlocks, 2, enrichmentMean)^2, na.rm = TRUE)
    sqrt(enrichmentVar)
}

# =============================================================================
# Genome build utilities
# =============================================================================

#' @title Validate Genome Build Consistency
#' @description Check that genome builds match between objects. Each object
#'   contributes a single genome build (from its \code{genome} slot). An
#'   object that records no build is skipped rather than treated as a
#'   mismatch: an LD reference often names no build (a genotype panel's
#'   variants rarely carry one), and an unknown build cannot contradict a
#'   known one.
#' @param ... Objects with a \code{genome} slot.
#' @return TRUE if all match, error otherwise.
#' @keywords internal
checkGenomeBuild <- function(...) {
    objects <- list(...)
    # Everything supported keeps its build in seqinfo() rather than in a slot
    # of its own (that redundancy is why the LdBlocks wrapper went away), so
    # one seqinfo read serves them all. Reject anything else up front --
    # notably a plain list(genome = ...), the pre-seqinfo shape this check
    # replaced -- instead of letting seqinfo() fail obscurely.
    supported <- c(
        "GenomicRanges",
        "GwasSumStats",
        "QtlSumStats",
        "LdStatistic",
        "AnnotationMatrix"
    )
    known <- map_lgl(map(map(objects, methods::is), is_in, supported), any)
    if (!all(known)) {
        abort("Unknown object type for genome build check")
    }
    # genome() is per-seqlevel, so flatten across objects before uniquing: one
    # object naming two builds and two objects disagreeing are the same error.
    # An object naming none contributes NA and is skipped, because an LD
    # reference often names no build and demanding one here would reject most
    # of them. The coordinate check in .h2CheckPositions() is the stronger
    # guard regardless.
    builds <- unname(discard(
        unique(list_c(map(objects, GenomeInfoDb::genome))),
        is.na
    ))
    if (length(builds) > 1L) {
        msg <- glue(
            "Genome build mismatch: {str_flatten(builds, ', ')}"
        )
        abort(msg)
    }
    invisible(TRUE)
}

# =============================================================================
# Gazal tau* standardization
# =============================================================================

#' @title Standardize Tau to Tau-Star (Gazal et al. 2017)
#' @description Compute the Gazal-standardized per-annotation effect
#'   \eqn{\tau^*_C = \tau_C \cdot sd_C \cdot M_{ref} / h^2_g}, with jackknife SE
#'   from block-level tau values.
#' @param tau Numeric vector of per-annotation regression coefficients.
#' @param annotationJackknifeCoefs Numeric matrix (nBlocks x nAnnotations) of
#'   block-level tau estimates from delete-one jackknife.
#' @param sdAnnot Numeric vector of per-annotation standard deviations, same
#'   length as \code{tau}.
#' @param MRef Scalar integer, total number of reference-panel SNPs.
#' @param h2g Numeric scalar, total estimated SNP heritability.
#' @return A list with:
#'   \describe{
#'     \item{tauStar}{Numeric vector of standardized tau values.}
#'     \item{tauStarSe}{Numeric vector of jackknife SE for tauStar.}
#'   }
#' @keywords internal
standardizeTauStar <- function(
    tau,
    annotationJackknifeCoefs,
    sdAnnot,
    MRef,
    h2g
) {
    if (length(tau) != length(sdAnnot)) {
        abort("standardizeTauStar: tau and sdAnnot must have the same length.")
    }
    if (h2g == 0) {
        abort("standardizeTauStar: h2g must be non-zero.")
    }

    # Gazal standardization: tau* = tau * sdAnnot * MRef / h2g
    coef <- sdAnnot * MRef / h2g
    tauStar <- tau * coef

    # Jackknife SE from block-level tau
    tauStarBlocks <- sweep(annotationJackknifeCoefs, 2L, coef, FUN = "*")
    nBlocks <- nrow(tauStarBlocks)
    jkVar <- apply(tauStarBlocks, 2L, var, na.rm = TRUE)
    tauStarSe <- sqrt((nBlocks - 1)^2 / nBlocks * jkVar)

    list(tauStar = tauStar, tauStarSe = tauStarSe)
}

# =============================================================================
# DerSimonian-Laird random-effects meta-analysis
# =============================================================================

# Degenerate meta-analysis inputs: k = 0 (all NA) or k = 1 (pass-through), plus
# positivity validation. Returns NULL when the full estimator should run.
# @noRd
.rmaMetaEarly <- function(means, ses, k) {
    if (k == 0L) {
        return(list(
            mean = NA_real_,
            se = NA_real_,
            tau2 = NA_real_,
            I2 = NA_real_,
            Q = NA_real_
        ))
    }
    if (k == 1L) {
        return(list(mean = means[1], se = ses[1], tau2 = 0, I2 = 0, Q = 0))
    }
    if (any(!is.finite(ses) | ses <= 0)) {
        abort(".rmaMeta: all ses must be positive and finite.")
    }
    NULL
}

#' @title Options for the Random-Effects Meta-Analysis Engine
#' @description Build a record of extra arguments for
#'   \code{metafor::rma()}, the random-effects meta-analysis engine behind
#'   \code{\link{metaSldscRandom}},
#'   \code{\link{metaAnalysisPerCondition}} and
#'   \code{\link{calculateFeatureScores}}.
#' @param ... Arguments for \code{metafor::rma()}: \code{test} (the
#'   null-distribution choice -- \code{"knha"} is the small-\eqn{k}
#'   correction, and these meta-analyses often run on few studies),
#'   \code{level}, \code{weighted}, \code{tau2} and \code{control}.
#'   \code{rma()}'s signature ends in \code{...}, so names cannot be
#'   checked here and a misspelling is passed through. \code{yi} and
#'   \code{sei} are the estimate/standard-error pair pecotmr assembles and
#'   are refused, as is \code{method}: that is the caller's
#'   \code{metaMethod}.
#'
#'   These settings apply to the first fit. pecotmr falls back to the
#'   closed-form DerSimonian-Laird estimator when an iterative one fails, and
#'   the fallback carries them too.
#' @return A \code{MethodOptions} record for the \code{metaArgs} argument.
#' @seealso \code{\link{metaSldscRandom}},
#'   \code{\link{metaAnalysisPerCondition}}
#' @examples
#' RmaOptions(test = "knha")
#' @export
RmaOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            yi = "the per-study estimates pecotmr assembles",
            sei = "the per-study standard errors pecotmr assembles",
            method = "the caller's `metaMethod`"
        ),
        "RmaOptions"
    )
    .newMethodOptions(
        "metafor::rma",
        defaults = list(),
        extra = extra,
        label = "RmaOptions",
        engine = "rma"
    )
}

# tryCatch handler for the iterative estimators (only reached when method !=
# "DL", so DL is always a real fallback here): they can fail on small /
# near-homogeneous inputs, so fall back to closed-form DerSimonian-Laird, which
# never iterates.
.rmaMetaFallbackToDL <- function(e, means, ses, method, metaArgs = list()) {
    msg <- glue(
        ".rmaMeta: metafor::rma(method = '{method}') failed ",
        "({conditionMessage(e)}); falling back to DL."
    )
    warn(msg)
    exec(
        metafor::rma,
        yi = means,
        sei = ses,
        method = "DL",
        !!!as.list(metaArgs)
    )
}

#' @title Random-Effects Meta-Analysis (metafor wrapper)
#' @description Thin convenience wrapper around \code{metafor::rma()}. pecotmr
#'   does not implement its own meta-analysis engine -- this helper only guards
#'   the degenerate 0- and 1-study cases (where \code{rma} is undefined) and
#'   reshapes the \code{metafor} fit into a small list. The between-study
#'   variance estimator is chosen by \code{method} and passed straight through.
#' @param means Numeric vector of study-level point estimates.
#' @param ses Numeric vector of study-level standard errors (must be positive
#'   and finite).
#' @param method Between-study variance estimator forwarded to
#'   \code{metafor::rma(method = )} (e.g. \code{"DL"} (default), \code{"REML"},
#'   \code{"ML"}, \code{"EB"}).
#' @param metaArgs Extra \code{metafor::rma()} arguments, an
#'   \code{\link{RmaOptions}} record.
#' @return A list with \code{mean}, \code{se}, \code{tau2}, \code{I2} (a [0, 1]
#'   proportion) and \code{Q} (Cochran's Q).
#' @importFrom metafor rma
#' @noRd
#' @importFrom rlang try_fetch
.rmaMeta <- function(means, ses, method = "DL", metaArgs = list()) {
    k <- length(means)
    if (k != length(ses)) {
        abort(".rmaMeta: means and ses must have the same length.")
    }
    early <- .rmaMetaEarly(means, ses, k)
    if (!is.null(early)) {
        return(early)
    }
    # metafor is a hard dependency; it reports I-squared as a percentage.
    # DL (closed-form DerSimonian-Laird) is the last-resort estimator: run it
    # directly and let any error propagate. Only the iterative estimators get a
    # tryCatch fallback to DL.
    fit <- if (identical(method, "DL")) {
        exec(
            metafor::rma,
            yi = means,
            sei = ses,
            method = "DL",
            !!!as.list(metaArgs)
        )
    } else {
        try_fetch(
            exec(
                metafor::rma,
                yi = means,
                sei = ses,
                method = method,
                !!!as.list(metaArgs)
            ),
            error = function(cnd) {
                .rmaMetaFallbackToDL(cnd, means, ses, method, metaArgs)
            }
        )
    }
    list(
        mean = as.numeric(fit$b),
        se = as.numeric(fit$se),
        tau2 = as.numeric(fit$tau2),
        I2 = as.numeric(fit$I2) / 100,
        Q = as.numeric(fit$QE)
    )
}


#' @title LDER: LD Eigenvalue Regression
#' @description Estimate heritability using LD eigenvalue regression (Song et
#'   al. 2022). Supports univariate global and local estimation, with optional
#'   annotation stratification.
#' @name pecotmr-h2-lder
#' @keywords internal
#' @references
#'   Song S, Jiang W, Zhang Y, Hou L, Zhao H (2022). Leveraging LD
#'   eigenvalue regression to improve the estimation of SNP heritability
#'   and confounding inflation. Am J Hum Genet, 109(5):802-811.
NULL

# =============================================================================
# Univariate LDER
# =============================================================================

#' @title Univariate LDER
#' @description Estimate SNP heritability by LD eigenvalue regression (Song et
#'   al. 2022). Within each LD block the z-scores are whitened in the
#'   eigenbasis,
#'   \eqn{x_i = (u_i'z)/\sqrt{\lambda_i}}, and the whitened \eqn{\chi^2 = x_i^2}
#'   is regressed on \eqn{\lambda_i} via IRWLS:
#'   \eqn{E[x_i^2] = (N h2/M)\lambda_i
#'   + (1 + N a)}, so \eqn{Cov(z) = (N h2/M) R^2 + (1 + N a) R}. The IRWLS
#'   weight
#'   uses an in-sample (\eqn{\min(1,\lambda)}) or out-of-sample (logistic) form
#'   selected by \code{inSample}, and a two-stage step estimates the
#'   confounding inflation \eqn{a} on a \eqn{\chi^2}-filtered subset. The
#'   stratified extension regresses on \eqn{\lambda_i \ell_{ia}}
#'   (\eqn{\ell_{ia} = \sum_j A_{ja} V_{ji}^2}) and reduces to the univariate
#'   model when there is a single genome-wide annotation.
#' @name pecotmr-h2-lder
#' @keywords internal
#' @references
#'   Song S, Jiang W, Zhang Y, et al. (2022). Leveraging LD eigenvalue
#'   regression to improve the estimation of SNP heritability and confounding
#'   inflation. Am J Hum Genet, 109:802-811.
#' @importFrom stats lm.fit quantile cor
NULL

# Baseline annotation matrix (NULL when absent or empty). Shared by LDER + HDL.
.h2BaselineMat <- function(annotations) {
    if (is.null(annotations)) {
        return(NULL)
    }
    baseline <- baselineAnnotations(annotations)
    if (ncol(assay(baseline, "annotations")) > 0) {
        assay(baseline, "annotations")
    } else {
        NULL
    }
}

# The baseline design all four estimators fit, NULL when unstratified: the
# user's baseline columns, with an all-ones "base" column put first unless
# they already span it (a partition of the genome does; one annotation
# alone, or a few that leave variants uncovered, do not). Without the base,
# the variants no baseline annotation covers would be modelled as carrying
# no heritability, and each annotation's coefficient would absorb the
# genome-wide signal instead of its excess over it. Adding it to a baseline
# that spans it would make the design singular. `userCols` records which
# columns are the user's.
# @noRd
.h2DesignMat <- function(annotations) {
    baselineMat <- .h2BaselineMat(annotations)
    if (is.null(baselineMat)) {
        return(NULL)
    }
    A <- if (.h2SpansOnes(baselineMat)) {
        baselineMat
    } else {
        cbind(base = 1, baselineMat)
    }
    `attr<-`(
        A,
        "userCols",
        match(colnames(baselineMat), colnames(A))
    )
}

# TRUE when the all-ones vector lies in the column span of `B`, by least
# squares on B'B (pivoted, so overlapping or duplicated columns do no harm),
# which keeps a genome-wide B from being copied into a QR.
# @noRd
.h2SpansOnes <- function(B) {
    decomp <- qr(crossprod(B))
    beta <- qr.coef(decomp, colSums(B))
    beta[is.na(beta)] <- 0
    max(abs(1 - as.vector(B %*% beta))) < 1e-6
}

# The partitioned heritabilities, enrichment and per-annotation coefficients
# of a stratified fit, reported for the user's baseline annotations alone,
# from its per-variant-effect coefficients `tau` over the design `A` and their
# delete-one-block values `tauBlocks` (folds x columns). Annotations overlap
# in general -- the base overlaps every one -- so an annotation's
# heritability is (A'A tau)_a, not tau_a M_a, and the total is 1'A tau.
# @noRd
.h2StratifiedReport <- function(A, tau, tauBlocks) {
    user <- attr(A, "userCols")
    AtA <- crossprod(A)
    size <- colSums(A)
    propSnps <- size[user] / nrow(A)
    propH2 <- as.vector(AtA %*% tau)[user] / sum(size * tau)
    enrichment <- propH2 / propSnps
    blockH2 <- (tauBlocks %*% AtA)[, user, drop = FALSE]
    blockEnrichment <- sweep(
        blockH2 / as.vector(tauBlocks %*% size),
        2L,
        propSnps,
        "/"
    )
    enrichmentSe <- apply(blockEnrichment, 2L, .jackknifeSeFromLoo)
    tauSe <- apply(tauBlocks, 2L, .jackknifeSeFromLoo)[user]
    list(
        tau = tau[user],
        tauSe = tauSe,
        annotationJackknifeCoefs = tauBlocks[, user, drop = FALSE],
        enrichment = tibble(
            annotation = colnames(A)[user],
            tau = unname(tau[user]),
            tauSe = unname(tauSe),
            enrichment = unname(enrichment),
            enrichmentSe = unname(enrichmentSe),
            enrichmentP = unname(
                pnorm(abs(enrichment - 1) / enrichmentSe, lower.tail = FALSE) *
                    2
            ),
            propH2 = unname(propH2),
            propSnps = unname(propSnps)
        )
    )
}

# The result of a stratified or unstratified fit: `tau` and `tauBlocks` over
# the design `A` (NULL when unstratified, where the heritability stands in
# for the coefficient).
# @noRd
.h2FitResult <- function(
    h2,
    h2Se,
    intercept,
    interceptSe,
    A,
    tau,
    tauBlocks,
    localDf,
    scoreStats
) {
    if (is.null(A)) {
        jkRes <- list(
            se = c(h2Se, interceptSe),
            tauSe = h2Se,
            annotationJackknifeCoefs = NULL
        )
        return(.h2Result(h2, jkRes, intercept, h2, localDf, NULL, scoreStats))
    }
    report <- .h2StratifiedReport(A, tau, tauBlocks)
    jkRes <- list(
        se = c(h2Se, interceptSe),
        tauSe = report$tauSe,
        annotationJackknifeCoefs = report$annotationJackknifeCoefs
    )
    .h2Result(
        h2,
        jkRes,
        intercept,
        report$tau,
        localDf,
        report$enrichment,
        scoreStats
    )
}

# Per-block whitened statistics x_i = (u_i'z)/sqrt(lam_i) and stratified
# eigenvalue scores l_{ia} = sum_j A_{ja} V_{ji}^2 (all-ones column when
# unstratified), plus the candidate annotations' scores, from the block's
# decomposition. Near-null directions (lam < 1e-6) are zeroed, as in the
# original LDER implementation, which keeps every direction of the block:
# the zero-eigenvalue directions a decomposition from a panel's factor does
# not list are added as rows with x = 0 and lam = 0.
.lderBlockStats <- function(block, eigenRef, z, baselineMat, candMat) {
    dec <- .ldEigenBlockAt(block, eigenRef)
    lam <- dec$values
    raw <- as.vector(
        crossprod(dec$vectors, z[dec$snpIdx]) / sqrt(pmax(lam, 1e-12))
    )
    nAll <- length(lam) + dec$nNull
    list(
        x = c(replace(raw, lam < 1e-6, 0), numeric(dec$nNull)),
        lam = c(lam, numeric(dec$nNull)),
        ldAnnot = .h2EigenAnnot(dec, baselineMat) %||% matrix(1, nAll, 1),
        ldCand = .h2EigenAnnot(dec, candMat)
    )
}

# One block's eigenvalue scores sum_k V_ki^2 a_k for each column of
# `annotMat`, with a zero row for each zero-eigenvalue direction the
# decomposition does not list: every use multiplies a direction's row by its
# eigenvalue. NULL without annotations.
.h2EigenAnnot <- function(dec, annotMat) {
    if (is.null(annotMat)) {
        return(NULL)
    }
    listed <- crossprod(
        dec$vectors^2,
        annotMat[dec$snpIdx, , drop = FALSE]
    )
    rbind(listed, matrix(0, dec$nNull, ncol(annotMat)))
}

# The candidate (tested) annotations, or NULL when there are none to score.
.h2CandidateMat <- function(annotations) {
    if (is.null(annotations)) {
        return(NULL)
    }
    candMat <- assay(testedAnnotations(annotations), "annotations")
    if (ncol(candMat) == 0L) NULL else candMat
}

# Concatenate per-block stats into the genome-wide regression design. M_a is
# the annotation SNP counts (univariate: total number of directions).
.lderDesign <- function(blockStats, baselineMat) {
    x <- list_c(map(blockStats, "x"))
    ldAnnotList <- map(blockStats, "ldAnnot")
    list(
        x = x,
        lam = list_c(map(blockStats, "lam")),
        ldAnnot = exec(rbind, !!!ldAnnotList),
        blockId = rep(
            seq_along(blockStats),
            lengths(map(blockStats, "lam"))
        ),
        M_a = if (is.null(baselineMat)) length(x) else colSums(baselineMat)
    )
}

# IRWLS weight w = 1/(1 + predGenetic)^2 * w2(lam); predGenetic uses the current
# (clamped) slopes, w2 = pmin(1,lam) (in-sample) or logistic (out-of-sample).
.lderWeights <- function(lam, ldAnnot, slopes, N, M, rough) {
    # A correlation matrix's eigenvalues cannot be negative, but a real,
    # rank-deficient LD block decomposes with a tail of tiny negative values
    # (order -1e-15) that are pure floating-point noise. Left alone they make
    # `pmin(1, lam)` negative, and the IRWLS step then takes sqrt() of a
    # negative weight -- NaN designs and an "NA/NaN/Inf in 'x'" abort from
    # lm.fit. Clamping at zero drops those null directions from the fit,
    # which is what a zero-information direction should contribute anyway.
    lam <- pmax(lam, 0)
    predG <- lam * as.vector(ldAnnot %*% pmin(pmax(slopes, 0), N / M))
    w2 <- if (rough) 1 / (1 + exp(-5 * (lam - 1))) else pmin(1, lam)
    1 / (1 + predG)^2 * w2
}

# 3-iteration IRWLS. The first nA design columns are the eigenvalue-score
# columns (coeff[1:nA] = slopes); a trailing intercept column may follow.
.lderIrwls <- function(design, y, lam, ldAnnot, N, M, rough, iter = 3L) {
    nA <- ncol(ldAnnot)
    w <- rep(1, length(y))
    for (i in seq_len(iter)) {
        coeff <- lm.fit(design * sqrt(w), y * sqrt(w))$coefficients
        w <- .lderWeights(lam, ldAnnot, coeff[seq_len(nA)], N, M, rough)
    }
    coeff
}

# Whitened-chi2 regression on lam*ldAnnot (+ intercept). Free intercept when a
# is NULL (estimates inflation a); fixed N*a+1 otherwise. tau_a = coeff_a / N.
.lderCalH2 <- function(x, lam, ldAnnot, N, M, M_a, a = NULL, rough = FALSE) {
    chi2 <- x^2
    nA <- ncol(ldAnnot)
    scores <- lam * ldAnnot
    if (is.null(a)) {
        coeff <- .lderIrwls(cbind(scores, 1), chi2, lam, ldAnnot, N, M, rough)
        a <- (coeff[nA + 1L] - 1) / N
    } else {
        coeff <- .lderIrwls(
            scores,
            chi2 - (N * a + 1),
            lam,
            ldAnnot,
            N,
            M,
            rough
        )
    }
    tau <- coeff[seq_len(nA)] / N
    h2a <- tau * M_a
    list(h2 = sum(h2a), h2a = h2a, tau = tau, a = a)
}

# Two-stage LDER (Song et al.): a rough GC estimate sets a chi2 outlier
# threshold; inflation a is estimated on the filtered subset, then h2 is fit
# with a fixed (>= 0). twostage = FALSE returns the free-intercept fit.
.lderGetRes <- function(x, lam, ldAnnot, N, M, M_a, rough, twostage) {
    if (!twostage) {
        return(.lderCalH2(x, lam, ldAnnot, N, M, M_a, a = NULL, rough = rough))
    }
    m <- length(x)
    newGC <- 1 + N * .lderCalH2(x, lam, ldAnnot, N, M, M_a, rough = FALSE)$a
    s2thld <- if (newGC < 1.05) {
        3.5 * sqrt(m / N) + 2.5
    } else {
        quantile(x^2, 1 - 5 / m)
    }
    keep <- x^2 <= s2thld
    aFit <- .lderCalH2(
        x[keep],
        lam[keep],
        ldAnnot[keep, , drop = FALSE],
        N,
        M,
        M_a,
        a = NULL,
        rough = rough
    )
    .lderCalH2(x, lam, ldAnnot, N, M, M_a, a = max(aFit$a, 0), rough = rough)
}

# Delete-one-block jackknife SE: sqrt((nB-1)/nB * sum((jk - mean)^2)). M_aFull
# is fixed across folds so h2_loo = sum(tau_loo * M_aFull) shares a scale.
.lderJackknife <- function(
    blockStats,
    baselineMat,
    N,
    M,
    M_aFull,
    rough,
    twostage
) {
    nB <- length(blockStats)
    loo <- map(
        seq_len(nB),
        .lderLooBlockFit,
        blockStats = blockStats,
        baselineMat = baselineMat,
        N = N,
        M = M,
        M_aFull = M_aFull,
        rough = rough,
        twostage = twostage
    )
    h2Loo <- map_dbl(loo, "h2")
    aLoo <- map_dbl(loo, "a")
    tauList <- map(loo, "tau")
    jkCoefs <- exec(rbind, !!!tauList)
    list(
        h2Se = .jackknifeSeFromLoo(h2Loo),
        aSe = .jackknifeSeFromLoo(aLoo),
        tauSe = apply(jkCoefs, 2, .jackknifeSeFromLoo),
        annotationJackknifeCoefs = jkCoefs
    )
}

# Delete-one-block jackknife standard error from the n leave-one-out
# estimates `v`: sqrt((n - 1) / n * sum((v - mean(v))^2)). The LDER, gLDSC,
# HDL and S-LDSC jackknifes all report this quantity.
# @noRd
.jackknifeSeFromLoo <- function(v) {
    n <- length(v)
    sqrt((n - 1) / n * sum((v - mean(v))^2))
}

# Local heritability for one block: single-slope fit with intercept fixed at
# the global inflation; h2Local = tau * (block directions).
.lderLocalBlock <- function(bs, b, N, a) {
    d <- bs$lam
    if (length(d) < 3) {
        return(tibble(
            blockId = b,
            h2Local = NA_real_,
            h2LocalSe = NA_real_
        ))
    }
    chi2 <- bs$x^2
    w <- 1 / (2 * pmax(chi2, 1)^2)
    slope <- sum(w * d * (chi2 - (N * a + 1))) / sum(w * d^2)
    info <- sum(w * (N * d)^2)
    tibble(
        blockId = b,
        h2Local = slope / N * length(d),
        h2LocalSe = length(d) / sqrt(max(info, 1e-10))
    )
}

# Per-block local heritability (NULL unless requested).
.lderLocal <- function(blockStats, N, a) {
    bind_rows(
        map(
            seq_along(blockStats),
            .lderLocalRow,
            blockStats = blockStats,
            N = N,
            a = a
        )
    )
}

# LDER candidate-annotation score statistics against the fitted baseline:
# the regression of chi2 = x^2 on mu = N lam (l_a tau) + N a + 1, scored for
# each candidate's column N lam l_c with the weights 1 / (2 mu^2) and the
# baseline and inflation as nuisance (.h2EfficientScoreBlocks()).
.lderScoreStats <- function(blockStats, annotations, design, fit, N) {
    candidate <- testedAnnotations(annotations)
    candMat <- assay(candidate, "annotations")
    if (ncol(candMat) == 0) {
        return(list(enrichment = NULL, scoreStats = NULL))
    }
    lam <- design$lam
    mu <- pmax(
        N * lam * as.vector(design$ldAnnot %*% fit$tau) + (N * fit$a + 1),
        1e-8
    )
    w <- 1 / (2 * mu^2)
    ldCand <- exec(rbind, !!!map(blockStats, "ldCand"))
    Ub <- .h2EfficientScoreBlocks(
        N * lam * ldCand,
        cbind(N * lam * design$ldAnnot, N),
        w,
        w * (design$x^2 - mu),
        design$blockId
    )
    sc <- .h2BlockScoreStats(Ub)
    .buildEnrichmentResult(
        SummarizedExperiment::colData(candidate),
        sc$scoreZ,
        sc$R
    )
}

lderUnivariate <- function(
    z,
    n,
    eigenRef,
    annotations = NULL,
    local = FALSE
) {
    rough <- !inSample(eigenRef)
    baselineMat <- .h2DesignMat(annotations)
    # One pass over the blocks, each decomposed only while its statistics are
    # taken, so a reference read from a panel holds one block at a time.
    blockStats <- map(
        .ldEigenBlocks(eigenRef),
        .lderBlockStats,
        eigenRef = eigenRef,
        z = z,
        baselineMat = baselineMat,
        candMat = .h2CandidateMat(annotations)
    )
    design <- .lderDesign(blockStats, baselineMat)
    M <- length(design$x)
    fit <- .lderGetRes(
        design$x,
        design$lam,
        design$ldAnnot,
        n,
        M,
        design$M_a,
        rough,
        TRUE
    )
    jk <- .lderJackknife(blockStats, baselineMat, n, M, design$M_a, rough, TRUE)
    localDf <- if (local) .lderLocal(blockStats, n, fit$a) else NULL
    .h2FitResult(
        h2 = fit$h2,
        h2Se = jk$h2Se,
        intercept = fit$a,
        interceptSe = jk$aSe,
        A = baselineMat,
        tau = fit$tau,
        tauBlocks = jk$annotationJackknifeCoefs,
        localDf = localDf,
        scoreStats = .lderCandidateScores(
            blockStats,
            annotations,
            design,
            fit,
            n
        )
    )
}

# LDER's score statistics for the candidate annotations; NULL when there are
# none to test.
# @noRd
.lderCandidateScores <- function(blockStats, annotations, design, fit, n) {
    if (is.null(annotations)) {
        return(NULL)
    }
    .lderScoreStats(blockStats, annotations, design, fit, n)$scoreStats
}


# Assemble the enrichment result shared by the h2 enrichment helpers: a
# per-annotation table (z + two-tailed p) plus the scoreStats bundle
# (z, R, annotation names) consumed downstream.
.buildEnrichmentResult <- function(candMeta, scoreZ, R) {
    list(
        enrichment = tibble(
            annotation = candMeta$name,
            scoreZ = unname(scoreZ),
            scoreP = unname(.zToPvalue(scoreZ))
        ),
        scoreStats = list(z = scoreZ, R = R, annotationNames = candMeta$name)
    )
}


# Assemble the shared H2Estimate-style result list (used by all methods).
.h2Result <- function(
    h2,
    jkRes,
    intercept,
    tau,
    localDf,
    enrichmentDf,
    scoreStats
) {
    list(
        h2 = h2,
        h2Se = jkRes$se[1],
        intercept = intercept,
        interceptSe = jkRes$se[2],
        tau = tau,
        tauSe = jkRes$tauSe,
        annotationJackknifeCoefs = jkRes$annotationJackknifeCoefs,
        local = localDf,
        enrichment = enrichmentDf,
        scoreStats = scoreStats
    )
}

#' @title g-LDSC: Generalized LD Score Regression
#' @description Estimate heritability and functional enrichment by feasible
#'   generalized least squares (FGLS) on LD scores (Xiong et al. 2024). The
#'   response \eqn{Z_j^2 - 1} is regressed on stratified LD scores
#'   \eqn{\ell_{ja} = \sum_k R^2_{jk} A_{ka}} (with a confounding intercept;
#'   the raw sample \eqn{R^2}, as the original gLDSC implementation forms
#'   them, not the unbiased \eqn{r^2} of the S-LDSC scores) using the full
#'   residual covariance
#'   \eqn{\Omega = (\hat{N\tau}\, R D_a R + R)^2} (elementwise), inverted per LD
#'   block via a 99\%-eigenvalue-mass truncation. Each block's \eqn{\Omega}
#'   is formed, decomposed and dropped in turn; for a reference read from a
#'   genotype panel it is formed from the panel's factor, never from the
#'   block's \eqn{R}. \eqn{\hat{N\tau} = \sum(Z^2-1)
#'   / \sum \ell} is a method-of-moments plug-in;
#'   \eqn{D_a = diag(\sum_a A_{ja})}.
#'   Coefficients are \eqn{\tau = (X'\Omega^{-1}X)^{-1} X'\Omega^{-1} y / N} and
#'   partitioned heritabilities \eqn{A'A\,\tau}. Standard errors come from a
#'   delete-one-block jackknife that subtracts each block's
#'   already-\eqn{\Omega}-
#'   weighted contribution from the accumulated normal equations.
#' @name pecotmr-h2-gldsc
#' @keywords internal
#' @references
#'   Xiong Z, Thach TQ, Zhang YD, Sham PC (2024). Improved estimation
#'   of functional enrichment in SNP heritability using feasible
#'   generalized least squares. HGG Advances, 5(2):100272.
#' @importFrom stats var pnorm
NULL

# The annotation matrix A the LD-score estimators regress on: the baseline
# design (.h2DesignMat()), or the all-ones base alone when unstratified.
.gldscAnnotMatrix <- function(annotations, M) {
    .h2DesignMat(annotations) %||%
        matrix(1, M, 1, dimnames = list(NULL, "base"))
}

# Per-block design: the block's annotations, its stratified LD scores
# l = (R o R) A, the candidate annotations' scores (R o R) C, and the
# response y = Z^2 - 1. The one LD score here that stays on the raw r^2: the
# original gLDSC implementation squares its LD matrix as is (middle_fun.R
# LDSCmatrix), and the GLS design is reproduced from it. A block read from a
# genotype panel gets its scores from the panel's factor a stripe of rows at
# a time (.ldSquaredProduct()).
.gldscBlockDesign <- function(block, panel, z, A, candMat) {
    idx <- block$idx
    Ab <- A[idx, , drop = FALSE]
    Cb <- if (!is.null(candMat)) candMat[idx, , drop = FALSE]
    scores <- .ldSquaredProduct(.ldStatBlockLd(block, panel), cbind(Ab, Cb))
    nA <- ncol(Ab)
    list(
        Ab = Ab,
        ldsc = scores[, seq_len(nA), drop = FALSE],
        ldCand = if (!is.null(Cb)) scores[, -seq_len(nA), drop = FALSE],
        y = z[idx]^2 - 1
    )
}

# Omega = S o S, with S = rawNtau * R Da R + R and Da = diag(rowSums(A)),
# reduced to its leading eigenpairs carrying `cutpoint` of the eigenvalue
# mass (eigen.cut in the original gLDSC implementation): Omega^-1 is applied
# as V diag(1 / values) V' and never formed.
.gldscOmegaEigen <- function(ld, Ab, rawNtau, cutpoint = 0.99) {
    S <- .gldscS(ld, rowSums(Ab), rawNtau)
    ei <- eigen(S^2, symmetric = TRUE)
    k <- which(cumsum(ei$values) / sum(ei$values) >= cutpoint)[1]
    list(
        vectors = ei$vectors[, seq_len(k), drop = FALSE],
        values = ei$values[seq_len(k)]
    )
}

# S = rawNtau * R Da R + R for a block's LD source, with Da = diag(d). From a
# stored R, R Da R = crossprod(R sqrt(Da)). From a factor F, with
# R = F'F / (m - 1), S = F' (rawNtau / (m - 1) * (F Da F') F + F) / (m - 1):
# only S is formed, which Omega needs whole, never the block's R.
.gldscS <- function(ld, d, rawNtau) {
    if (!is.null(ld$R)) {
        return(rawNtau * crossprod(ld$R * sqrt(d)) + ld$R)
    }
    ldFactor <- ld$F
    m1 <- nrow(ldFactor) - 1
    gram <- tcrossprod(sweep(ldFactor, 2L, sqrt(d), "*"))
    crossprod(ldFactor, rawNtau / m1 * (gram %*% ldFactor) + ldFactor) / m1
}

# One block's contribution to the GLS normal equations, L = X'Omega^-1 X and
# R = X'Omega^-1 y with X = [1 | l], plus its candidate-score pieces. Omega
# is formed for the block, decomposed and dropped, so only one block's is
# ever in memory.
.gldscBlockContrib <- function(block, design, panel, rawNtau) {
    om <- .gldscOmegaEigen(.ldStatBlockLd(block, panel), design$Ab, rawNtau)
    X <- cbind(1, design$ldsc)
    P <- crossprod(om$vectors, X)
    Py <- crossprod(om$vectors, design$y)
    Pw <- P / om$values
    list(
        L = crossprod(Pw, P),
        R = as.vector(crossprod(Pw, Py)),
        p = length(block$idx),
        cand = .gldscCandidatePieces(om, design$ldCand, P, Py)
    )
}

# A block's candidate-score pieces. The score is taken against the FITTED
# baseline residual, U_c = l_c' Omega^-1 (y - X coef), which is linear in
# coef: keeping l_c' Omega^-1 y and l_c' Omega^-1 X lets the block's score be
# formed once the fit is known, without revisiting the block. NULL without
# candidates.
.gldscCandidatePieces <- function(om, ldCand, P, Py) {
    if (is.null(ldCand)) {
        return(NULL)
    }
    Qw <- crossprod(om$vectors, ldCand) / om$values
    list(
        Uy = as.vector(crossprod(Qw, Py)),
        UX = crossprod(Qw, P)
    )
}

# GLS estimate from accumulated normal equations: tau = solve(L, R)/N (the
# leading entry is the confounding intercept), partitioned h2 = A'A tau.
.gldscEstimate <- function(left, right, A, N) {
    ridge <- 1e-8 * mean(abs(diag(left)))
    coef <- as.vector(solve(left + diag(ridge, nrow(left)), right))
    tauAll <- coef / N
    tau <- tauAll[-1]
    list(
        tau = tau,
        coef = coef,
        h2 = sum(colSums(A) * tau),
        intercept = N * tauAll[1] + 1
    )
}

# Delete-one-block jackknife: subtract each block's Omega-weighted contribution
# from the accumulated normal equations (so folds stay GLS-weighted). Jackknife
# variance = var * (nB-1)^2 / nB (Xiong et al.).
.gldscJackknife <- function(contrib, left, right, A, N) {
    loo <- map(
        contrib,
        .gldscLooEstimate,
        left = left,
        right = right,
        A = A,
        N = N
    )
    jkCoefs <- exec(rbind, !!!map(loo, "tau"))
    list(
        h2Se = .jackknifeSeFromLoo(map_dbl(loo, "h2")),
        intSe = .jackknifeSeFromLoo(map_dbl(loo, "intercept")),
        annotationJackknifeCoefs = jkCoefs
    )
}

# Local heritability for one LD block: a per-block GLS on [1 | base LD score],
# whose normal equations are the leading 2 x 2 of the block's own;
# h2Local = tau_base * (block SNPs), SE from the GLS coefficient covariance.
.gldscLocalBlock <- function(cb, b, N) {
    L <- cb$L[1:2, 1:2]
    coef <- solve(L, cb$R[1:2])
    covTau <- (solve(L) / N^2)[2, 2]
    tibble(
        blockId = b,
        h2Local = coef[2] / N * cb$p,
        h2LocalSe = sqrt(max(covTau, 0)) * cb$p
    )
}

# Per-block local heritability (NULL unless requested).
.gldscLocal <- function(contrib, N) {
    bind_rows(
        map(
            seq_along(contrib),
            .gldscLocalRow,
            contrib = contrib,
            N = N
        )
    )
}

# One block's score for each candidate, against the fitted baseline GLS
# residual r = y - X coef, with the candidate's LD scores l_c first
# residualized on the design: U_c = (l_c - X beta_c)' Omega^-1 r, beta_c =
# (X' Omega^-1 X)^-1 X' Omega^-1 l_c. Summed over blocks the adjustment
# vanishes (X' Omega^-1 r = 0 at the fit); block by block it does not, and
# it is what makes the blocks' scores the efficient ones the jackknife
# variance is taken over.
.gldscCandidateBlock <- function(cb, coef, beta) {
    U <- cb$cand$Uy - as.vector(cb$cand$UX %*% coef)
    Xr <- cb$R - as.vector(cb$L %*% coef)
    U - as.vector(crossprod(beta, Xr))
}

# Candidate score statistics at the real fitted coefficients (not a refit),
# calibrated by the blocks' own spread (.h2BlockScoreStats()) rather than by
# the model information l_c' Omega^-1 l_c, which ignored the estimated
# baseline and trusted the plug-in Omega.
.gldscScoreStats <- function(contrib, coef) {
    left <- reduce(map(contrib, "L"), `+`)
    lX <- reduce(map(contrib, function(cb) t(cb$cand$UX)), `+`)
    ridge <- 1e-8 * mean(abs(diag(left)))
    beta <- solve(left + diag(ridge, nrow(left)), lX)
    Ub <- exec(
        rbind,
        !!!map(contrib, .gldscCandidateBlock, coef = coef, beta = beta)
    )
    .h2BlockScoreStats(Ub)
}

# Each LD block's efficient score for every candidate (a blocks x candidates
# matrix). An observation i -- a variant, or an eigen-direction -- adds
# gTilde_ic e_i, where e_i is its score residual (weighted residual, or the
# likelihood's), gC[i, c] the derivative of its modelled mean or variance in
# the candidate's parameter, and gTilde that derivative residualized on the
# nuisance derivatives G under the information weights `info`:
# gTilde = gC - G (G' info G)^-1 G' info gC. Residualizing is what makes the
# test account for the baseline and intercept having been estimated -- a
# candidate overlapping the baseline otherwise picks up the baseline's fitted
# signal -- and it needs no refit. `blockOf` gives each observation's block.
# @noRd
.h2EfficientScoreBlocks <- function(gC, G, info, e, blockOf) {
    gC <- as.matrix(gC)
    gtg <- crossprod(G, G * info)
    ridge <- 1e-8 * mean(abs(diag(gtg)))
    beta <- solve(gtg + diag(ridge, ncol(G)), crossprod(G, gC * info))
    rowsum((gC - G %*% beta) * e, blockOf)
}

# Candidate score statistics from each LD block's score (a blocks x
# candidates matrix Ub): z_c = U_c / sqrt(V_c) with U_c = sum_b Ub[b, c] and
# V_c = nB / (nB - 1) sum_b (Ub[b, c] - U_c / nB)^2, the delete-one-block
# jackknife variance of the sum. Taking the variance from the blocks makes the
# test valid however wrong the regression weights are -- correlated
# residuals, a misspecified Omega -- as the heritability's own jackknife SE
# is. The candidates' correlation is that of their block scores.
# @noRd
.h2BlockScoreStats <- function(Ub) {
    Ub <- as.matrix(Ub)
    nB <- nrow(Ub)
    U <- colSums(Ub)
    centered <- sweep(Ub, 2L, U / nB)
    V <- nB / (nB - 1) * colSums(centered^2)
    R <- if (ncol(Ub) == 1L) {
        matrix(1, 1, 1)
    } else {
        Rc <- stats::cov2cor(crossprod(centered))
        replace(Rc, is.na(Rc), 0)
    }
    list(scoreZ = unname(U / sqrt(pmax(V, .Machine$double.xmin))), R = R)
}

gldscUnivariate <- function(
    z,
    n,
    ldRef,
    annotations = NULL,
    local = FALSE
) {
    if (!.ldScoreHasLd(ldRef)) {
        abort(c(
            "g-LDSC needs each LD block's LD, and this LdScore carries none.",
            i = str_c(
                "Build it from a genotype panel, or from LdMatrix blocks ",
                "with `keepLdMatrices = TRUE`."
            )
        ))
    }
    M <- length(ldRef)
    A <- .gldscAnnotMatrix(annotations, M)
    panel <- ldReference(ldRef)
    blocks <- .ldScoreBlocks(ldRef)
    # Two passes over the blocks: Omega depends on rawNtau, which needs every
    # block's design first. A block read from a panel is read once per pass.
    designs <- map(
        blocks,
        .gldscBlockDesign,
        panel = panel,
        z = z,
        A = A,
        candMat = .h2CandidateMat(annotations)
    )
    rawNtau <- sum(list_c(map(designs, "y"))) /
        sum(list_c(map(designs, "ldsc")))
    contrib <- map2(
        blocks,
        designs,
        .gldscBlockContrib,
        panel = panel,
        rawNtau = rawNtau
    )
    left <- reduce(map(contrib, "L"), `+`)
    right <- reduce(map(contrib, "R"), `+`)
    fit <- .gldscEstimate(left, right, A, n)
    jk <- .gldscJackknife(contrib, left, right, A, n)
    .gldscResult(fit, jk, annotations, contrib, local, n, M, A)
}

# The candidate-annotation score statistics, when candidates were supplied.
# NULL otherwise.
# @noRd
.gldscCandidateScores <- function(annotations, contrib, fit) {
    if (is.null(.h2CandidateMat(annotations))) {
        return(NULL)
    }
    sc <- .gldscScoreStats(contrib, fit$coef)
    list(
        z = sc$scoreZ,
        R = sc$R,
        annotationNames = SummarizedExperiment::colData(
            testedAnnotations(annotations)
        )$name
    )
}

# Assemble the gldscUnivariate return (enrichment, local, candidate scores).
.gldscResult <- function(fit, jk, annotations, contrib, local, n, M, A) {
    .h2FitResult(
        h2 = fit$h2,
        h2Se = jk$h2Se,
        intercept = fit$intercept,
        interceptSe = jk$intSe,
        A = if (is.null(.h2BaselineMat(annotations))) NULL else A,
        tau = fit$tau,
        tauBlocks = jk$annotationJackknifeCoefs,
        localDf = if (local) .gldscLocal(contrib, n) else NULL,
        scoreStats = .gldscCandidateScores(annotations, contrib, fit)
    )
}


#' @title HDL/sHDL: High-Definition Likelihood
#' @description Estimate heritability by maximizing the eigenvalue-based
#'   Gaussian likelihood of GWAS z-scores (Ning et al. 2020). Within each LD
#'   block the z-scores are rotated into the eigenbasis; rotated coordinate
#'   \eqn{i} (in \eqn{b^* = V'z/\sqrt{N}} units) has modeled variance
#'   \deqn{h2/M \cdot \lambda_i^2 - h2 \cdot \lambda_i / N_{ref} +
#'     a \cdot \lambda_i / N,}
#'   where the middle term corrects for the finite LD reference panel of size
#'   \eqn{N_{ref}}. The stratified extension replaces the constant genetic
#'   density \eqn{h2/M} with a per-eigenvalue density
#'   \eqn{s_i = \sum_a \tau_a \ell_{ia}} using the stratified eigenvalue scores
#'   \eqn{\ell_{ia} = \sum_j A_{ja} V_{ji}^2}; it reduces exactly to the
#'   univariate model when there is a single genome-wide annotation.
#' @name pecotmr-h2-hdl
#' @keywords internal
#' @references
#'   Ning Z, Pawitan Y, Shen X (2020). High-definition likelihood
#'   inference of genetic correlations across human complex traits.
#'   Nat Genet, 52:859-864.
#' @importFrom stats optim optimize cor
NULL

# =============================================================================
# HDL likelihood + design
# =============================================================================

# Unified HDL negative log-likelihood in b* = V'(z/sqrt(N)) units. Modeled
# variance of rotated coordinate i is
#   lamh2_i = s_i lam_i^2 - lam_i h2Total / nRef + int lam_i / n,
# s_i = (ldAnnotScaled %*% h2a)_i and h2Total = sum(h2a). With a single 1/M
# score column this is exactly llfun in the original HDL implementation
# (Ning et al. 2020). param = c(h2_1, .., h2_nTau, int).
# `lim` is the original HDL's floor on the modelled per-direction variance
# (llfun.R uses exp(-10)). It has to match: with 40% of a real reference's
# eigen-directions near zero the floor is active for most of them, and the
# `bstar^2 / lamh2` term it feeds dominates the likelihood -- a smaller floor
# changes the objective by orders of magnitude and moves the optimum.
# @noRd
.hdlNll <- function(
    param,
    n,
    nSamples,
    lam,
    bstar,
    ldAnnotScaled,
    lim = exp(-10)
) {
    nTau <- ncol(ldAnnotScaled)
    h2a <- param[seq_len(nTau)]
    int <- param[nTau + 1L]
    lamh2 <- pmax(
        as.vector(ldAnnotScaled %*% h2a) *
            lam^2 -
            lam * sum(h2a) / nSamples +
            int * lam / n,
        lim
    )
    sum(log(lamh2)) + sum(bstar^2 / lamh2)
}

# Per-block HDL quantities: eigenvalues, rotated N-normalized effects
# bstar = V'(z/sqrt(n)), stratified eigenvalue scores l_{ia} (NULL when
# unstratified) and the candidate annotations' scores, from the block's
# decomposition. The zero-eigenvalue directions a decomposition from a
# panel's factor does not list are added with lam = 0.
.hdlBlockData <- function(block, eigenRef, z, n, baselineMat, candMat) {
    dec <- .ldEigenBlockAt(block, eigenRef)
    zb <- z[dec$snpIdx] / sqrt(n)
    bstar <- as.vector(crossprod(dec$vectors, zb))
    list(
        lam = c(dec$values, numeric(dec$nNull)),
        bstar = c(bstar, .hdlNullBstar(zb, bstar, dec$nNull)),
        ldAnnot = .h2EigenAnnot(dec, baselineMat),
        ldCand = .h2EigenAnnot(dec, candMat)
    )
}

# The rotated effects along a block's unlisted zero-eigenvalue directions.
# HDL gives every such direction the same modelled variance, its floor, so
# the likelihood sees only their total sum of squares: the part of z outside
# the listed directions, spread evenly over them, leaves every estimate as a
# full decomposition would give it.
.hdlNullBstar <- function(zb, bstar, nNull) {
    if (nNull == 0L) {
        return(numeric())
    }
    rep(sqrt(max(sum(zb^2) - sum(bstar^2), 0) / nNull), nNull)
}

# Genome-wide design: eigenvalues, rotated effects, block ids, the M_a-scaled
# eigenvalue-score matrix, and M_a. Univariate uses a single all-ones score
# column (l_{i,base} = 1) scaled by total M.
.hdlDesign <- function(blockData, M, baselineMat) {
    lam <- list_c(map(blockData, "lam"))
    bstar <- list_c(map(blockData, "bstar"))
    blockId <- rep(
        seq_along(blockData),
        lengths(map(blockData, "lam"))
    )
    if (is.null(baselineMat)) {
        return(list(
            lam = lam,
            bstar = bstar,
            blockId = blockId,
            ldAnnotScaled = matrix(1 / M, length(lam), 1),
            M_a = M
        ))
    }
    M_a <- colSums(baselineMat)
    ldAnnotList <- map(blockData, "ldAnnot")
    ldAnnot <- exec(rbind, !!!ldAnnotList)
    list(
        lam = lam,
        bstar = bstar,
        blockId = blockId,
        ldAnnotScaled = sweep(ldAnnot, 2, M_a, "/"),
        M_a = M_a,
        baseAdded = !is_in(1L, attr(baselineMat, "userCols"))
    )
}

# Fit HDL by L-BFGS-B over (h2_1..h2_nTau, int); bounds h2_a in [0,1], int in
# [0,10] (Ning et al. 2020). Parametrizing by per-annotation h2 (scale ~1)
# keeps the L-BFGS-B finite differences well conditioned.
.hdlFit <- function(design, n, nSamples, warnOnBound = TRUE) {
    nTau <- ncol(design$ldAnnotScaled)
    # Beside an added base column an annotation's parameter is its excess
    # over the genome-wide rate, which a depleted annotation makes negative;
    # only the base, and the columns of a baseline that spans the genome on
    # its own, are heritabilities bounded below by 0.
    lower <- c(
        if (isTRUE(design$baseAdded)) {
            c(0, rep(-1, nTau - 1L))
        } else {
            rep(0, nTau)
        },
        0
    )
    upper <- c(rep(1, nTau), 10)
    opt <- optim(
        c(rep(0.1, nTau), 1),
        .hdlNll,
        n = n,
        nSamples = nSamples,
        lam = design$lam,
        bstar = design$bstar,
        ldAnnotScaled = design$ldAnnotScaled,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper
    )
    if (warnOnBound) {
        .hdlWarnIfAtBound(opt$par, upper, n, nSamples)
    }
    h2a <- opt$par[seq_len(nTau)]
    list(h2a = h2a, int = opt$par[nTau + 1L], h2 = sum(h2a))
}

# HDL's likelihood carries a finite-reference term, -lam * h2 / nRef, that
# subtracts from the modelled variance. Once the reference is small enough
# relative to the GWAS, that term outgrows the signal term, the modelled
# variance is driven negative across most directions, and the optimiser walks
# the intercept out to its upper bound. The heritability that comes back then
# reflects the boundary rather than the data -- on the bundled reference,
# nRef = 1000 against N = 1e5 returns 0.08 for a true 0.4 with the intercept
# pinned at 10.
#
# It is a silent failure otherwise, and a detectable one: say so.
#
# Only the UPPER bounds are diagnostic. An intercept pinned at 0 simply means
# no inflation was detected and leaves the heritability sound -- it happens
# routinely on a well-sized reference -- and a heritability of 0 is a
# legitimate estimate for a null trait. Warning on those would fire on
# healthy fits and train the reader to ignore the message.
# @noRd
.hdlWarnIfAtBound <- function(par, upper, n, nSamples) {
    tol <- 1e-6
    atUpper <- par >= upper - tol
    if (!any(atUpper)) {
        return(invisible(FALSE))
    }
    which <- if (atUpper[[length(atUpper)]]) {
        "intercept"
    } else {
        "heritability"
    }
    warn(glue(
        "HDL: the fitted {which} sits at the top of its range, so the ",
        "estimate reflects that bound rather than the data. This usually ",
        "means the LD reference is too small for the GWAS -- the ",
        "finite-reference correction scales as 1/nSamples, and here ",
        "nSamples = {nSamples} against N = {round(n)}. Compare the other ",
        "estimators, or ",
        "use a larger reference panel."
    ))
    invisible(TRUE)
}

# Baseline modeled variance lamh2_i at the fitted parameters (for local h2 and
# the sHDL candidate score test).
.hdlBaselineVar <- function(design, ft, n, nSamples, lim = exp(-10)) {
    s <- as.vector(design$ldAnnotScaled %*% ft$h2a)
    pmax(
        s *
            design$lam^2 -
            design$lam * ft$h2 / nSamples +
            ft$int * design$lam / n,
        lim
    )
}

# =============================================================================
# HDL standard errors (delete-one-block jackknife) + enrichment inputs
# =============================================================================

# Delete-one-block jackknife (Ning et al. 2020): refit leaving each block out,
# SE = sqrt(mean((jack - mean)^2) * (nBlocks - 1)). Returns per-annotation and
# total-h2 SEs, the intercept SE, and the LOO per-annotation blocks.
.hdlJackknife <- function(blockData, M, baselineMat, n, nSamples) {
    nBlocks <- length(blockData)
    loo <- map(
        seq_len(nBlocks),
        .hdlLooFit,
        blockData = blockData,
        M = M,
        baselineMat = baselineMat,
        n = n,
        nSamples = nSamples
    )
    h2aList <- map(loo, "h2a")
    h2aBlocks <- exec(rbind, !!!h2aList)
    intLoo <- map_dbl(loo, "int")
    list(
        h2aSe = apply(h2aBlocks, 2, .jackknifeSeFromLoo),
        h2Se = .jackknifeSeFromLoo(rowSums(h2aBlocks)),
        intSe = .jackknifeSeFromLoo(intLoo),
        h2aBlocks = h2aBlocks
    )
}

# =============================================================================
# HDL local heritability + sHDL candidate score test
# =============================================================================

# Local h2 for one block: MLE of the block's genome-wide-scaled h2 with the
# intercept fixed at the global fit; SE from the Gaussian-variance Fisher
# information sum((dlamh2/dh2)^2 / (2 lamh2^2)).
.hdlLocalBlock <- function(bd, b, n, nSamples, int, M, lim = exp(-10)) {
    # `intercept`, not `int`: optimize() would partially match `int` to its
    # own `interval` argument.
    h2 <- optimize(
        .hdlLocalNll,
        c(0, 1),
        bd = bd,
        M = M,
        nSamples = nSamples,
        intercept = int,
        n = n,
        lim = lim
    )$minimum
    lamh2 <- .hdlLocalVar(h2, bd, M, nSamples, int, n, lim)
    g <- bd$lam^2 / M - bd$lam / nSamples
    info <- sum(g^2 / (2 * lamh2^2))
    tibble(
        blockId = b,
        h2Local = h2,
        h2LocalSe = 1 / sqrt(max(info, 1e-10))
    )
}

# Modeled variance of one block's bstar^2 at local heritability `h2`, the
# intercept fixed at the global fit; floored at `lim`.
# @noRd
.hdlLocalVar <- function(h2, bd, M, nSamples, intercept, n, lim) {
    pmax(
        h2 / M * bd$lam^2 - h2 * bd$lam / nSamples + intercept * bd$lam / n,
        lim
    )
}

# Gaussian negative log-likelihood (up to a constant) of one block's bstar at
# local heritability `h2` -- what .hdlLocalBlock() minimizes.
# @noRd
.hdlLocalNll <- function(h2, bd, M, nSamples, intercept, n, lim) {
    lamh2 <- .hdlLocalVar(h2, bd, M, nSamples, intercept, n, lim)
    sum(log(lamh2) + bd$bstar^2 / lamh2)
}

# Per-block local heritability (NULL unless requested).
.hdlLocal <- function(blockData, n, nSamples, int, M) {
    bind_rows(
        map(
            seq_along(blockData),
            .hdlLocalRow,
            blockData = blockData,
            n = n,
            nSamples = nSamples,
            int = int,
            M = M
        )
    )
}

# sHDL candidate-annotation score statistics. Each eigen-direction's rotated
# effect has modelled variance v_i = lamh2_i, and the likelihood's score for
# a candidate's heritability is sum_i g_ic (bstar_i^2 / v_i^2 - 1 / v_i),
# with g_ic = dv_i/dh2_c = l_ic lam_i^2 / M_c - lam_i / nPanel and
# information weight 2 / v_i^2. The baseline heritabilities and the intercept
# are the nuisance (.h2EfficientScoreBlocks()), with derivatives
# l_ia lam_i^2 / M_a - lam_i / nPanel and lam_i / n.
.shdlStratified <- function(blockData, nPanel, annotations, design, lamh2, n) {
    candidate <- testedAnnotations(annotations)
    candMat <- assay(candidate, "annotations")
    if (ncol(candMat) == 0) {
        return(list(enrichment = NULL, scoreStats = NULL))
    }
    lam <- design$lam
    gradient <- function(scaled) scaled * lam^2 - lam / nPanel
    ldCand <- exec(rbind, !!!map(blockData, "ldCand"))
    gC <- gradient(sweep(ldCand, 2L, colSums(candMat), "/"))
    G <- cbind(gradient(design$ldAnnotScaled), lam / n)
    Ub <- .h2EfficientScoreBlocks(
        gC,
        G,
        2 / lamh2^2,
        design$bstar^2 / lamh2^2 - 1 / lamh2,
        design$blockId
    )
    sc <- .h2BlockScoreStats(Ub)
    .buildEnrichmentResult(
        SummarizedExperiment::colData(candidate),
        sc$scoreZ,
        sc$R
    )
}

# =============================================================================
# hdlUnivariate -- orchestrator
# =============================================================================

hdlUnivariate <- function(
    z,
    n,
    eigenRef,
    annotations = NULL,
    local = FALSE
) {
    M <- length(eigenRef)
    nPanel <- ldPanelNSamples(eigenRef)
    baselineMat <- .h2DesignMat(annotations)
    # One pass over the blocks, each decomposed only while its data are taken.
    blockData <- map(
        .ldEigenBlocks(eigenRef),
        .hdlBlockData,
        eigenRef = eigenRef,
        z = z,
        n = n,
        baselineMat = baselineMat,
        candMat = .h2CandidateMat(annotations)
    )
    design <- .hdlDesign(blockData, M, baselineMat)
    ft <- .hdlFit(design, n, nPanel)
    jk <- .hdlJackknife(blockData, M, baselineMat, n, nPanel)
    localDf <- if (local) .hdlLocal(blockData, n, nPanel, ft$int, M) else NULL
    .h2FitResult(
        h2 = ft$h2,
        h2Se = jk$h2Se,
        intercept = ft$int,
        interceptSe = jk$intSe,
        A = baselineMat,
        tau = ft$h2a / design$M_a,
        tauBlocks = sweep(jk$h2aBlocks, 2L, design$M_a, "/"),
        localDf = localDf,
        scoreStats = .hdlCandidateScores(
            blockData,
            nPanel,
            annotations,
            design,
            ft,
            n
        )
    )
}

# sHDL's score statistics for the candidate annotations, under the variances
# the baseline fit implies; NULL when there are none to test.
# @noRd
.hdlCandidateScores <- function(blockData, nPanel, annotations, design, ft, n) {
    if (is.null(annotations)) {
        return(NULL)
    }
    lamh2 <- .hdlBaselineVar(design, ft, n, nPanel)
    .shdlStratified(blockData, nPanel, annotations, design, lamh2, n)$scoreStats
}


#' @title Heritability Estimation Entry Points and Converters
#' @description Top-level entry point for heritability estimation, LD score
#'   computation methods, H2Estimate accessors, and a converter to bridge
#'   H2Estimate into the sldscWrapper.R postprocessing pipeline.
#' @name pecotmr-h2-wrappers
#' @keywords internal
#' @include AllGenerics.R
#' @importFrom stats median
NULL

# =============================================================================
# S-LDSC (Finucane et al. 2015) -- stratified LD score regression
# -----------------------------------------------------------------------------
# The original method the other three descend from. It regresses per-variant
# chi-square on annotation-stratified LD scores:
#
#   E[chi2_j] = N * sum_C tau_C * l(j, C) + N*a + 1
#   with l(j, C) = sum_{k in C} r2_jk
#
# Unlike g-LDSC it needs only the per-variant scores, not the per-block LD
# matrices, so an unstratified run works off a scores-only LdScore. The
# stratified run needs the matrices, since the per-annotation scores are
# formed from them exactly as g-LDSC forms its own.
# =============================================================================

#' @title S-LDSC: Stratified LD Score Regression
#' @description Estimate heritability and functional enrichment by weighted
#'   least squares on LD scores (Finucane et al. 2015) -- the original method
#'   the other three estimators descend from. Per-variant \eqn{\chi^2_j} is
#'   regressed on annotation-stratified LD scores
#'   \eqn{\ell_{jC} = \sum_{k \in C} r^2_{jk}}, each \eqn{r^2} being ldsc's
#'   unbiased estimate \code{r^2 - (1 - r^2) / (n - 2)} at the reference
#'   panel size \code{n} (see \code{\link{buildLdScore}}), with a
#'   confounding intercept, fitting
#'   \eqn{E[\chi^2_j] = N \sum_C \tau_C \ell_{jC} + Na + 1}.
#'
#'   Weights are LDSC's two multiplied factors, both floored at 1: the
#'   heteroskedasticity term \eqn{1 / (2 E[\chi^2_j]^2)}, iterated on the
#'   fitted values, and the overcounting term \eqn{1 / \ell_j} that stops
#'   variants in strong LD from being counted many times. Standard errors come
#'   from a delete-one-block jackknife over the LD blocks.
#'
#'   Two routes are taken, matching \code{ldsc.py} -- the split is on
#'   annotation count, not on preference:
#'   \describe{
#'     \item{unstratified}{the two-step estimator: the intercept is fit on
#'       the variants with \eqn{\chi^2 < 30}, then held fixed while the
#'       slope is fit over all of them.}
#'     \item{stratified}{no two-step -- the original \code{ldsc} refuses
#'       one on a partitioned design. Variants with \eqn{\chi^2 >}
#'       \code{max(0.001
#'       N, 80)} are dropped, the design is weighted once, and a single
#'       weighted least-squares fit follows.}
#'   }
#'
#'   The intercept is reported as estimated and is not floored at 1: the
#'   model implies at least 1, but a value below it is informative rather
#'   than an error to clamp away, and clamping would move the heritability
#'   too.
#'
#'   Estimates agree with the original \code{ldsc} implementation to
#'   roughly 1e-6 relative; the residual is a small ridge in the weighted
#'   solve, which the original has no equivalent of and which exists so an
#'   exactly singular design returns a value rather than \code{NA}.
#'
#'   Unlike \code{\link{pecotmr-h2-gldsc}} this needs only the per-variant
#'   scores when unstratified, so it runs on a scores-only
#'   \code{\link{LdScore}}; the stratified fit forms per-annotation scores
#'   from each block's LD, so it needs an \code{LdScore} that carries it:
#'   one built from a genotype panel, or from \code{LdMatrix} blocks with
#'   their matrices kept.
#'
#'   \strong{Scale.} LD score regression is a genome-scale method -- real
#'   analyses regress over roughly a million variants. On a reference of a
#'   few thousand the intercept is not identifiable (an unweighted OLS on
#'   such data returns an intercept several times its true value), and the
#'   heritability estimate is attenuated as a result. Treat small-reference
#'   estimates as illustrative. \code{\link{pecotmr-h2-lder}} is much less
#'   sensitive to reference size and is the better choice there.
#' @name pecotmr-h2-sldsc
#' @keywords internal
#' @references
#'   Finucane HK, Bulik-Sullivan B, Gusev A, et al. (2015). Partitioning
#'   heritability by functional annotation using genome-wide association
#'   summary statistics. Nature Genetics 47(11):1228-1235.
NULL

# @noRd
sldscUnivariate <- function(
    z,
    n,
    ldRef,
    annotations = NULL,
    local = FALSE,
    nIter = 2L
) {
    M <- length(ldRef)
    chi2 <- z^2
    # Same design g-LDSC uses (.h2DesignMat()), so `crossprod(A)` is the
    # annotation overlap matrix.
    A <- .gldscAnnotMatrix(annotations, M)
    scores <- .sldscScoreMatrix(ldRef, A)
    # The weights read the stored base score whatever the design holds.
    baseScore <- as.vector(mcols(ldRef)$ldScores[, 1L])
    blockIdx <- .sldscBlockIndex(ldRef)
    fit <- .sldscFit(chi2, scores, baseScore, n, A, nIter)
    jk <- .sldscJackknife(chi2, scores, baseScore, n, A, nIter, blockIdx)
    scoreStats <- .sldscCandidateScores(fit, ldRef, annotations, blockIdx)
    .h2FitResult(
        h2 = fit$h2,
        h2Se = jk$h2Se,
        intercept = fit$intercept,
        interceptSe = jk$intSe,
        A = if (is.null(.h2BaselineMat(annotations))) NULL else A,
        tau = fit$tau,
        tauBlocks = jk$annotationJackknifeCoefs,
        localDf = if (local) .sldscLocal(chi2, baseScore, fit) else NULL,
        scoreStats = scoreStats
    )
}

# The candidate annotations' score statistics, NULL without candidates: each
# candidate's stratified LD score is scored against the fitted residual in
# the fit's final weights, with the fitted design as nuisance
# (.h2EfficientScoreBlocks()).
# @noRd
.sldscCandidateScores <- function(fit, ldRef, annotations, blockIdx) {
    candMat <- .h2CandidateMat(annotations)
    if (is.null(candMat)) {
        return(NULL)
    }
    if (!.ldScoreHasLd(ldRef)) {
        abort(str_c(
            "S-LDSC's candidate-annotation score test needs each LD block's ",
            "LD to form the candidates' LD scores. Build the reference from ",
            "a genotype panel, or with `buildLdScore(keepLdMatrices = TRUE)`."
        ))
    }
    parts <- attr(fit, "scoreParts")
    ldCand <- .h2StratifiedScores(ldRef, candMat)[parts$keep, , drop = FALSE]
    blockOf <- rep(seq_along(blockIdx), lengths(blockIdx))[
        order(unlist(blockIdx, use.names = FALSE))
    ]
    Ub <- .h2EfficientScoreBlocks(
        ldCand,
        parts$design,
        parts$w,
        parts$w * parts$resid,
        blockOf[parts$keep]
    )
    sc <- .h2BlockScoreStats(Ub)
    list(
        z = sc$scoreZ,
        R = sc$R,
        annotationNames = SummarizedExperiment::colData(
            testedAnnotations(annotations)
        )$name
    )
}

# The design's LD-score columns, one per column of the annotation matrix `A`
# (.gldscAnnotMatrix()): the stored base score for the base column, and a
# stratified score for each baseline annotation. Stratified scores come from
# each block's LD, the same source g-LDSC uses, and take ldsc's unbiased r^2
# as the stored base score does: ldsc applies it to every annotation column.
# @noRd
.sldscScoreMatrix <- function(ldRef, A) {
    base <- matrix(
        as.vector(mcols(ldRef)$ldScores[, 1L]),
        ncol = 1L,
        dimnames = list(NULL, "base_l2")
    )
    hasBase <- identical(colnames(A)[[1L]], "base")
    baselineMat <- if (hasBase) A[, -1L, drop = FALSE] else A
    if (ncol(baselineMat) == 0L) {
        return(base)
    }
    if (!.ldScoreHasLd(ldRef)) {
        abort(str_c(
            "stratified S-LDSC needs each LD block's LD to form the ",
            "per-annotation LD scores. Build the reference from a genotype ",
            "panel, or with `buildLdScore(keepLdMatrices = TRUE)`."
        ))
    }
    stratified <- .h2StratifiedScores(ldRef, baselineMat)
    if (hasBase) cbind(base, stratified) else stratified
}

# Annotation-stratified LD scores sum_k r^2_jk a_k under ldsc's unbiased r^2,
# one column per annotation, from each block's LD: stored, or read from the
# genotype panel a block at a time.
# @noRd
.h2StratifiedScores <- function(ldRef, annotMat) {
    .h2ScatterBlockRows(
        matrix(0, length(ldRef), ncol(annotMat)),
        map(
            .ldScoreBlocks(ldRef),
            .h2BlockAnnotScores,
            panel = ldReference(ldRef),
            annotMat = annotMat,
            nRef = .ldscPanelN(ldRef)
        )
    )
}

# One block's per-annotation LD scores from its LD, summing ldsc's unbiased
# r^2 by its summed form.
# @noRd
.h2BlockAnnotScores <- function(block, panel, annotMat, nRef) {
    idx <- block$idx
    Ab <- annotMat[idx, , drop = FALSE]
    raw <- .ldSquaredProduct(.ldStatBlockLd(block, panel), Ab)
    list(idx = idx, values = .ldscUnbiasedSum(raw, colSums(Ab), nRef))
}

# Per-block variant indices: the LD blocks the estimators walk, which the
# jackknife leaves out one at a time.
# @noRd
.sldscBlockIndex <- function(ldRef) {
    map(.ldScoreBlocks(ldRef), "idx")
}

# `Hsq.weights` from the original ldsc, verbatim (ldscore/regressions.py):
#
#   hsq   = clamp(hsq, 0, 1)
#   ld    = max(ld, 1);  w_ld = max(w_ld, 1)
#   c     = hsq * N / M
#   het_w = 1 / (2 (intercept + c*ld)^2)
#   w     = het_w / w_ld
#
# `ld` is the LD score entering the model and `w_ld` the one computed over
# just the regression variants. pecotmr carries a single score, so both read
# the same column -- also the usual case with the original ldsc, where
# --w-ld-chr names the same file as --ref-ld-chr.
# @noRd
.sldscWeights <- function(ld, wLd, n, M, hsq, intercept) {
    hsq <- min(max(hsq, 0), 1)
    ld <- pmax(ld, 1)
    wLd <- pmax(wLd, 1)
    hetW <- 1 / (2 * (intercept + (hsq * n / M) * ld)^2)
    hetW / wLd
}

# LDSC's `aggregate`: the method-of-moments h2 that seeds the first weights.
#   h2 = M (mean(chi2) - intercept) / mean(N * ld)
# @noRd
.sldscAggregateH2 <- function(chi2, baseScore, n, M, intercept = 1) {
    M * (mean(chi2) - intercept) / mean(n * baseScore)
}

# The original ldsc takes two different routes and the split is on annotation
# count, not on preference (ldscore/sumstats.py):
#
#   univariate  -- two-step estimator, chi2 < 30 selects step 1, IRWLS twice
#   partitioned -- NO two-step; it raises rather than attempt one. Variants
#                  with chi2 > max(0.001*N, 80) are dropped outright, the
#                  design is weighted once by the initial weights, and a plain
#                  least-squares jackknife follows (its `old_weights` path).
#
# Both are reproduced so the estimates line up with `ldsc.py`.
# @noRd
.sldscFit <- function(chi2, scores, baseScore, n, A, nIter) {
    if (ncol(scores) == 1L) {
        return(.sldscFitUnivariate(chi2, scores, baseScore, n, A, nIter))
    }
    .sldscFitPartitioned(chi2, scores, baseScore, n, A)
}

# @noRd
.sldscFitUnivariate <- function(chi2, scores, baseScore, n, A, nIter) {
    M <- nrow(A)
    # Too few variants survive the chi2 cap to fit on: keep everything rather
    # than fit on a handful.
    capped <- chi2 < 30
    keep <- if (sum(capped) < 4L) rep(TRUE, length(chi2)) else capped
    # One set of initial weights, formed on the full data, reused by both
    # steps -- the original ldsc computes `initial_w` once and subsets it for
    # step 1.
    xTot <- rowSums(scores)
    aggH2 <- .sldscAggregateH2(chi2, xTot, n, M)
    w0 <- .sldscWeights(xTot, baseScore, n, M, aggH2, 1)
    step1 <- .sldscIrwls(
        chi2[keep],
        scores[keep, , drop = FALSE],
        baseScore[keep],
        n,
        M,
        nIter,
        freeIntercept = TRUE,
        w0 = w0[keep]
    )
    # Not floored at 1. The model implies an intercept of at least 1, but the
    # original ldsc reports whatever step 1 estimates -- a value below 1 is
    # informative (it says the data are inconsistent with the model, often
    # deflation) and clamping it silently would change the h2 as well.
    intercept <- step1$intercept
    step2 <- .sldscIrwls(
        chi2 - intercept,
        scores,
        baseScore,
        n,
        M,
        nIter,
        freeIntercept = FALSE,
        intercept = intercept,
        w0 = w0
    )
    fit <- .sldscCoefToFit(c(step2$coef, intercept), n, A)
    # What the candidate score test reads (.sldscCandidateScores()): the
    # step-2 regression, with the intercept step 1 estimated as a column so
    # the candidates are residualized on it too.
    `attr<-`(
        fit,
        "scoreParts",
        list(
            keep = rep(TRUE, length(chi2)),
            design = cbind(scores, 1),
            w = step2$w,
            resid = as.vector(chi2 - intercept - scores %*% step2$coef)
        )
    )
}

# @noRd
.sldscFitPartitioned <- function(chi2, scores, baseScore, n, A) {
    M <- nrow(A)
    # Too few variants survive the chi2 cap for the design's width: keep
    # everything rather than fit an underdetermined system.
    capped <- chi2 < max(0.001 * n, 80)
    keep <- if (sum(capped) < ncol(scores) + 2L) {
        rep(TRUE, length(chi2))
    } else {
        capped
    }
    chi2 <- chi2[keep]
    scores <- scores[keep, , drop = FALSE]
    baseScore <- baseScore[keep]
    # `x_tot` in the original ldsc: the row-sum of the design's annotation
    # columns, which coincides with the base score only when there is a
    # single column.
    xTot <- rowSums(scores)
    hsq <- .sldscAggregateH2(chi2, xTot, n, M)
    w <- .sldscWeights(xTot, baseScore, n, M, hsq, 1)
    design <- cbind(scores, 1)
    coef <- .sldscWlsCoef(design, chi2, w)
    `attr<-`(
        .sldscCoefToFit(coef, n, A),
        "scoreParts",
        list(
            keep = keep,
            design = design,
            w = w,
            resid = as.vector(chi2 - design %*% coef)
        )
    )
}

# Two IRWLS passes, matching the original ldsc's `for i in range(2)`.
# @noRd
.sldscIrwls <- function(
    y,
    scores,
    baseScore,
    n,
    M,
    nIter,
    freeIntercept,
    intercept = 1,
    w0 = NULL
) {
    design <- if (freeIntercept) cbind(scores, 1) else scores
    xTot <- rowSums(scores)
    if (is.null(w0)) {
        hsq <- .sldscAggregateH2(y, xTot, n, M)
        w0 <- .sldscWeights(xTot, baseScore, n, M, hsq, intercept)
    }
    w <- w0
    coef <- NULL
    for (i in seq_len(max(1L, nIter))) {
        coef <- .sldscWlsCoef(design, y, w)
        # The original ldsc builds the weight from the FIRST coefficient and
        # `x_tot` alone -- a univariate approximation even for a partitioned
        # design -- rather than from the full fitted values.
        hsq <- M * coef[[1L]] / n
        if (freeIntercept) {
            intercept <- max(coef[[length(coef)]], 0)
        }
        w <- .sldscWeights(xTot, baseScore, n, M, hsq, intercept)
    }
    # The original ldsc's loop updates the weights `nIter` times and then
    # solves once more with the converged weights -- three solves for its two
    # passes, not two. Stopping at the loop's last solve leaves the final
    # weight update unused and shifts the estimate by around a percent.
    coef <- .sldscWlsCoef(design, y, w)
    list(
        coef = coef[seq_len(ncol(scores))],
        intercept = if (freeIntercept) coef[[length(coef)]] else intercept,
        w = w
    )
}


# Solved through ridge-regularised normal equations rather than lm.fit. The
# base LD-score column is the sum of the stratified ones whenever the baseline
# annotations partition the genome, which makes the design exactly singular;
# lm.fit then returns NA coefficients that poison the next IRWLS weight. The
# ridge is the same negligible one g-LDSC uses for the same reason.
# @noRd
.sldscWlsCoef <- function(design, y, w) {
    sw <- sqrt(pmax(w, 0))
    x <- design * sw
    xtx <- crossprod(x)
    ridge <- 1e-8 * mean(abs(diag(xtx)))
    as.vector(solve(xtx + diag(ridge, ncol(xtx)), crossprod(x, y * sw)))
}

# Regressing chi2 on raw LD scores gives coef_C = N * tau_C, so tau = coef/N.
# Partitioned heritabilities are A'A tau, not tau * M_C: annotations overlap
# in general (the all-ones base column overlaps every one of them), and the
# overlap matrix is what turns per-annotation effects into per-annotation
# heritability. Its leading entry is the total, since A[, 1] is all ones.
# @noRd
.sldscCoefToFit <- function(coef, n, A) {
    nS <- ncol(A)
    tau <- coef[seq_len(nS)] / n
    list(
        tau = tau,
        h2 = sum(colSums(A) * tau),
        intercept = coef[[nS + 1L]],
        coef = coef
    )
}

# Delete-one-block jackknife over the LD blocks.
# @noRd
.sldscJackknife <- function(chi2, scores, baseScore, n, A, nIter, blockIdx) {
    loo <- map(
        blockIdx,
        .sldscLooFit,
        chi2 = chi2,
        scores = scores,
        baseScore = baseScore,
        n = n,
        A = A,
        nIter = nIter
    )
    list(
        h2Se = .jackknifeSeFromLoo(map_dbl(loo, "h2")),
        intSe = .jackknifeSeFromLoo(map_dbl(loo, "intercept")),
        annotationJackknifeCoefs = exec(rbind, !!!map(loo, "tau"))
    )
}

# The left-out block is dropped from the regression rows, but `A` stays whole:
# the annotation sizes are genome-wide constants, not per-fold quantities.
# @noRd
.sldscLooFit <- function(idx, chi2, scores, baseScore, n, A, nIter) {
    keep <- setdiff(seq_along(chi2), idx)
    .sldscFit(
        chi2[keep],
        scores[keep, , drop = FALSE],
        baseScore[keep],
        n,
        A,
        nIter
    )
}

# Per-variant local h2 contribution: the fitted non-intercept signal spread
# over the variant's own LD score.
# @noRd
.sldscLocal <- function(chi2, baseScore, fit) {
    tibble(
        variantIdx = seq_along(chi2),
        chi2 = chi2,
        ldScore = baseScore,
        h2Local = fit$tau[[1L]] * pmax(baseScore, 0)
    )
}

# =============================================================================
# estimateH2 -- main dispatch
# =============================================================================

# Resolve which study's stats to operate on; defaults to the single entry.
.estimateH2ResolveStudy <- function(sumstats, studyName) {
    if (!is.null(studyName)) {
        return(studyName)
    }
    # Studies, not rows: a genome-wide study is one row per chromosome (or
    # per LD block), and its accessors stitch those rows back together.
    studies <- unique(as.character(sumstats$study))
    if (length(studies) != 1L) {
        msg <- glue(
            "`studyName` is required when the GwasSumStats has ",
            "{length(studies)} studies."
        )
        abort(msg)
    }
    studies
}

# Dispatch to the method-specific univariate estimator.
.estimateH2Dispatch <- function(
    method,
    z,
    n,
    ldRef,
    annotations,
    local,
    estimatorArgs = list()
) {
    # `estimatorArgs` rather than `...`: the estimators take different
    # trailing arguments (`nIter` for sldsc; the others take none today), so
    # an unknown name should fail here rather than at whichever estimator the
    # method token happens to select. lder/gldsc/hdl used to advertise a
    # ridge `lambda` here that nothing read -- the original LDER and HDL
    # implementations have no such penalty, and the function behind it had no
    # callers.
    base <- list(z, n, ldRef, annotations, local)
    fn <- switch(
        method,
        "lder" = lderUnivariate,
        "gldsc" = gldscUnivariate,
        "sldsc" = sldscUnivariate,
        "hdl" = hdlUnivariate
    )
    exec(fn, !!!base, !!!estimatorArgs)
}

# The estimators read `z` and the annotation rows positionally, by each LD
# block's own snpIdx -- nothing is matched on variant id along the way. A
# sumstats or annotation set that is not the reference's variant list, in the
# reference's order, therefore yields a silently wrong h2 rather than an
# error. This is the one place holding all three objects, so check here.
# @noRd
.h2ValidateInputs <- function(sumstats, studyName, ldRef, annotations) {
    checkGenomeBuild(sumstats, ldRef)
    .h2CheckBlockCount(ldRef)
    entry <- sumStats(sumstats, studyName = studyName)
    .h2CheckVariantCount(
        length(entry),
        length(ldRef),
        glue("study '{studyName}'")
    )
    .h2CheckPositions(entry, ldRef, glue("study '{studyName}'"))
    if (is.null(annotations)) {
        return(invisible(TRUE))
    }
    checkGenomeBuild(ldRef, annotations)
    .h2CheckVariantCount(nrow(annotations), length(ldRef), "`annotations`")
    .h2CheckPositions(rowRanges(annotations), ldRef, "`annotations`")
    invisible(TRUE)
}

# Every estimator takes its standard error from a delete-one-block jackknife,
# which needs at least two blocks. With one, LDER and g-LDSC fail deep inside
# a fold fit ("'x' must be a matrix") and HDL reports an SE of exactly 0, so
# say so here instead.
# @noRd
.h2CheckBlockCount <- function(ldRef) {
    nBlocks <- length(ldBlocks(ldRef))
    if (nBlocks >= 2L) {
        return(invisible(TRUE))
    }
    msg <- glue(
        "estimateH2: the LD reference carries {nBlocks} LD block(s). The ",
        "estimators take their standard error from a delete-one-block ",
        "jackknife, so at least two are required, and a usable SE needs ",
        "many more. Build the reference from per-block LD -- see ",
        "`buildLdEigen()` / `buildLdScore()`."
    )
    abort(msg)
}

# @noRd
.h2CheckVariantCount <- function(n, M, what) {
    if (n == M) {
        return(invisible(TRUE))
    }
    msg <- glue(
        "estimateH2: {what} carries {n} variant(s) but the LD reference ",
        "carries {M}. The estimators index by the reference's variant ",
        "order, so both must hold the same variants in the same order."
    )
    abort(msg)
}

# Element-wise coordinate check of `gr` (a study's variants, or the
# annotation rows) against the reference, `what` naming it. Coordinates are
# compared rather than variant ids because the two sides routinely name the
# same variant differently (rsid vs chr:pos:A2:A1), while the same variant
# always sits at the same coordinate.
#' @importFrom SummarizedExperiment rowRanges
#' @noRd
.h2CheckPositions <- function(gr, ldRef, what) {
    same <- as.character(seqnames(gr)) == as.character(seqnames(ldRef)) &
        start(gr) == start(ldRef)
    if (all(same)) {
        return(invisible(TRUE))
    }
    nBad <- sum(!same)
    firstBad <- which(!same)[[1L]]
    msg <- glue(
        "estimateH2: {what} and the LD reference disagree at ",
        "{nBad} of {length(same)} position(s), the first at index ",
        "{firstBad}. The estimators index by the reference's variant ",
        "order, so both must hold the same variants in the same order."
    )
    abort(msg)
}

# Wrap a univariate estimator result list into an H2Estimate S4 object.
.h2EstimateFromResult <- function(result, method, M, studyName) {
    new(
        "H2Estimate",
        estimate = result$h2,
        estimateSe = result$h2Se,
        intercept = result$intercept %||% NA_real_,
        interceptSe = result$interceptSe %||% NA_real_,
        localBlocks = result$local,
        enrichment = result$enrichment,
        annotationJackknifeCoefs = result$annotationJackknifeCoefs,
        scoreStats = result$scoreStats,
        method = method,
        nSnps = as.integer(M),
        traitName = studyName
    )
}

#' @rdname estimateH2
#' @export
setMethod(
    "estimateH2",
    signature(sumstats = "GwasSumStats", ldRef = "LdStatistic"),
    function(
        sumstats,
        ldRef,
        method = "lder",
        annotations = NULL,
        local = FALSE,
        studyName = NULL,
        estimatorArgs = list()
    ) {
        method <- arg_match(method, c("lder", "gldsc", "sldsc", "hdl"))
        .validateMethodRef(method, ldRef)
        studyName <- .estimateH2ResolveStudy(sumstats, studyName)
        .h2ValidateInputs(sumstats, studyName, ldRef, annotations)
        z <- z(sumstats, studyName = studyName)
        rawN <- median(sampleSize(sumstats, studyName = studyName))
        M <- nSnps(sumstats, studyName = studyName)
        # Legacy heritability-wrapper correction (separate from the SuSiE RSS
        # binaryTraitModel handling in the fine-mapping pipeline).
        varY <- varY(sumstats, studyName = studyName)
        n <- if (is.null(varY)) rawN else rawN / varY
        result <- .estimateH2Dispatch(
            method,
            z,
            n,
            ldRef,
            annotations,
            local,
            estimatorArgs = estimatorArgs
        )
        .h2EstimateFromResult(result, method, M, studyName)
    }
)

#' @keywords internal
.validateMethodRef <- function(method, ldRef) {
    if (is_in(method, c("lder", "hdl")) && !is(ldRef, "LdEigen")) {
        msg <- glue(
            "Method '{method}' requires an LdEigen object, ",
            "got {class(ldRef)}"
        )
        abort(msg)
    }
    if (is_in(method, c("gldsc", "sldsc")) && !is(ldRef, "LdScore")) {
        msg <- glue(
            "Method '{method}' requires an LdScore object, ",
            "got {class(ldRef)}"
        )
        abort(msg)
    }
    invisible(TRUE)
}

# =============================================================================
# computeLdScores -- LD score computation
# =============================================================================

# Write each block's rows into `base`. Deliberate scatter: the score matrices
# are variants-long and the blocks partition the rows, so every block writes
# its own slice exactly once and rows in no block keep their initial value.
# Rebuilding by rbind would require the blocks to cover every variant.
# @noRd
.h2ScatterBlockRows <- function(base, blocks) {
    for (b in blocks) {
        base[b$idx, ] <- b$values
    }
    base
}

# One block's base LD scores. sum_k r^2_jk is sum_i (V[j,i] d[i])^2 because
# R = V D V' with orthonormal V, so the raw sum needs no R; the summed form of
# ldsc's unbiased r^2 then corrects it exactly, the block size being the
# weight sum.
# @noRd
.h2BlockBaseScores <- function(block, nRef) {
    Vd <- sweep(block$vectors, 2, block$values, "*")
    list(
        idx = block$snpIdx,
        values = .ldscUnbiasedSum(
            cbind(rowSums(Vd^2)),
            length(block$snpIdx),
            nRef
        )
    )
}

# Base LD scores: ldsc's unbiased r^2 summed over each variant's block.
.computeLdScoresBase <- function(ldRef, nRef) {
    .h2ScatterBlockRows(
        matrix(
            0,
            nrow = length(ldRef),
            ncol = 1,
            dimnames = list(NULL, "base_l2")
        ),
        map(
            .ldEigenBlocks(ldRef),
            .h2BlockBaseScoresAt,
            ldRef = ldRef,
            nRef = nRef
        )
    )
}

# @noRd
.h2BlockBaseScoresAt <- function(block, ldRef, nRef) {
    .h2BlockBaseScores(.ldEigenBlockAt(block, ldRef), nRef)
}

# Base + annotation-stratified LD scores. Stratified column a:
# l2_a[j] = sum_k r^2_jk annot[k, a], under ldsc's unbiased r^2.
.computeLdScoresStratified <- function(ldRef, annotations, nRef) {
    annotMat <- assay(annotations, "annotations")
    annotNames <- SummarizedExperiment::colData(annotations)$name
    .h2ScatterBlockRows(
        matrix(
            0,
            nrow = length(ldRef),
            ncol = 1 + ncol(annotMat),
            dimnames = list(NULL, c("base_l2", annotNames))
        ),
        map(
            .ldEigenBlocks(ldRef),
            .h2BlockStratScores,
            ldRef = ldRef,
            annotMat = annotMat,
            nRef = nRef
        )
    )
}

# One block's base + annotation-stratified scores. The weighted sums need R's
# entries: the eigen shortcut sum_i (V[j,i] d[i])^2 * sum_k V[k,i]^2
# annot[k,a] drops the cross terms between eigenvectors, which cancel only
# when every annot[k, a] is 1, so it is not the stratified score. R's rows
# are rebuilt from the decomposition a stripe at a time, and the summed form
# of ldsc's unbiased r^2 corrects every column at once.
# @noRd
.h2BlockStratScores <- function(block, ldRef, annotMat, nRef) {
    dec <- .ldEigenBlockAt(block, ldRef)
    idx <- dec$snpIdx
    Ab <- annotMat[idx, , drop = FALSE]
    raw <- .ldSquaredProduct(dec[c("vectors", "values")], cbind(1, Ab))
    list(
        idx = idx,
        values = .ldscUnbiasedSum(raw, c(length(idx), colSums(Ab)), nRef)
    )
}

#' @rdname computeLdScores
#' @export
setMethod(
    "computeLdScores",
    signature(ldRef = "LdEigen"),
    function(ldRef, annotations = NULL) {
        # LD scores from the eigendecompositions, under ldsc's unbiased r^2,
        # so they match buildLdScore() on the same reference.
        nRef <- .ldscPanelN(ldRef)
        if (is.null(annotations)) {
            return(.computeLdScoresBase(ldRef, nRef))
        }
        .computeLdScoresStratified(ldRef, annotations, nRef)
    }
)

#' @rdname computeLdScores
#' @export
setMethod(
    "computeLdScores",
    signature(ldRef = "LdScore"),
    function(ldRef, annotations = NULL) {
        if (is.null(annotations)) {
            return(mcols(ldRef)$ldScores)
        }
        if (!.ldScoreHasLd(ldRef)) {
            abort(str_c(
                "Annotation-stratified LD scores need each LD block's LD, ",
                "and this LdScore carries none. Build it from a genotype ",
                "panel, or from LdMatrix blocks with `keepLdMatrices = TRUE`."
            ))
        }
        # The stratified columns l2_a[j] = sum_k r^2_jk annot[k, a], under
        # ldsc's unbiased r^2, are scattered per block; the base column comes
        # from the reference and is joined on afterwards, so the scatter
        # never touches it.
        strat <- .h2StratifiedScores(
            ldRef,
            assay(annotations, "annotations")
        )
        `colnames<-`(
            cbind(mcols(ldRef)$ldScores[, 1], strat),
            c("base_l2", SummarizedExperiment::colData(annotations)$name)
        )
    }
)

# H2Estimate accessor methods (localH2Blocks /
# stratifiedHeritabilityEnrichment / scoreStats) live
# in R/h2Estimate.R alongside the class definition.

# =============================================================================
# Converter: H2Estimate -> sldsc_wrapper list format
# =============================================================================

#' @importFrom checkmate assertClass
#' @title Convert H2Estimate to S-LDSC Trait Format
#' @description Convert an \code{H2Estimate} object into the list format
#'   expected by \code{\link{standardizeSldscTrait}} and
#'   \code{\link{metaSldscRandom}}. This bridges the h2 estimation methods
#'   (LDER, gLDSC, HDL) into the sldscWrapper.R postprocessing pipeline.
#' @param h2Est An \code{H2Estimate} object with enrichment and
#'   \code{annotationJackknifeCoefs}.
#' @return A named list matching the format of \code{\link{readSldscTrait}}:
#'   \describe{
#'     \item{categories}{Character vector of annotation names}
#'     \item{tau}{Named numeric vector of per-annotation coefficients}
#'     \item{tauSe}{Named numeric vector of tau standard errors}
#'     \item{enrichment}{Named numeric vector of enrichment ratios}
#'     \item{enrichmentSe}{Named numeric vector of enrichment SEs}
#'     \item{enrichmentP}{Named numeric vector of enrichment p-values}
#'     \item{propH2}{Named numeric vector of proportion of h2}
#'     \item{propSnps}{Named numeric vector of proportion of SNPs}
#'     \item{h2g}{Numeric scalar, global h2 estimate}
#'     \item{annotationJackknifeCoefs}{Matrix (nBlocks x nCategories) for
#'   jackknife} \item{nBlocks}{Integer, number of jackknife blocks}
#'   }
#' @examples
#' data(h2EstimateExample)
#' h2EstimateToSldscTrait(h2EstimateExample)
#' @export
h2EstimateToSldscTrait <- function(h2Est) {
    assertClass(h2Est, "H2Estimate")

    enrichDf <- stratifiedHeritabilityEnrichment(h2Est)
    if (is.null(enrichDf)) {
        msg <- glue(
            "H2Estimate has no enrichment results. ",
            "Run estimateH2 with annotations to get enrichment estimates."
        )
        abort(msg)
    }

    cats <- as.character(enrichDf$annotation)
    nCats <- length(cats)

    rawTau <- annotationJackknifeCoefs(h2Est)
    # Without per-block estimates the point estimates stand in as one block.
    jkCoefs <- `colnames<-`(
        rawTau %||% matrix(enrichDf$tau, nrow = 1),
        colnames(rawTau) %||% cats
    )
    nBlocks <- nrow(jkCoefs)

    list(
        categories = cats,
        tau = set_names(enrichDf$tau, cats),
        tauSe = set_names(enrichDf$tauSe, cats),
        enrichment = set_names(enrichDf$enrichment, cats),
        enrichmentSe = set_names(enrichDf$enrichmentSe, cats),
        enrichmentP = set_names(enrichDf$enrichmentP, cats),
        propH2 = set_names(enrichDf$propH2, cats),
        propSnps = set_names(enrichDf$propSnps, cats),
        h2g = heritabilityEstimate(h2Est),
        annotationJackknifeCoefs = jkCoefs,
        nBlocks = nBlocks
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One LDER delete-one-block jackknife fold: design and fit with block `b` left
# out.
# @noRd
.lderLooBlockFit <- function(
    b,
    blockStats,
    baselineMat,
    N,
    M,
    M_aFull,
    rough,
    twostage
) {
    d <- .lderDesign(blockStats[-b], baselineMat)
    .lderGetRes(d$x, d$lam, d$ldAnnot, N, M, M_aFull, rough, twostage)
}

# One LDER local-heritability row for block `b`.
# @noRd
.lderLocalRow <- function(b, blockStats, N, a) {
    .lderLocalBlock(blockStats[[b]], b, N, a)
}

# One g-LDSC delete-one-block jackknife fold: re-solve the normal equations with
# block `cb`'s Omega-weighted contribution removed.
# @noRd
.gldscLooEstimate <- function(cb, left, right, A, N) {
    .gldscEstimate(left - cb$L, right - cb$R, A, N)
}

# One g-LDSC local-heritability row for block `b`.
# @noRd
.gldscLocalRow <- function(b, contrib, N) {
    .gldscLocalBlock(contrib[[b]], b, N)
}

# g-LDSC leave-one-block-out standardized candidate score vector for block `b`.
# @noRd
# One HDL delete-one-block jackknife fold: refit leaving block `b` out.
# @noRd
.hdlLooFit <- function(b, blockData, M, baselineMat, n, nSamples) {
    # The full-data fit already reported any bound it hit; repeating that
    # once per left-out block would bury it in noise.
    .hdlFit(
        .hdlDesign(blockData[-b], M, baselineMat),
        n,
        nSamples,
        warnOnBound = FALSE
    )
}

# One HDL local-heritability row for block `b`.
# @noRd
.hdlLocalRow <- function(b, blockData, n, nSamples, int, M) {
    .hdlLocalBlock(blockData[[b]], b, n, nSamples, int, M)
}
