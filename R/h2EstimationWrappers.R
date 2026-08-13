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
# Plus shared utilities (block ops, WLS, jackknife SE, enrichment, LD
# shrinkage, genome-build checks, tauStar standardization, meta-analysis)
# and the converter that bridges H2Estimate results into the sLDSC
# postprocessing pipeline.
#
# The H2Estimate result class itself lives in R/h2Estimate.R.
# =============================================================================

#' @title Shared Utilities for Heritability Estimation
#' @description Internal helper functions for block operations, regression,
#'   jackknife SE, enrichment computation, and meta-analysis.
#' @name pecotmr-h2-utils
#' @keywords internal
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel bplapply bpparam
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
#' @param ldBlocks An \code{LdBlocks} object.
#' @return A list of integer vectors, one per block.
#' @keywords internal
snpsPerBlock <- function(snpInfo, ldBlocks) {
    blocksGr <- getBlocks(ldBlocks)
    snpGr <- GRanges(
        seqnames = snpInfo$CHR,
        ranges = IRanges(start = snpInfo$BP, width = 1L)
    )
    hits <- findOverlaps(snpGr, blocksGr)
    split(queryHits(hits), subjectHits(hits))
}

#' @title Apply Function Per Block with BiocParallel
#' @description Apply a function to each LD block in parallel.
#' @param blockIndices List of SNP index vectors per block.
#' @param FUN Function to apply to each block's indices.
#' @param BPPARAM BiocParallel parameter object.
#' @param ... Additional arguments passed to FUN.
#' @return A list of results, one per block.
#' @keywords internal
bplapplyBlocks <- function(blockIndices, FUN, BPPARAM = NULL, ...) {
    if (is.null(BPPARAM)) {
        BPPARAM <- bpparam()
    }
    bplapply(blockIndices, FUN, BPPARAM = BPPARAM, ...)
}

# =============================================================================
# Regression utilities
# =============================================================================

#' @title Weighted Least Squares
#' @description Compute WLS estimate with standard errors.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, predictors.
#' @param w Numeric vector, weights (inverse variance).
#' @return A list with coefficients, SE, residuals, fitted values.
#' @keywords internal
weightedLs <- function(y, X, w) {
    if (is.null(dim(X))) {
        X <- matrix(X, ncol = 1)
    }
    W <- diag(sqrt(w))
    Xw <- W %*% X
    yw <- W %*% y
    XtX <- crossprod(Xw)
    Xty <- crossprod(Xw, yw)
    coef <- solve(XtX, Xty)
    fitted <- X %*% coef
    resid <- y - fitted
    # Heteroskedasticity-robust SE (HC0)
    meat <- crossprod(Xw * as.vector(resid))
    bread <- solve(XtX)
    vcov <- bread %*% meat %*% bread
    se <- sqrt(diag(vcov))
    list(
        coef = as.vector(coef),
        se = se,
        residuals = as.vector(resid),
        fitted = as.vector(fitted),
        vcov = vcov
    )
}

#' @title Jackknife Standard Errors by Block
#' @description Compute jackknife SE estimates using leave-one-block-out.
#' @param estimatesFull Numeric vector, full-sample parameter estimates.
#' @param estimatesLoo A matrix (nBlocks x nParams), leave-one-out estimates.
#' @return Numeric vector of jackknife SEs.
#' @keywords internal
jackknifeSe <- function(estimatesFull, estimatesLoo) {
    nBlocks <- nrow(estimatesLoo)
    pseudoVals <- nBlocks *
        matrix(
            estimatesFull,
            nrow = nBlocks,
            ncol = length(estimatesFull),
            byrow = TRUE
        ) -
        (nBlocks - 1) * estimatesLoo
    jkVar <- apply(pseudoVals, 2, var) / nBlocks
    sqrt(jkVar)
}

# =============================================================================
# Ridge-regularized WLS
# =============================================================================

#' @title Ridge-Regularized Weighted Least Squares
#' @description WLS with optional L2 penalty on coefficients.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, predictors.
#' @param w Numeric vector, weights (inverse variance).
#' @param lambda Numeric, ridge penalty. 0 = no penalty (delegates to
#'   \code{weightedLs}).
#' @param penalizeIntercept Logical. If FALSE (default), the last column of X
#'   (assumed to be the intercept) is not penalized.
#' @return Same structure as \code{weightedLs}: coef, se, residuals, fitted,
#'   vcov.
#' @keywords internal
weightedLsRidge <- function(y, X, w, lambda = 0, penalizeIntercept = FALSE) {
    if (lambda == 0) {
        return(weightedLs(y, X, w))
    }
    if (is.null(dim(X))) {
        X <- matrix(X, ncol = 1)
    }
    p <- ncol(X)
    W <- diag(sqrt(w))
    Xw <- W %*% X
    yw <- W %*% y
    XtX <- crossprod(Xw)
    Xty <- crossprod(Xw, yw)
    # Ridge penalty matrix (don't penalize intercept by default)
    penalty <- diag(lambda, p)
    if (!penalizeIntercept && p > 1) {
        penalty[p, p] <- 0
    }
    coef <- solve(XtX + penalty, Xty)
    fitted <- X %*% coef
    resid <- y - fitted
    # Sandwich SE accounting for ridge shrinkage
    bread <- solve(XtX + penalty)
    meat <- crossprod(Xw * as.vector(resid))
    vcov <- bread %*% meat %*% bread
    se <- sqrt(pmax(diag(vcov), 0))
    list(
        coef = as.vector(coef),
        se = se,
        residuals = as.vector(resid),
        fitted = as.vector(fitted),
        vcov = vcov
    )
}

# =============================================================================
# Baseline enrichment computation
# =============================================================================

#' @title Compute Baseline Annotation Enrichment Quantities
#' @description Given tau coefficients and a baseline annotation matrix, compute
#'   the full set of enrichment quantities: propH2, propSnps, enrichment ratio,
#'   enrichment SE (from jackknife or delta method), and p-value.
#' @param tau Numeric vector of per-annotation regression coefficients.
#' @param tauSe Numeric vector of SE for tau.
#' @param tauBlocks Numeric matrix (nBlocks x nAnnotations) of jackknife
#'   block-level tau values, or NULL.
#' @param baselineMat Numeric matrix (nSnps x nAnnotations).
#' @param annotNames Character vector of annotation names.
#' @param h2 Numeric scalar, total estimated h2.
#' @return A data.frame with columns: annotation, tau, tauSe, enrichment,
#'   enrichmentSe, enrichmentP, propH2, propSnps.
#' @keywords internal
# Enrichment SE: jackknife over per-block enrichment (preferred) or the
# delta-method fallback tauSe * M / |h2| when no blocks are available.
.baselineEnrichmentSe <- function(tauBlocks, tauSe, M_a, M, h2) {
    if (is.null(tauBlocks)) {
        return(tauSe * M / abs(h2))
    }
    nBlocks <- nrow(tauBlocks)
    h2Blocks <- as.vector(tauBlocks %*% M_a)
    h2Blocks[h2Blocks == 0] <- NA
    enrichmentBlocks <- sweep(tauBlocks, 1, h2Blocks, FUN = "/") * M
    enrichmentMean <- colMeans(enrichmentBlocks, na.rm = TRUE)
    enrichmentVar <- (nBlocks - 1) /
        nBlocks *
        colSums(sweep(enrichmentBlocks, 2, enrichmentMean)^2, na.rm = TRUE)
    sqrt(enrichmentVar)
}

computeBaselineEnrichment <- function(
    tau,
    tauSe,
    tauBlocks,
    baselineMat,
    annotNames,
    h2
) {
    M <- nrow(baselineMat)
    M_a <- colSums(baselineMat)
    propSnps <- M_a / M
    propH2 <- (tau * M_a) / h2
    # Enrichment ratio (propH2 / propSnps) = tau * M / h2.
    enrichment <- tau * M / h2
    enrichmentSe <- .baselineEnrichmentSe(tauBlocks, tauSe, M_a, M, h2)
    enrichmentP <- .zToPvalue(enrichment / enrichmentSe)
    data.frame(
        annotation = annotNames,
        tau = tau,
        tauSe = tauSe,
        enrichment = enrichment,
        enrichmentSe = enrichmentSe,
        enrichmentP = enrichmentP,
        propH2 = propH2,
        propSnps = propSnps,
        stringsAsFactors = FALSE
    )
}

# =============================================================================
# LD shrinkage
# =============================================================================

#' @title Apply LD Shrinkage
#' @description Apply shrinkage to sample LD matrix to reduce noise from finite
#'   reference panel size, following Wen & Stephens (2010).
#' @param R Numeric matrix, sample LD correlation matrix.
#' @param nRef Integer, reference panel sample size.
#' @param shrinkageType Character, one of "wen_stephens", "constant".
#' @param geneticMap Numeric vector, genetic map positions for SNPs in R.
#' @return Shrunk LD correlation matrix.
#' @keywords internal
shrinkLd <- function(
    R,
    nRef,
    shrinkageType = "wen_stephens",
    geneticMap = NULL
) {
    if (shrinkageType == "wen_stephens" && !is.null(geneticMap)) {
        # Wen & Stephens (2010) shrinkage based on genetic distance
        p <- nrow(R)
        theta <- 2 * nRef / (22 * nRef + 16) # effective recombination
        distCm <- abs(outer(geneticMap, geneticMap, "-"))
        shrinkFactor <- exp(-4 * nRef * distCm / (100 * (2 * nRef + 16)))
        RShrunk <- R * shrinkFactor
        diag(RShrunk) <- 1
    } else {
        # Simple constant shrinkage
        lambda <- 1 / sqrt(nRef)
        RShrunk <- (1 - lambda) * R + lambda * diag(nrow(R))
    }
    RShrunk
}

# =============================================================================
# Genome build utilities
# =============================================================================

#' @title Validate Genome Build Consistency
#' @description Check that genome builds match between objects. Each object
#'   contributes a single genome build (from its \code{genome} slot).
#' @param ... Objects with a \code{genome} slot.
#' @return TRUE if all match, error otherwise.
#' @keywords internal
checkGenomeBuild <- function(...) {
    objects <- list(...)
    genomes <- vapply(
        objects,
        function(x) {
            if (
                is(x, "GwasSumStats") ||
                    is(x, "QtlSumStats") ||
                    is(x, "LdStatistic") ||
                    is(x, "AnnotationMatrix") ||
                    is(x, "LdBlocks")
            ) {
                getGenome(x)
            } else {
                stop("Unknown object type for genome build check")
            }
        },
        character(1)
    )
    if (length(unique(genomes)) > 1L) {
        stop("Genome build mismatch: ", paste(genomes, collapse = ", "))
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
#' @param tauBlocks Numeric matrix (nBlocks x nAnnotations) of block-level tau
#'   estimates from delete-one jackknife.
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
standardizeTauStar <- function(tau, tauBlocks, sdAnnot, MRef, h2g) {
    if (length(tau) != length(sdAnnot)) {
        stop("standardizeTauStar: tau and sdAnnot must have the same length.")
    }
    if (h2g == 0) {
        stop("standardizeTauStar: h2g must be non-zero.")
    }

    # Gazal standardization: tau* = tau * sdAnnot * MRef / h2g
    coef <- sdAnnot * MRef / h2g
    tauStar <- tau * coef

    # Jackknife SE from block-level tau
    tauStarBlocks <- sweep(tauBlocks, 2L, coef, FUN = "*")
    nBlocks <- nrow(tauStarBlocks)
    jkVar <- apply(tauStarBlocks, 2L, function(x) var(x, na.rm = TRUE))
    tauStarSe <- sqrt((nBlocks - 1)^2 / nBlocks * jkVar)

    list(tauStar = tauStar, tauStarSe = tauStarSe)
}

# =============================================================================
# DerSimonian-Laird random-effects meta-analysis
# =============================================================================

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
#' @return A list with \code{mean}, \code{se}, \code{tau2}, \code{I2} (a [0, 1]
#'   proportion) and \code{Q} (Cochran's Q).
#' @importFrom metafor rma
#' @noRd
# Degenerate meta-analysis inputs: k = 0 (all NA) or k = 1 (pass-through), plus
# positivity validation. Returns NULL when the full estimator should run.
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
        stop(".rmaMeta: all ses must be positive and finite.")
    }
    NULL
}

# tryCatch handler: iterative estimators can fail on small/near-homogeneous
# inputs; fall back to closed-form DerSimonian-Laird (never iterates).
.rmaMetaFallback <- function(e, means, ses, method) {
    if (identical(method, "DL")) {
        stop(e)
    }
    warning(
        ".rmaMeta: metafor::rma(method = '",
        method,
        "') failed (",
        conditionMessage(e),
        "); falling back to DL.",
        call. = FALSE
    )
    metafor::rma(yi = means, sei = ses, method = "DL")
}

.rmaMeta <- function(means, ses, method = "DL") {
    k <- length(means)
    if (k != length(ses)) {
        stop(".rmaMeta: means and ses must have the same length.")
    }
    early <- .rmaMetaEarly(means, ses, k)
    if (!is.null(early)) {
        return(early)
    }
    # metafor is a hard dependency; it reports I-squared as a percentage.
    fit <- tryCatch(
        metafor::rma(yi = means, sei = ses, method = method),
        error = function(e) .rmaMetaFallback(e, means, ses, method)
    )
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
#'   al. 2022). Within each LD block the z-scores are whitened in the eigenbasis,
#'   \eqn{x_i = (u_i'z)/\sqrt{\lambda_i}}, and the whitened \eqn{\chi^2 = x_i^2}
#'   is regressed on \eqn{\lambda_i} via IRWLS: \eqn{E[x_i^2] = (N h2/M)\lambda_i
#'   + (1 + N a)}, so \eqn{Cov(z) = (N h2/M) R^2 + (1 + N a) R}. The IRWLS weight
#'   uses an in-sample (\eqn{\min(1,\lambda)}) or out-of-sample (logistic) form
#'   selected by \code{getInSample}, and a two-stage step estimates the
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
    baseline <- getBaseline(annotations)
    if (ncol(getAnnotations(baseline)) > 0) {
        getAnnotations(baseline)
    } else {
        NULL
    }
}

# Per-block whitened statistics x_i = (u_i'z)/sqrt(lam_i) and stratified
# eigenvalue scores l_{ia} = sum_j A_{ja} V_{ji}^2 (all-ones column when
# unstratified). Near-null directions (lam < 1e-6) are zeroed (as upstream).
.lderBlockStats <- function(block, z, baselineMat) {
    V <- block$vectors
    lam <- block$values
    x <- as.vector(crossprod(V, z[block$snpIdx]) / sqrt(pmax(lam, 1e-12)))
    x[lam < 1e-6] <- 0
    ldAnnot <- if (is.null(baselineMat)) {
        matrix(1, length(lam), 1)
    } else {
        crossprod(V^2, baselineMat[block$snpIdx, , drop = FALSE])
    }
    list(x = x, lam = lam, ldAnnot = ldAnnot)
}

# Concatenate per-block stats into the genome-wide regression design. M_a is
# the annotation SNP counts (univariate: total number of directions).
.lderDesign <- function(blockStats, baselineMat) {
    x <- unlist(lapply(blockStats, `[[`, "x"))
    list(
        x = x,
        lam = unlist(lapply(blockStats, `[[`, "lam")),
        ldAnnot = do.call(rbind, lapply(blockStats, `[[`, "ldAnnot")),
        blockId = rep(
            seq_along(blockStats),
            lengths(lapply(blockStats, `[[`, "lam"))
        ),
        M_a = if (is.null(baselineMat)) length(x) else colSums(baselineMat)
    )
}

# IRWLS weight w = 1/(1 + predGenetic)^2 * w2(lam); predGenetic uses the current
# (clamped) slopes, w2 = pmin(1,lam) (in-sample) or logistic (out-of-sample).
.lderWeights <- function(lam, ldAnnot, slopes, N, M, rough) {
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
    loo <- lapply(seq_len(nB), function(b) {
        d <- .lderDesign(blockStats[-b], baselineMat)
        .lderGetRes(d$x, d$lam, d$ldAnnot, N, M, M_aFull, rough, twostage)
    })
    h2Loo <- unlist(lapply(loo, `[[`, "h2"))
    aLoo <- unlist(lapply(loo, `[[`, "a"))
    tauBlocks <- do.call(rbind, lapply(loo, `[[`, "tau"))
    jkSe <- function(v) sqrt((nB - 1) / nB * sum((v - mean(v))^2))
    list(
        h2Se = jkSe(h2Loo),
        aSe = jkSe(aLoo),
        tauSe = apply(tauBlocks, 2, jkSe),
        tauBlocks = tauBlocks
    )
}

# tau / tauSe / tauBlocks for the enrichment framework (univariate: tau = h2).
.lderEnrichParts <- function(fit, jk, baselineMat) {
    if (is.null(baselineMat)) {
        return(list(tau = fit$h2, tauSe = jk$h2Se, tauBlocks = NULL))
    }
    list(tau = fit$tau, tauSe = jk$tauSe, tauBlocks = jk$tauBlocks)
}

# Local heritability for one block: single-slope fit with intercept fixed at
# the global inflation; h2Local = tau * (block directions).
.lderLocalBlock <- function(bs, b, N, a) {
    d <- bs$lam
    if (length(d) < 3) {
        return(data.frame(blockId = b, h2Local = NA, h2LocalSe = NA))
    }
    chi2 <- bs$x^2
    w <- 1 / (2 * pmax(chi2, 1)^2)
    slope <- sum(w * d * (chi2 - (N * a + 1))) / sum(w * d^2)
    info <- sum(w * (N * d)^2)
    data.frame(
        blockId = b,
        h2Local = slope / N * length(d),
        h2LocalSe = length(d) / sqrt(max(info, 1e-10))
    )
}

# Per-block local heritability (NULL unless requested).
.lderLocal <- function(blockStats, N, a) {
    do.call(
        rbind,
        lapply(seq_along(blockStats), function(b) {
            .lderLocalBlock(blockStats[[b]], b, N, a)
        })
    )
}

# Weighted score z for one candidate: standardized weighted covariance of the
# baseline residual with the candidate design column N*lam*l_ic.
.lderCandidateScore <- function(ldCandCol, lam, w, resid, N, keep) {
    d <- N * lam[keep] * ldCandCol[keep]
    sum(w[keep] * resid[keep] * d) / sqrt(sum(w[keep] * d^2))
}

# Jackknife (leave-one-block-out) correlation of the candidate score z's.
.lderScoreCor <- function(ldCand, lam, w, resid, N, blockId, nCand) {
    if (nCand == 1) {
        return(matrix(1, 1, 1))
    }
    looZ <- do.call(
        rbind,
        lapply(sort(unique(blockId)), function(b) {
            keep <- blockId != b
            map_dbl(seq_len(nCand), function(c) {
                .lderCandidateScore(ldCand[, c], lam, w, resid, N, keep)
            })
        })
    )
    R <- cor(looZ)
    R[is.na(R)] <- 0
    R
}

# LDER candidate-annotation score statistics against the fitted baseline.
.lderScoreStats <- function(eigenRef, annotations, design, fit, N) {
    candidate <- getCandidates(annotations)
    candMat <- getAnnotations(candidate)
    nCand <- ncol(candMat)
    if (nCand == 0) {
        return(list(enrichment = NULL, scoreStats = NULL))
    }
    lam <- design$lam
    mu <- pmax(
        N * lam * as.vector(design$ldAnnot %*% fit$tau) + (N * fit$a + 1),
        1e-8
    )
    w <- 1 / (2 * mu^2)
    resid <- design$x^2 - mu
    ldCand <- do.call(
        rbind,
        lapply(getEigenList(eigenRef), function(bk) {
            crossprod(bk$vectors^2, candMat[bk$snpIdx, , drop = FALSE])
        })
    )
    keepAll <- rep(TRUE, length(lam))
    scoreZ <- map_dbl(seq_len(nCand), function(c) {
        .lderCandidateScore(ldCand[, c], lam, w, resid, N, keepAll)
    })
    R <- .lderScoreCor(ldCand, lam, w, resid, N, design$blockId, nCand)
    .buildEnrichmentResult(getAnnotationMeta(candidate), scoreZ, R)
}

lderUnivariate <- function(
    z,
    n,
    eigenRef,
    annotations = NULL,
    local = FALSE,
    lambda = 0
) {
    rough <- !getInSample(eigenRef)
    baselineMat <- .h2BaselineMat(annotations)
    blockStats <- lapply(
        getEigenList(eigenRef),
        .lderBlockStats,
        z = z,
        baselineMat = baselineMat
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
    ep <- .lderEnrichParts(fit, jk, baselineMat)
    localDf <- if (local) .lderLocal(blockStats, n, fit$a) else NULL
    scoreStats <- if (is.null(annotations)) {
        NULL
    } else {
        .lderScoreStats(eigenRef, annotations, design, fit, n)$scoreStats
    }
    jkRes <- list(
        se = c(jk$h2Se, jk$aSe),
        tauSe = ep$tauSe,
        tauBlocks = ep$tauBlocks
    )
    enrichmentDf <- .h2Enrichment(
        baselineMat,
        ep$tau,
        ep$tauSe,
        ep$tauBlocks,
        fit$h2
    )
    .h2Result(fit$h2, jkRes, fit$a, ep$tau, localDf, enrichmentDf, scoreStats)
}


# Assemble the enrichment result shared by the h2 enrichment helpers: a
# per-annotation table (z + two-tailed p) plus the scoreStats bundle
# (z, R, annotation names) consumed downstream.
.buildEnrichmentResult <- function(candMeta, scoreZ, R) {
    list(
        enrichment = data.frame(
            annotation = candMeta$name,
            scoreZ = scoreZ,
            scoreP = .zToPvalue(scoreZ),
            stringsAsFactors = FALSE
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
        tauBlocks = jkRes$tauBlocks,
        local = localDf,
        enrichment = enrichmentDf,
        scoreStats = scoreStats
    )
}

# Per-annotation enrichment table via the shared tau-based framework (used by
# HDL and LDER; g-LDSC builds its own from partitioned h2). NULL when
# unstratified.
.h2Enrichment <- function(baselineMat, tau, tauSe, tauBlocks, h2) {
    if (is.null(baselineMat)) {
        return(NULL)
    }
    annotNames <- if (!is.null(colnames(baselineMat))) {
        colnames(baselineMat)
    } else {
        paste0("annot_", seq_len(ncol(baselineMat)))
    }
    computeBaselineEnrichment(
        tau,
        tauSe,
        tauBlocks,
        baselineMat,
        annotNames,
        h2
    )
}

#' @title g-LDSC: Generalized LD Score Regression
#' @description Estimate heritability and functional enrichment by feasible
#'   generalized least squares (FGLS) on LD scores (Xiong et al. 2024). The
#'   response \eqn{Z_j^2 - 1} is regressed on stratified LD scores
#'   \eqn{\ell_{ja} = \sum_k R^2_{jk} A_{ka}} (with a confounding intercept)
#'   using the full residual covariance
#'   \eqn{\Omega = (\hat{N\tau}\, R D_a R + R)^2} (elementwise), inverted per LD
#'   block via a 99\%-eigenvalue-mass truncation. \eqn{\hat{N\tau} = \sum(Z^2-1)
#'   / \sum \ell} is a method-of-moments plug-in; \eqn{D_a = diag(\sum_a A_{ja})}.
#'   Coefficients are \eqn{\tau = (X'\Omega^{-1}X)^{-1} X'\Omega^{-1} y / N} and
#'   partitioned heritabilities \eqn{A'A\,\tau}. Standard errors come from a
#'   delete-one-block jackknife that subtracts each block's already-\eqn{\Omega}-
#'   weighted contribution from the accumulated normal equations.
#' @name pecotmr-h2-gldsc
#' @keywords internal
#' @references
#'   Xiong Z, Thach TQ, Zhang YD, Sham PC (2024). Improved estimation
#'   of functional enrichment in SNP heritability using feasible
#'   generalized least squares. HGG Advances, 5(2):100272.
#' @importFrom stats var pnorm
NULL

# Full annotation matrix A = [base (all-ones) | baseline columns]; the base
# column carries genome-wide h2 (est.h[1] = total).
.gldscAnnotMatrix <- function(annotations, M) {
    baselineMat <- .h2BaselineMat(annotations)
    if (is.null(baselineMat)) matrix(1, M, 1) else cbind(base = 1, baselineMat)
}

# Per-block quantities: LD matrix R, block annotations, stratified LD scores
# l = (R^2) A, and response y = Z^2 - 1.
.gldscBlockPrep <- function(block, z, A) {
    R <- block$R
    idx <- block$snpIdx
    Ab <- A[idx, , drop = FALSE]
    list(R = R, Ab = Ab, ldsc = (R^2) %*% Ab, y = z[idx]^2 - 1, snpIdx = idx)
}

# 99%-mass truncated inverse of Omega = (rawNtau * R Da R + R)^2 (upstream
# eigen.cut). Da = diag(rowSums(A)) so R Da R = crossprod(R sqrt(Da)).
.gldscOmegaInv <- function(R, Ab, rawNtau, cutpoint = 0.99) {
    Omega <- (rawNtau * crossprod(R * sqrt(rowSums(Ab))) + R)^2
    ei <- eigen(Omega, symmetric = TRUE)
    k <- which(cumsum(ei$values) / sum(ei$values) >= cutpoint)[1]
    Vc <- ei$vectors[, seq_len(k), drop = FALSE]
    Vc %*% (t(Vc) / ei$values[seq_len(k)])
}

# Block contribution to the GLS normal equations: L = X'Omega^-1 X,
# R = X'Omega^-1 y, with X = [1 | l].
.gldscBlockGls <- function(prep, rawNtau) {
    OmInv <- .gldscOmegaInv(prep$R, prep$Ab, rawNtau)
    X <- cbind(1, prep$ldsc)
    XtOi <- crossprod(X, OmInv)
    list(L = XtOi %*% X, R = as.vector(XtOi %*% prep$y))
}

# GLS estimate from accumulated normal equations: tau = solve(L, R)/N (the
# leading entry is the confounding intercept), partitioned h2 = A'A tau.
.gldscEstimate <- function(left, right, A, N) {
    ridge <- 1e-8 * mean(abs(diag(left)))
    coef <- as.vector(solve(left + diag(ridge, nrow(left)), right))
    tauAll <- coef / N
    tau <- tauAll[-1]
    estH <- as.vector(crossprod(A) %*% tau)
    list(
        tau = tau,
        coef = coef,
        estH = estH,
        h2 = estH[1],
        intercept = N * tauAll[1] + 1
    )
}

# Delete-one-block jackknife: subtract each block's Omega-weighted contribution
# from the accumulated normal equations (so folds stay GLS-weighted). Jackknife
# variance = var * (nB-1)^2 / nB (Xiong et al.).
.gldscJackknife <- function(contrib, left, right, A, N) {
    nB <- length(contrib)
    loo <- lapply(contrib, function(cb) {
        .gldscEstimate(left - cb$L, right - cb$R, A, N)
    })
    estHBlocks <- do.call(rbind, lapply(loo, `[[`, "estH"))
    tauBlocks <- do.call(rbind, lapply(loo, `[[`, "tau"))
    intLoo <- unlist(lapply(loo, `[[`, "intercept"))
    jkSe <- function(v) sqrt(var(v) * (nB - 1)^2 / nB)
    list(
        h2Se = jkSe(estHBlocks[, 1]),
        intSe = jkSe(intLoo),
        tauSe = apply(tauBlocks, 2, jkSe),
        tauBlocks = tauBlocks,
        estHBlocks = estHBlocks
    )
}

# Enrichment for baseline (non-base) annotations from partitioned h2 est.h:
# propH2 = est.h_a / total, propSnps = M_a / M, enrichment = propH2 / propSnps.
.gldscEnrichmentDf <- function(fit, jk, baselineMat, annotNames, M) {
    total <- fit$estH[1]
    M_a <- colSums(baselineMat)
    propSnps <- M_a / M
    propH2 <- fit$estH[-1] / total
    enrichment <- propH2 / propSnps
    eBlocks <- sweep(
        jk$estHBlocks[, -1, drop = FALSE] / jk$estHBlocks[, 1],
        2,
        propSnps,
        "/"
    )
    nB <- nrow(eBlocks)
    eSe <- sqrt(apply(eBlocks, 2, var) * (nB - 1)^2 / nB)
    data.frame(
        annotation = annotNames,
        tau = fit$tau[-1],
        tauSe = jk$tauSe[-1],
        enrichment = enrichment,
        enrichmentSe = eSe,
        enrichmentP = pnorm(abs(enrichment - 1) / eSe, lower.tail = FALSE) * 2,
        propH2 = propH2,
        propSnps = propSnps,
        stringsAsFactors = FALSE
    )
}

# Local heritability for one LD block: a per-block GLS on [1 | base LD score];
# h2Local = tau_base * (block SNPs), SE from the GLS coefficient covariance.
.gldscLocalBlock <- function(prep, b, rawNtau, N) {
    OmInv <- .gldscOmegaInv(prep$R, prep$Ab, rawNtau)
    X <- cbind(1, prep$ldsc[, 1])
    XtOi <- crossprod(X, OmInv)
    L <- XtOi %*% X
    coef <- solve(L, as.vector(XtOi %*% prep$y))
    p <- nrow(prep$R)
    covTau <- (solve(L) / N^2)[2, 2]
    data.frame(
        blockId = b,
        h2Local = coef[2] / N * p,
        h2LocalSe = sqrt(max(covTau, 0)) * p
    )
}

# Per-block local heritability (NULL unless requested).
.gldscLocal <- function(preps, rawNtau, N) {
    do.call(
        rbind,
        lapply(seq_along(preps), function(b) {
            .gldscLocalBlock(preps[[b]], b, rawNtau, N)
        })
    )
}

# Candidate-annotation score against the FITTED baseline GLS residual
# r = y - X coef: U_c = sum_b l_c' Omega^-1 r, I_c = sum_b l_c' Omega^-1 l_c,
# scoreZ = U_c / sqrt(I_c). Uses the real fitted coefficients (not a refit).
.gldscCandidateBlock <- function(prep, candMat, rawNtau, coef) {
    OmInv <- .gldscOmegaInv(prep$R, prep$Ab, rawNtau)
    r <- prep$y - as.vector(cbind(1, prep$ldsc) %*% coef)
    ldCand <- (prep$R^2) %*% candMat[prep$snpIdx, , drop = FALSE]
    ct <- crossprod(ldCand, OmInv)
    list(U = as.vector(ct %*% r), I = diag(ct %*% ldCand))
}

# sHDL-style candidate score statistics + jackknife score correlation.
.gldscScoreStats <- function(preps, candMat, rawNtau, coef) {
    nCand <- ncol(candMat)
    if (nCand == 0) {
        return(NULL)
    }
    per <- lapply(
        preps,
        .gldscCandidateBlock,
        candMat = candMat,
        rawNtau = rawNtau,
        coef = coef
    )
    U <- Reduce(`+`, lapply(per, `[[`, "U"))
    I <- Reduce(`+`, lapply(per, `[[`, "I"))
    scoreZ <- U / sqrt(pmax(I, 1e-10))
    R <- if (nCand == 1) {
        matrix(1, 1, 1)
    } else {
        looZ <- do.call(
            rbind,
            lapply(seq_along(per), function(b) {
                (U - per[[b]]$U) / sqrt(pmax(I - per[[b]]$I, 1e-10))
            })
        )
        Rc <- cor(looZ)
        Rc[is.na(Rc)] <- 0
        Rc
    }
    list(scoreZ = scoreZ, R = R)
}

gldscUnivariate <- function(
    z,
    n,
    ldRef,
    annotations = NULL,
    local = FALSE,
    lambda = 0
) {
    ldMatrixList <- getLdMatrixList(ldRef)
    if (length(ldMatrixList) == 0L) {
        stop(
            "g-LDSC requires full per-block LD matrices (ldMatrixList). ",
            "Read the LD reference with the LD matrices retained."
        )
    }
    M <- nrow(getSnpInfo(ldRef))
    A <- .gldscAnnotMatrix(annotations, M)
    preps <- lapply(ldMatrixList, .gldscBlockPrep, z = z, A = A)
    rawNtau <- sum(unlist(lapply(preps, function(p) sum(p$y)))) /
        sum(unlist(lapply(preps, function(p) sum(p$ldsc))))
    contrib <- lapply(preps, .gldscBlockGls, rawNtau = rawNtau)
    left <- Reduce(`+`, lapply(contrib, `[[`, "L"))
    right <- Reduce(`+`, lapply(contrib, `[[`, "R"))
    fit <- .gldscEstimate(left, right, A, n)
    jk <- .gldscJackknife(contrib, left, right, A, n)
    .gldscResult(fit, jk, annotations, preps, rawNtau, local, n, M)
}

# Assemble the gldscUnivariate return (enrichment, local, candidate scores).
.gldscResult <- function(fit, jk, annotations, preps, rawNtau, local, n, M) {
    baselineMat <- .h2BaselineMat(annotations)
    localDf <- if (local) .gldscLocal(preps, rawNtau, n) else NULL
    enrichmentDf <- NULL
    scoreStats <- NULL
    if (!is.null(baselineMat)) {
        nm <- getAnnotationMeta(getBaseline(annotations))$name
        enrichmentDf <- .gldscEnrichmentDf(fit, jk, baselineMat, nm, M)
    }
    if (!is.null(annotations)) {
        candMat <- getAnnotations(getCandidates(annotations))
        sc <- .gldscScoreStats(preps, candMat, rawNtau, fit$coef)
        if (!is.null(sc)) {
            cn <- getAnnotationMeta(getCandidates(annotations))$name
            scoreStats <- list(z = sc$scoreZ, R = sc$R, annotationNames = cn)
        }
    }
    tau <- if (is.null(baselineMat)) fit$h2 else fit$tau[-1]
    tauSe <- if (is.null(baselineMat)) jk$h2Se else jk$tauSe[-1]
    jkRes <- list(
        se = c(jk$h2Se, jk$intSe),
        tauSe = tauSe,
        tauBlocks = if (is.null(baselineMat)) {
            NULL
        } else {
            jk$tauBlocks[, -1, drop = FALSE]
        }
    )
    .h2Result(
        fit$h2,
        jkRes,
        fit$intercept,
        tau,
        localDf,
        enrichmentDf,
        scoreStats
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
# score column this is exactly upstream HDL's llfun (Ning et al. 2020).
# param = c(h2_1, .., h2_nTau, int).
.hdlNll <- function(param, n, nRef, lam, bstar, ldAnnotScaled, lim = exp(-18)) {
    nTau <- ncol(ldAnnotScaled)
    h2a <- param[seq_len(nTau)]
    int <- param[nTau + 1L]
    lamh2 <- as.vector(ldAnnotScaled %*% h2a) *
        lam^2 -
        lam * sum(h2a) / nRef +
        int * lam / n
    lamh2 <- pmax(lamh2, lim)
    sum(log(lamh2)) + sum(bstar^2 / lamh2)
}

# Per-block HDL quantities: eigenvalues, rotated N-normalized effects
# bstar = V'(z/sqrt(n)), and stratified eigenvalue scores l_{ia} (NULL when
# unstratified).
.hdlBlockData <- function(block, z, n, baselineMat) {
    idx <- block$snpIdx
    V <- block$vectors
    ldAnnot <- if (is.null(baselineMat)) {
        NULL
    } else {
        crossprod(V^2, baselineMat[idx, , drop = FALSE])
    }
    list(
        lam = block$values,
        bstar = as.vector(crossprod(V, z[idx] / sqrt(n))),
        ldAnnot = ldAnnot
    )
}

# Genome-wide design: eigenvalues, rotated effects, block ids, the M_a-scaled
# eigenvalue-score matrix, and M_a. Univariate uses a single all-ones score
# column (l_{i,base} = 1) scaled by total M.
.hdlDesign <- function(blockData, M, baselineMat) {
    lam <- unlist(lapply(blockData, `[[`, "lam"))
    bstar <- unlist(lapply(blockData, `[[`, "bstar"))
    blockId <- rep(
        seq_along(blockData),
        lengths(lapply(blockData, `[[`, "lam"))
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
    ldAnnot <- do.call(rbind, lapply(blockData, `[[`, "ldAnnot"))
    list(
        lam = lam,
        bstar = bstar,
        blockId = blockId,
        ldAnnotScaled = sweep(ldAnnot, 2, M_a, "/"),
        M_a = M_a
    )
}

# Fit HDL by L-BFGS-B over (h2_1..h2_nTau, int); bounds h2_a in [0,1], int in
# [0,10] (Ning et al. 2020). Parametrizing by per-annotation h2 (scale ~1)
# keeps the L-BFGS-B finite differences well conditioned.
.hdlFit <- function(design, n, nRef) {
    nTau <- ncol(design$ldAnnotScaled)
    opt <- optim(
        c(rep(0.1, nTau), 1),
        .hdlNll,
        n = n,
        nRef = nRef,
        lam = design$lam,
        bstar = design$bstar,
        ldAnnotScaled = design$ldAnnotScaled,
        method = "L-BFGS-B",
        lower = c(rep(0, nTau), 0),
        upper = c(rep(1, nTau), 10)
    )
    h2a <- opt$par[seq_len(nTau)]
    list(h2a = h2a, int = opt$par[nTau + 1L], h2 = sum(h2a))
}

# Baseline modeled variance lamh2_i at the fitted parameters (for local h2 and
# the sHDL candidate score test).
.hdlBaselineVar <- function(design, ft, n, nRef, lim = exp(-18)) {
    s <- as.vector(design$ldAnnotScaled %*% ft$h2a)
    pmax(
        s * design$lam^2 - design$lam * ft$h2 / nRef + ft$int * design$lam / n,
        lim
    )
}

# =============================================================================
# HDL standard errors (delete-one-block jackknife) + enrichment inputs
# =============================================================================

# Delete-one-block jackknife (Ning et al. 2020): refit leaving each block out,
# SE = sqrt(mean((jack - mean)^2) * (nBlocks - 1)). Returns per-annotation and
# total-h2 SEs, the intercept SE, and the LOO per-annotation blocks.
.hdlJackknife <- function(blockData, M, baselineMat, n, nRef) {
    nBlocks <- length(blockData)
    loo <- lapply(seq_len(nBlocks), function(b) {
        .hdlFit(.hdlDesign(blockData[-b], M, baselineMat), n, nRef)
    })
    h2aBlocks <- do.call(rbind, lapply(loo, `[[`, "h2a"))
    intLoo <- unlist(lapply(loo, `[[`, "int"))
    jkSe <- function(x) sqrt(mean((x - mean(x))^2) * (nBlocks - 1))
    list(
        h2aSe = apply(h2aBlocks, 2, jkSe),
        h2Se = jkSe(rowSums(h2aBlocks)),
        intSe = jkSe(intLoo),
        h2aBlocks = h2aBlocks
    )
}

# Convert fitted per-annotation h2 to per-SNP tau = h2_a / M_a for the
# enrichment framework (univariate: tau = h2 scalar, no blocks).
.hdlEnrichParts <- function(ft, jk, design, baselineMat) {
    if (is.null(baselineMat)) {
        return(list(tau = ft$h2, tauSe = jk$h2Se, tauBlocks = NULL))
    }
    list(
        tau = ft$h2a / design$M_a,
        tauSe = jk$h2aSe / design$M_a,
        tauBlocks = sweep(jk$h2aBlocks, 2, design$M_a, "/")
    )
}

# =============================================================================
# HDL local heritability + sHDL candidate score test
# =============================================================================

# Local h2 for one block: MLE of the block's genome-wide-scaled h2 with the
# intercept fixed at the global fit; SE from the Gaussian-variance Fisher
# information sum((dlamh2/dh2)^2 / (2 lamh2^2)).
.hdlLocalBlock <- function(bd, b, n, nRef, int, M, lim = exp(-18)) {
    varOf <- function(h2) {
        pmax(h2 / M * bd$lam^2 - h2 * bd$lam / nRef + int * bd$lam / n, lim)
    }
    nll <- function(h2) sum(log(varOf(h2)) + bd$bstar^2 / varOf(h2))
    h2 <- optimize(nll, c(0, 1))$minimum
    lamh2 <- varOf(h2)
    g <- bd$lam^2 / M - bd$lam / nRef
    info <- sum(g^2 / (2 * lamh2^2))
    data.frame(
        blockId = b,
        h2Local = h2,
        h2LocalSe = 1 / sqrt(max(info, 1e-10))
    )
}

# Per-block local heritability (NULL unless requested).
.hdlLocal <- function(blockData, n, nRef, int, M) {
    do.call(
        rbind,
        lapply(seq_along(blockData), function(b) {
            .hdlLocalBlock(blockData[[b]], b, n, nRef, int, M)
        })
    )
}

# sHDL score z for one candidate annotation: score U_c = sum_i g_ic *
# (bstar_i^2/lamh2_i^2 - 1/lamh2_i) with g_ic = dlamh2_i/dh2_c, standardized by
# sqrt(Var(U_c)) = sqrt(2 sum_i g_ic^2 / lamh2_i^2).
.shdlCandidateScore <- function(
    ldCandCol,
    M_c,
    lam,
    lamh2,
    scoreResid,
    nRef,
    keep
) {
    g <- ldCandCol[keep] * lam[keep]^2 / M_c - lam[keep] / nRef
    lh <- lamh2[keep]
    sum(g * scoreResid[keep]) / sqrt(2 * sum(g^2 / lh^2))
}

# sHDL candidate-annotation score statistics + jackknife score correlation.
.shdlStratified <- function(eigenRef, annotations, design, lamh2) {
    candidate <- getCandidates(annotations)
    candMat <- getAnnotations(candidate)
    nCand <- ncol(candMat)
    if (nCand == 0) {
        return(list(enrichment = NULL, scoreStats = NULL))
    }
    nRef <- getNRef(eigenRef)
    M_c <- colSums(candMat)
    ldCand <- do.call(
        rbind,
        lapply(getEigenList(eigenRef), function(bk) {
            crossprod(bk$vectors^2, candMat[bk$snpIdx, , drop = FALSE])
        })
    )
    scoreResid <- design$bstar^2 / lamh2^2 - 1 / lamh2
    keepAll <- rep(TRUE, length(lamh2))
    scoreZ <- map_dbl(seq_len(nCand), function(c) {
        .shdlCandidateScore(
            ldCand[, c],
            M_c[c],
            design$lam,
            lamh2,
            scoreResid,
            nRef,
            keepAll
        )
    })
    R <- .shdlScoreCor(ldCand, M_c, design, lamh2, scoreResid, nRef, nCand)
    .buildEnrichmentResult(getAnnotationMeta(candidate), scoreZ, R)
}

# Jackknife (leave-one-block-out) correlation of the candidate score z's.
.shdlScoreCor <- function(ldCand, M_c, design, lamh2, scoreResid, nRef, nCand) {
    if (nCand == 1) {
        return(matrix(1, 1, 1))
    }
    blocks <- sort(unique(design$blockId))
    looZ <- do.call(
        rbind,
        lapply(blocks, function(b) {
            keep <- design$blockId != b
            map_dbl(seq_len(nCand), function(c) {
                .shdlCandidateScore(
                    ldCand[, c],
                    M_c[c],
                    design$lam,
                    lamh2,
                    scoreResid,
                    nRef,
                    keep
                )
            })
        })
    )
    R <- cor(looZ)
    R[is.na(R)] <- 0
    R
}

# =============================================================================
# hdlUnivariate -- orchestrator
# =============================================================================

hdlUnivariate <- function(
    z,
    n,
    eigenRef,
    annotations = NULL,
    local = FALSE,
    lambda = 0
) {
    eigenList <- getEigenList(eigenRef)
    M <- nrow(getSnpInfo(eigenRef))
    nRef <- getNRef(eigenRef)
    baselineMat <- .h2BaselineMat(annotations)
    blockData <- lapply(
        eigenList,
        .hdlBlockData,
        z = z,
        n = n,
        baselineMat = baselineMat
    )
    design <- .hdlDesign(blockData, M, baselineMat)
    ft <- .hdlFit(design, n, nRef)
    jk <- .hdlJackknife(blockData, M, baselineMat, n, nRef)
    ep <- .hdlEnrichParts(ft, jk, design, baselineMat)
    localDf <- if (local) .hdlLocal(blockData, n, nRef, ft$int, M) else NULL
    scoreStats <- if (is.null(annotations)) {
        NULL
    } else {
        lamh2 <- .hdlBaselineVar(design, ft, n, nRef)
        .shdlStratified(eigenRef, annotations, design, lamh2)$scoreStats
    }
    jkRes <- list(
        se = c(jk$h2Se, jk$intSe),
        tauSe = ep$tauSe,
        tauBlocks = ep$tauBlocks
    )
    enrichmentDf <- .h2Enrichment(
        baselineMat,
        ep$tau,
        ep$tauSe,
        ep$tauBlocks,
        ft$h2
    )
    .h2Result(ft$h2, jkRes, ft$int, ep$tau, localDf, enrichmentDf, scoreStats)
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
# estimateH2 -- main dispatch
# =============================================================================

# Resolve which study's stats to operate on; defaults to the single entry.
.estimateH2ResolveStudy <- function(sumstats, study) {
    if (!is.null(study)) {
        return(study)
    }
    if (nrow(sumstats) != 1L) {
        stop(
            "`study` is required when the GwasSumStats has ",
            nrow(sumstats),
            " entries."
        )
    }
    as.character(sumstats$study[[1L]])
}

# Dispatch to the method-specific univariate estimator.
.estimateH2Dispatch <- function(
    method,
    z,
    n,
    ldRef,
    annotations,
    local,
    ...
) {
    switch(
        method,
        "lder" = lderUnivariate(z, n, ldRef, annotations, local, ...),
        "gldsc" = gldscUnivariate(z, n, ldRef, annotations, local, ...),
        "hdl" = hdlUnivariate(z, n, ldRef, annotations, local, ...)
    )
}

# Wrap a univariate estimator result list into an H2Estimate S4 object.
.h2EstimateFromResult <- function(result, method, M, study) {
    new(
        "H2Estimate",
        h2 = result$h2,
        h2Se = result$h2Se,
        intercept = result$intercept %||% NA_real_,
        interceptSe = result$interceptSe %||% NA_real_,
        local = result$local,
        enrichment = result$enrichment,
        tauBlocks = result$tauBlocks,
        scoreStats = result$scoreStats,
        method = method,
        nSnps = as.integer(M),
        traitName = study
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
        study = NULL,
        ...
    ) {
        method <- match.arg(method, c("lder", "gldsc", "hdl"))
        .validateMethodRef(method, ldRef)
        study <- .estimateH2ResolveStudy(sumstats, study)
        z <- getZ(sumstats, study = study)
        n <- median(getN(sumstats, study = study))
        M <- nSnps(sumstats, study = study)
        # Legacy heritability-wrapper correction (separate from the SuSiE RSS
        # binaryTraitModel handling in the fine-mapping pipeline).
        varY <- getVarY(sumstats, study = study)
        if (!is.null(varY)) {
            n <- n / varY
        }
        result <- .estimateH2Dispatch(
            method,
            z,
            n,
            ldRef,
            annotations,
            local,
            ...
        )
        .h2EstimateFromResult(result, method, M, study)
    }
)

#' @keywords internal
.validateMethodRef <- function(method, ldRef) {
    if (method %in% c("lder", "hdl") && !is(ldRef, "LdEigen")) {
        stop(
            "Method '",
            method,
            "' requires an LdEigen object, ",
            "got ",
            class(ldRef)
        )
    }
    if (method == "gldsc" && !is(ldRef, "LdScore")) {
        stop(
            "Method 'gldsc' requires an LdScore object, ",
            "got ",
            class(ldRef)
        )
    }
    invisible(TRUE)
}

# =============================================================================
# computeLdScores -- LD score computation
# =============================================================================

# Base LD scores l2[j] = sum_i (V[j,i] * d[i])^2 (since R = V D V').
.computeLdScoresBase <- function(eigenList, nSnps) {
    l2 <- numeric(nSnps)
    for (b in seq_along(eigenList)) {
        block <- eigenList[[b]]
        Vd <- sweep(block$vectors, 2, block$values, "*")
        l2[block$snpIdx] <- rowSums(Vd^2)
    }
    matrix(l2, ncol = 1, dimnames = list(NULL, "base_l2"))
}

# Base + annotation-stratified LD scores. Stratified column a:
# l2_a[j] = sum_i (V[j,i] d[i])^2 * (sum_k V[k,i]^2 annot[k,a]).
.computeLdScoresStratified <- function(eigenList, annotations, nSnps) {
    annotMat <- getAnnotations(annotations)
    nAnnot <- ncol(annotMat)
    l2Strat <- matrix(0, nrow = nSnps, ncol = 1 + nAnnot)
    for (b in seq_along(eigenList)) {
        block <- eigenList[[b]]
        idx <- block$snpIdx
        V <- block$vectors
        Vd2 <- sweep(V, 2, block$values, "*")^2
        l2Strat[idx, 1] <- rowSums(Vd2)
        for (a in seq_len(nAnnot)) {
            annotWeights <- as.vector(crossprod(V^2, annotMat[idx, a]))
            l2Strat[idx, 1 + a] <- as.vector(Vd2 %*% annotWeights)
        }
    }
    colnames(l2Strat) <- c("base_l2", getAnnotationMeta(annotations)$name)
    l2Strat
}

#' @rdname computeLdScores
#' @export
setMethod(
    "computeLdScores",
    signature(ldRef = "LdEigen"),
    function(ldRef, annotations = NULL, ...) {
        # Reconstruct LD scores from eigendecompositions
        # l2[j] = sum_k r^2_{jk} = sum_b sum_{eigenvalues in b} V[j,.]^2 * d
        nSnps <- nrow(getSnpInfo(ldRef))
        eigenList <- getEigenList(ldRef)
        if (is.null(annotations)) {
            return(.computeLdScoresBase(eigenList, nSnps))
        }
        .computeLdScoresStratified(eigenList, annotations, nSnps)
    }
)

#' @rdname computeLdScores
#' @export
setMethod(
    "computeLdScores",
    signature(ldRef = "LdScore"),
    function(ldRef, annotations = NULL, ...) {
        if (is.null(annotations)) {
            return(getLdScores(ldRef))
        }

        # Compute annotation-stratified LD scores using LD matrices
        ldMatrixList <- getLdMatrixList(ldRef)
        if (length(ldMatrixList) == 0) {
            stop(
                "Annotation-stratified LD scores require ldMatrixList in ",
                "LdScore. ",
                "Recompute the LD reference with full LD matrices."
            )
        }

        nSnps <- nrow(getSnpInfo(ldRef))
        annotMat <- getAnnotations(annotations)
        nAnnot <- ncol(annotMat)

        # Base L2 + annotation-stratified columns
        l2Strat <- matrix(0, nrow = nSnps, ncol = 1 + nAnnot)
        l2Strat[, 1] <- getLdScores(ldRef)[, 1]

        for (b in seq_along(ldMatrixList)) {
            block <- ldMatrixList[[b]]
            R <- block$R
            idx <- block$snpIdx
            R2 <- R^2
            for (a in seq_len(nAnnot)) {
                # l2_a[j] = sum_k R^2_{jk} * annot[k, a]
                l2Strat[idx, 1 + a] <- as.vector(R2 %*% annotMat[idx, a])
            }
        }

        colNames <- c("base_l2", getAnnotationMeta(annotations)$name)
        colnames(l2Strat) <- colNames
        l2Strat
    }
)

# H2Estimate accessor methods (getLocal/getEnrichment/getScoreStats) live
# in R/h2Estimate.R alongside the class definition.

# =============================================================================
# Converter: H2Estimate -> sldsc_wrapper list format
# =============================================================================

#' @title Convert H2Estimate to S-LDSC Trait Format
#' @description Convert an \code{H2Estimate} object into the list format
#'   expected by \code{\link{standardizeSldscTrait}} and
#'   \code{\link{metaSldscRandom}}. This bridges the h2 estimation methods
#'   (LDER, gLDSC, HDL) into the sldscWrapper.R postprocessing pipeline.
#' @param h2Est An \code{H2Estimate} object with enrichment and tauBlocks.
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
#'     \item{tauBlocks}{Matrix (nBlocks x nCategories) for jackknife}
#'     \item{nBlocks}{Integer, number of jackknife blocks}
#'   }
#' @examples
#' data(h2EstimateExample)
#' h2EstimateToSldscTrait(h2EstimateExample)
#' @export
h2EstimateToSldscTrait <- function(h2Est) {
    if (!is(h2Est, "H2Estimate")) {
        stop("h2Est must be an H2Estimate object")
    }

    enrichDf <- getEnrichment(h2Est)
    if (is.null(enrichDf)) {
        stop(
            "H2Estimate has no enrichment results. ",
            "Run estimateH2 with annotations to get enrichment estimates."
        )
    }

    cats <- as.character(enrichDf$annotation)
    nCats <- length(cats)

    tauBlocks <- getTauBlocks(h2Est)
    if (is.null(tauBlocks)) {
        # Create a dummy single-block matrix from the point estimates
        tauBlocks <- matrix(enrichDf$tau, nrow = 1)
        colnames(tauBlocks) <- cats
        nBlocks <- 1L
    } else {
        nBlocks <- nrow(tauBlocks)
        if (is.null(colnames(tauBlocks))) {
            colnames(tauBlocks) <- cats
        }
    }

    list(
        categories = cats,
        tau = setNames(enrichDf$tau, cats),
        tauSe = setNames(enrichDf$tauSe, cats),
        enrichment = setNames(enrichDf$enrichment, cats),
        enrichmentSe = setNames(enrichDf$enrichmentSe, cats),
        enrichmentP = setNames(enrichDf$enrichmentP, cats),
        propH2 = setNames(enrichDf$propH2, cats),
        propSnps = setNames(enrichDf$propSnps, cats),
        h2g = getH2(h2Est),
        tauBlocks = tauBlocks,
        nBlocks = nBlocks
    )
}
