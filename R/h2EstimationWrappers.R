# =============================================================================
# Heritability Estimation Wrappers
# -----------------------------------------------------------------------------
# Consolidated entry points for univariate heritability estimation: both
# `estimateH2` methods, with the estimators they call.
#
# From GWAS summary statistics, one genome-wide estimate per study:
#   * gLDSC (Generalized LD Score Regression, Xiong et al. 2024)
#   * HDL/sHDL (High-Definition Likelihood, Ning et al. 2020 + Zhao 2023)
#   * LDER (LD Eigenvalue Regression, Song et al. 2022)
#   * S-LDSC (stratified LD score regression, Finucane et al. 2015)
#
# From QTL summary statistics, one cis estimate per gene:
#   * LAVA (Werme et al. 2022)
#   * HDL-L (Li et al. 2025)
#   * FastHer (Svishcheva 2025; HEELS, Li et al. 2023, maximizes the same
#     likelihood)
#
# The estimators are pecotmr's own implementations rather than an external
# engine, so the entry points live beside them instead of in a *Pipeline.R
# file.
#
# Plus shared utilities (block ops, WLS, jackknife SE, enrichment,
# genome-build checks, tauStar standardization, meta-analysis) and the
# converter that bridges H2Estimate results into the sLDSC postprocessing
# pipeline.
#
# The H2Estimate result class itself lives in R/H2Estimate.R.
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
    user <- as.list(metaArgs)
    exec(
        metafor::rma,
        yi = means,
        sei = ses,
        method = "DL",
        !!!user
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
    user <- as.list(metaArgs)
    fit <- if (identical(method, "DL")) {
        exec(
            metafor::rma,
            yi = means,
            sei = ses,
            method = "DL",
            !!!user
        )
    } else {
        try_fetch(
            exec(
                metafor::rma,
                yi = means,
                sei = ses,
                method = method,
                !!!user
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
#'   al. 2022): one genome-wide estimate, with optional annotation
#'   stratification.
#' @name pecotmr-h2-lder
#' @keywords internal
#' @references
#'   Song S, Jiang W, Zhang Y, Hou L, Zhao H (2022). Leveraging LD
#'   eigenvalue regression to improve the estimation of SNP heritability
#'   and confounding inflation. Am J Hum Genet, 109(5):802-811.
NULL

# =============================================================================
# The eigen regression, genome-wide
# =============================================================================

#' @title Univariate LDER
#' @description Estimate SNP heritability by LD eigenvalue regression (Song et
#'   al. 2022). Within each LD block the z-scores are whitened in the
#'   eigenbasis,
#'   \eqn{x_i = (u_i'z)/\sqrt{\lambda_i}}, and the whitened \eqn{\chi^2 = x_i^2}
#'   is regressed on \eqn{\lambda_i} via IRWLS:
#'   \eqn{E[x_i^2] = (N h2/M)\lambda_i
#'   + (1 + N a)}, so \eqn{Cov(z) = (N h2/M) R^2 + (1 + N a) R}. Under
#'   LDER's settings the IRWLS weight takes an in-sample
#'   (\eqn{\min(1,\lambda)}) or out-of-sample (logistic) eigenvalue taper,
#'   as the reference's \code{inSample} sets, and a two-stage step estimates
#'   the confounding inflation \eqn{a} on a \eqn{\chi^2}-filtered subset. The
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

# Baseline annotation matrix (NULL when absent or empty). Shared by the eigen
# regression and the eigen likelihood.
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

# Whether a fit on `annotations` is stratified: its baseline design
# (.h2DesignMat()) has more than the one column an unstratified fit regresses
# on. An all-ones baseline alone does not.
# @noRd
.h2Stratified <- function(annotations) {
    design <- .h2DesignMat(annotations)
    !is.null(design) && ncol(design) > 1L
}

# TRUE when the all-ones vector lies in the column span of `B`.
# @noRd
.h2SpansOnes <- function(B) {
    max(abs(1 - as.vector(B %*% .h2OnesCoef(B)))) < 1e-6
}

# The least-squares coefficients of the all-ones vector on `B`, from B'B
# (pivoted, so overlapping or duplicated columns do no harm), which keeps a
# genome-wide B from being copied into a QR.
# @noRd
.h2OnesCoef <- function(B) {
    unname(coalesce(qr.coef(qr(crossprod(B)), colSums(B)), 0))
}

# The partitioned heritabilities, enrichment and per-annotation coefficients
# of a stratified fit, reported for the user's baseline annotations alone,
# from its per-variant-effect coefficients `tau` over the design `A` and their
# delete-one-block values `tauBlocks` (folds x columns), under the frequency
# weights `w`.
# @noRd
.h2StratifiedReport <- function(A, tau, tauBlocks, w = NULL) {
    user <- attr(A, "userCols")
    tauSe <- apply(tauBlocks, 2L, .jackknifeSeFromLoo)[user]
    list(
        tau = tau[user],
        tauSe = tauSe,
        annotationJackknifeCoefs = tauBlocks[, user, drop = FALSE],
        enrichment = bind_cols(
            tibble(
                annotation = colnames(A)[user],
                tau = unname(tau[user]),
                tauSe = unname(tauSe)
            ),
            .h2EnrichmentColumns(A, tau, tauBlocks, w)
        )
    )
}

# The enrichment of each of the user's annotations and the two shares it is
# the ratio of. Annotations overlap in general -- the base overlaps every one
# -- so an annotation's heritability is (A'A tau)_a, not tau_a M_a, and the
# total is 1'A tau; under frequency weights `w`, which scale each variant's
# variance, they are (A'WA tau)_a and w'A tau. These are shares of variants
# and of heritability only for values in [0, 1], the domain Dey et al.
# (2020) define enrichment on; outside it they depend on where the column's
# zero happens to be and are NA, while tau stays reported.
# @noRd
.h2EnrichmentColumns <- function(A, tau, tauBlocks, w = NULL) {
    user <- attr(A, "userCols")
    weighted <- .h2WeightedDesign(A, w)
    AtA <- as.matrix(if (is.null(w)) crossprod(A) else crossprod(A, weighted))
    size <- colSums(A)
    total <- colSums(weighted)
    propSnps <- size[user] / nrow(A)
    propH2 <- as.vector(AtA %*% tau)[user] / sum(total * tau)
    blockH2 <- (tauBlocks %*% AtA)[, user, drop = FALSE]
    blockEnrichment <- sweep(
        blockH2 / as.vector(tauBlocks %*% total),
        2L,
        propSnps,
        "/"
    )
    userCols <- A[, user, drop = FALSE]
    outside <- which(colSums(userCols < 0 | userCols > 1) > 0)
    columns <- list(
        enrichment = unname(propH2 / propSnps),
        enrichmentSe = unname(apply(blockEnrichment, 2L, .jackknifeSeFromLoo)),
        enrichmentP = .h2EnrichmentP(
            AtA,
            size,
            nrow(A),
            tau,
            tauBlocks,
            total
        )[user],
        propH2 = unname(propH2),
        propSnps = unname(propSnps)
    )
    as_tibble(map(columns, replace, list = outside, values = NA_real_))
}

# ldsc's test of enrichment (regressions.py, _overlap_output): a t test, on as
# many degrees of freedom as blocks, of the per-variant heritability inside each
# column against outside it, h2_a / M_a - (h2 - h2_a) / (M - M_a). That is
# linear in tau and 0 exactly when the enrichment is 1, and unlike the
# enrichment it has no denominator to take a block's noise from. NA where it is
# undefined: a column covering every variant has no outside. `total` is what a
# unit of each coefficient adds to the total heritability: the column's size
# without frequency weights.
# @noRd
#' @importFrom stats pt
.h2EnrichmentP <- function(AtA, size, M, tau, tauBlocks, total = size) {
    difference <- AtA / size - sweep(-AtA, 2L, total, "+") / (M - size)
    estimate <- as.vector(difference %*% tau)
    se <- apply(tauBlocks %*% t(difference), 2L, .jackknifeSeFromLoo)
    p <- 2 * pt(abs(estimate / se), nrow(tauBlocks), lower.tail = FALSE)
    unname(if_else(se > 0, p, NA_real_))
}

# The result of a stratified or unstratified fit: `tau` and `tauBlocks` over
# the design `A` (NULL when unstratified, where the heritability stands in
# for the coefficient), whose variants' variances carry the frequency weights
# `w` (none: every variant alike).
# @noRd
.h2FitResult <- function(
    h2,
    h2Se,
    intercept,
    interceptSe,
    A,
    tau,
    tauBlocks,
    scoreStats,
    w = NULL
) {
    if (is.null(A)) {
        jkRes <- list(
            se = c(h2Se, interceptSe),
            tauSe = h2Se,
            annotationJackknifeCoefs = NULL
        )
        return(.h2Result(h2, jkRes, intercept, h2, NULL, scoreStats))
    }
    report <- .h2StratifiedReport(A, tau, tauBlocks, w)
    jkRes <- list(
        se = c(h2Se, interceptSe),
        tauSe = report$tauSe,
        annotationJackknifeCoefs = report$annotationJackknifeCoefs
    )
    .h2Result(h2, jkRes, intercept, report$tau, report$enrichment, scoreStats)
}

# Per-block whitened statistics x_i = (u_i'z)/sqrt(lam_i) and stratified
# eigenvalue scores l_{ia} = sum_j A_{ja} V_{ji}^2 (all-ones column when
# unstratified), plus the candidate annotations' scores, from the block's
# decomposition. Near-null directions (lam < 1e-6) are zeroed, as in the
# original LDER implementation, which keeps every direction of the block:
# the zero-eigenvalue directions a decomposition from a panel's factor does
# not list are added as rows with x = 0 and lam = 0. Under frequency weights
# `w` the block's pseudo-eigenpairs stand in for its eigenpairs
# (.h2WeightedBlock()), and the annotation matrices are the weighted ones.
.h2EigenRegressionBlockStats <- function(
    block,
    eigenRef,
    z,
    baselineMat,
    candMat,
    w = NULL
) {
    dec <- .h2WeightedBlock(.ldEigenBlockAt(block, eigenRef), w)
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

# Rows `idx` of an annotation matrix, dense: the matrix may be stored sparse
# (.annotStored()), and each estimator computes on its block as a base matrix.
# Its column and row sums read either form (Matrix's colSums()/rowSums()).
# @noRd
#' @importFrom Matrix colSums rowSums
.denseBlock <- function(A, idx) {
    methods::as(A[idx, , drop = FALSE], "matrix")
}

# One block's eigenvalue scores sum_k V_ki^2 a_k for each column of
# `annotMat`, with a zero row for each zero-eigenvalue direction the
# decomposition does not list: every use multiplies a direction's row by its
# eigenvalue. A weighted decomposition's are taken on its loadings
# (.h2WeightedBlock()), `annotMat` then weighted too. NULL without
# annotations.
.h2EigenAnnot <- function(dec, annotMat) {
    if (is.null(annotMat)) {
        return(NULL)
    }
    listed <- crossprod(
        (dec$loadings %||% dec$vectors)^2,
        .denseBlock(annotMat, dec$snpIdx)
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

# Delete-one-block jackknife standard error from the n leave-one-out
# estimates `v`: sqrt((n - 1) / n * sum((v - mean(v))^2)). Every genome-wide
# estimator's jackknife reports this quantity.
# @noRd
.jackknifeSeFromLoo <- function(v) {
    n <- length(v)
    sqrt((n - 1) / n * sum((v - mean(v))^2))
}

# The eigen regression's candidate-annotation score statistics against the
# fitted baseline: the regression of x^2 on mu, the rows' fitted mean, scored
# for each candidate's column N lam l_c with the weights 1 / (2 mu^2) and the
# baseline and inflation (the rows' columns) as nuisance
# (.h2EfficientScoreBlocks()).
.h2EigenRegressionScoreStats <- function(
    rows,
    fit,
    blockStats,
    annotations,
    N
) {
    candidate <- testedAnnotations(annotations)
    candMat <- assay(candidate, "annotations")
    if (ncol(candMat) == 0) {
        return(list(enrichment = NULL, scoreStats = NULL))
    }
    mu <- pmax(as.vector(rows$X %*% fit$coef), 1e-8)
    w <- 1 / (2 * mu^2)
    blockCands <- map(blockStats, "ldCand")
    ldCand <- exec(rbind, !!!blockCands)
    Ub <- .h2EfficientScoreBlocks(
        N * rows$lam * ldCand,
        rows$X,
        w,
        w * (rows$y - mu),
        rows$block
    )
    sc <- .h2BlockScoreStats(Ub)
    .buildEnrichmentResult(
        SummarizedExperiment::colData(candidate),
        sc$scoreZ,
        sc$R
    )
}

# The genome-wide eigen regression, stratified by `annotations` where given,
# under the settings `param` (the lder preset's by default:
# .h2GenomeWideSettings()): each block's whitened statistics as rows of the
# regression every eigen regression fits (.h2EigenRegressionRows(),
# .h2Regression()), with the delete-one-block jackknife.
eigenRegressionUnivariate <- function(
    z,
    n,
    eigenRef,
    annotations = NULL,
    param = .h2GenomeWideSettings(
        "lder",
        .h2Stratified(annotations),
        inSample(eigenRef)
    ),
    weights = .h2RefWeights(param, eigenRef)
) {
    baselineMat <- .h2DesignMat(annotations)
    modelMat <- .h2WeightedDesign(baselineMat, weights)
    blockStats <- .h2EigenRegressionAllBlockStats(
        eigenRef,
        z,
        modelMat,
        annotations,
        weights
    )
    rows <- .h2EigenRegressionRows(
        blockStats,
        modelMat,
        n,
        length(eigenRef),
        param
    )
    fit <- .h2Regression(param, rows)
    .h2RowsResult(
        fit,
        .h2RegressionLooCoefs(param, rows),
        rows,
        baselineMat,
        .h2EigenRegressionCandidateScores(
            rows,
            fit,
            blockStats,
            annotations,
            n
        ),
        weights
    )
}

# The eigen regression's rows: each direction's whitened chi-square x^2 on its
# eigenvalue scores, column a scaled by n / M_a (.h2ColumnSizes()) so that its
# coefficient times its share (.h2ColumnShares()) is annotation a's
# heritability, and the noise column; with the eigenvalue taper, the weights'
# cap on each column (M_a / M, a per-variant rate of at most 1 / M), the
# unstratified column for the first weights, and each direction's block for the
# jackknife.
# @noRd
.h2EigenRegressionRows <- function(blockStats, modelMat, n, M, param) {
    x <- list_c(map(blockStats, "x"))
    lam <- list_c(map(blockStats, "lam"))
    ldAnnot <- exec(rbind, !!!map(blockStats, "ldAnnot"))
    size <- .h2ColumnSizes(modelMat, M)
    list(
        y = x^2,
        X = cbind(sweep(n * lam * ldAnnot, 2L, size, "/"), 1),
        chi2 = x^2,
        n = n,
        taper = .h2Taper(param@procedures@eigenvalueTaper, pmax(lam, 0)),
        redundancy = rep(1, length(x)),
        caps = size / M,
        size = size,
        share = .h2ColumnShares(modelMat, size),
        base = n * lam / M,
        lam = lam,
        block = rep(seq_along(blockStats), lengths(map(blockStats, "lam"))),
        weights = "variance"
    )
}

# A genome-wide regression's result from its fit and its leave-one-block-out
# coefficients `loo` (blocks x columns): the heritability, the intercept
# (the chi-square regression's, about 1 without inflation), each column's
# coefficient over its size, and their jackknife.
# @noRd
.h2RowsResult <- function(fit, loo, rows, A, scoreStats, w) {
    k <- ncol(rows$X) - 1L
    genetic <- seq_len(k)
    .h2FitResult(
        h2 = fit$h2,
        h2Se = .jackknifeSeFromLoo(.h2RowsH2(loo, rows)),
        intercept = fit$intercept,
        interceptSe = .jackknifeSeFromLoo(loo[, k + 1L]),
        A = A,
        tau = fit$coef[genetic] / rows$size,
        tauBlocks = sweep(loo[, genetic, drop = FALSE], 2L, rows$size, "/"),
        scoreStats = scoreStats,
        w = w
    )
}

# Each LD block's eigen-regression statistics, in one pass over the blocks, each
# decomposed only while its statistics are taken, so a reference read from a
# panel holds one block at a time; under frequency weights `w`, the weighted
# baseline `baselineMat` and candidates.
# @noRd
.h2EigenRegressionAllBlockStats <- function(
    eigenRef,
    z,
    baselineMat,
    annotations,
    w = NULL
) {
    map(
        .ldEigenBlocks(eigenRef),
        .h2EigenRegressionBlockStats,
        eigenRef = eigenRef,
        z = z,
        baselineMat = baselineMat,
        candMat = .h2WeightedDesign(.h2CandidateMat(annotations), w),
        w = w
    )
}

# The eigen regression's score statistics for the candidate annotations; NULL
# when there are none to test.
# @noRd
.h2EigenRegressionCandidateScores <- function(
    rows,
    fit,
    blockStats,
    annotations,
    n
) {
    if (is.null(annotations)) {
        return(NULL)
    }
    .h2EigenRegressionScoreStats(
        rows,
        fit,
        blockStats,
        annotations,
        n
    )$scoreStats
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
.h2Result <- function(h2, jkRes, intercept, tau, enrichmentDf, scoreStats) {
    list(
        h2 = h2,
        h2Se = jkRes$se[1],
        intercept = intercept,
        interceptSe = jkRes$se[2],
        tau = tau,
        tauSe = jkRes$tauSe,
        annotationJackknifeCoefs = jkRes$annotationJackknifeCoefs,
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
.h2VariantAnnotMatrix <- function(annotations, M) {
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
.h2GlsBlockDesign <- function(block, panel, z, A, candMat) {
    idx <- block$idx
    Ab <- .denseBlock(A, idx)
    Cb <- if (!is.null(candMat)) .denseBlock(candMat, idx)
    scores <- .ldSquaredProduct(.ldStatBlockLd(block, panel), cbind(Ab, Cb))
    nA <- ncol(Ab)
    list(
        Ab = Ab,
        ldScores = scores[, seq_len(nA), drop = FALSE],
        ldCand = if (!is.null(Cb)) scores[, -seq_len(nA), drop = FALSE],
        y = z[idx]^2 - 1
    )
}

# Omega = S o S, with S = rawNtau * R Da R + R and Da = diag(rowSums(A)),
# reduced to its leading eigenpairs carrying `cutpoint` of the eigenvalue
# mass (.h2TruncatedEigen()). A variant's entry of Da is its plug-in
# variance, which a column taking negative values can drive below 0 and the
# square root in .h2GlsS() cannot take, so it is floored there; Omega only
# weights the fit, whose mean model this leaves alone.
.h2GlsOmegaEigen <- function(ld, Ab, rawNtau, cutpoint = 0.99) {
    S <- .h2GlsS(ld, pmax(rowSums(Ab), 0), rawNtau)
    .h2TruncatedEigen(S^2, cutpoint)
}

# The leading eigenpairs of a symmetric weight carrying `cutpoint` of its
# eigenvalue mass (eigen.cut in the original gLDSC implementation), as every
# GLS fit takes its inverse: applied as V diag(1 / values) V', never formed.
# @noRd
.h2TruncatedEigen <- function(Omega, cutpoint) {
    ei <- eigen(Omega, symmetric = TRUE)
    k <- which(cumsum(ei$values) / sum(ei$values) >= cutpoint)[[1L]]
    list(
        vectors = ei$vectors[, seq_len(k), drop = FALSE],
        values = ei$values[seq_len(k)]
    )
}

# The GLS of `y` on `X` under the truncated weight `om`
# (.h2TruncatedEigen()), with the relative `ridge` on its normal equations
# (.h2RidgeSolve()): its coefficients, the map `A` from the response to them,
# and the inverse of its ridged normal equations.
# @noRd
.h2GlsFit <- function(y, X, om, ridge) {
    P <- crossprod(om$vectors, X)
    Pw <- P / om$values
    inverse <- .h2RidgeSolve(crossprod(Pw, P), diag(ncol(X)), ridge)
    A <- inverse %*% t(om$vectors %*% Pw)
    list(coef = as.vector(A %*% y), A = A, inverse = inverse)
}

# A ridge-regularised solve of the normal equations `left` against `right`,
# as every regression and GLS fit solves them: each equation's ridge is
# `ridge` of its own diagonal entry (of the diagonal's mean for one with next
# to none, as an empty column's), so that it does not depend on the columns'
# scales. It keeps a design whose columns are dependent solvable: the base
# LD-score column is the sum of the stratified ones whenever the baseline
# annotations partition the genome, and a jackknife fold can empty a column.
# @noRd
.h2RidgeSolve <- function(left, right, ridge) {
    d <- abs(diag(left))
    scale <- if_else(d > 1e-8 * mean(d), d, mean(d))
    solve(left + diag(ridge * scale, nrow(left)), right)
}

# S = rawNtau * R Da R + noise * R for a block's LD source, with Da =
# diag(d). From a stored R, R Da R = crossprod(R sqrt(Da)). From a factor F,
# with R = F'F / (m - 1), S = F' (rawNtau / (m - 1) * (F Da F') F + noise F) /
# (m - 1): only S is formed, which Omega needs whole, never the block's R.
# Any other source's R is formed first.
.h2GlsS <- function(ld, d, rawNtau, noise = 1) {
    if (is.null(ld$R) && is.null(ld$ldFactor)) {
        return(.h2GlsS(list(R = .ldSourceMatrix(ld)), d, rawNtau, noise))
    }
    if (!is.null(ld$R)) {
        return(rawNtau * crossprod(ld$R * sqrt(d)) + noise * ld$R)
    }
    ldFactor <- ld$ldFactor
    m1 <- nrow(ldFactor) - 1
    gram <- tcrossprod(sweep(ldFactor, 2L, sqrt(d), "*"))
    crossprod(
        ldFactor,
        rawNtau / m1 * (gram %*% ldFactor) + noise * ldFactor
    ) /
        m1
}

# One block's contribution to the GLS normal equations, L = X'Omega^-1 X and
# R = X'Omega^-1 y with X = [1 | l], plus its candidate-score pieces. Omega
# is formed for the block, decomposed and dropped, so only one block's is
# ever in memory; its inverse keeps `cutpoint` of the eigenvalue mass.
.h2GlsBlockContrib <- function(block, design, panel, rawNtau, cutpoint) {
    om <- .h2GlsOmegaEigen(
        .ldStatBlockLd(block, panel),
        design$Ab,
        rawNtau,
        cutpoint
    )
    X <- cbind(1, design$ldScores)
    P <- crossprod(om$vectors, X)
    Py <- crossprod(om$vectors, design$y)
    Pw <- P / om$values
    list(
        L = crossprod(Pw, P),
        R = as.vector(crossprod(Pw, Py)),
        cand = .h2GlsCandidatePieces(om, design$ldCand, P, Py)
    )
}

# A block's candidate-score pieces. The score is taken against the FITTED
# baseline residual, U_c = l_c' Omega^-1 (y - X coef), which is linear in
# coef: keeping l_c' Omega^-1 y and l_c' Omega^-1 X lets the block's score be
# formed once the fit is known, without revisiting the block. NULL without
# candidates.
.h2GlsCandidatePieces <- function(om, ldCand, P, Py) {
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
# leading entry is the confounding intercept), partitioned h2 = A'A tau, the
# relative `ridge` on the normal equations (.h2RidgeSolve()).
.h2GlsEstimate <- function(left, right, A, N, ridge) {
    coef <- as.vector(.h2RidgeSolve(left, right, ridge))
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
# from the accumulated normal equations (so folds stay GLS-weighted), each
# fold's Omega the full data's, as the original gldsc's oneout() keeps it.
# Jackknife variance = var * (nB-1)^2 / nB (Xiong et al.).
.h2GlsJackknife <- function(contrib, left, right, A, N, ridge) {
    loo <- map(
        contrib,
        .h2GlsLooEstimate,
        left = left,
        right = right,
        A = A,
        N = N,
        ridge = ridge
    )
    taus <- map(loo, "tau")
    list(
        h2Se = .jackknifeSeFromLoo(map_dbl(loo, "h2")),
        intSe = .jackknifeSeFromLoo(map_dbl(loo, "intercept")),
        annotationJackknifeCoefs = exec(rbind, !!!taus)
    )
}

# One block's score for each candidate, against the fitted baseline GLS
# residual r = y - X coef, with the candidate's LD scores l_c first
# residualized on the design: U_c = (l_c - X beta_c)' Omega^-1 r, beta_c =
# (X' Omega^-1 X)^-1 X' Omega^-1 l_c. Summed over blocks the adjustment
# vanishes (X' Omega^-1 r = 0 at the fit); block by block it does not, and
# it is what makes the blocks' scores the efficient ones the jackknife
# variance is taken over.
.h2GlsCandidateBlock <- function(cb, coef, beta) {
    U <- cb$cand$Uy - as.vector(cb$cand$UX %*% coef)
    Xr <- cb$R - as.vector(cb$L %*% coef)
    U - as.vector(crossprod(beta, Xr))
}

# Candidate score statistics at the real fitted coefficients (not a refit),
# calibrated by the blocks' own spread (.h2BlockScoreStats()) rather than by
# the model information l_c' Omega^-1 l_c, which ignored the estimated
# baseline and trusted the plug-in Omega.
.h2GlsScoreStats <- function(contrib, coef, ridge) {
    left <- reduce(map(contrib, "L"), `+`)
    lX <- reduce(map(contrib, .h2GlsCandidateCross), `+`)
    beta <- .h2RidgeSolve(left, lX, ridge)
    blockScores <- map(
        contrib,
        .h2GlsCandidateBlock,
        coef = coef,
        beta = beta
    )
    .h2BlockScoreStats(exec(rbind, !!!blockScores))
}

# One block's candidate cross-products X' Omega^-1 l_c, design by candidate.
# @noRd
.h2GlsCandidateCross <- function(cb) {
    t(cb$cand$UX)
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
    grad <- as.matrix(gC)
    gtg <- crossprod(G, G * info)
    beta <- .h2RidgeSolve(gtg, crossprod(G, grad * info), 1e-8)
    rowsum((grad - G %*% beta) * e, blockOf)
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
    blockScores <- as.matrix(Ub)
    nB <- nrow(blockScores)
    U <- colSums(blockScores)
    centered <- sweep(blockScores, 2L, U / nB)
    V <- nB / (nB - 1) * colSums(centered^2)
    R <- if (ncol(blockScores) == 1L) {
        matrix(1, 1, 1)
    } else {
        Rc <- stats::cov2cor(crossprod(centered))
        replace(Rc, is.na(Rc), 0)
    }
    list(scoreZ = unname(U / sqrt(pmax(V, .Machine$double.xmin))), R = R)
}

# The genome-wide variant GLS, stratified by `annotations` where given, under
# the settings `param` (the gldsc preset's by default: .h2GenomeWideSettings()).
variantGlsUnivariate <- function(
    z,
    n,
    ldRef,
    annotations = NULL,
    param = .h2GenomeWideSettings("gldsc", .h2Stratified(annotations), FALSE),
    weights = .h2RefWeights(param, ldRef)
) {
    proc <- param@procedures
    .h2GlsAssertLd(ldRef)
    M <- length(ldRef)
    A <- .h2VariantAnnotMatrix(annotations, M)
    # The model's columns: each variant's variance under frequency weights
    # is its weight times A's; A itself is what the annotations cover.
    modelMat <- .h2WeightedDesign(A, weights)
    panel <- ldReference(ldRef)
    blocks <- .ldScoreBlocks(ldRef)
    # Two passes over the blocks: Omega depends on rawNtau, which needs every
    # block's design first. A block read from a panel is read once per pass.
    designs <- map(
        blocks,
        .h2GlsBlockDesign,
        panel = panel,
        z = z,
        A = modelMat,
        candMat = .h2WeightedDesign(.h2CandidateMat(annotations), weights)
    )
    rawNtau <- sum(list_c(map(designs, "y"))) /
        sum(list_c(map(designs, "ldScores")))
    contrib <- map2(
        blocks,
        designs,
        .h2GlsBlockContrib,
        panel = panel,
        rawNtau = rawNtau,
        cutpoint = proc@omegaFraction
    )
    left <- reduce(map(contrib, "L"), `+`)
    right <- reduce(map(contrib, "R"), `+`)
    fit <- .h2GlsEstimate(left, right, modelMat, n, proc@ridge)
    jk <- .h2GlsJackknife(contrib, left, right, modelMat, n, proc@ridge)
    .h2GlsResult(fit, jk, annotations, contrib, A, proc@ridge, weights)
}

# The candidate-annotation score statistics, when candidates were supplied.
# NULL otherwise.
# @noRd
.h2GlsCandidateScores <- function(annotations, contrib, fit, ridge) {
    if (is.null(.h2CandidateMat(annotations))) {
        return(NULL)
    }
    sc <- .h2GlsScoreStats(contrib, fit$coef, ridge)
    list(
        z = sc$scoreZ,
        R = sc$R,
        annotationNames = SummarizedExperiment::colData(
            testedAnnotations(annotations)
        )$name
    )
}

# The variant GLS forms each block's Omega from its LD, which a scores-only
# LdScore does not carry.
# @noRd
.h2GlsAssertLd <- function(ldRef) {
    if (.ldScoreHasLd(ldRef)) {
        return(invisible(NULL))
    }
    abort(c(
        str_c(
            "The variant GLS needs each LD block's LD, and this LdScore ",
            "carries none."
        ),
        i = str_c(
            "Build it from a genotype panel, or from LdMatrix blocks ",
            "with `keepLdMatrices = TRUE`."
        )
    ))
}

# Assemble the variantGlsUnivariate return (enrichment, candidate scores).
.h2GlsResult <- function(fit, jk, annotations, contrib, A, ridge, w = NULL) {
    .h2FitResult(
        h2 = fit$h2,
        h2Se = jk$h2Se,
        intercept = fit$intercept,
        interceptSe = jk$intSe,
        A = if (is.null(.h2BaselineMat(annotations))) NULL else A,
        tau = fit$tau,
        tauBlocks = jk$annotationJackknifeCoefs,
        scoreStats = .h2GlsCandidateScores(annotations, contrib, fit, ridge),
        w = w
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
# Every eigen likelihood's maximum: projected Fisher scoring
# =============================================================================

# The maximum of the log-likelihood `model` gives `data` from each of `starts`
# (a list of parameter vectors) within [lower, upper], by projected Fisher
# scoring in compiled code (eigenLikelihoodFitCpp(), src/eigen_likelihood.cpp):
# dlib's box-constrained minimizer along a projected Newton direction whose
# Hessian is the expected information. `model` is one of "variance" (b ~ N(0,
# offset + G theta), floored), "tied" (the noise tied to 1 - h2, with a scale
# where fitted), "nullSpace" (the tied model over all n samples), "cross" (a
# pair's second statistics given the first's), "reduced" (a pair on two
# references as independent pairs) and "joint" (a group, on one reference or
# several); `data` holds what it reads. With `free`, the fit is in those
# parameters alone, the rest held at `full`'s. Each iteration is one of
# dlib's, its line search from the whole step. A run stops where the next
# step's predicted rise is below 1e-14 of the log-likelihood's size, or below
# 1e-10 where its value alone is read (`valueOnly`: a null's or a profile's),
# closer than an interval's root or a p-value reads it (at a value of 200,
# 2e-8 moves a root by under 1e-8), or where an iteration does not raise it;
# it has not converged where it reaches `maxit` iterations first. Each
# start's fit: list(par, value, converged, information).
# @noRd
.h2Scoring <- function(
    model,
    data,
    starts,
    lower,
    upper,
    maxit,
    valueOnly = FALSE,
    full = NULL,
    free = NULL
) {
    eigenLikelihoodFitCpp(
        model,
        data,
        map(starts, as.double),
        as.double(lower),
        as.double(upper),
        as.double(full %||% numeric()),
        as.integer(free %||% integer()),
        as.integer(maxit),
        valueOnly
    )
}

# The log-likelihood `model` gives `data` at theta (.h2Scoring()'s models)
# and, where `scoring` asks, its score and expected information
# (eigenLikelihoodCpp()).
# @noRd
.h2Likelihood <- function(model, data, theta, scoring = FALSE) {
    eigenLikelihoodCpp(model, data, as.double(theta), scoring)
}

# The start among `starts`, each kept within [lower, upper], where `model`'s
# log-likelihood is highest, the first where none is finite: a start moved
# far from its fit (one extrapolated past a jump between modes) can sit where
# the conditional variances are floored, where scoring crawls.
# @noRd
.h2BestStart <- function(model, data, starts, lower, upper) {
    boxed <- map(starts, .h2Boxed, lower = lower, upper = upper)
    values <- map_dbl(boxed, .h2LikelihoodValue, model = model, data = data)
    boxed[[which.max(replace(values, !is.finite(values), -Inf))]]
}

# @noRd
.h2Boxed <- function(x, lower, upper) {
    pmin(pmax(x, lower), upper)
}

# @noRd
.h2LikelihoodValue <- function(theta, model, data) {
    .h2Likelihood(model, data, theta)$value
}

# The best of `fits` that converged with a value no worse than `nullValue`,
# within rounding (a fit at the null's boundary is the null); NULL when none
# does.
# @noRd
.h2ScoringBest <- function(fits, nullValue = -Inf) {
    values <- map_dbl(fits, "value")
    floor <- nullValue - 1e-10 * max(abs(nullValue), 1)
    accepted <- map_lgl(fits, "converged") & is.finite(values) & values >= floor
    if (!any(accepted)) {
        return(NULL)
    }
    fits[[which.max(replace(values, !accepted, -Inf))]]
}

# A parameter's standard error, the first's, from the expected information at
# a fit; NA unless the settings ask for it.
# @noRd
.h2ScoringSe <- function(param, fit) {
    if (!identical(param@inference, "information")) {
        return(NA_real_)
    }
    sqrt(solve(fit$information)[[1L, 1L]])
}

# The 95% interval of a coordinate from inverting the likelihood ratio: the
# values whose profile log-likelihood `profile(value, ...)` lies within half
# the chi-square's 95% point on one degree of freedom of `top`, its value at
# the fit `estimate`; each end root-found between the fit and its side of
# `bounds`, the bound itself where the profile never falls that far.
# @noRd
.h2LrInterval <- function(estimate, bounds, top, profile, ...) {
    c(
        .h2LrSide(bounds[[1L]], estimate, top, profile, ...),
        .h2LrSide(bounds[[2L]], estimate, top, profile, ...)
    )
}

# One end of the interval, toward `bound`.
# @noRd
#' @importFrom stats uniroot
.h2LrSide <- function(bound, estimate, top, profile, ...) {
    if (bound == estimate || .h2LrExcess(bound, top, profile, ...) <= 0) {
        return(bound)
    }
    uniroot(
        .h2LrExcess,
        sort(c(estimate, bound)),
        top = top,
        profile = profile,
        ...,
        tol = 1e-8
    )$root
}

# How far the profile at `value` lies below `top`, past the cutoff, capped
# where no fit there converges.
# @noRd
#' @importFrom stats qchisq
.h2LrExcess <- function(value, top, profile, ...) {
    min(top - profile(value, ...) - qchisq(0.95, 1) / 2, 1e10)
}

# The starts of a variance model's fit in (h2 per annotation..., intercept):
# the eigen regression's solution where the procedures take it, then every
# pair of their h2 and intercept starts, each annotation's h2 at the pair's.
# @noRd
.h2VarianceStarts <- function(proc, data, nTau, lower, upper) {
    grid <- expand.grid(intercept = proc@interceptStarts, h2 = proc@h2Starts)
    c(
        if (proc@regressionStart) {
            list(.h2VarianceRegressionStart(data, lower, upper))
        },
        map2(grid$h2, grid$intercept, .h2VarianceGridStart, nTau = nTau)
    )
}

# @noRd
.h2VarianceGridStart <- function(h2, intercept, nTau) {
    c(rep(h2, nTau), intercept)
}

# The eigen regression's solution as a start: the squared statistics less the
# offset on G by least squares, then weighted at the variance that fit
# implies, kept in the box.
# @noRd
.h2VarianceRegressionStart <- function(data, lower, upper) {
    y <- data$b2 - data$offset
    first <- .h2LeastSquares(data$G, y)
    v <- pmax(data$offset + as.vector(data$G %*% first), data$floor)
    pmin(pmax(.h2LeastSquares(data$G / v, y / v), lower), upper)
}

# Least-squares coefficients of y on X, pivoted, a column that cannot be
# fitted at 0.
# @noRd
.h2LeastSquares <- function(X, y) {
    coalesce(qr.coef(qr(X), y), 0)
}

# =============================================================================
# The eigen likelihood: the genome-wide design and fit
# =============================================================================

# Per-block quantities of the eigen likelihood: eigenvalues, rotated
# N-normalized effects bstar = V'(z/sqrt(n)), stratified eigenvalue scores
# l_{ia} (NULL when unstratified) and the candidate annotations' scores, from
# the block's decomposition. The zero-eigenvalue directions a decomposition from
# a panel's factor does not list are added with lam = 0. Under frequency weights
# `w` the pseudo-eigenpairs stand in (.h2WeightedBlock()), with the weighted
# annotation matrices.
.h2EigenBlockData <- function(
    block,
    eigenRef,
    z,
    n,
    baselineMat,
    candMat,
    w = NULL
) {
    dec <- .h2WeightedBlock(.ldEigenBlockAt(block, eigenRef), w)
    zb <- z[dec$snpIdx] / sqrt(n)
    bstar <- as.vector(crossprod(dec$vectors, zb))
    list(
        lam = c(dec$values, numeric(dec$nNull)),
        bstar = c(bstar, .h2UnlistedBstar(zb, bstar, dec$nNull)),
        ldAnnot = .h2EigenAnnot(dec, baselineMat),
        ldCand = .h2EigenAnnot(dec, candMat)
    )
}

# The rotated effects along a block's unlisted zero-eigenvalue directions. The
# likelihood gives every such direction the same modelled variance, its floor,
# so the likelihood sees only their total sum of squares: the part of z outside
# the listed directions, spread evenly over them, leaves every estimate as a
# full decomposition would give it. With frequency weights, whose
# pseudo-eigenvectors are not orthonormal, the share is not that part, but no
# estimate reads it: no parameter reaches the floor.
.h2UnlistedBstar <- function(zb, bstar, nNull) {
    if (nNull == 0L) {
        return(numeric())
    }
    rep(sqrt(max(sum(zb^2) - sum(bstar^2), 0) / nNull), nNull)
}

# Genome-wide design: eigenvalues, rotated effects, block ids and the scaled
# eigenvalue-score matrix, with the parametrization `basis`
# (.h2AnnotationBasis()). Univariate uses a single all-ones score column
# (l_{i,base} = 1) scaled by total M, h2 within the procedures' bounds.
.h2EigenDesign <- function(blockData, M, basis, proc) {
    lam <- list_c(map(blockData, "lam"))
    fit <- basis %||%
        list(
            transform = matrix(1),
            size = M,
            share = 1,
            lower = proc@h2Bounds[[1L]],
            upper = proc@h2Bounds[[2L]]
        )
    ldAnnot <- if (is.null(basis)) {
        matrix(1, length(lam), 1)
    } else {
        ldAnnotList <- map(blockData, "ldAnnot")
        exec(rbind, !!!ldAnnotList)
    }
    list_assign(
        fit,
        lam = lam,
        bstar = list_c(map(blockData, "bstar")),
        blockId = rep(seq_along(blockData), lengths(map(blockData, "lam"))),
        ldAnnotScaled = sweep(ldAnnot %*% fit$transform, 2L, fit$size, "/")
    )
}

# How a stratified fit parametrizes the design `A`, whose columns need not read
# as heritabilities: the variance of a direction depends on a column only
# through tau_a l_ia, and a column whose zero its values never reach, like a
# z-score or x + 3, trades its coefficient against the base's without bound.
# Such a column is centered, so the base absorbs its location, except where it
# is part of what spans the all-ones vector: A* = A T with T = I - b c', for
# Ab = 1 and c the columns' means. Then each parameter is tau*_a size_a, with
# size_a the column's total absolute value, so it is on a scale near 1 even at
# a mean of 0; it adds tau*_a sum(A*_a) = parameter * share_a to the total, and
# tau = T tau* recovers the coefficients of `A`. Columns reaching 0 and never
# going below it, binary ones and PIPs among them, stay as they are: the
# parameter is the column's heritability, share 1 and T = I.
#
# Bounds, from the procedures' h2 `bounds` ([0, 1] by default): a parameter
# is within them only if some variant carries that column alone among those
# left uncentered, since then tau_a is that variant's own rate, as for the
# base or a partition's levels. Otherwise it is an excess over another
# column's rate, within plus or minus the bounds' span. A centered column's
# parameter is within twice that: its positive and negative parts cancel,
# and with non-negative variant rates each is at most the total.
# @noRd
.h2AnnotationBasis <- function(A, w = NULL, bounds = c(0, 1)) {
    weighted <- .h2WeightedDesign(A, w)
    ones <- .h2OnesCoef(A)
    anchored <- colSums(A != 0) < nrow(A) & colSums(A < 0) == 0
    centered <- !anchored & abs(ones) < 1e-8
    means <- colSums(weighted) / nrow(A)
    shift <- replace(numeric(ncol(A)), centered, means[centered])
    size <- replace(
        colSums(abs(weighted)),
        centered,
        .h2CenteredSize(A[, centered, drop = FALSE], means[centered], w)
    )
    carried <- .h2CarriedColumns(A[, !centered, drop = FALSE])
    alone <- replace(logical(ncol(A)), !centered, carried)
    span <- bounds[[2L]] - bounds[[1L]]
    list(
        transform = diag(ncol(A)) - tcrossprod(ones, shift),
        size = unname(size),
        share = unname(replace(colSums(weighted) / size, centered, 0)),
        lower = unname(
            if_else(centered, -2 * span, if_else(alone, bounds[[1L]], -span))
        ),
        upper = unname(
            if_else(centered, 2 * span, if_else(alone, bounds[[2L]], span))
        )
    )
}

# The total absolute value of each column of `A` about its mean, each
# variant's term times its frequency weight in `w` where there are weights.
# @noRd
.h2CenteredSize <- function(A, means, w = NULL) {
    colSums(abs(as.matrix(A) - rep(means, each = nrow(A))) * (w %||% 1))
}

# Whether some variant is non-zero on each column of `A` and on no other.
# @noRd
.h2CarriedColumns <- function(A) {
    nonZero <- A != 0
    alone <- rowSums(nonZero) == 1L
    colSums(nonZero[alone, , drop = FALSE]) > 0
}

# The coefficients of the design's columns from its parameters, one set per row
# of `h2a`.
# @noRd
.h2AnnotationTau <- function(h2a, design) {
    sweep(h2a, 2L, design$size, "/") %*% t(design$transform)
}

# Fit the genome-wide eigen likelihood over (h2_1..h2_nTau, int) by Fisher
# scoring (.h2Scoring()), with each h2_a in the design's bounds
# (.h2AnnotationBasis(); the procedures' h2 bounds unstratified) and int in
# theirs ([0, 10] in Ning et al. 2020), from the procedures' starts
# (.h2VarianceStarts()), keeping the best that converges; NA where none does.
# Parametrizing by per-annotation h2 (scale ~1) keeps the information well
# conditioned.
.h2VarianceFit <- function(design, n, nSamples, proc, warnOnBound = TRUE) {
    nTau <- ncol(design$ldAnnotScaled)
    data <- .h2VarianceData(design, n, nSamples, proc)
    lower <- c(design$lower, proc@interceptBounds[[1L]])
    upper <- c(design$upper, proc@interceptBounds[[2L]])
    best <- .h2ScoringBest(.h2Scoring(
        "variance",
        data,
        .h2VarianceStarts(proc, data, nTau, lower, upper),
        lower,
        upper,
        proc@maxit
    ))
    if (is.null(best)) {
        return(list(h2a = rep(NA_real_, nTau), int = NA_real_, h2 = NA_real_))
    }
    if (warnOnBound) {
        .h2WarnIfAtBound(best$par, upper, n, nSamples)
    }
    h2a <- best$par[seq_len(nTau)]
    list(h2a = h2a, int = best$par[[nTau + 1L]], h2 = sum(h2a * design$share))
}

# What the "variance" likelihood reads of a genome-wide design: each
# direction's variance in (h2_1..h2_nTau, int), G = [its scaled eigenvalue
# scores times lam^2, less lam share / nRef; lam / n], its squared rotated
# effect, and the procedures' floor.
# @noRd
.h2VarianceData <- function(design, n, nSamples, proc) {
    list(
        G = cbind(
            design$ldAnnotScaled *
                design$lam^2 -
                tcrossprod(design$lam, design$share) / nSamples,
            design$lam / n
        ),
        b2 = design$bstar^2,
        offset = 0,
        floor = proc@varianceFloor
    )
}

# The eigen likelihood carries a finite-reference term, -lam * h2 / nRef, that
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
.h2WarnIfAtBound <- function(par, upper, n, nSamples) {
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
        "In the eigen likelihood's fit, the fitted {which} sits at the top ",
        "of its range, so the estimate reflects that bound rather than the ",
        "data. This usually ",
        "means the LD reference is too small for the GWAS -- the ",
        "finite-reference correction scales as 1/nSamples, and here ",
        "nSamples = {nSamples} against N = {round(n)}. Compare the other ",
        "estimators, or ",
        "use a larger reference panel."
    ))
    invisible(TRUE)
}

# Baseline modeled variance lamh2_i at the fitted parameters (for the
# candidate annotations' score test).
.h2BaselineVariance <- function(design, ft, n, nSamples, lim = exp(-10)) {
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
# The eigen likelihood's standard errors (jackknife) + enrichment inputs
# =============================================================================

# Delete-one-block jackknife (Ning et al. 2020): refit leaving each block out,
# SE = sqrt(mean((jack - mean)^2) * (nBlocks - 1)). Returns per-annotation and
# total-h2 SEs, the intercept SE, and the LOO per-annotation blocks.
.h2VarianceJackknife <- function(blockData, M, basis, n, nSamples, proc) {
    nBlocks <- length(blockData)
    loo <- map(
        seq_len(nBlocks),
        .h2VarianceLooFit,
        blockData = blockData,
        M = M,
        basis = basis,
        n = n,
        nSamples = nSamples,
        proc = proc
    )
    h2aList <- map(loo, "h2a")
    h2aBlocks <- exec(rbind, !!!h2aList)
    intLoo <- map_dbl(loo, "int")
    list(
        h2aSe = apply(h2aBlocks, 2, .jackknifeSeFromLoo),
        h2Se = .jackknifeSeFromLoo(map_dbl(loo, "h2")),
        intSe = .jackknifeSeFromLoo(intLoo),
        h2aBlocks = h2aBlocks
    )
}

# =============================================================================
# The local engine: one window's heritability under an H2EstimateParam
# =============================================================================
#
# Every per-trait setting fits one model of a window's marginal correlations
# r. Rotated onto the LD eigenvectors, u = Q'r, each direction's variance is
#   s (h2 lam^2 / M - e h2 lam / nRef + c lam / n),
# for the window's M variants, a study of n and a reference of nRef (e = 1
# with the eigenvalue term, 0 without it), the noise term c (tied: 1 - h2;
# free: an intercept) and a scale s on the phenotype's variance (1 unless
# profiled). Each variant's squared correlation has mean h2 l / M + c / n on
# its LD score l in the window. The statistic, the criterion and the
# procedures say how that is fitted: LAVA, HDL-L, FastHer and GREML are
# settings of it (.h2LocalSettings()), and the eigen regression, the variant
# GLS and the variant regression take the genome-wide estimators' criteria
# into one window.

# One window's estimate under `param`, from the window's z (t statistics),
# per-variant sample sizes and LD source (.rssLdSource(), carrying the
# reference's sample size and the study's share in it).
# @noRd
.h2LocalFit <- function(param, z, nVariant, ld) {
    fit <- switch(
        str_c(param@statistic, "/", param@criterion),
        "eigen/moments" = .h2EigenMoments,
        "eigen/likelihood" = .h2EigenLikelihood,
        "eigen/regression" = .h2EigenRegression,
        "variant/regression" = .h2VariantRegression,
        "variant/gls" = .h2VariantGls
    )
    if (is.null(fit)) {
        abort(c(
            "The phenotype statistic needs individual-level data.",
            i = "Estimate from a QtlDataset, or choose another statistic."
        ))
    }
    fit(param, z, nVariant, ld)
}

# The settings a local estimate cannot be made under: a supplied noise term,
# which the local functions take no value for, a jackknife with nothing to
# leave out, and, on summary statistics, the phenotype.
# @noRd
.h2RefuseLocalSetting <- function(param, individual, label) {
    problem <- if (identical(param@noise, "supplied")) {
        "supplied noise needs a value, which a local estimate takes none of"
    } else if (
        identical(param@inference, "jackknife") &&
            !(identical(param@statistic, "variant") &&
                identical(param@criterion, "regression"))
    ) {
        str_c(
            "in one window, only the variant regression has blocks to leave ",
            "out for a jackknife"
        )
    } else if (!individual && identical(param@statistic, "phenotype")) {
        "the phenotype statistic needs individual-level data (a QtlDataset)"
    }
    if (is.null(problem)) {
        return(invisible(NULL))
    }
    abort(glue("{label}: {problem}."))
}

# One window's result; what the settings do not report stays NA.
# @noRd
.cisH2Result <- function(
    nComponents,
    h2 = NA_real_,
    h2Se = NA_real_,
    p = NA_real_,
    intercept = NA_real_,
    scale = NA_real_,
    smallReference = FALSE,
    undefined = FALSE,
    h2Interval = c(NA_real_, NA_real_)
) {
    list(
        h2 = h2,
        h2Se = h2Se,
        h2Lower = h2Interval[[1L]],
        h2Upper = h2Interval[[2L]],
        p = p,
        nComponents = as.integer(nComponents),
        intercept = intercept,
        scale = scale,
        smallReference = smallReference,
        undefined = undefined
    )
}

# The per-variant sample sizes summarized as the settings say.
# @noRd
.h2SampleSize <- function(how, nVariant) {
    switch(
        how,
        mean = mean(nVariant),
        median = median(nVariant),
        harmonic = 1 / mean(1 / nVariant)
    )
}

# Whether the settings read the LD reference's sample size.
# @noRd
.h2ReadsReferenceSize <- function(param) {
    is_in(param@referenceNoise, c("overlapFactor", "eigenTerm", "unbiasedR2"))
}

# ---- the eigen statistic -----------------------------------------------------

# The window's LD eigenpairs, largest first: from the thin SVD of a factor,
# or the matrix's own decomposition.
# @noRd
.h2Eigen <- function(ld) {
    if (!is.null(ld$ldFactor)) {
        return(.ldFactorEigen(ld$ldFactor))
    }
    eigen(ld$R, symmetric = TRUE)
}

# The directions `truncation` keeps of `dec` (largest first) for a study of
# n and a window of M variants: those that first reach its fraction of the
# positive eigenvalues' sum or of M (every one when they never do, and at a
# fraction of 1 every direction, null ones included), those tied with the
# last if asked, and at most maxPerSample n of them, never fewer than
# `minimum` by that cap.
# @noRd
.h2Truncate <- function(dec, truncation, n, M) {
    values <- dec$values
    candidates <- if (
        identical(truncation@of, "eigenvalues") &&
            truncation@fraction < 1
    ) {
        which(values > 0)
    } else {
        seq_along(values)
    }
    K <- .h2TruncationCount(values[candidates], truncation, M)
    reached <- candidates[seq_len(K)]
    tied <- if (truncation@keepTies && K > 0L) {
        candidates[values[candidates] >= values[[reached[[K]]]]]
    } else {
        reached
    }
    cap <- if (is.null(truncation@maxPerSample)) {
        length(tied)
    } else {
        max(floor(truncation@maxPerSample * n), truncation@minimum)
    }
    keep <- tied[seq_len(min(length(tied), cap))]
    list(values = values[keep], vectors = dec$vectors[, keep, drop = FALSE])
}

# How many of the candidate eigenvalues `values` reach the fraction.
# @noRd
.h2TruncationCount <- function(values, truncation, M) {
    if (truncation@fraction >= 1) {
        return(length(values))
    }
    total <- if (identical(truncation@of, "eigenvalues")) sum(values) else M
    reached <- which(cumsum(values) / total >= truncation@fraction)
    if (length(reached) == 0L) length(values) else reached[[1L]]
}

# ---- the frequency exponent --------------------------------------------------

# The variants' weights under the frequency exponent alpha: Var(beta_j) =
# h2 w_j / M with w_j = [2 p_j (1 - p_j)]^(1 + alpha) scaled to mean 1 over
# the window, none for a variant the reference holds monomorphic; NULL at
# alpha = -1, where every variant weighs alike (the per-SD model).
# @noRd
.h2FrequencyWeights <- function(het, alpha) {
    if (alpha == -1) {
        return(NULL)
    }
    raw <- if_else(het > 0, het^(1 + alpha), 0)
    raw / mean(raw)
}

# A window's weights under the settings, from its LD source's
# heterozygosities 2p(1 - p), which the readers add where alpha asks for
# them (.h2WithHeterozygosity(), a QtlDataset's block frequencies). A bare
# genotype matrix, as a weight function's in-sample estimate reads, has
# none: it may be residualized, and its frequencies are not its own.
# @noRd
.h2SourceWeights <- function(param, ld) {
    if (param@alpha != -1 && is.null(ld$het)) {
        abort(c(
            glue(
                "alpha = {param@alpha} weighs each variant by its allele ",
                "frequency, which a genotype matrix alone does not record."
            ),
            i = str_c(
                "Estimate the heritability with estimateLocalH2() on a ",
                "QtlDataset and pass it as `h2Estimate`."
            )
        ))
    }
    .h2FrequencyWeights(ld$het, param@alpha)
}

# The kept eigenpairs under frequency weights `w` (none: as they are). With
# Var(r) = h2 R W R / M + c R, the whitened statistics Lambda^-1/2 Q'r have
# covariance h2 B + c I, B = Lambda^1/2 Q'WQ Lambda^1/2 / M, so in B's
# eigenvectors U each direction is independent again, of variance h2 b + c.
# The pseudo-eigenpairs lambda = M b, with vectors Q Lambda^-1/2 U
# lambda^1/2, give those statistics the per-SD form h2 lambda^2 / M + c
# lambda that every eigen setting fits; a direction B leaves empty is
# dropped. The plain pairs, the rotation and the weights are kept for two
# references' cross shapes.
# @noRd
.h2Weighted <- function(dec, w) {
    if (is.null(w)) {
        return(dec)
    }
    M <- length(w)
    roots <- sqrt(dec$values)
    B <- crossprod(dec$vectors * sqrt(w)) * tcrossprod(roots) / M
    eb <- eigen(B, symmetric = TRUE)
    keep <- eb$values > max(eb$values) * 1e-10
    values <- M * eb$values[keep]
    U <- eb$vectors[, keep, drop = FALSE]
    list(
        values = values,
        vectors = sweep(dec$vectors %*% (U / roots), 2L, sqrt(values), `*`),
        plain = dec,
        rotation = U,
        weights = w
    )
}

# The window's kept eigenpairs under the settings: its LD's, truncated
# (.h2Truncate()), under the frequency weights alpha gives (.h2Weighted()).
# @noRd
.h2Basis <- function(param, ld, n, M) {
    .h2Weighted(
        .h2Truncate(.h2Eigen(ld), param@truncation, n, M),
        .h2SourceWeights(param, ld)
    )
}

# A genome-wide block's decomposition under the reference's frequency weights
# `w` (none: as it is): its directions' pseudo-eigenpairs (.h2Weighted())
# under the block's share of `w`, with the loadings Q Lambda^1/2 U b^-1/2 on
# which a direction's annotation scores sum_k L_ki^2 w_k a_k are taken
# (.h2EigenAnnot()); under them a direction's genetic variance is its
# scores' sum times b_i^2. The weights are scaled over the whole reference,
# so its blocks' heritabilities add up. A direction with no eigenvalue above
# LDER's 1e-6, or one the weights leave empty, joins the null ones: no
# genetic variance reaches it.
# @noRd
.h2WeightedBlock <- function(dec, w) {
    if (is.null(w)) {
        return(dec)
    }
    positive <- dec$values > 1e-6
    plain <- list(
        values = dec$values[positive],
        vectors = dec$vectors[, positive, drop = FALSE]
    )
    pseudo <- .h2Weighted(plain, w[dec$snpIdx])
    loadings <- plain$vectors %*% (pseudo$rotation * sqrt(plain$values))
    list(
        values = pseudo$values,
        vectors = pseudo$vectors,
        loadings = sweep(loadings, 2L, sqrt(pseudo$values), "/"),
        snpIdx = dec$snpIdx,
        nNull = dec$nNull + length(dec$values) - length(pseudo$values)
    )
}

# A design's or annotation matrix's columns under the frequency weights `w`,
# each variant's row times its weight: the model's per-variant variance is
# w_j sum_a tau_a A_ja. As they are without weights.
# @noRd
.h2WeightedDesign <- function(A, w) {
    if (is.null(A) || is.null(w)) {
        return(A)
    }
    A * w
}

# A genome-wide reference's variants' weights under the settings' frequency
# exponent (.h2FrequencyWeights()), scaled to mean 1 over the whole
# reference; NULL at alpha = -1.
# @noRd
.h2RefWeights <- function(param, ldRef) {
    if (param@alpha == -1) {
        return(NULL)
    }
    .h2FrequencyWeights(.h2RefHeterozygosity(ldRef), param@alpha)
}

# Each of a genome-wide reference's variants' heterozygosity 2p(1 - p): from
# the MAF it records (the AF of the LdMatrix blocks or the sketch it was built
# from), which serves as well as the AF, or else from its genotype panel's:
# the frequencies the panel records, or its dosages read a block at a time.
# @noRd
.h2RefHeterozygosity <- function(ldRef) {
    maf <- S4Vectors::mcols(ldRef)[["MAF"]]
    if (!is.null(maf)) {
        return(2 * maf * (1 - maf))
    }
    panel <- ldReference(ldRef)
    if (is.null(panel)) {
        abort(c(
            str_c(
                "A frequency exponent (alpha) weighs each variant by its ",
                "allele frequency, and this reference records none."
            ),
            i = str_c(
                "Build it from LdMatrix blocks with an AF column, or from a ",
                "genotype panel."
            )
        ))
    }
    af <- .ldRefAf(panel)
    p <- if (is.null(af)) {
        .h2PanelFrequencies(ldRef, panel)
    } else {
        af[S4Vectors::mcols(ldRef)[["panelRow"]]]
    }
    2 * p * (1 - p)
}

# The A1 frequency of each of a reference's variants from its panel's
# dosages, a block at a time.
# @noRd
.h2PanelFrequencies <- function(ldRef, panel) {
    .h2ScatterBlockRows(
        matrix(0, length(ldRef), 1L),
        map(.ldStatLazyBlocks(ldRef), .h2PanelBlockFrequencies, panel = panel)
    )[, 1L]
}

# @noRd
.h2PanelBlockFrequencies <- function(block, panel) {
    dosage <- .ldRefDosage(panel, block$rows, meanImpute = TRUE)
    list(idx = block$idx, values = colMeans(dosage) / 2)
}

# ---- the moments -------------------------------------------------------------

# The moments (LAVA's univariate estimator, written from Werme et al. 2022,
# Methods, and LAVA's documentation: LAVA's code is not licensed for reuse,
# and none of it is used here). The correlations r are projected onto the K
# kept standardized directions, delta = Lambda_K^-1/2 Q_K' r, for the
# summarized sample size n, with no estimate from fewer than `minimum`
# directions or with no degrees of freedom left. With phi how far the
# reference's LD lies from the study's own (.h2OverlapFactor()), the noise per
# direction is (1 - h2 + h2 phi) / (n - 1), which gives
#   h2 = (delta'delta - m) / v, m = K / (n - 1), v = 1 + m (phi - 1).
# At phi = 0 that is LAVA's adjusted R^2 exactly; for a reference sharing no
# samples it matches, to first order, the correction (nRef - K - 1) / (nRef -
# 1) that LAVA's release notes describe. LAVA's own noise term, (1 -
# delta'delta) / (n - K - 1), holds at phi = 0 alone: past it, it leaves h2 m
# in the estimate, and it has no answer once delta'delta reaches 1, as a
# reference's components can with K / n large. Without the overlap factor,
# phi is 1: a reference with no sampling noise of its own. p is .h2MomentsP()'s.
# LAVA withholds h2 where K / n exceeds 0.1; it is reported here at any K,
# since a cis window's K usually does.
# @noRd
.h2EigenMoments <- function(param, z, nVariant, ld) {
    n <- .h2SampleSize(param@sampleSize, nVariant)
    dec <- .h2Basis(param, ld, n, length(z))
    K <- length(dec$values)
    if (K < param@truncation@minimum || n <= K + 1L) {
        return(.cisH2Result(K))
    }
    delta <- crossprod(dec$vectors, .zToCor(z, nVariant)) / sqrt(dec$values)
    dtd <- sum(delta^2)
    m <- K / (n - 1)
    factor <- identical(param@referenceNoise, "overlapFactor")
    overlap <- ld$referenceOverlap %||% 0
    phi <- if (factor) .h2OverlapFactor(n, ld$nSamples, overlap) else 1
    v <- 1 + m * (phi - 1)
    .cisH2Result(
        K,
        h2 = (dtd - m) / v,
        p = .h2MomentsP(dtd, m, v, n),
        smallReference = factor && overlap == 0 && ld$nSamples <= 2 * n
    )
}

# How far a reference's LD lies from the study's own, in units of its sampling
# noise: n Var(r_study - r_reference), for n study samples, a reference of nRef
# and the share `overlap` of the study's samples in it. 0 when the reference is
# the study, about 1 for a panel of thousands, whatever the two share.
# Divergence between the populations adds to it, unmodelled: measuring that
# takes the study's own LD.
# @noRd
.h2OverlapFactor <- function(n, nRef, overlap) {
    1 + (1 - 2 * overlap) * n / nRef
}

# The test of delta = 0. Under it delta'delta = y'Ay / y'y, a ratio of quadratic
# forms in the trait, with mean m and variance about 2 m v / (n + 1): a Beta
# with those moments, at phi = 0 exactly the F test of the adjusted R^2.
# NA where they fit no Beta, at a very small n.
# @noRd
#' @importFrom stats pbeta
.h2MomentsP <- function(dtd, m, v, n) {
    s <- (1 - m) * (n + 1) / (2 * v) - 1
    if (s <= 0) {
        return(NA_real_)
    }
    pbeta(dtd, m * s, (1 - m) * s, lower.tail = FALSE)
}

# ---- the likelihood ----------------------------------------------------------

# The Gaussian likelihood of the rotated correlations bstar = Q'r over the
# kept directions, with the noise free (an intercept) or tied (1 - h2): the
# data each fit reads.
# @noRd
.h2EigenLikelihood <- function(param, z, nVariant, ld) {
    n <- .h2SampleSize(param@sampleSize, nVariant)
    M <- length(z)
    dec <- .h2Basis(param, ld, n, M)
    if (length(dec$values) < param@truncation@minimum) {
        return(.cisH2Result(length(dec$values)))
    }
    r <- .zToCor(z, nVariant)
    eigenTerm <- identical(param@referenceNoise, "eigenTerm")
    data <- list(
        lam = dec$values,
        bstar = as.vector(crossprod(dec$vectors, r)),
        n = n,
        nRef = if (eigenTerm) ld$nSamples else Inf,
        M = M,
        floor = param@procedures@varianceFloor
    )
    if (identical(param@noise, "free")) {
        return(.h2LikelihoodFree(param, data))
    }
    if (param@procedures@nullSpace) {
        return(.h2NullSpaceFit(param, data))
    }
    .h2TiedFit(param, data)
}

# The likelihood-ratio p of a fit against no heritability, halved where the
# settings refer it to the boundary mixture; NA unless they ask for it.
# @noRd
#' @importFrom stats pchisq
.h2LikelihoodP <- function(param, value, nullValue) {
    if (!identical(param@inference, "lrt")) {
        return(NA_real_)
    }
    p <- pchisq(2 * (value - nullValue), 1, lower.tail = FALSE)
    if (param@procedures@halfBoundaryP) p / 2 else p
}

# Free noise: h2 and the intercept together
# by Fisher scoring (.h2Scoring()) from every start the procedures give
# (.h2VarianceStarts()), keeping the best fit that converges and, where the
# procedures say, beats the intercept-only fit (.h2WindowNull()); every
# estimate is NA when none does.
# @noRd
.h2LikelihoodFree <- function(param, data) {
    proc <- param@procedures
    K <- length(data$lam)
    full <- .h2WindowData(data)
    lower <- c(proc@h2Bounds[[1L]], proc@interceptBounds[[1L]])
    upper <- c(proc@h2Bounds[[2L]], proc@interceptBounds[[2L]])
    null <- .h2WindowNull(data, proc)
    fits <- .h2Scoring(
        "variance",
        full,
        .h2VarianceStarts(proc, full, 1L, lower, upper),
        lower,
        upper,
        proc@maxit
    )
    beat <- if (proc@mustBeatNull) null$value %||% Inf else -Inf
    best <- .h2ScoringBest(fits, beat)
    if (is.null(best) || is.null(null) || .h2AllFloored(best$par, full)) {
        return(.cisH2Result(K))
    }
    .cisH2Result(
        K,
        h2 = best$par[[1L]],
        h2Se = .h2ScoringSe(param, best),
        p = .h2LikelihoodP(param, best$value, null$value),
        intercept = best$par[[2L]],
        h2Interval = .h2ProfileInterval(param, data, best$par[[1L]])
    )
}

# What the "variance" likelihood reads of one window's rotated statistics: their
# variance in (h2, intercept), G = [lam^2 / M - lam / nRef, lam / n]; with h2
# held at `h2`, the intercept's column alone; with the intercept held at
# `intercept`, h2's alone.
# @noRd
.h2WindowData <- function(data, h2 = NULL, intercept = NULL) {
    genetic <- data$lam^2 / data$M - data$lam / data$nRef
    noise <- data$lam / data$n
    G <- if (!is.null(h2)) {
        cbind(noise)
    } else if (!is.null(intercept)) {
        cbind(genetic)
    } else {
        cbind(genetic, noise)
    }
    list(
        G = G,
        b2 = data$bstar^2,
        offset = (h2 %||% 0) * genetic + (intercept %||% 0) * noise,
        floor = data$floor
    )
}

# Whether every direction's variance at theta sits at the floor: a fit with
# no direction left to read, as where every statistic is 0.
# @noRd
.h2AllFloored <- function(theta, data) {
    all(data$offset + as.vector(data$G %*% theta) <= data$floor)
}

# The fit at h2 = `h2`, the intercept refitted from the procedures' intercept
# starts, the best that converges; NULL when none does. At 0, the
# intercept-only null.
# @noRd
.h2WindowNull <- function(data, proc, h2 = 0) {
    .h2ScoringBest(.h2Scoring(
        "variance",
        .h2WindowData(data, h2 = h2),
        as.list(unique(proc@interceptStarts)),
        proc@interceptBounds[[1L]],
        proc@interceptBounds[[2L]],
        proc@maxit,
        valueOnly = TRUE
    ))
}

# The 95% interval for h2 from inverting the likelihood ratio, its profile
# the intercept refitted within its bounds (.h2WindowProfile(),
# .h2LrInterval()). NA unless the settings ask for the likelihood-ratio test.
# @noRd
.h2ProfileInterval <- function(param, data, estimate) {
    if (!identical(param@inference, "lrt")) {
        return(c(NA_real_, NA_real_))
    }
    proc <- param@procedures
    .h2LrInterval(
        estimate,
        proc@h2Bounds,
        .h2WindowProfile(estimate, data, proc),
        .h2WindowProfile,
        data = data,
        proc = proc
    )
}

# The log-likelihood at h2 with the intercept refitted within its bounds.
# @noRd
.h2WindowProfile <- function(h2, data, proc) {
    .h2WindowNull(data, proc, h2)$value %||% -Inf
}


# The modelled variance of the rotated correlations at heritability `h2` and
# `intercept`: h2 / M lam^2 - h2 lam / nRef + intercept lam / n, floored.
# @noRd
.h2WindowVariance <- function(h2, data, intercept) {
    pmax(
        h2 /
            data$M *
            data$lam^2 -
            h2 * data$lam / data$nRef +
            intercept * data$lam / data$n,
        data$floor
    )
}

# Tied noise without the null space: h2 alone, or with a scale on the
# phenotypic variance, by Fisher scoring (.h2Scoring()) within the
# procedures' bounds from their h2 starts (the scale from 1); the
# likelihood ratio against h2 = 0, the scale refitted.
# @noRd
.h2TiedFit <- function(param, data) {
    proc <- param@procedures
    scaled <- proc@profileScale
    model <- list(
        offset = data$lam / data$n,
        g = data$lam^2 / data$M - data$lam / data$nRef - data$lam / data$n,
        b2 = data$bstar^2,
        floor = data$floor,
        scaled = scaled
    )
    fit <- .h2TiedScoring("tied", model, proc)
    if (is.null(fit$best)) {
        return(.cisH2Result(length(data$lam)))
    }
    .cisH2Result(
        length(data$lam),
        h2 = fit$best$par[[1L]],
        h2Se = .h2ScoringSe(param, fit$best),
        p = .h2LikelihoodP(param, fit$best$value, fit$null),
        scale = if (scaled) fit$best$par[[2L]] else NA_real_
    )
}

# A tied fit of the likelihood `likelihood` ("tied" or "nullSpace") by
# Fisher scoring in (h2, and the scale where `model$scaled`) from the
# procedures' h2 starts, the scale from 1, with the log-likelihood at h2 = 0,
# the scale refitted there: list(best, null).
# @noRd
.h2TiedScoring <- function(likelihood, model, proc) {
    scaled <- model$scaled
    lower <- c(proc@h2Bounds[[1L]], if (scaled) 1e-8)
    upper <- c(proc@h2Bounds[[2L]], if (scaled) Inf)
    starts <- map(unique(proc@h2Starts), c, if (scaled) 1)
    best <- .h2ScoringBest(
        .h2Scoring(likelihood, model, starts, lower, upper, proc@maxit)
    )
    null <- if (scaled) {
        .h2ScoringBest(.h2Scoring(
            likelihood,
            model,
            list(1),
            1e-8,
            Inf,
            proc@maxit,
            valueOnly = TRUE,
            full = c(0, 1),
            free = 2L
        ))$value
    } else {
        .h2Likelihood(likelihood, model, 0)$value
    }
    list(best = best, null = null %||% NA_real_)
}

# Tied noise with the null space: the n study samples' y ~ N(0, d (h2 X X' /
# m + (1 - h2) I)) over the window's m variants (FastHer's likelihood,
# Svishcheva 2025, which HEELS, Li, Mazumder and Lin 2023, reaches by
# iteration), following FastHer's R code (MIT; see inst/COPYRIGHTS) for the
# likelihood, with d fitted beside h2 by Fisher scoring (.h2Scoring()) or
# held at 1. In the eigenbasis of the LD, X'X = n R, X'y = n r and y'y = n;
# with a_k = n (Q'r)_k^2 / lambda_k and S their sum, the likelihood needs r <
# n and S < n, which LD from the study's own samples guarantees; the estimate
# is NA otherwise, as it can be with LD from outside the study. The standard
# error is FastHer's, from the expected information in (h2, d), and so holds
# under the infinitesimal model that likelihood assumes.
# @noRd
.h2NullSpaceFit <- function(param, data) {
    n <- data$n
    a <- n * data$bstar^2 / data$lam
    r <- length(a)
    if (r >= n || sum(a) >= n) {
        return(.cisH2Result(r, undefined = TRUE))
    }
    proc <- param@procedures
    model <- list(
        a = a,
        ck = n * data$lam / data$M - 1,
        n = n,
        scaled = proc@profileScale
    )
    fit <- .h2TiedScoring("nullSpace", model, proc)
    if (is.null(fit$best)) {
        return(.cisH2Result(r))
    }
    .cisH2Result(
        r,
        h2 = fit$best$par[[1L]],
        h2Se = .h2ScoringSe(param, fit$best),
        p = .h2LikelihoodP(param, fit$best$value, fit$null),
        scale = if (model$scaled) fit$best$par[[2L]] else NA_real_
    )
}

# ---- the regressions ---------------------------------------------------------

# The eigen regression (LDER's model, Song et al. 2022, fitted in one
# window): each kept direction's u^2 / lambda on lambda / M (less 1 / nRef
# with the eigenvalue term) and 1 / n, by weighted least squares
# (.h2Regression()), its coefficients h2 and the noise term.
# @noRd
.h2EigenRegression <- function(param, z, nVariant, ld) {
    n <- .h2SampleSize(param@sampleSize, nVariant)
    M <- length(z)
    dec <- .h2Basis(param, ld, n, M)
    lam <- dec$values
    if (length(lam) < max(param@truncation@minimum, 3L)) {
        return(.cisH2Result(length(lam)))
    }
    nRef <- if (identical(param@referenceNoise, "eigenTerm")) ld$nSamples
    y <- as.vector(crossprod(dec$vectors, .zToCor(z, nVariant)))^2 / lam
    fit <- .h2Regression(
        param,
        list(
            y = y,
            X = cbind(lam / M - 1 / (nRef %||% Inf), 1 / n),
            chi2 = n * y,
            n = n,
            taper = .h2Taper(param@procedures@eigenvalueTaper, lam),
            redundancy = rep(1, length(y)),
            lam = lam,
            weights = "variance"
        )
    )
    .cisH2Result(
        length(lam),
        h2 = fit$h2,
        h2Se = fit$h2Se,
        intercept = fit$intercept
    )
}

# The variant regression (S-LDSC's model, as MESC fits it per gene): each
# variant's n r^2 on n l / M and 1, its LD score l in the window from the
# unbiased r^2 where the settings say, with a jackknife over contiguous blocks
# of the window where they ask for one.
# @noRd
.h2VariantRegression <- function(param, z, nVariant, ld) {
    n <- .h2SampleSize(param@sampleSize, nVariant)
    M <- length(z)
    if (M < 3L) {
        return(.cisH2Result(NA_integer_))
    }
    l <- .h2WindowLdScores(param, ld, M)
    y <- n * .zToCor(z, nVariant)^2
    wLd <- .h2RedundancyScores(param, ld, l)
    rows <- list(
        y = y,
        X = cbind(n * l / M, 1),
        chi2 = y,
        n = n,
        taper = rep(1, M),
        redundancy = 1 / pmax(wLd, 1),
        ldTotal = l,
        wLd = wLd,
        M = M,
        hsqScale = 1,
        weights = "ldScore",
        index = seq_len(M),
        ldModel = list(ld = ld, w = .h2SourceWeights(param, ld), n = n, M = M)
    )
    fit <- .h2Regression(param, rows)
    .cisH2Result(
        NA_integer_,
        h2 = fit$h2,
        h2Se = if (identical(param@inference, "jackknife")) {
            .h2RegressionJackknife(param, rows)
        } else {
            fit$h2Se
        },
        intercept = fit$intercept
    )
}

# A window's LD scores sum_k w_k r_jk^2, w the settings' frequency weights
# (each 1 without them) unless `weights` gives others, the unbiased r^2's
# where the settings say.
# @noRd
.h2WindowLdScores <- function(
    param,
    ld,
    M,
    weights = .h2SourceWeights(param, ld) %||% rep(1, M)
) {
    raw <- .ldSquaredProduct(ld, matrix(weights, M, 1L))
    if (!identical(param@referenceNoise, "unbiasedR2")) {
        return(as.vector(raw))
    }
    as.vector(.ldScoreUnbiasedSum(
        raw,
        sum(weights),
        .ldScoreAssertN(ld$nSamples)
    ))
}

# The LD scores the LD-redundancy weight divides by, where the procedures
# take it (each 1 otherwise): the window's with every variant counted once,
# as ldsc's w_ld counts how many of the regression's variants tag each one,
# whatever frequency weights the regression's own scores `l` carry (they are
# those at alpha = -1).
# @noRd
.h2RedundancyScores <- function(param, ld, l) {
    M <- length(l)
    if (!param@procedures@ldRedundancy) {
        return(rep(1, M))
    }
    if (param@alpha == -1) {
        return(l)
    }
    .h2WindowLdScores(param, ld, M, rep(1, M))
}

# The eigenvalue taper on each direction's weight.
# @noRd
.h2Taper <- function(taper, lam) {
    switch(
        taper,
        none = rep(1, length(lam)),
        linear = pmin(1, pmax(lam, 0)),
        logistic = 1 / (1 + exp(-5 * (lam - 1)))
    )
}

# A weighted regression of `rows$y` on `rows$X` = [the genetic columns, the
# noise column] as the settings say, on the rows the chi-square filter keeps
# (.h2Chi2Kept()). `rows` also carries each row's chi-square, the sample
# size, the taper and the LD-redundancy factor, and how its weights read the
# fit (.h2RegressionWeights()). The heritability is the genetic
# coefficients' sum: each column is scaled so that its coefficient is its
# annotation's heritability.
# @noRd
.h2Regression <- function(param, rows) {
    proc <- param@procedures
    kept <- .h2SubsetRows(
        rows,
        .h2Chi2Kept(proc, rows$chi2, rows$n, ncol(rows$X))
    )
    fit <- .h2RegressionFit(param, kept)
    k <- ncol(rows$X)
    list(
        coef = fit$coef,
        h2 = .h2RowsH2(rbind(fit$coef), rows),
        h2Se = .h2RegressionSe(param, kept, fit),
        intercept = if (identical(param@noise, "tied")) {
            NA_real_
        } else {
            fit$coef[[k]]
        },
        fit = fit,
        kept = kept
    )
}

# The fit the noise settings call for, from the first weights on these rows
# (.h2FirstWeights()): in two stages, the noise from the first held in the
# second; free; tied to 1 - h2; or floored, held at the floor where the free
# fit falls below it.
# @noRd
.h2RegressionFit <- function(param, rows) {
    proc <- param@procedures
    w0 <- .h2FirstWeights(param, rows)
    if (proc@twoStage) {
        noise <- .h2FlooredNoise(param, .h2StageNoise(param, rows, w0))
        return(.h2SolveRegression(param, rows, noise, w0))
    }
    free <- .h2SolveRegression(param, rows, NULL, w0)
    noiseAt <- ncol(rows$X)
    if (
        !identical(param@noise, "floored") ||
            free$coef[[noiseAt]] >= proc@interceptFloor
    ) {
        return(free)
    }
    .h2SolveRegression(param, rows, proc@interceptFloor, w0)
}

.h2SubsetRows <- function(rows, keep) {
    list_modify(
        rows,
        y = rows$y[keep],
        X = rows$X[keep, , drop = FALSE],
        chi2 = rows$chi2[keep],
        taper = rows$taper[keep],
        redundancy = rows$redundancy[keep],
        index = rows$index[keep],
        base = rows$base[keep],
        lam = rows$lam[keep],
        ldTotal = rows$ldTotal[keep],
        wLd = rows$wLd[keep],
        block = rows$block[keep]
    )
}

# The rows the chi-square filter keeps, every row without one or where too
# few survive for a design `width` columns wide.
# @noRd
.h2Chi2Kept <- function(proc, chi2, n, width = 2L) {
    if (is.null(proc@chi2Filter)) {
        return(rep(TRUE, length(chi2)))
    }
    kept <- chi2 <= max(proc@chi2Filter[[1L]] * n, proc@chi2Filter[[2L]])
    if (sum(kept) < width + 2L) rep(TRUE, length(chi2)) else kept
}

# The noise term from the first of two stages: the free fit on the rows at or
# below the chi-square threshold (every row where fewer than four are),
# .h2Chi2Threshold()'s where none is set, from the first weights `w0`.
# @noRd
.h2StageNoise <- function(param, rows, w0) {
    threshold <- param@procedures@chi2Threshold %||%
        .h2Chi2Threshold(param, rows)
    below <- rows$chi2 <= threshold
    stage <- if (sum(below) < 4L) rep(TRUE, length(below)) else below
    fit <- .h2SolveRegression(
        param,
        .h2SubsetRows(rows, stage),
        NULL,
        w0[stage]
    )
    fit$coef[[ncol(rows$X)]]
}

# The default chi-square threshold, LDER's: 3.5 sqrt(m / N) + 2.5, or the
# chi-square's top 5 / m where a rough free fit's inflation reaches 1.05; the
# rough fit takes the linear eigenvalue taper on eigen rows, as LDER's does.
# @noRd
#' @importFrom stats quantile
.h2Chi2Threshold <- function(param, rows) {
    m <- length(rows$y)
    rough <- if (is.null(rows$lam)) {
        rows
    } else {
        list_assign(rows, taper = .h2Taper("linear", pmax(rows$lam, 0)))
    }
    inflation <- .h2SolveRegression(
        param,
        rough,
        NULL,
        .h2FirstWeights(param, rough)
    )$coef[[ncol(rows$X)]]
    if (inflation < 1.05) {
        return(3.5 * sqrt(m / rows$n) + 2.5)
    }
    quantile(rows$chi2, 1 - 5 / m)
}

# The noise term held where the settings floor it.
# @noRd
.h2FlooredNoise <- function(param, noise) {
    if (!identical(param@noise, "floored")) {
        return(noise)
    }
    max(noise, param@procedures@interceptFloor)
}

# The weighted fit of `rows` with the noise term `noise` held (NULL: fitted,
# or tied to 1 - h2 with tied noise): the coefficients (h2, noise), the final
# weights and the design they were fitted on, its first column h2's.
# @noRd
.h2SolveRegression <- function(
    param,
    rows,
    noise,
    w0 = NULL,
    tied = identical(param@noise, "tied") && is.null(noise)
) {
    design <- .h2RegressionDesign(rows, noise, tied)
    fit <- .h2Irwls(param, rows, design, noise, tied, w0)
    list(
        coef = .h2DesignCoef(fit$coef, noise, tied),
        w = fit$w,
        design = design
    )
}

# Iteratively reweighted least squares: the first weights, then
# .h2WeightUpdates() updates, each followed by a solve, stopping early on
# convergence when the settings iterate to it. The weights evolve with the
# fit, so this stays a loop.
# @noRd
.h2Irwls <- function(param, rows, design, noise, tied, w0 = NULL) {
    proc <- param@procedures
    w <- w0 %||% .h2FirstWeights(param, rows)
    coef <- .h2WlsCoef(design$X, design$y, w)
    for (i in seq_len(.h2WeightUpdates(proc))) {
        w <- .h2RegressionWeights(
            proc,
            rows,
            .h2DesignCoef(coef, noise, tied)
        )
        previous <- coef
        coef <- .h2IrwlsStep(
            proc,
            rows,
            list(previous = previous, noise = noise, tied = tied),
            .h2WlsCoef(design$X, design$y, w)
        )
        if (.h2Converged(proc, previous, coef)) {
            break
        }
    }
    list(coef = coef, w = w)
}

# One weight update's coefficients. A fixed count of updates (the published
# procedures) takes each new fit whole. A fit to convergence takes it only
# where it lowers the regression's deviance, and otherwise halves the step
# toward it until it does (none lowers it: converged), as step-halving IRLS
# does: plain updates can cycle between two fits when a window's signal sits
# on a few large eigenvalues, and the halved steps reach the fixed point the
# likelihood shares. `at` holds the update's start and the noise settings.
# @noRd
.h2IrwlsStep <- function(proc, rows, at, candidate) {
    if (!is.null(proc@iterations) || identical(rows$weights, "cross")) {
        return(candidate)
    }
    start <- .h2Deviance(
        proc,
        rows,
        .h2DesignCoef(at$previous, at$noise, at$tied)
    )
    halving <- detect(
        0:30,
        .h2StepLowers,
        proc = proc,
        rows = rows,
        at = at,
        candidate = candidate,
        start = start
    )
    if (is.null(halving)) {
        return(at$previous)
    }
    .h2HalvedStep(at$previous, candidate, halving)
}

# The coefficients `halving` halvings of the step from `previous` toward
# `candidate`.
# @noRd
.h2HalvedStep <- function(previous, candidate, halving) {
    previous + (candidate - previous) / 2^halving
}

# Whether the step halved `halving` times lowers the deviance below `start`.
# @noRd
.h2StepLowers <- function(halving, proc, rows, at, candidate, start) {
    step <- .h2HalvedStep(at$previous, candidate, halving)
    .h2Deviance(proc, rows, .h2DesignCoef(step, at$noise, at$tied)) <
        start - 1e-12 * abs(start)
}

# The regression's deviance at (h2, noise), on the mean its weights read:
# each row's y / mu + log(mu), the scaled chi-squares' -2 log-likelihood up to
# a constant, times the taper and the LD-redundancy factor.
# @noRd
.h2Deviance <- function(proc, rows, coef) {
    mu <- .h2WeightMean(proc, rows, coef)
    sum(rows$taper * rows$redundancy * (rows$y / mu + log(mu)))
}

# The design and response with the noise held or tied: (X, y) as given; the
# genetic columns alone with y less the held noise's part; or, tied, the
# noise the held total (1 for a heritability's) less the genetic part: y less
# that total's part on the genetic columns less the noise column.
# @noRd
.h2RegressionDesign <- function(rows, noise, tied) {
    k <- ncol(rows$X)
    genetic <- rows$X[, -k, drop = FALSE]
    if (tied) {
        return(list(
            X = genetic - rows$X[, k],
            y = rows$y - (noise %||% 1) * rows$X[, k]
        ))
    }
    if (!is.null(noise)) {
        return(list(X = genetic, y = rows$y - noise * rows$X[, k]))
    }
    list(X = rows$X, y = rows$y)
}

# (genetic coefficients, noise) from a fit's coefficients.
# @noRd
.h2DesignCoef <- function(coef, noise, tied) {
    if (tied) {
        return(c(coef, (noise %||% 1) - sum(coef)))
    }
    if (!is.null(noise)) {
        return(c(coef, noise))
    }
    coef
}

# How many weight updates: none for fixed weights, the procedures' count, or
# up to 100 until convergence.
# @noRd
.h2WeightUpdates <- function(proc) {
    if (identical(proc@weights, "fixed")) {
        return(0L)
    }
    as.integer(proc@iterations %||% 100L)
}

# Whether a fit to convergence has stopped moving.
# @noRd
.h2Converged <- function(proc, previous, coef) {
    is.null(proc@iterations) &&
        max(abs(coef - previous) / pmax(abs(previous), 1e-8)) < 1e-8
}

# The first weights: flat, all equal, as LDER's; or, for cross-products, at
# no covariance; or the weights at the moments' heritability and no
# inflation, as S-LDSC's: its aggregate
# estimate, M (mean(chi2) - 1) / mean(N l) on the LD scores' total where the
# rows take ldsc's weights, the mean excess over the noise column on the
# unstratified column otherwise.
# @noRd
.h2FirstWeights <- function(param, rows) {
    proc <- param@procedures
    if (identical(proc@initialWeights, "flat")) {
        return(rep(1, length(rows$y)))
    }
    k <- ncol(rows$X)
    if (identical(rows$weights, "cross")) {
        return(rows$taper * rows$redundancy / (rows$mu11 * rows$mu22))
    }
    if (identical(rows$weights, "ldScore")) {
        hsq <- .h2AggregateH2(rows$y, rows$ldTotal, rows$n, rows$M)
        return(
            rows$taper *
                .h2LdScoreWeights(
                    rows$ldTotal,
                    rows$wLd,
                    rows$n,
                    rows$M,
                    hsq,
                    1
                )
        )
    }
    base <- rows$base %||% rows$X[, 1L]
    h2 <- min(max(mean(rows$y - rows$X[, k]) / mean(base), 0), 1)
    mu <- pmax(h2 * base + rows$X[, k], 1e-12)
    rows$taper * rows$redundancy / (2 * mu^2)
}

# The weights at (genetic coefficients, noise): a cross-product's
# (.h2CrossWeights()) where the rows are two conditions' cross-products;
# ldsc's (.h2LdScoreWeights()), from the first coefficient and the LD scores'
# total, where the rows take them; otherwise the inverse variance 1 / (2
# mu^2) for the fitted mean mu (.h2WeightMean()); each times the taper and
# the LD-redundancy factor.
# @noRd
.h2RegressionWeights <- function(proc, rows, coef) {
    if (identical(rows$weights, "cross")) {
        return(.h2CrossWeights(rows, coef))
    }
    if (!identical(rows$weights, "ldScore")) {
        return(
            rows$taper *
                rows$redundancy /
                (2 * .h2WeightMean(proc, rows, coef)^2)
        )
    }
    k <- ncol(rows$X)
    noise <- if (proc@fitWeightIntercept) max(coef[[k]], 0) else 1
    rows$taper *
        .h2LdScoreWeights(
            rows$ldTotal,
            rows$wLd,
            rows$n,
            rows$M,
            rows$hsqScale * coef[[1L]],
            noise
        )
}

# Cross-products' weights at (covariance, shared noise): the inverse of a
# product's variance mu11 mu22 + mu12^2, mu12 its fitted mean kept within the
# root of mu11 mu22, the two conditions' own means at their fits.
# @noRd
.h2CrossWeights <- function(rows, coef) {
    bound <- sqrt(rows$mu11 * rows$mu22)
    mu12 <- pmin(pmax(as.vector(rows$X %*% coef), -bound), bound)
    rows$taper * rows$redundancy / (rows$mu11 * rows$mu22 + mu12^2)
}

# The mean the weights read at (genetic coefficients, noise): each genetic
# coefficient held in [0, its cap] (1 unless the rows set one, h2's range),
# the noise at no less than 0, or at 1 where the procedures do not fit it,
# floored above 0.
# @noRd
.h2WeightMean <- function(proc, rows, coef) {
    k <- ncol(rows$X)
    caps <- rows$caps %||% rep(1, k - 1L)
    genetic <- pmin(pmax(coef[-k], 0), caps)
    noise <- if (proc@fitWeightIntercept) max(coef[[k]], 0) else 1
    fitted <- as.vector(rows$X[, -k, drop = FALSE] %*% genetic)
    pmax(fitted + noise * rows$X[, k], 1e-12)
}

# h2's standard error from the information: a sandwich on the final weights
# whose middle is the responses' covariance under the model at the fit
# (.h2RegressionMeat()). NA unless the settings ask.
# @noRd
.h2RegressionSe <- function(param, rows, fit) {
    if (!identical(param@inference, "information")) {
        return(NA_real_)
    }
    X <- fit$design$X
    bread <- solve(crossprod(X, X * fit$w))
    sqrt((bread %*% .h2RegressionMeat(rows, fit) %*% bread)[[1L, 1L]])
}

# The sandwich's middle, the weighted scores' covariance under the model at
# the fit: each row's variance 2 mu^2 where the rows are independent (an
# eigen regression's directions); where they are variants, which LD couples,
# Cov(y_j, y_k) = 2 S_jk^2 over every pair of them (.h2RowsModelS()).
# @noRd
.h2RegressionMeat <- function(rows, fit) {
    X <- fit$design$X
    if (is.null(rows$ldModel)) {
        mu <- pmax(as.vector(rows$X %*% fit$coef), 1e-12)
        return(crossprod(X, X * (fit$w^2 * 2 * mu^2)))
    }
    Xw <- X * fit$w
    crossprod(Xw, (2 * .h2RowsModelS(rows, fit$coef)^2) %*% Xw)
}

# A variant regression's statistics sqrt(n) r's covariance under the model at
# `coef` = (h2, noise), over its rows' variants: (n h2 / M) R W R + noise R,
# W the frequency weights, with h2 and the noise held at no less than 0 so
# that it is a covariance.
# @noRd
.h2RowsModelS <- function(rows, coef) {
    model <- rows$ldModel
    S <- .h2GlsS(
        model$ld,
        model$w %||% rep(1, model$M),
        model$n * max(coef[[1L]], 0) / model$M,
        max(coef[[2L]], 0)
    )
    S[rows$index, rows$index, drop = FALSE]
}

# h2's delete-one-block jackknife standard error over `jackknifeBlocks`
# contiguous blocks of the window's rows (.h2RegressionLooCoefs()).
# @noRd
.h2RegressionJackknife <- function(param, rows) {
    .jackknifeSeFromLoo(.h2RowsH2(.h2RegressionLooCoefs(param, rows), rows))
}

# The heritability of each row of coefficients `coef` (fits x columns): the
# genetic coefficients times each column's share of its size (1 unless the
# rows set one), summed.
# @noRd
.h2RowsH2 <- function(coef, rows) {
    genetic <- coef[, -ncol(coef), drop = FALSE]
    as.vector(genetic %*% (rows$share %||% rep(1, ncol(genetic))))
}

# A design's column sizes: each column's total absolute value (`M` without
# a design), on which a column's coefficient is a rate's scale.
# @noRd
.h2ColumnSizes <- function(modelMat, M) {
    if (is.null(modelMat)) {
        return(M)
    }
    size <- colSums(abs(modelMat))
    replace(size, size == 0, 1)
}

# Each column's sum over its size: 1 for a column that never goes below 0,
# what a unit of its coefficient adds to the heritability.
# @noRd
.h2ColumnShares <- function(modelMat, size) {
    if (is.null(modelMat)) 1 else colSums(modelMat) / size
}

# Each block's leave-one-out coefficients (blocks x columns): the rows' own
# blocks (`rows$block`, a reference's LD blocks), or `jackknifeBlocks`
# contiguous blocks of them.
# @noRd
.h2RegressionLooCoefs <- function(param, rows) {
    blockOf <- rows$block %||% .h2ContiguousBlocks(param, length(rows$y))
    exec(
        rbind,
        !!!map(
            unique(blockOf),
            .h2RegressionLoo,
            blockOf = blockOf,
            param = param,
            rows = rows
        )
    )
}

# `n` rows cut into the procedures' `jackknifeBlocks` contiguous blocks.
# @noRd
.h2ContiguousBlocks <- function(param, n) {
    nBlocks <- param@procedures@jackknifeBlocks
    if (is.null(nBlocks)) {
        abort(str_c(
            "In one window, a jackknife needs `jackknifeBlocks` ",
            "(H2RegressionParam()): there are no LD blocks to leave out."
        ))
    }
    cut(seq_len(n), min(nBlocks, n), labels = FALSE)
}

# The coefficients with block `b` left out.
# @noRd
.h2RegressionLoo <- function(b, blockOf, param, rows) {
    .h2Regression(param, .h2SubsetRows(rows, blockOf != b))$coef
}

# ---- the GLS -----------------------------------------------------------------

# The variant GLS (gLDSC's, Xiong et al. 2024, in one window): each variant's
# n r^2 - 1 on [1, l] by generalized least squares, with Omega = S o S for S
# = rawNtau R + R (no annotations), its rawNtau from a first pass of the
# same statistics in the window, inverted over the eigenpairs that hold the
# procedures' share of its mass (.h2GlsOmegaEigen()). With Var(z^2) = 2 S^2,
# the coefficients' covariance is 2 (X' Omega^-1 X)^-1.
# @noRd
.h2VariantGls <- function(param, z, nVariant, ld) {
    proc <- param@procedures
    n <- .h2SampleSize(param@sampleSize, nVariant)
    M <- length(z)
    l <- .h2WindowLdScores(param, ld, M)
    y <- n * .zToCor(z, nVariant)^2 - 1
    om <- .h2GlsOmegaEigen(
        ld,
        matrix(.h2SourceWeights(param, ld) %||% 1, M, 1L),
        sum(y) / sum(l),
        cutpoint = proc@omegaFraction
    )
    fit <- .h2GlsFit(y, cbind(1, l), om, proc@ridge)
    coef <- fit$coef
    covariance <- 2 * fit$inverse
    .cisH2Result(
        NA_integer_,
        h2 = coef[[2L]] * M / n,
        h2Se = if (identical(param@inference, "information")) {
            sqrt(covariance[[2L, 2L]]) * M / n
        } else {
            NA_real_
        },
        intercept = coef[[1L]] + 1
    )
}

# Candidate-annotation score statistics of the stratified eigen likelihood. Each
# eigen-direction's rotated effect has modelled variance v_i = lamh2_i, and the
# likelihood's score for a candidate's heritability is sum_i g_ic (bstar_i^2 /
# v_i^2 - 1 / v_i), with g_ic = dv_i/dh2_c = l_ic lam_i^2 / size_c - share_c
# lam_i / nPanel, for a column's total absolute value size_c and share_c = M_c /
# size_c, and information weight 2 / v_i^2. The baseline heritabilities and the
# intercept are the nuisance (.h2EfficientScoreBlocks()), with derivatives from
# the fit's own parameters (.h2AnnotationBasis()) and lam_i / n. Under frequency
# weights `w` a candidate's column, and so its size and share, is the weighted
# one.
.h2VarianceScoreStats <- function(
    blockData,
    nPanel,
    annotations,
    design,
    lamh2,
    n,
    w = NULL
) {
    candidate <- testedAnnotations(annotations)
    candMat <- assay(candidate, "annotations")
    if (ncol(candMat) == 0) {
        return(list(enrichment = NULL, scoreStats = NULL))
    }
    lam <- design$lam
    blockCands <- map(blockData, "ldCand")
    ldCand <- exec(rbind, !!!blockCands)
    weighted <- .h2WeightedDesign(candMat, w)
    candSize <- colSums(abs(weighted))
    gC <- .h2VarianceGradient(
        sweep(ldCand, 2L, candSize, "/"),
        lam,
        nPanel,
        colSums(weighted) / candSize
    )
    G <- cbind(
        .h2VarianceGradient(design$ldAnnotScaled, lam, nPanel, design$share),
        lam / n
    )
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

# The derivative of each eigen-direction's modelled variance in a column's
# parameter, from the column's LD scores scaled by its size: l_ia lam_i^2 /
# size_a - share_a lam_i / nPanel.
# @noRd
.h2VarianceGradient <- function(scaled, lam, nPanel, share) {
    scaled * lam^2 - tcrossprod(lam, share) / nPanel
}

# =============================================================================
# eigenLikelihoodUnivariate -- orchestrator
# =============================================================================

# The genome-wide eigen likelihood, stratified by `annotations` where given,
# under the settings `param` (the hdl preset's by default:
# .h2GenomeWideSettings()). Without the eigenvalue term the reference is taken
# to have no sampling noise.
eigenLikelihoodUnivariate <- function(
    z,
    n,
    eigenRef,
    annotations = NULL,
    param = .h2GenomeWideSettings("hdl", .h2Stratified(annotations), FALSE),
    weights = .h2RefWeights(param, eigenRef)
) {
    proc <- param@procedures
    M <- length(eigenRef)
    nPanel <- .h2PanelSize(param, eigenRef)
    baselineMat <- .h2DesignMat(annotations)
    # One pass over the blocks, each decomposed only while its data are taken.
    blockData <- map(
        .ldEigenBlocks(eigenRef),
        .h2EigenBlockData,
        eigenRef = eigenRef,
        z = z,
        n = n,
        baselineMat = .h2WeightedDesign(baselineMat, weights),
        candMat = .h2WeightedDesign(.h2CandidateMat(annotations), weights),
        w = weights
    )
    basis <- if (!is.null(baselineMat)) {
        .h2AnnotationBasis(baselineMat, weights, proc@h2Bounds)
    }
    design <- .h2EigenDesign(blockData, M, basis, proc)
    ft <- .h2VarianceFit(design, n, nPanel, proc)
    jk <- .h2VarianceJackknife(blockData, M, basis, n, nPanel, proc)
    .h2FitResult(
        h2 = ft$h2,
        h2Se = jk$h2Se,
        intercept = ft$int,
        interceptSe = jk$intSe,
        A = baselineMat,
        tau = as.vector(.h2AnnotationTau(matrix(ft$h2a, 1L), design)),
        tauBlocks = .h2AnnotationTau(jk$h2aBlocks, design),
        scoreStats = .h2CandidateScores(
            blockData,
            nPanel,
            annotations,
            design,
            ft,
            n,
            proc,
            weights
        ),
        w = weights
    )
}

# The reference size the eigenvalue term reads, or none (infinite) where the
# settings take the reference to have no sampling noise.
# @noRd
.h2PanelSize <- function(param, eigenRef) {
    if (identical(param@referenceNoise, "eigenTerm")) {
        ldPanelNSamples(eigenRef)
    } else {
        Inf
    }
}

# The score statistics for the candidate annotations, under the variances
# the baseline fit implies; NULL when there are none to test.
# @noRd
.h2CandidateScores <- function(
    blockData,
    nPanel,
    annotations,
    design,
    ft,
    n,
    proc,
    w = NULL
) {
    if (is.null(annotations)) {
        return(NULL)
    }
    lamh2 <- .h2BaselineVariance(design, ft, n, nPanel, proc@varianceFloor)
    .h2VarianceScoreStats(
        blockData,
        nPanel,
        annotations,
        design,
        lamh2,
        n,
        w
    )$scoreStats
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
# The variant regression: stratified LD score regression (Finucane et al. 2015)
# -----------------------------------------------------------------------------
# The original method the other three descend from. It regresses per-variant
# chi-square on annotation-stratified LD scores:
#
#   E[chi2_j] = N * sum_C tau_C * l(j, C) + N*a + 1
#   with l(j, C) = sum_{k in C} r2_jk
#
# Unlike the variant GLS it needs only the per-variant scores, not the per-block
# LD matrices, so an unstratified run works off a scores-only LdScore. The
# stratified run needs the matrices, since the per-annotation scores are formed
# from them exactly as the variant GLS forms its own.
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
#'   Its regression settings (\code{\link{H2RegressionParam}}) choose the
#'   route. S-LDSC's own take the two of \code{ldsc.py}, whose split is on
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
#'   The intercept is reported as estimated and is not floored at 1 unless
#'   the settings floor the noise: the model implies at least 1, but a value
#'   below it is informative rather than an error to clamp away, and
#'   clamping would move the heritability too.
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

# The genome-wide variant regression, stratified by `annotations` where given,
# under the settings `param` (the sldsc preset's by default:
# .h2GenomeWideSettings()): each variant's chi-square on its LD scores as rows
# of the regression every variant regression fits (.h2VariantRegressionRows(),
# .h2Regression()), with the delete-one-block jackknife.
# @noRd
variantRegressionUnivariate <- function(
    z,
    n,
    ldRef,
    annotations = NULL,
    param = .h2GenomeWideSettings("sldsc", .h2Stratified(annotations), FALSE),
    weights = .h2RefWeights(param, ldRef)
) {
    M <- length(ldRef)
    # The design the variant GLS uses (.h2DesignMat()), so `crossprod(A)` is the
    # annotation overlap matrix; the model's columns are its weighted ones.
    A <- .h2VariantAnnotMatrix(annotations, M)
    modelMat <- .h2WeightedDesign(A, weights)
    rows <- .h2VariantRegressionRows(
        z^2,
        .h2VariantScoreMatrix(ldRef, A, weights),
        modelMat,
        as.vector(mcols(ldRef)$ldScores[, 1L]),
        n,
        param,
        .h2VariantBlockIndex(ldRef)
    )
    fit <- .h2Regression(param, rows)
    .h2RowsResult(
        fit,
        .h2RegressionLooCoefs(param, rows),
        rows,
        if (is.null(.h2BaselineMat(annotations))) NULL else A,
        .h2VariantRegressionCandidateScores(fit, ldRef, annotations, weights),
        weights
    )
}

# The variant regression's rows: each variant's chi-square on its LD scores,
# column a scaled by n / M_a (.h2ColumnSizes()) so that its coefficient times
# its share is annotation a's heritability, and the noise column; with ldsc's
# weights (.h2LdScoreWeights()), which read the LD scores' total, the stored
# base score the weights divide by and the first coefficient's ldsc scale, and
# each variant's block for the jackknife.
# @noRd
.h2VariantRegressionRows <- function(
    chi2,
    scores,
    modelMat,
    baseScore,
    n,
    param,
    blockIdx
) {
    M <- length(chi2)
    size <- .h2ColumnSizes(modelMat, M)
    share <- .h2ColumnShares(modelMat, size)
    blockOf <- rep(seq_along(blockIdx), lengths(blockIdx))[
        order(list_c(blockIdx))
    ]
    list(
        y = chi2,
        X = cbind(sweep(n * scores, 2L, size, "/"), 1),
        chi2 = chi2,
        n = n,
        taper = rep(1, M),
        redundancy = 1 /
            pmax(.h2RedundancyScore(param@procedures, baseScore), 1),
        size = size,
        share = share,
        ldTotal = rowSums(scores),
        wLd = .h2RedundancyScore(param@procedures, baseScore),
        M = M,
        hsqScale = M / size[[1L]],
        index = seq_len(M),
        block = blockOf,
        weights = "ldScore"
    )
}

# The candidate annotations' score statistics, NULL without candidates: each
# candidate's stratified LD score (on its weighted column, under frequency
# weights `w`) is scored against the fitted residual in the fit's final
# weights, with the fitted design as nuisance (.h2EfficientScoreBlocks()).
# @noRd
.h2VariantRegressionCandidateScores <- function(
    fit,
    ldRef,
    annotations,
    w = NULL
) {
    candMat <- .h2CandidateMat(annotations)
    if (is.null(candMat)) {
        return(NULL)
    }
    if (!.ldScoreHasLd(ldRef)) {
        abort(str_c(
            "The variant regression's candidate-annotation score test needs ",
            "each LD block's LD to form the candidates' LD scores. Build the ",
            "reference from ",
            "a genotype panel, or with `buildLdScore(keepLdMatrices = TRUE)`."
        ))
    }
    kept <- fit$kept
    ldCand <- .h2StratifiedScores(
        ldRef,
        .h2WeightedDesign(candMat, w)
    )[kept$index, , drop = FALSE]
    regression <- fit$fit$w
    Ub <- .h2EfficientScoreBlocks(
        ldCand,
        kept$X,
        regression,
        regression * as.vector(kept$y - kept$X %*% fit$coef),
        kept$block
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
# (.h2VariantAnnotMatrix()): the stored base score for the base column, and a
# stratified score for each baseline annotation. Stratified scores come from
# each block's LD, the same source the variant GLS uses, and take ldsc's
# unbiased r^2 as the stored base score does: ldsc applies it to every
# annotation column. Under frequency weights `w` every column, the base's too,
# is the weighted one's score sum_k r^2_jk w_k A_ka, which the stored scores are
# not.
# @noRd
.h2VariantScoreMatrix <- function(ldRef, A, w = NULL) {
    if (!is.null(w)) {
        return(.h2VariantWeightedScores(ldRef, .h2WeightedDesign(A, w)))
    }
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
            "The stratified variant regression needs each LD block's LD to ",
            "form the ",
            "per-annotation LD scores. Build the reference from a genotype ",
            "panel, or with `buildLdScore(keepLdMatrices = TRUE)`."
        ))
    }
    stratified <- .h2StratifiedScores(ldRef, baselineMat)
    if (hasBase) cbind(base, stratified) else stratified
}

# The weighted design's LD scores, from each block's LD, which a scores-only
# LdScore does not carry.
# @noRd
.h2VariantWeightedScores <- function(ldRef, weighted) {
    if (!.ldScoreHasLd(ldRef)) {
        abort(str_c(
            "The variant regression with a frequency exponent (alpha) needs ",
            "each LD block's ",
            "LD to form the weighted LD scores. Build the reference from a ",
            "genotype panel, or with `buildLdScore(keepLdMatrices = TRUE)`."
        ))
    }
    .h2StratifiedScores(ldRef, weighted)
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
            nRef = .ldScorePanelN(ldRef)
        )
    )
}

# One block's per-annotation LD scores from its LD, summing ldsc's unbiased
# r^2 by its summed form.
# @noRd
.h2BlockAnnotScores <- function(block, panel, annotMat, nRef) {
    idx <- block$idx
    Ab <- .denseBlock(annotMat, idx)
    raw <- .ldSquaredProduct(.ldStatBlockLd(block, panel), Ab)
    list(idx = idx, values = .ldScoreUnbiasedSum(raw, colSums(Ab), nRef))
}

# Per-block variant indices: the LD blocks the estimators walk, which the
# jackknife leaves out one at a time.
# @noRd
.h2VariantBlockIndex <- function(ldRef) {
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
.h2LdScoreWeights <- function(ld, wLd, n, M, hsq, intercept) {
    hsq <- min(max(hsq, 0), 1)
    ld <- pmax(ld, 1)
    wLd <- pmax(wLd, 1)
    hetW <- 1 / (2 * (intercept + (hsq * n / M) * ld)^2)
    hetW / wLd
}

# LDSC's `aggregate`: the method-of-moments h2 that seeds the first weights.
#   h2 = M (mean(chi2) - intercept) / mean(N * ld)
# @noRd
.h2AggregateH2 <- function(chi2, baseScore, n, M, intercept = 1) {
    M * (mean(chi2) - intercept) / mean(n * baseScore)
}

# The LD score the weights divide by: the base score with the LD-redundancy
# correction, 1 (no division) without it.
# @noRd
.h2RedundancyScore <- function(proc, baseScore) {
    if (proc@ldRedundancy) baseScore else rep(1, length(baseScore))
}

# The weighted least-squares solve of every regression (.h2RidgeSolve()).
# @noRd
.h2WlsCoef <- function(design, y, w) {
    sw <- sqrt(pmax(w, 0))
    x <- design * sw
    as.vector(.h2RidgeSolve(crossprod(x), crossprod(x, y * sw), 1e-8))
}

# =============================================================================
# The genome-wide engine: one study's estimate under an H2EstimateParam
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

# A study's genome-wide estimate under `param`, stratified by `annotations`
# where given: the settings and reference checked, the inputs lined up, the
# fit run (.h2GenomeWideFit()) and its result recorded with the settings, on
# the liability scale too where the study records its prevalence
# (.h2LiabilityFactor()). The fit reads z alone, so the trait's own variance
# (varY), whatever its scale, does not enter it.
# @noRd
.h2GenomeWideEstimate <- function(
    data,
    ldRef,
    param,
    annotations,
    studyName,
    label
) {
    .h2RefuseGenomeWideSetting(param, label)
    .h2AssertSettingsRef(param, ldRef, label)
    study <- .estimateH2ResolveStudy(data, studyName)
    .h2ValidateInputs(data, study, ldRef, annotations, label)
    result <- .h2GenomeWideFit(
        param,
        z(data, studyName = study),
        .h2SampleSize(param@sampleSize, sampleSize(data, studyName = study)),
        ldRef,
        annotations
    )
    .h2EstimateFromResult(
        result,
        param,
        nSnps(data, studyName = study),
        study,
        ldRef,
        .h2LiabilityFactor(data, study)
    )
}

# What takes a study's heritability to the liability scale (Lee et al. 2011,
# equation 23): K^2 (1 - K)^2 / (P (1 - P) phi(t)^2) at its population
# prevalence K, its liability threshold t = qnorm(1 - K) and the share P of
# cases among the samples its N counts (.h2CaseShare()); NA for a study that
# records no prevalence, as a quantitative trait's does not.
# @noRd
#' @importFrom stats dnorm qnorm
.h2LiabilityFactor <- function(data, study) {
    K <- .h2StudyValue(data, "prevalence", study)
    if (is.na(K)) {
        return(NA_real_)
    }
    P <- .h2CaseShare(data, study)
    K^2 * (1 - K)^2 / (P * (1 - P) * dnorm(qnorm(1 - K))^2)
}

# The share of cases among the samples a study's N counts: one half where
# summaryStatsQc() set N to the effective sample size (qcInfo()'s nSource),
# which counts a balanced study's samples (Grotzinger et al. 2023);
# otherwise the share in its case/control counts, its variants' where it
# records them per variant, as that QC reads them first.
# @noRd
.h2CaseShare <- function(data, study) {
    if (identical(.h2StudyNSource(data, study), "effective")) {
        return(0.5)
    }
    md <- mcols(sumStats(data, studyName = study))
    cases <- as.numeric(md[["N_CASE"]])
    controls <- as.numeric(md[["N_CONTROL"]])
    counted <- is.finite(cases) & is.finite(controls)
    if (any(counted)) {
        return(sum(cases[counted]) / sum(cases[counted] + controls[counted]))
    }
    nCase <- .h2StudyValue(data, "nCase", study)
    nCase / (nCase + .h2StudyValue(data, "nControl", study))
}

# How summaryStatsQc() set a study's N (qcInfo()'s nSource), the one way its
# rows record; NA where it records none, as a collection the QC has not run
# on does not.
# @noRd
.h2StudyNSource <- function(data, study) {
    info <- qcInfo(data)
    if (!is_in("nSource", colnames(info))) {
        return(NA_character_)
    }
    source <- unique(discard(info[["nSource"]][info$study == study], is.na))
    if (length(source) > 1L) {
        abort(glue(
            "The rows of study '{study}' disagree on how the QC set N ",
            "({str_flatten(source, ', ')}; qcInfo()'s nSource), so no one ",
            "share of cases takes its heritability to the liability scale."
        ))
    }
    if (length(source) == 0L) NA_character_ else source
}

# A study's one value of a per-study column, NA where the collection does not
# record it; its rows disagreeing is an error, as there is no one value.
# @noRd
.h2StudyValue <- function(data, column, study) {
    if (!is_in(column, .tupleColumnNames(data))) {
        return(NA_real_)
    }
    rows <- as.character(data$study) == study
    value <- unique(as.numeric(.tupleColumn(data, column))[rows])
    if (length(value) > 1L) {
        abort(glue(
            "The rows of study '{study}' disagree on {column} ",
            "({str_flatten(value, ', ')})."
        ))
    }
    value
}

# The genome-wide fit each setting runs, over the reference's LD blocks with a
# delete-one-block jackknife: the variant regression (S-LDSC's), the variant
# GLS (g-LDSC's), the eigen regression (LDER's) and the eigen likelihood
# (HDL's, sHDL's stratified). Each reads a variant's correlation as z /
# sqrt(N), as the published implementations do; a local fit reads the exact
# z / sqrt(z^2 + n - 2) (.zToCor()), which a strong cis signal needs and
# which genome-wide statistics, each small, do not.
# @noRd
.h2GenomeWideFit <- function(param, z, n, ldRef, annotations) {
    fit <- switch(
        str_c(param@statistic, "/", param@criterion),
        "variant/regression" = variantRegressionUnivariate,
        "variant/gls" = variantGlsUnivariate,
        "eigen/regression" = eigenRegressionUnivariate,
        "eigen/likelihood" = eigenLikelihoodUnivariate
    )
    fit(z, n, ldRef, annotations, param = param)
}

# The settings a genome-wide estimate cannot be made under, refused with what
# it can be made under.
# @noRd
.h2RefuseGenomeWideSetting <- function(param, label) {
    problem <- .h2GenomeWideProblem(param)
    if (is.null(problem)) {
        return(invisible(NULL))
    }
    abort(glue("{label}: {problem}."))
}

# @noRd
.h2GenomeWideProblem <- function(param) {
    pair <- str_c(param@statistic, "/", param@criterion)
    fits <- c(
        "variant/regression",
        "variant/gls",
        "eigen/regression",
        "eigen/likelihood"
    )
    if (!is_in(pair, fits)) {
        return(str_c(
            "a genome-wide estimate fits the variant regression or GLS, or ",
            "the eigen regression or likelihood"
        ))
    }
    if (!is_in(param@noise, c("free", "floored"))) {
        return("a genome-wide estimate fits a free or floored noise term")
    }
    if (!identical(param@inference, "jackknife")) {
        return(str_c(
            "a genome-wide estimate reports the jackknife over the ",
            "reference's LD blocks"
        ))
    }
    .h2GenomeWideReferenceProblem(param, pair) %||%
        .h2GenomeWideProcedureProblem(param, pair)
}

# The reference terms each genome-wide fit reads: the LD score's unbiased r^2
# for the variant regression (an LdScore stores those), the raw r^2 for the
# variant GLS, the eigenvalue term or none for the likelihood, and the
# weights' or none for the eigen regression.
# @noRd
.h2GenomeWideReferenceProblem <- function(param, pair) {
    allowed <- switch(
        pair,
        "variant/regression" = "unbiasedR2",
        "variant/gls" = "none",
        "eigen/regression" = c("weightsOnly", "none"),
        "eigen/likelihood" = c("eigenTerm", "none")
    )
    if (is_in(param@referenceNoise, allowed)) {
        return(NULL)
    }
    as.character(glue(
        "the genome-wide {param@statistic} {param@criterion} reads ",
        "referenceNoise = {str_flatten(dQuote(allowed, FALSE), ' or ')}"
    ))
}

# What the genome-wide fits do not do: truncate or cap an LdEigen's
# directions or ask a block for a least number of them, cut their own
# jackknife blocks, or start the likelihood from a regression, test it
# against the null, or fit the null space or a scale.
# @noRd
.h2GenomeWideProcedureProblem <- function(param, pair) {
    proc <- param@procedures
    truncation <- param@truncation
    if (
        !is.null(truncation) &&
            (truncation@fraction < 1 ||
                !is.null(truncation@maxPerSample) ||
                truncation@minimum != 1)
    ) {
        return(str_c(
            "a genome-wide estimate keeps every direction its LdEigen holds: ",
            "use H2TruncationParam(fraction = 1, maxPerSample = NULL, ",
            "minimum = 1)"
        ))
    }
    if (is(proc, "H2RegressionParam") && !is.null(proc@jackknifeBlocks)) {
        return(str_c(
            "a genome-wide jackknife leaves out the reference's LD blocks: ",
            "leave `jackknifeBlocks` unset"
        ))
    }
    if (!is(proc, "H2LikelihoodParam")) {
        return(NULL)
    }
    unsupported <- c(
        mustBeatNull = proc@mustBeatNull,
        nullSpace = proc@nullSpace,
        profileScale = proc@profileScale
    )
    if (!any(unsupported)) {
        return(NULL)
    }
    as.character(glue(
        "the genome-wide likelihood does not take ",
        "{str_flatten(names(unsupported)[unsupported], ', ')}: set ",
        "{if (sum(unsupported) == 1L) 'it' else 'them'} to FALSE"
    ))
}

# The LD statistic each fit reads: an LdEigen for the eigen statistic, an
# LdScore for the variants.
# @noRd
.h2AssertSettingsRef <- function(param, ldRef, label) {
    want <- if (identical(param@statistic, "eigen")) "LdEigen" else "LdScore"
    if (is(ldRef, want)) {
        return(invisible(TRUE))
    }
    abort(glue(
        "{label}: the {param@statistic} statistic requires an {want} ",
        "object, got {class(ldRef)[[1L]]}."
    ))
}

# The estimators read `z` and the annotation rows positionally, by each LD
# block's own snpIdx -- nothing is matched on variant id along the way. A
# sumstats or annotation set that is not the reference's variant list, in the
# reference's order, therefore yields a silently wrong h2 rather than an
# error. This is the one place holding all three objects, so check here.
# @noRd
.h2ValidateInputs <- function(sumstats, studyName, ldRef, annotations, label) {
    checkGenomeBuild(sumstats, ldRef)
    .h2CheckBlockCount(ldRef, label)
    entry <- sumStats(sumstats, studyName = studyName)
    .h2CheckVariantCount(
        length(entry),
        length(ldRef),
        glue("study '{studyName}'"),
        label
    )
    .h2CheckPositions(entry, ldRef, glue("study '{studyName}'"), label)
    if (is.null(annotations)) {
        return(invisible(TRUE))
    }
    checkGenomeBuild(ldRef, annotations)
    .h2CheckVariantCount(
        nrow(annotations),
        length(ldRef),
        "`annotations`",
        label
    )
    .h2CheckPositions(rowRanges(annotations), ldRef, "`annotations`", label)
    invisible(TRUE)
}

# Every estimator takes its standard error from a delete-one-block jackknife,
# which needs at least two blocks. With one, LDER and g-LDSC fail deep inside
# a fold fit ("'x' must be a matrix") and HDL reports an SE of exactly 0, so
# say so here instead.
# @noRd
.h2CheckBlockCount <- function(ldRef, label) {
    nBlocks <- length(ldBlocks(ldRef))
    if (nBlocks >= 2L) {
        return(invisible(TRUE))
    }
    msg <- glue(
        "{label}: the LD reference carries {nBlocks} LD block(s). The ",
        "estimators take their standard error from a delete-one-block ",
        "jackknife, so at least two are required, and a usable SE needs ",
        "many more. Build the reference from per-block LD -- see ",
        "`buildLdEigen()` / `buildLdScore()`."
    )
    abort(msg)
}

# @noRd
.h2CheckVariantCount <- function(n, M, what, label) {
    if (n == M) {
        return(invisible(TRUE))
    }
    msg <- glue(
        "{label}: {what} carries {n} variant(s) but the LD reference ",
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
.h2CheckPositions <- function(gr, ldRef, what, label) {
    same <- as.character(seqnames(gr)) == as.character(seqnames(ldRef)) &
        start(gr) == start(ldRef)
    if (all(same)) {
        return(invisible(TRUE))
    }
    nBad <- sum(!same)
    firstBad <- which(!same)[[1L]]
    msg <- glue(
        "{label}: {what} and the LD reference disagree at ",
        "{nBad} of {length(same)} position(s), the first at index ",
        "{firstBad}. The estimators index by the reference's variant ",
        "order, so both must hold the same variants in the same order."
    )
    abort(msg)
}

# A univariate estimator's result as an H2Estimate of one genome-wide
# estimate, one row for the study, under the settings `param`, on the
# liability scale too by the factor `liability` (NA for none). The estimate
# covers the whole reference, which spans more chromosomes than an element
# may, so its element holds no ranges; the row takes the reference's genome.
# @noRd
.h2EstimateFromResult <- function(
    result,
    param,
    M,
    studyName,
    ldRef,
    liability
) {
    .h2EstimateNew(
        list(GenomicRanges::GRanges()),
        list(
            study = studyName,
            context = NA_character_,
            trait = NA_character_,
            h2 = unname(result$h2),
            h2_se = unname(result$h2Se),
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = liability * unname(result$h2),
            h2_liability_se = liability * unname(result$h2Se),
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = unname(result$intercept %||% NA_real_),
            intercept_se = unname(result$interceptSe %||% NA_real_),
            p_value = NA_real_,
            n_snps = as.integer(M),
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_
        ),
        list(
            enrichment = list(result$enrichment),
            jackknife = list(result$annotationJackknifeCoefs),
            score_stats = list(result$scoreStats)
        ),
        ldRef,
        param
    )
}

#' @rdname estimateH2
#' @export
setMethod(
    "estimateH2",
    signature(data = "GwasSumStats", ldRef = "LdStatistic"),
    function(data, ldRef, param = NULL, studyName = NULL) {
        .assertMethodParam(param %||% list(), "H2EstimateParam", "param")
        # An empty list carries no settings, as everywhere a Param is taken.
        # The default is the eigen likelihood (the hdl preset), as every
        # estimate's is.
        settings <- if (length(param) == 0L) {
            .h2GenomeWideSettings("hdl", FALSE, inSample(ldRef))
        } else {
            param
        }
        .h2GenomeWideEstimate(
            data,
            ldRef,
            settings,
            NULL,
            studyName,
            "estimateH2"
        )
    }
)

#' @rdname estimateStratifiedH2
#' @export
setMethod(
    "estimateStratifiedH2",
    signature(data = "GwasSumStats", ldRef = "LdStatistic"),
    function(data, ldRef, annotations, param = NULL, studyName = NULL) {
        assertClass(annotations, "AnnotationMatrix", .var.name = "annotations")
        .assertMethodParam(param %||% list(), "H2EstimateParam", "param")
        # An empty list carries no settings, as everywhere a Param is taken.
        # The default is the eigen likelihood (the hdl preset), stratified, as
        # every estimate's is.
        settings <- if (length(param) == 0L) {
            .h2GenomeWideSettings(
                "hdl",
                .h2Stratified(annotations),
                inSample(ldRef)
            )
        } else {
            param
        }
        .h2GenomeWideEstimate(
            data,
            ldRef,
            settings,
            annotations,
            studyName,
            "estimateStratifiedH2"
        )
    }
)

# =============================================================================
# estimateLocalH2 on QTL summary statistics -- one cis estimate per gene
# =============================================================================

#' @rdname estimateLocalH2
#' @param contexts Optional context selection, as
#'   \code{\link{twasWeightsPipeline}} takes it; \code{NULL} (the default)
#'   uses every context.
#' @param traitId Optional trait selection, as
#'   \code{\link{twasWeightsPipeline}} takes it; \code{NULL} (the default)
#'   uses every trait.
#' @param cisWindow Optionally, each selected trait is analysed over its
#'   \code{traitPos} \eqn{\pm} \code{cisWindow}, within its \code{region}.
#'   For a \code{QtlSumStats}, every selected row then needs a
#'   \code{traitPos}.
#' @param BPPARAM For a \code{QtlSumStats}: a \code{BiocParallelParam} over
#'   the selected rows, each row's estimate on one worker. Default
#'   \code{SerialParam()}.
#' @param referenceOverlap For a \code{QtlSumStats} (and a
#'   \code{MultiStudyQtlDataset}'s summary statistics): the share of the
#'   study's samples that are in the LD reference, from 0 to 1, one value or
#'   one per selected row. Default 0. Usually unknown, and with a reference of
#'   thousands it hardly matters (see Details): set it when the LD came from
#'   the study's own samples or a small cohort holding them, 1 if it holds
#'   them all. Only the overlap factor (\code{referenceNoise =
#'   "overlapFactor"}) reads it, and its estimates record it in
#'   \code{reference_overlap}.
#' @export
#' @importFrom BiocParallel SerialParam
setMethod(
    "estimateLocalH2",
    "QtlSumStats",
    function(
        data,
        param = NULL,
        contexts = NULL,
        traitId = NULL,
        cisWindow = NULL,
        BPPARAM = SerialParam(),
        referenceOverlap = 0
    ) {
        settings <- .h2LocalParam(param, individual = FALSE)
        .cisH2AssertQc(data)
        selection <- .qssSelection(
            data,
            contexts,
            traitId,
            cisWindow,
            "estimateLocalH2"
        )
        rows <- .cisH2SelectedRows(selection)
        fits <- .bpMap(
            rows,
            .cisH2Row,
            data = selection$data,
            param = settings,
            referenceOverlap = .cisH2RowOverlaps(
                referenceOverlap,
                rows,
                length(selection$data)
            ),
            BPPARAM = BPPARAM
        )
        .cisH2Estimate(fits, data, settings)
    }
)

# The settings a local estimate is made under: `param`, or its defaults,
# refused where a local estimate cannot be made under them.
# @noRd
.h2LocalParam <- function(param, individual) {
    .assertMethodParam(param %||% list(), "H2EstimateParam", "param")
    # An empty list carries no settings, as everywhere a Param is taken.
    settings <- if (length(param) == 0L) H2EstimateParam() else param
    .h2RefuseLocalSetting(settings, individual, "estimateLocalH2")
    settings
}

# The rows' fits as an H2Estimate on the data's genome, with one warning for
# each kind of gene that has no estimate or rests on an assumption.
# @noRd
.cisH2Estimate <- function(fits, data, settings) {
    out <- .h2EstimateNew(
        map(fits, "window"),
        as.list(list_rbind(map(fits, "values"))),
        .h2NoObjects(length(fits)),
        data,
        settings
    )
    .cisH2WarnEmpty(out)
    .cisH2WarnUndefined(out, map_lgl(fits, "undefined"))
    .cisH2WarnSmallReference(out, map_lgl(fits, "smallReference"))
    out
}

# A cis estimate needs QC'd QTL summary statistics; `label` names the caller.
# @noRd
.cisH2AssertQc <- function(data, label = "estimateLocalH2") {
    if (length(qcOptions(data)) > 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: the QtlSumStats has no QC record (qcOptions() ",
        "is empty). Call summaryStatsQc() first and pass the QC-applied ",
        "result."
    ))
}

# The rows the selection keeps, refusing an empty selection.
# @noRd
.cisH2SelectedRows <- function(selection, label = "estimateLocalH2") {
    rows <- which(.sumStatsRowsSelected(
        selection$data,
        selection$contexts,
        selection$traitId
    ))
    if (length(rows) == 0L) {
        abort(glue(
            "{label}: no rows matched the supplied contexts / traitId ",
            "filters."
        ))
    }
    rows
}

# Each row of the data's reference overlap, from one value or one per selected
# row.
# @noRd
.cisH2RowOverlaps <- function(
    referenceOverlap,
    rows,
    nRows,
    label = "estimateLocalH2"
) {
    assertNumeric(
        referenceOverlap,
        lower = 0,
        upper = 1,
        any.missing = FALSE,
        min.len = 1L
    )
    if (!is_in(length(referenceOverlap), c(1L, length(rows)))) {
        abort(glue(
            "{label}: `referenceOverlap` has {length(referenceOverlap)} ",
            "values for {length(rows)} selected row(s); give one, or one per ",
            "selected row."
        ))
    }
    replace(numeric(nRows), rows, referenceOverlap)
}

# One row's estimate, from its z and per-variant sample sizes against its LD
# reference's block for those variants: the window it covers, its values, and
# whether it is undefined or rested on assuming no overlap with a small
# reference.
# @noRd
.cisH2Row <- function(i, data, param, referenceOverlap) {
    st <- as.character(data$study[[i]])
    ctx <- as.character(data$context[[i]])
    tr <- as.character(data$trait[[i]])
    label <- glue(
        "estimateLocalH2: study='{st}', context='{ctx}', trait='{tr}'"
    )
    df <- .sumStatsDf(
        data,
        studyName = st,
        context = ctx,
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    fit <- .cisH2Fit(
        df,
        .rowLdRef(data, i),
        param,
        label,
        referenceOverlap[[i]]
    )
    list(
        window = .rowRegionAt(data, i),
        smallReference = fit$smallReference,
        undefined = fit$undefined,
        values = .cisH2Values(
            st,
            ctx,
            tr,
            param,
            fit,
            nrow(df),
            referenceOverlap[[i]]
        )
    )
}

# One estimate's H2Estimate values, from a per-trait fit (.cisH2Result()) on
# `nSnps` variants. Only the overlap factor reads the reference overlap, so
# only its estimates record it.
# @noRd
.cisH2Values <- function(st, ctx, tr, param, fit, nSnps, overlap) {
    factor <- identical(param@referenceNoise, "overlapFactor")
    tibble(
        study = st,
        context = ctx,
        trait = tr,
        h2 = fit$h2,
        h2_se = fit$h2Se,
        h2_lower = fit$h2Lower,
        h2_upper = fit$h2Upper,
        h2_liability = NA_real_,
        h2_liability_se = NA_real_,
        h2_liability_lower = NA_real_,
        h2_liability_upper = NA_real_,
        intercept = fit$intercept,
        intercept_se = NA_real_,
        p_value = fit$p,
        n_snps = as.integer(nSnps),
        n_components = fit$nComponents,
        scale = fit$scale,
        reference_overlap = if (factor) overlap else NA_real_
    )
}

# One row's fit under `param`; NA when its window holds no variants, as
# there is nothing to estimate from.
# @noRd
.cisH2Fit <- function(df, ldReference, param, label, referenceOverlap) {
    if (nrow(df) == 0L) {
        return(.cisH2Result(NA_integer_))
    }
    ld <- .cisH2LdSource(df$variant_id, ldReference, param, label)
    if (identical(param@referenceNoise, "overlapFactor")) {
        .cisH2AssertOverlap(mean(df$N), ld$nSamples, referenceOverlap, label)
    }
    .h2LocalFit(
        param,
        z = df$z,
        nVariant = df$N,
        ld = list_assign(ld, referenceOverlap = referenceOverlap)
    )
}

# A reference cannot hold more of the study's samples than it has samples.
# @noRd
.cisH2AssertOverlap <- function(n, nRef, overlap, label) {
    if (overlap * n <= nRef) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: `referenceOverlap` = {overlap} puts {overlap * n} of the ",
        "{n} study samples in an LD reference of {nRef}."
    ))
}

# The row's LD source (.rssLdSource()) for its variants, refusing a missing
# reference, and for settings that read the reference's sample size, a
# reference that does not record it.
# @noRd
.cisH2LdSource <- function(ids, ldReference, param, label) {
    if (is.null(ldReference)) {
        abort(glue("{label}: the row has no LD reference."))
    }
    ld <- .rssLdSource(list(variantNames = ids), ldReference, label)
    .h2AssertReferenceSize(param, ld$nSamples, label)
    .h2WithHeterozygosity(ld, ldReference, ids, param, label)
}

# An LD source with its variants' heterozygosities 2p(1 - p) in the
# reference (.h2Heterozygosity()) where the settings' frequency exponent
# reads them, refusing a reference that records no frequencies.
# @noRd
.h2WithHeterozygosity <- function(ld, ldReference, ids, param, label) {
    if (param@alpha == -1) {
        return(ld)
    }
    het <- .h2Heterozygosity(ldReference, ids, label)
    if (is.null(het)) {
        abort(c(
            glue(
                "{label}: alpha = {param@alpha} weighs each variant by its ",
                "allele frequency in the LD reference, which records none."
            ),
            i = "Give an LdMatrix with an AF column, or a genotype panel."
        ))
    }
    list_assign(ld, het = het)
}

# Each of `ids`' heterozygosity 2p(1 - p) in an LD reference: from the A1
# frequencies it records (an LdMatrix's AF, a sketch's, a mixture's), or from
# a panel's own dosages; NULL where it has neither. A swap of alleles leaves
# it as it is.
# @noRd
.h2Heterozygosity <- function(ldReference, ids, label) {
    matched <- .ldFromRefMatch(ldReference, ids, label, "error")
    af <- .ldRefAf(ldReference)
    p <- if (!is.null(af)) {
        af[matched$idx]
    } else if (.ldIsPanel(ldReference)) {
        colMeans(.ldRefDosage(ldReference, matched$idx, meanImpute = TRUE)) / 2
    }
    if (is.null(p)) NULL else 2 * p * (1 - p)
}

# The overlap factor, the eigenvalue term and the unbiased r^2 read the LD
# reference's sample size, which an LdMatrix records only when built with it.
# @noRd
.h2AssertReferenceSize <- function(param, nSamples, label) {
    if (!.h2ReadsReferenceSize(param) || !is.na(nSamples)) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "{label}: referenceNoise = \"{param@referenceNoise}\" reads the ",
            "LD reference's sample size."
        ),
        i = "Build the LdMatrix with `nSamples`."
    ))
}

# A gene whose window holds no variants has no estimate; one warning lists
# them.
# @noRd
.cisH2WarnEmpty <- function(out) {
    empty <- out$n_snps == 0L
    if (!any(empty)) {
        return(invisible(NULL))
    }
    warn(glue(
        "estimateLocalH2: {sum(empty)} of {length(empty)} gene(s) hold no ",
        "variants in their window; their h2 is NA: ",
        "{str_flatten(out$trait[empty], ', ')}."
    ))
}

# The likelihood with the null space is undefined for some genes (see
# estimateLocalH2()); one warning lists them.
# @noRd
.cisH2WarnUndefined <- function(out, undefined) {
    if (!any(undefined)) {
        return(invisible(NULL))
    }
    warn(glue(
        "estimateLocalH2: the likelihood with the null space is undefined ",
        "for {sum(undefined)} of {nrow(out)} gene(s), whose LD reference's ",
        "rank or statistics reach the sample size, as they can when the LD ",
        "is not the study's own; their h2 is NA: ",
        "{str_flatten(out$trait[undefined], ', ')}."
    ))
}

# Where the reference is no bigger than about the study, the overlap factor's
# estimate turns on the samples they share, and it assumed none: LD from the
# study's own samples comes out far too low. One warning lists the genes.
# @noRd
.cisH2WarnSmallReference <- function(out, small) {
    if (!any(small)) {
        return(invisible(NULL))
    }
    warn(c(
        glue(
            "estimateLocalH2: the LD reference of {sum(small)} of ",
            "{length(small)} gene(s) has at most twice the study's samples, ",
            "so the overlap factor's estimate depends on how many they share, ",
            "and `referenceOverlap = 0` assumes none: ",
            "{str_flatten(out$trait[small], ', ')}."
        ),
        i = str_c(
            "If the reference holds the study's samples, set ",
            "`referenceOverlap`: 1 if it holds them all."
        )
    ))
}

# =============================================================================
# estimateLocalH2 on a QtlDataset -- one in-sample estimate per trait
# =============================================================================

#' @rdname estimateLocalH2
#' @param region For a \code{QtlDataset} (each study of a
#'   \code{MultiStudyQtlDataset}): optionally, one region that every
#'   selected trait is analysed over, as \code{\link{twasWeightsPipeline}}
#'   takes it; otherwise each trait's own window (see \code{cisWindow}).
#' @param genotypeFilterParam,residualizationParam,phenotypeFilterParam For a
#'   \code{QtlDataset} (each study of a \code{MultiStudyQtlDataset}): the
#'   genotype filters, the covariate residualization
#'   and the phenotype filters, as \code{\link{twasWeightsPipeline}} takes
#'   them, so that each estimate is on the genotypes and phenotype that the
#'   trait's TWAS fits are on.
#' @export
setMethod(
    "estimateLocalH2",
    "QtlDataset",
    function(
        data,
        param = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        genotypeFilterParam = GenotypeFilterParam(),
        residualizationParam = ResidualizationParam(),
        phenotypeFilterParam = PhenotypeFilterParam()
    ) {
        settings <- .h2LocalParam(param, individual = TRUE)
        .h2AssertOneRegion(region, cisWindow)
        filter <- .qtlResolveFilter(genotypeFilterParam)
        qd <- .qtlRestrictSamples(data, filter$keepSamples)
        work <- .h2QtlDatasetWork(
            qd,
            settings,
            parseContexts(contexts, qd),
            parseTraitIds(traitId, qd),
            region,
            list(
                cisWindow = cisWindow,
                residualizationParam = residualizationParam,
                phenotypeFilter = .qtlResolvePhenotypeFilter(
                    phenotypeFilterParam,
                    naDefault = "drop",
                    naAllowed = c("drop", "impute"),
                    caller = "estimateLocalH2"
                )
            )
        )
        .h2IndividualEstimate(
            .qtlRunWork(
                work,
                qd,
                filter,
                residualizationParam,
                withAf = settings@alpha != -1
            ),
            settings,
            qd
        )
    }
)

# One region at most: a trait analysed over several would have several
# estimates, and an H2Estimate (or an RgEstimate, per pair) holds one per
# (study, context, trait).
# @noRd
.h2AssertOneRegion <- function(
    region,
    cisWindow,
    label = "estimateLocalH2(QtlDataset)"
) {
    .rejectRegionWithCisWindow(region, cisWindow, label)
    if (length(region) <= 1L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: `region` has {length(region)} ranges; ",
        "estimate one region at a time."
    ))
}

# The plan of a QtlDataset's estimates: one per selected (context, trait),
# enumerated and built as twasWeightsPipeline() builds its univariate fits,
# so that each is on the same window, samples, filters and residualized,
# scaled genotypes and phenotype.
# @noRd
.h2QtlDatasetWork <- function(qd, param, contexts, traitId, region, args) {
    xRegion <- .makeXRegions(region, FALSE)[[1L]]
    useCtx <- .twasQdsResolveContexts(qd, contexts)
    scope <- .twasQdsScope(
        studyName(qd),
        useCtx,
        .twasQdsResolveTraits(
            qd,
            useCtx,
            traitId,
            xRegion,
            "estimateLocalH2(QtlDataset)"
        ),
        if (.hasNestedTraitSelection(traitId)) traitId
    )
    regionArgs <- c(args, list(region = xRegion))
    units <- .enumUnivariateIndividual(qd, scope, regionArgs)
    .qtlWork(
        units,
        map(units, .unitSpec, data = qd, args = regionArgs),
        .h2UnitFit,
        list(param = param, selection = scope$traitSelection),
        .h2FlattenRows
    )
}

# One unit's estimates, one per (context, trait) group it builds on `block`.
# @noRd
.h2UnitFit <- function(unit, param, selection, block = NULL) {
    groups <- .unitGroups(unit, block) |>
        map(.jointGroupSelect, traitId = selection) |>
        compact()
    map(groups, .h2GroupRow, param = param, block = block)
}

# A group's genotypes' heterozygosities 2p(1 - p), from its block's allele
# frequencies; NULL where the block carries none.
# @noRd
.h2GroupHet <- function(g, block) {
    if (is.null(block$af)) {
        return(NULL)
    }
    af <- block$af[match(colnames(.jgX(g)), colnames(block$X))]
    2 * af * (1 - af)
}

# One group's estimate: its window and its H2Estimate values.
# @noRd
.h2GroupRow <- function(g, param, block) {
    cond <- .jgConditions(g)
    fit <- .h2IndividualFitParam(
        param,
        .jgX(g),
        .jgY(g)[, 1L],
        .h2GroupHet(g, block)
    )
    list(
        window = block$window,
        undefined = fit$undefined,
        values = .cisH2Values(
            cond$study[[1L]],
            cond$context[[1L]],
            cond$trait[[1L]],
            param,
            fit,
            fit$nSnps,
            1
        )
    )
}

# The estimates of a plan's units as one list, the units with none left out.
# @noRd
.h2FlattenRows <- function(results) {
    list_flatten(compact(results))
}

# A trait's estimate under `param` from its genotypes `X` (samples x
# variants) and phenotype `y`, on the variants that vary, with how many there
# are. The summary-statistic settings read each variant's marginal t
# statistic against the trait, with the genotypes as the LD and the reference
# the study's own samples: the moments' estimate is then the adjusted R^2 of
# the trait on the genotypes' components. `het` holds each column's
# heterozygosity 2p(1 - p) where the settings' frequency exponent reads it.
# @noRd
#' @importFrom matrixStats colSds
#' @importFrom stats cor
.h2IndividualFitParam <- function(param, X, y, het = NULL) {
    sds <- colSds(X)
    varies <- !is.na(sds) & sds > 0
    Xs <- scale(X[, varies, drop = FALSE])
    ld <- list(
        ldFactor = Xs,
        nSamples = nrow(Xs),
        referenceOverlap = 1,
        het = het[varies]
    )
    fit <- if (ncol(Xs) == 0L) {
        .cisH2Result(NA_integer_)
    } else if (identical(param@criterion, "reml")) {
        .h2Reml(.h2WeightedColumns(Xs, param, ld), y, param@procedures)
    } else {
        n <- nrow(Xs)
        .h2LocalFit(
            param,
            z = .corToZ(as.vector(cor(Xs, y)), n),
            nVariant = rep(n, ncol(Xs)),
            ld = ld
        )
    }
    c(fit, list(nSnps = ncol(Xs)))
}

# The standardized genotypes with each column times the root of its
# frequency weight, so that their GRM is X W X' / M (as they are without
# weights).
# @noRd
.h2WeightedColumns <- function(Xs, param, ld) {
    w <- .h2SourceWeights(param, ld)
    if (is.null(w)) Xs else sweep(Xs, 2L, sqrt(w), `*`)
}

# The t statistic of a marginal regression whose correlation is `r`, on `n`
# samples: the inverse of .zToCor().
# @noRd
.corToZ <- function(r, n) {
    r * sqrt((n - 2) / (1 - r^2))
}

# REML: the trait's genetic and residual variances by qgg's REML on the
# relationship matrix of the standardized genotypes `Xs`, with the further
# qgg::greml() arguments `options` (GremlOptions()); h2 the genetic share,
# its standard error the delta method's from the variances' inverse average
# information.
# @noRd
.h2Reml <- function(Xs, y, options) {
    if (!requireNamespace("qgg", quietly = TRUE)) {
        abort("Package 'qgg' is required for the reml criterion.")
    }
    fit <- exec(
        qgg::greml,
        y = as.numeric(y),
        GRM = list(tcrossprod(Xs) / ncol(Xs)),
        !!!as.list(options)
    )
    theta <- as.numeric(fit$theta)
    total <- sum(theta)
    gradient <- c(theta[[2L]], -theta[[1L]]) / total^2
    .cisH2Result(
        NA_integer_,
        h2 = theta[[1L]] / total,
        h2Se = sqrt(drop(crossprod(gradient, fit$asd %*% gradient)))
    )
}

# A QtlDataset's estimates as an H2Estimate, on the genotypes' genome.
# @noRd
.h2IndividualEstimate <- function(rows, param, qd) {
    if (length(rows) == 0L) {
        abort(str_c(
            "estimateLocalH2(QtlDataset): no selected trait shares two or ",
            "more samples with its genotypes."
        ))
    }
    out <- .h2EstimateNew(
        map(rows, "window"),
        as.list(list_rbind(map(rows, "values"))),
        .h2NoObjects(length(rows)),
        .qtlGenotypeSe(qd),
        param
    )
    .cisH2WarnUndefined(out, map_lgl(rows, "undefined"))
    out
}

# =============================================================================
# estimateLocalH2 on a MultiStudyQtlDataset -- each study as its own kind
# =============================================================================

#' @rdname estimateLocalH2
#' @export
setMethod(
    "estimateLocalH2",
    "MultiStudyQtlDataset",
    function(
        data,
        param = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        genotypeFilterParam = GenotypeFilterParam(),
        residualizationParam = ResidualizationParam(),
        phenotypeFilterParam = PhenotypeFilterParam(),
        referenceOverlap = 0
    ) {
        label <- "estimateLocalH2(MultiStudyQtlDataset)"
        selection <- list(contexts = contexts, traitId = traitId)
        summary <- .h2MsSummarySelection(sumStats(data), selection)
        settings <- .h2LocalParam(param, individual = is.null(summary))
        .h2AssertOneRegion(region, cisWindow, label)
        studies <- map(
            as.list(data),
            .h2MsStudyEstimate,
            settings = settings,
            selection = selection,
            args = list(
                region = region,
                cisWindow = cisWindow,
                genotypeFilterParam = genotypeFilterParam,
                residualizationParam = residualizationParam,
                phenotypeFilterParam = phenotypeFilterParam
            )
        )
        parts <- compact(c(
            unname(studies),
            list(.h2MsSummaryEstimate(
                sumStats(data),
                summary,
                settings,
                cisWindow,
                referenceOverlap
            ))
        ))
        if (length(parts) == 0L) {
            abort(glue("{label}: no study holds a selected context and trait."))
        }
        exec(c, !!!parts)
    }
)

# The contexts and traits of a selection that a component holds, a selection
# not given by name left to it; NULL where it holds none of them.
# @noRd
.h2MsHeld <- function(selection, contexts, traits) {
    held <- list(
        contexts = .rgMsHeld(selection$contexts, contexts),
        traitId = .rgMsHeld(selection$traitId, traits)
    )
    none <- identical(held$contexts, character(0)) ||
        identical(held$traitId, character(0))
    if (none) NULL else held
}

# One study's estimates, as a QtlDataset's are, under its share of the
# selection; NULL where it holds none of it.
# @noRd
.h2MsStudyEstimate <- function(qd, settings, selection, args) {
    held <- .h2MsHeld(selection, contexts(qd), .spListTraits(qd))
    if (is.null(held)) {
        return(NULL)
    }
    exec(
        estimateLocalH2,
        qd,
        settings,
        contexts = held$contexts,
        traitId = held$traitId,
        !!!args
    )
}

# The selection's share of the summary statistics' rows; NULL where there
# are none or it holds none of them.
# @noRd
.h2MsSummarySelection <- function(ss, selection) {
    if (is.null(ss)) {
        return(NULL)
    }
    .h2MsHeld(
        selection,
        unique(as.character(ss$context)),
        unique(as.character(ss$trait))
    )
}

# The summary statistics' estimates, as a QtlSumStats's are, under their
# share of the selection; NULL without one.
# @noRd
.h2MsSummaryEstimate <- function(ss, held, settings, cisWindow, overlap) {
    if (is.null(held)) {
        return(NULL)
    }
    estimateLocalH2(
        ss,
        settings,
        contexts = held$contexts,
        traitId = held$traitId,
        cisWindow = cisWindow,
        referenceOverlap = overlap
    )
}

# =============================================================================
# A trait's heritability for the methods that start from one
# =============================================================================
#
# The Bayesian alphabet, SuSiE-RSS and LDpred2 can start from a trait's cis
# heritability: a value (`h2Estimate`), or an estimate under the settings
# `h2EstimateParam` on the method's own data. A weight function takes either;
# a pipeline takes a value for every trait, an H2Estimate to read each
# trait's from, or the settings to estimate each trait's under.

# The heritability a method starts from: NULL for neither, which leaves each
# engine its own default; `h2Estimate`, a number or a one-row H2Estimate; or
# the estimate under `h2EstimateParam`, made by `fit(settings, !!!fitArgs)`,
# which only runs then (`individual` where the method's data are a
# QtlDataset's genotypes and phenotype). Settled by .h2Settle().
# @noRd
.h2Resolve <- function(
    h2Estimate,
    h2EstimateParam,
    label,
    fit,
    fitArgs,
    individual = FALSE
) {
    .h2AssertOneSource(h2Estimate, h2EstimateParam, label)
    if (!is.null(h2Estimate)) {
        return(.h2Settle(.h2SuppliedValue(h2Estimate, label), label))
    }
    if (is.null(h2EstimateParam)) {
        return(NULL)
    }
    settings <- .h2EstimateSettings(h2EstimateParam, individual, label)
    .h2Settle(exec(fit, settings, !!!fitArgs)$h2, label)
}

# A heritability is supplied or estimated, not both.
# @noRd
.h2AssertOneSource <- function(h2Estimate, h2EstimateParam, label) {
    if (is.null(h2Estimate) || is.null(h2EstimateParam)) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: give the heritability (`h2Estimate`) or the settings to ",
        "estimate it under (`h2EstimateParam`), not both."
    ))
}

# The settings a heritability is estimated under, refused where the method's
# data cannot be fitted under them.
# @noRd
.h2EstimateSettings <- function(param, individual, label) {
    .assertMethodParam(param, "H2EstimateParam", "h2EstimateParam")
    .h2RefuseLocalSetting(param, individual, label)
    param
}

# A supplied heritability: a number, or the estimate of a one-row
# H2Estimate.
# @noRd
.h2SuppliedValue <- function(h2Estimate, label) {
    if (!is(h2Estimate, "H2Estimate")) {
        return(.h2AssertValue(h2Estimate, label))
    }
    if (length(h2Estimate) != 1L) {
        abort(glue(
            "{label}: an `h2Estimate` H2Estimate gives one heritability, so ",
            "it holds one row, not {length(h2Estimate)}."
        ))
    }
    as.numeric(h2Estimate$h2)
}

# A supplied heritability: one number in [0, 1).
# @noRd
#' @importFrom checkmate assertNumber
.h2AssertValue <- function(h2, label) {
    assertNumber(h2, lower = 0, upper = 1, .var.name = "h2Estimate")
    if (h2 >= 1) {
        abort(glue("{label}: `h2Estimate` must be below 1, not {h2}."))
    }
    h2
}

# A heritability as the methods use it: floored at 0.001 and capped at 0.99,
# which keeps every engine's starting variances and SuSiE-RSS's residual
# variance positive. An estimate that could not be made (NA) gives NULL, so
# the method keeps its own default, with a warning naming `label`.
# @noRd
.h2Settle <- function(value, label) {
    if (is.na(value)) {
        warn(glue(
            "{label}: no heritability estimate could be made, so the method ",
            "keeps its own default."
        ))
        return(NULL)
    }
    min(max(value, 0.001), 0.99)
}

# A summary-statistic heritability under `param` on a method's own z,
# per-variant sample sizes and LD source (.rssLdSource()), the reference
# sharing no samples with the study, as estimateLocalH2() assumes by default.
# Where that assumption decides the overlap factor's estimate, as
# estimateLocalH2() warns, so does this. A source without heterozygosities
# takes its variants' (`ids`) from `ldReference` where the frequency
# exponent reads them (.h2WithHeterozygosity()).
# @noRd
.h2SumStatsFit <- function(
    param,
    z,
    n,
    ld,
    label,
    ldReference = NULL,
    ids = NULL
) {
    .h2AssertReferenceSize(param, ld$nSamples, label)
    weighted <- if (is.null(ld$het) && !is.null(ldReference)) {
        .h2WithHeterozygosity(ld, ldReference, ids, param, label)
    } else {
        ld
    }
    fit <- .h2LocalFit(
        param,
        z = z,
        nVariant = n,
        ld = list_assign(weighted, referenceOverlap = 0)
    )
    if (isTRUE(fit$smallReference)) {
        warn(c(
            glue(
                "{label}: the LD reference has at most twice the study's ",
                "samples, so the overlap factor's estimate depends on how ",
                "many they share, and this assumes none."
            ),
            i = str_c(
                "If the reference holds the study's samples, pass ",
                "estimateLocalH2(..., referenceOverlap = ) as `h2Estimate` ",
                "instead."
            )
        ))
    }
    fit
}

# .h2SumStatsFit() on a summary-statistic weight function's own `stat` and
# LD reference, whose LD source is formed only here, when an estimate is
# asked for.
# @noRd
.h2StatFit <- function(param, stat, ldReference, label) {
    ld <- .h2WithHeterozygosity(
        .rssLdSource(stat, ldReference, label),
        ldReference,
        .rssStatVariantIds(stat, label),
        param,
        label
    )
    p <- .ldSourceSize(ld)
    n <- median(stat$n)
    .h2SumStatsFit(
        param,
        .corToZ(.rssMarginalCor(stat, n = n, p = p), n),
        rep_len(stat$n, p),
        ld,
        label
    )
}

# A pipeline's heritability, checked once at entry: NULL for none; one number
# for every trait; an H2Estimate, which must hold exactly one row for each of
# the (study, context, trait) `tuples` the run selected; or the settings to
# estimate each trait's under (`individual` for a QtlDataset's traits).
# @noRd
.h2PipelineSpec <- function(
    h2Estimate,
    h2EstimateParam,
    tuples,
    label,
    individual = FALSE
) {
    .h2AssertOneSource(h2Estimate, h2EstimateParam, label)
    if (!is.null(h2EstimateParam)) {
        settings <- .h2EstimateSettings(h2EstimateParam, individual, label)
        return(list(param = settings))
    }
    if (is.null(h2Estimate)) {
        return(NULL)
    }
    if (is(h2Estimate, "H2Estimate")) {
        return(list(values = .h2EstimateValues(h2Estimate, tuples, label)))
    }
    list(value = .h2AssertValue(h2Estimate, label))
}

# The (study, context, trait) of a summary-statistics collection's `rows`.
# @noRd
.h2RowTuples <- function(data, rows) {
    tibble(
        study = as.character(data$study[rows]),
        context = as.character(data$context[rows]),
        trait = as.character(data$trait[rows])
    )
}

# Row `i`'s LD reference records its sample size, where the settings `param`
# read it.
# @noRd
.h2AssertRowReference <- function(i, data, param, label) {
    ref <- .rowLdRef(data, i)
    if (is.null(param) || is.null(ref)) {
        return(invisible(NULL))
    }
    .h2AssertReferenceSize(param, .ldRefNoiseN(ref), label)
}

# The heritability of each of `tuples` in an H2Estimate, keyed by
# .h2TupleKey(), refusing a tuple it has no row for. An H2Estimate holds one
# row per tuple.
# @noRd
.h2EstimateValues <- function(est, tuples, label) {
    have <- .h2TupleKey(est$study, est$context, est$trait)
    want <- unique(.h2TupleKey(tuples$study, tuples$context, tuples$trait))
    .h2AssertTuples(setdiff(want, have), "holds no row for", label)
    set_names(as.numeric(est$h2), have)[want]
}

# One (study, context, trait) as a lookup key.
# @noRd
.h2TupleKey <- function(study, context, trait) {
    str_c(study, context, trait, sep = "\r")
}

# Refuse the tuples `keys`, which an H2Estimate `problem`, naming them.
# @noRd
.h2AssertTuples <- function(keys, problem, label) {
    if (length(keys) == 0L) {
        return(invisible(NULL))
    }
    shown <- str_flatten(
        str_replace_all(utils::head(keys, 5L), "\r", "/"),
        ", "
    )
    abort(glue(
        "{label}: the `h2Estimate` H2Estimate {problem} {length(keys)} ",
        "selected (study, context, trait): {shown}."
    ))
}

# One trait's heritability under a pipeline's settled `spec`: none, the run's
# one value, its H2Estimate row's (`key`), or its estimate under the settings
# by `fit(settings, !!!fitArgs)`. A condition the run derives itself, such as
# a top principal component, has no row and keeps its default. The spec is
# read by `[[`, which matches exactly: `$value` would take `values`.
# @noRd
.h2ForTrait <- function(spec, key, label, fit, fitArgs) {
    if (is.null(spec)) {
        return(NULL)
    }
    if (!is.null(spec[["value"]])) {
        return(.h2Settle(spec[["value"]], label))
    }
    values <- spec[["values"]]
    if (!is.null(values)) {
        if (!is_in(key, names(values))) {
            return(NULL)
        }
        return(.h2Settle(values[[key]], label))
    }
    .h2Settle(exec(fit, spec[["param"]], !!!fitArgs)$h2, label)
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
        values = .ldScoreUnbiasedSum(
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
    Ab <- .denseBlock(annotMat, idx)
    raw <- .ldSquaredProduct(dec[c("vectors", "values")], cbind(1, Ab))
    list(
        idx = idx,
        values = .ldScoreUnbiasedSum(raw, c(length(idx), colSums(Ab)), nRef)
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
        nRef <- .ldScorePanelN(ldRef)
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

# H2Estimate accessor methods (stratifiedHeritabilityEnrichment /
# scoreStats) live in R/H2Estimate.R alongside the class definition.

# =============================================================================
# Converter: H2Estimate -> sldsc_wrapper list format
# =============================================================================

#' @importFrom checkmate assertClass
#' @title Convert H2Estimate to S-LDSC Trait Format
#' @description Convert an \code{H2Estimate} object into the list format
#'   expected by \code{\link{standardizeSldscTrait}} and
#'   \code{\link{metaSldscRandom}}. This bridges the h2 estimation methods
#'   (LDER, gLDSC, HDL) into the sldscWrapper.R postprocessing pipeline.
#' @param h2Est An \code{H2Estimate} holding one stratified estimate, with
#'   its enrichment table and \code{annotationJackknifeCoefs}. Select one
#'   from several by its identity columns, e.g.
#'   \code{x[x$study == "trait1"]}.
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
    .h2AssertOneEstimate(h2Est)
    md <- mcols(h2Est, use.names = FALSE)
    enrichDf <- md$enrichment[[1L]]
    if (is.null(enrichDf)) {
        msg <- glue(
            "H2Estimate has no enrichment results. ",
            "Run estimateStratifiedH2() to get enrichment estimates."
        )
        abort(msg)
    }

    cats <- as.character(enrichDf$annotation)

    rawTau <- md$jackknife[[1L]]
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
        h2g = md$h2[[1L]],
        annotationJackknifeCoefs = jkCoefs,
        nBlocks = nBlocks
    )
}

# An error unless the H2Estimate holds one estimate.
# @noRd
.h2AssertOneEstimate <- function(h2Est) {
    if (length(h2Est) == 1L) {
        return(invisible(h2Est))
    }
    abort(c(
        glue(
            "h2EstimateToSldscTrait: the H2Estimate holds {length(h2Est)} ",
            "estimates; it converts one."
        ),
        i = "Select one by its identity, e.g. `x[x$study == \"trait1\"]`."
    ))
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One delete-one-block jackknife fold of the variant GLS: re-solve the normal
# equations with block `cb`'s Omega-weighted contribution removed.
# @noRd
.h2GlsLooEstimate <- function(cb, left, right, A, N, ridge) {
    .h2GlsEstimate(left - cb$L, right - cb$R, A, N, ridge)
}

# One delete-one-block jackknife fold of the genome-wide eigen likelihood: refit
# leaving block `b` out.
# @noRd
.h2VarianceLooFit <- function(b, blockData, M, basis, n, nSamples, proc) {
    # The full-data fit already reported any bound it hit; repeating that
    # once per left-out block would bury it in noise.
    .h2VarianceFit(
        .h2EigenDesign(blockData[-b], M, basis, proc),
        n,
        nSamples,
        proc,
        warnOnBound = FALSE
    )
}
