# Fine-mapping rows and fits shared by the fine-mapping, collection, TWAS and
# result-writer tests. The fitted fixtures are promises: each is built once,
# on first use, so a test run that never reads one never pays for it.

.makeAdjustEntry <- function(vids, maxNumSingleEffects = 2L) {
    p <- length(vids)
    set.seed(11L)
    lbf <- matrix(
        rnorm(maxNumSingleEffects * p),
        nrow = maxNumSingleEffects,
        ncol = p
    )
    colnames(lbf) <- vids
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    FineMappingRow(
        variantIds = vids,
        susieFit = list(
            pip = pip,
            alpha = alpha,
            lbf_variable = lbf,
            mu = matrix(0, maxNumSingleEffects, p),
            X_column_scale_factors = rep(1, p)
        ),
        topLoci = data.frame(
            variant_id = vids,
            pip = pip,
            betahat = rep(0, p),
            sebetahat = rep(1, p),
            stringsAsFactors = FALSE
        )
    )
}

# n = 100, p = 12 and L = 3 still give both credible sets (variants 3 and
# 10), bands and affected regions of either sign, at a fraction of the
# fitting time of a larger simulation.
.fsa_makeFit <- function(post_processing = "none") {
    set.seed(1)
    n <- 100L
    p <- 12L
    J <- 16L
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            paste0("chr1:", seq_len(p) * 100, ":A:G")
        )
    )
    b1 <- sin(seq(0, 2 * pi, length.out = J))
    b2 <- cos(seq(0, pi, length.out = J))
    Y <- X[, 3] %o% b1 + X[, 10] %o% b2 + matrix(rnorm(n * J, sd = 0.3), n, J)
    colnames(Y) <- paste0("f", seq_len(J))
    suppressWarnings(fsusieR::susiF(
        X = X,
        Y = Y,
        pos = seq_len(J),
        L = 3,
        post_processing = post_processing,
        verbose = FALSE
    ))
}

.fsa_entry <- function(fit, purityVal = 0.85) {
    vids <- names(fit$csd_X)
    pipv <- if (!is.null(fit$pip)) {
        as.numeric(fit$pip)
    } else {
        rep(0.1, length(vids))
    }
    # Canonical CS labels + purity from the fit's CS membership, mirroring what
    # The fsusie post-processing writes into a real FMR's topLoci
    # (cs_95 / cs_95_purity).
    cs95 <- rep("fsusie_0", length(vids))
    for (l in seq_along(fit$cs)) {
        idx <- fit$cs[[l]]
        if (is.numeric(idx)) {
            idx <- as.integer(idx)
        }
        cs95[idx] <- paste0("fsusie_", l)
    }
    tl <- data.frame(
        variant_id = vids,
        pip = pipv,
        cs_95 = cs95,
        cs_95_purity = ifelse(cs95 != "fsusie_0", purityVal, 0),
        stringsAsFactors = FALSE
    )
    FineMappingRow(variantIds = vids, susieFit = fit, topLoci = tl)
}

# Each fit is made once and shared: susiF() is the slow step, it is
# deterministic, and the tests only read the fits. NULL without fsusieR.
delayedAssign(
    ".fsa_fits",
    if (
        requireNamespace("fsusieR", quietly = TRUE) &&
            requireNamespace("wavethresh", quietly = TRUE)
    ) {
        list(none = .fsa_makeFit("none"), TI = .fsa_makeFit("TI"))
    }
)

.twcv_foldFits <- function(seed = 11, fold = 3) {
    set.seed(seed)
    # Loaded here, not into the global environment: this runs as a helper's
    # promise, whose environment does not see what data() puts there.
    data("eqtlRegionExample", package = "pecotmr", envir = environment())
    X <- eqtlRegionExample$X[, 1:40]
    y <- eqtlRegionExample$yRes
    Y <- matrix(y, ncol = 1, dimnames = list(rownames(X), "t1"))
    cv <- pecotmr:::.fmWeightsCv(
        X,
        Y,
        tokens = "susie",
        methodArgs = list(),
        fold = fold,
        verbose = 0,
        seed = 1
    )
    list(X = X, Y = Y, cv = cv, slice = pecotmr:::.fmSliceCv(cv, "susie"))
}

# Built once: each build runs a susie fit per fold, and the tests only read
# it. susieR is an Import, so the build needs no skip.
delayedAssign(".twcvFixture", suppressMessages(.twcv_foldFits()))

.thetaFit <- function(nV = 12, theta = NULL) {
    set.seed(1)
    alpha <- matrix(runif(2 * nV), nrow = 2)
    alpha <- alpha / rowSums(alpha)
    mu <- matrix(rnorm(2 * nV), nrow = 2)
    fit <- list(
        alpha = alpha,
        mu = mu,
        mu2 = mu^2 + 0.1,
        pip = runif(nV),
        X_column_scale_factors = rep(1, nV)
    )
    if (!is.null(theta)) {
        fit$theta <- theta
    }
    fit
}
