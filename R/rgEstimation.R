# =============================================================================
# Local genetic correlation: estimateLocalRg() and its fits
# -----------------------------------------------------------------------------
# One window's genetic covariance and correlation between two conditions, under
# the settings an H2EstimateParam gives: the bivariate forms of the per-trait
# estimators estimateLocalH2() fits. Within a window the two conditions'
# marginal correlations share the LD-weighted genetic covariance G12 lambda^2
# / M along each LD eigendirection, plus the noise their statistics share
# through shared samples, c12 lambda. Each pair is fitted in two steps: each
# condition's own heritability and noise, as estimateLocalH2() fits them, then
# the covariance with those held.
#
#   eigen / likelihood   the second condition's rotated statistics given the
#                        first's, the shared noise free (the window's own),
#                        held or zero
#   eigen / regression   the cross-products u1 u2 / lambda regressed on
#                        lambda / M and the shared-noise column
#   eigen / moments      LAVA's Omega = Delta'Delta / K - Sigma for every
#                        condition of a group at once, Sigma's off-diagonal
#                        supplied or zero; LAVA's simulation p-value
#   variant / regression cross-trait LD score regression in the window
#   variant / gls        the GLS of the cross-products with the cross-trait
#                        weight
#
# The covariance is on the heritability's scale; rg = G12 / sqrt(h1 h2), NA
# where either heritability is not positive.
# =============================================================================

#' @include AllGenerics.R RgEstimate.R h2EstimationWrappers.R
NULL

# ---- a pair's result ---------------------------------------------------------

# One pair's values, NA where the fit reports none, the members `i` and `j`
# of its unit naming the two conditions; the liability scale's are filled in
# where the pair's row is built (.rgPairRow()).
# @noRd
.rgPairResult <- function(
    i = 1L,
    j = 2L,
    nSnps = NA_integer_,
    nComponents = NA_integer_,
    rg = NA_real_,
    rgSe = NA_real_,
    rgLower = NA_real_,
    rgUpper = NA_real_,
    covariance = NA_real_,
    covarianceSe = NA_real_,
    h2 = c(NA_real_, NA_real_),
    noiseCorrelation = NA_real_,
    noiseSource = NA_character_,
    p = NA_real_,
    r2 = c(NA_real_, NA_real_),
    r2Lower = c(NA_real_, NA_real_),
    r2Upper = c(NA_real_, NA_real_),
    curve = NULL
) {
    list(
        i = i,
        j = j,
        values = c(
            list(
                rg = rg,
                rg_se = rgSe,
                rg_lower = rgLower,
                rg_upper = rgUpper,
                covariance = covariance,
                covariance_se = covarianceSe,
                h2_1 = h2[[1L]],
                h2_2 = h2[[2L]],
                covariance_liability = NA_real_,
                covariance_liability_se = NA_real_,
                h2_1_liability = NA_real_,
                h2_2_liability = NA_real_,
                noise_correlation = noiseCorrelation,
                noise_source = noiseSource,
                p_value = p
            ),
            .rgR2Columns(r2, r2Lower, r2Upper),
            list(
                n_snps = as.integer(nSnps),
                n_components = as.integer(nComponents),
                rg_loglik = list(curve)
            )
        )
    )
}

# The two conditions' r2 columns with their intervals.
# @noRd
.rgR2Columns <- function(r2, r2Lower, r2Upper) {
    list(
        r2_1 = r2[[1L]],
        r2_1_lower = r2Lower[[1L]],
        r2_1_upper = r2Upper[[1L]],
        r2_2 = r2[[2L]],
        r2_2_lower = r2Lower[[2L]],
        r2_2_upper = r2Upper[[2L]]
    )
}

# The genetic correlation of a covariance and two heritabilities, NA unless
# both heritabilities are positive.
# @noRd
.rgCorrelation <- function(covariance, h2) {
    if (anyNA(h2) || is.na(covariance) || any(h2 <= 0)) {
        return(NA_real_)
    }
    covariance / sqrt(h2[[1L]] * h2[[2L]])
}

# ---- the eigen statistics a unit's fits share -------------------------------

# The unit's eigen statistics: the LD's eigenpairs kept under the settings'
# truncation at the smallest sample size, as LAVA keeps one K for all its
# phenotypes, and each condition's rotated correlations u = Q'r on them. NULL
# when fewer directions are kept than the truncation's minimum.
# @noRd
.rgEigenData <- function(param, ud) {
    n <- map_dbl(ud$nVariant, .h2SampleSize, how = param@sampleSize)
    M <- length(ud$z[[1L]])
    dec <- .h2Basis(param, ud$ld, min(n), M)
    if (length(dec$values) < param@truncation@minimum) {
        return(NULL)
    }
    r <- map2(ud$z, ud$nVariant, .zToCor)
    list(
        lam = dec$values,
        M = M,
        n = n,
        r = r,
        u = map(r, .rgRotate, vectors = dec$vectors),
        nRef = if (identical(param@referenceNoise, "eigenTerm")) {
            ud$ld$nSamples
        } else {
            Inf
        }
    )
}

# The correlations `r` rotated onto the eigenvectors.
# @noRd
.rgRotate <- function(r, vectors) {
    as.vector(crossprod(vectors, r))
}

# The shared noise a pair's settings hold, on the scale of the fits' shared-
# noise term (whose noise correlation is it over `root`, the root of the
# product of the two conditions' own): the supplied correlation times `root`,
# or 0 for none. Tied noise holds the two phenotypes' correlation on their
# shared donors, Ph12, from which each fit takes the covariance, so that the
# shared noise is Ph12 - G12. NA when the settings fit it in the window.
# @noRd
.rgHeldCross <- function(param, ud, root) {
    switch(
        param@crossNoise,
        supplied = ud$noiseCorrelation[[1L, 2L]] * root,
        zero = 0,
        tied = ud$phenotypicCorrelation[[1L, 2L]],
        NA_real_
    )
}

# The noise correlation a pair reports: a window or tied fit's shared noise
# over the root of the product of the two conditions' own (`root`), NA where
# that is not positive; the held value otherwise, 0 or the supplied
# correlation.
# @noRd
.rgReportedNoise <- function(mode, shared, root, ud) {
    switch(
        mode,
        window = ,
        tied = if (root > 0) shared / root else NA_real_,
        zero = 0,
        ud$noiseCorrelation[[1L, 2L]]
    )
}

# ---- the likelihood's two steps ----------------------------------------------

# A pair's fit by the Gaussian likelihood of the two conditions' statistics
# rotated onto the kept LD eigenvectors: each condition's heritability and
# intercept from its own likelihood (step 1), then the covariance by the
# likelihood of the second condition's statistics given the first's, with those
# held (step 2). The p-value is the likelihood ratio against no covariance; rg
# is the fitted covariance over the root of the product of the heritabilities,
# its interval the likelihood ratio's inverted along the profile
# (.rgLikelihoodResult()); neither has a standard error.
# @noRd
.rgEigenLikelihood <- function(param, ud) {
    d <- .rgEigenData(param, ud)
    if (is.null(d)) {
        return(list(.rgPairResult(nSnps = length(ud$z[[1L]]))))
    }
    proc <- param@procedures
    members <- map(
        seq_len(2L),
        .rgLikelihoodMember,
        d = d,
        floor = proc@varianceFloor
    )
    steps <- map2(
        members,
        ud$heldIntercept,
        .rgLikelihoodStep1,
        param = param
    )
    h2 <- map_dbl(steps, "h2")
    intercepts <- map_dbl(steps, "intercept")
    base <- list(nSnps = d$M, nComponents = length(d$lam), h2 = h2)
    if (anyNA(h2) || anyNA(intercepts)) {
        return(list(exec(.rgPairResult, !!!base)))
    }
    # A heritability at its bound of 0 bounds the covariance at 0 too.
    if (any(h2 <= 0)) {
        return(list(exec(.rgPairResult, !!!base, covariance = 0)))
    }
    cross <- .rgCrossData(
        d,
        members,
        h2,
        intercepts,
        proc@varianceFloor,
        identical(param@crossNoise, "tied")
    )
    list(.rgLikelihoodCross(param, ud, d, cross, base))
}

# Condition `k`'s data for its own likelihood (.h2LikelihoodFree()).
# @noRd
.rgLikelihoodMember <- function(k, d, floor) {
    list(
        lam = d$lam,
        bstar = d$u[[k]],
        n = d$n[[k]],
        nRef = d$nRef,
        M = d$M,
        floor = floor
    )
}

# Step 1 for one condition: its own free-noise fit, or, where its intercept is
# held (a GWAS side's genome-wide intercept), its heritability alone.
# @noRd
.rgLikelihoodStep1 <- function(data, heldIntercept, param) {
    if (is.na(heldIntercept)) {
        return(.h2LikelihoodFree(param, data))
    }
    .rgLikelihoodHeld(param, data, heldIntercept)
}

# The heritability with the intercept held, by Fisher scoring from the
# eigen regression's solution where the procedures take it and their h2
# starts, keeping the best fit that converges and, where the procedures
# say, beats no heritability; NA when none does.
# @noRd
.rgLikelihoodHeld <- function(param, data, intercept) {
    proc <- param@procedures
    held <- .h2WindowData(data, intercept = intercept)
    lower <- proc@h2Bounds[[1L]]
    upper <- proc@h2Bounds[[2L]]
    starts <- c(
        if (proc@regressionStart) {
            list(.h2VarianceRegressionStart(held, lower, upper))
        },
        as.list(unique(proc@h2Starts))
    )
    nullValue <- .h2Likelihood("variance", held, 0)$value
    best <- .h2ScoringBest(
        .h2Scoring(
            "variance",
            held,
            starts,
            lower,
            upper,
            proc@maxit
        ),
        if (proc@mustBeatNull) nullValue else -Inf
    )
    if (is.null(best)) {
        return(.cisH2Result(length(data$lam), intercept = intercept))
    }
    .cisH2Result(length(data$lam), h2 = best$par[[1L]], intercept = intercept)
}

# What step 2 reads: the two conditions' rotated statistics and their own
# modelled variances at step 1's fits, and the scale of the shared noise, so
# that the cross intercept reads as the noise correlation times the root of
# the product of the two intercepts; `tied` takes the covariance off a held
# cross intercept (tied noise, Ph12 - h12).
# @noRd
.rgCrossData <- function(d, members, h2, intercepts, floor, tied = FALSE) {
    list(
        lam = d$lam,
        M = d$M,
        bstar1 = d$u[[1L]],
        bstar2 = d$u[[2L]],
        v11 = .h2WindowVariance(h2[[1L]], members[[1L]], intercepts[[1L]]),
        v22 = .h2WindowVariance(h2[[2L]], members[[2L]], intercepts[[2L]]),
        crossScale = 1 / sqrt(d$n[[1L]] * d$n[[2L]]),
        bound = sqrt(h2[[1L]] * h2[[2L]]),
        intercepts = intercepts,
        floor = floor,
        tied = as.numeric(tied)
    )
}

# What the "cross" likelihood reads of step 2's data: the two conditions'
# rotated statistics and own modelled variances, and their covariance v12 =
# offset + G theta in the parameters the fit leaves free: (h12, int12) where
# the window fits the shared noise; h12 alone with it held at `held` (less
# the covariance for tied noise); or int12 alone with the covariance held at
# `h12` (a profile, the null at 0).
# @noRd
.rgCrossScoringData <- function(cross, held = NA_real_, h12 = NULL) {
    genetic <- cross$lam^2 / cross$M
    noise <- cross$lam * cross$crossScale
    window <- is.na(held)
    G <- if (!is.null(h12)) {
        cbind(noise)
    } else if (window) {
        cbind(genetic, noise)
    } else {
        cbind(genetic - cross$tied * noise)
    }
    list(
        b1 = cross$bstar1,
        b2 = cross$bstar2,
        v11 = cross$v11,
        v22 = cross$v22,
        G = G,
        offset = (h12 %||% 0) * genetic + (if (window) 0 else held) * noise,
        floor = cross$floor
    )
}

# Step 2: the covariance, and the cross intercept where the settings fit it in
# the window, by Fisher scoring (.rgCrossBest()), keeping the best fit that
# converges and, where the procedures say, beats no covariance. A window fit
# that never does falls back to the supplied noise correlation where one is
# given; otherwise the pair has no estimate.
# @noRd
.rgLikelihoodCross <- function(param, ud, d, cross, base) {
    fit <- .rgLikelihoodCrossFit(param, ud, d, cross)
    if (is.null(fit)) {
        return(exec(.rgPairResult, !!!base))
    }
    .rgLikelihoodResult(param, ud, cross, fit, base)
}

# The best step-2 fit under the settings' shared noise, or under the supplied
# noise correlation where a window fit has none; NULL when there is none.
# @noRd
.rgLikelihoodCrossFit <- function(param, ud, d, cross) {
    mode <- param@crossNoise
    root <- sqrt(prod(cross$intercepts))
    held <- .rgHeldCross(param, ud, root)
    fit <- .rgCrossBest(param, ud, d, cross, mode, held)
    if (!is.null(fit) || !identical(mode, "window")) {
        return(fit)
    }
    if (is.null(ud$noiseCorrelation)) {
        return(NULL)
    }
    fallback <- ud$noiseCorrelation[[1L, 2L]] * root
    .rgCrossBest(param, ud, d, cross, "fallback", fallback)
}

# The best accepted step-2 fit under `mode` ("window", or a held cross
# intercept `held`: "supplied", "fallback" or "zero"), by Fisher scoring
# (.h2Scoring()) from .rgCrossStarts() and, where the window fits the cross
# intercept, from no covariance's own fit (which nests in the full model, so
# a run from it ends no lower), with the covariance within plus or minus the
# root of the product of the heritabilities and a fitted cross intercept
# within its bounds: the best that converges and, where the procedures say,
# beats no covariance (.rgCrossNullFit()); NULL when none does.
# @noRd
.rgCrossBest <- function(param, ud, d, cross, mode, held) {
    proc <- param@procedures
    window <- identical(mode, "window")
    data <- .rgCrossScoringData(cross, if (window) NA_real_ else held)
    starts <- .rgCrossStarts(proc, ud, cross, data, window)
    null <- .rgCrossNullFit(cross, proc, window, held, starts)
    if (is.null(null)) {
        return(NULL)
    }
    bounds <- proc@crossInterceptBounds
    best <- .h2ScoringBest(
        .h2Scoring(
            "cross",
            data,
            c(starts, if (window) list(c(0, null$par))),
            c(-cross$bound, if (window) bounds[[1L]]),
            c(cross$bound, if (window) bounds[[2L]]),
            proc@maxit
        ),
        if (proc@mustBeatNull) null$value else -Inf
    )
    if (is.null(best)) {
        return(NULL)
    }
    list_assign(best, nullValue = null$value, mode = mode, held = held)
}

# Step 2's starts: the cross-products' regression on the covariance's columns
# where the procedures take a regression start, then 0 and half the
# covariance's bound either way; where the window fits the cross intercept,
# each with it at the supplied noise correlation (when given), the
# correlation of the two statistics over |z| < 8, 1 and 0.
# @noRd
.rgCrossStarts <- function(proc, ud, cross, data, window) {
    regression <- if (proc@regressionStart) {
        .h2LeastSquares(data$G, data$b1 * data$b2 - data$offset)
    }
    h12 <- unique(c(regression[1L], 0, -cross$bound / 2, cross$bound / 2))
    if (!window) {
        return(as.list(pmin(pmax(h12, -cross$bound), cross$bound)))
    }
    supplied <- if (!is.null(ud$noiseCorrelation)) {
        ud$noiseCorrelation[[1L, 2L]] * sqrt(prod(cross$intercepts))
    }
    int12 <- c(regression[2L], supplied, .rgStatisticCorrelation(ud$z), 1, 0)
    grid <- expand.grid(int12 = unique(int12), h12 = h12)
    unique(map2(grid$h12, grid$int12, c))
}

# The no-covariance fit beside step 2: the cross intercept refitted from the
# starts' where the window fits it, the held one otherwise; NULL where no
# refit converges.
# @noRd
.rgCrossNullFit <- function(cross, proc, window, held, starts) {
    if (!window) {
        data <- .rgCrossScoringData(cross, held)
        return(list(value = .h2Likelihood("cross", data, 0)$value))
    }
    .h2ScoringBest(.h2Scoring(
        "cross",
        .rgCrossScoringData(cross, h12 = 0),
        as.list(unique(map_dbl(starts, 2L))),
        proc@crossInterceptBounds[[1L]],
        proc@crossInterceptBounds[[2L]],
        proc@maxit,
        valueOnly = TRUE
    ))
}

# The step-2 log-likelihood at covariance h12, the cross intercept refitted
# from `start` where the window fits it (`held` NA), held otherwise: a point
# of its profile (.rgProfile()), list(value, par) with the refitted intercept
# its par; NULL where the refit does not converge.
# @noRd
.rgCrossRefit <- function(h12, start, cross, proc, held) {
    if (!is.na(held)) {
        data <- .rgCrossScoringData(cross, held)
        value <- .h2Likelihood("cross", data, h12)$value
        return(list(value = value, par = numeric(0)))
    }
    .h2ScoringBest(.h2Scoring(
        "cross",
        .rgCrossScoringData(cross, h12 = h12),
        list(start),
        proc@crossInterceptBounds[[1L]],
        proc@crossInterceptBounds[[2L]],
        proc@maxit,
        valueOnly = TRUE
    ))
}

# The correlation of the two conditions' statistics over the variants where both
# are below 8 in size, a start for the cross intercept; 0 where there is none to
# take.
# @noRd
#' @importFrom stats cor
.rgStatisticCorrelation <- function(z) {
    small <- abs(z[[1L]]) < 8 & abs(z[[2L]]) < 8
    rho <- if (sum(small) > 2L) {
        suppressWarnings(cor(z[[1L]][small], z[[2L]][small]))
    } else {
        NA_real_
    }
    if (is.na(rho)) 0 else rho
}

# The pair's values from its step-2 fit: the covariance it fits and rg its
# ratio to the root of the product of the heritabilities, the likelihood
# ratio's p-value on one degree of freedom, and rg's 95% interval and curve
# from the profile (.rgCrossRefit(), .rgProfile()).
# @noRd
#' @importFrom stats pchisq
.rgLikelihoodResult <- function(param, ud, cross, fit, base) {
    window <- identical(fit$mode, "window")
    shared <- if (window) fit$par[[2L]] else fit$held
    int12 <- shared - cross$tied * fit$par[[1L]]
    bound <- cross$bound
    profile <- .rgProfile(
        .RG_CURVE_GRID * bound,
        list(estimate = fit$par[[1L]], par = fit$par[-1L], value = fit$value),
        .rgCrossRefit,
        cross = cross,
        proc = param@procedures,
        held = if (window) NA_real_ else fit$held
    )
    lr <- 2 * (fit$value - fit$nullValue)
    exec(
        .rgPairResult,
        !!!base,
        curve = profile$curve,
        rg = fit$par[[1L]] / bound,
        rgLower = profile$interval[[1L]] / bound,
        rgUpper = profile$interval[[2L]] / bound,
        covariance = fit$par[[1L]],
        noiseCorrelation = .rgReportedNoise(
            fit$mode,
            int12,
            sqrt(prod(cross$intercepts)),
            ud
        ),
        noiseSource = if (identical(fit$mode, "fallback")) {
            "fallback"
        } else {
            param@crossNoise
        },
        p = if (identical(param@inference, "lrt")) {
            pchisq(max(lr, 0), 1, lower.tail = FALSE)
        } else {
            NA_real_
        }
    )
}

# ---- the profile along rg ----------------------------------------------------

# rg's profile log-likelihood on the RgEstimate's rg grid (`grid`, on the
# fit's scale) and its 95% interval from inverting the likelihood ratio,
# about `fit` (its covariance `estimate`, its other parameters `par` and its
# `value`); `refit(value, start, ...)` refits the other parameters with the
# covariance held at `value` from `start`, list(value, par), NULL where that
# does not converge. The grid is walked out from the fit (.rgProfileWalk()),
# and each end of the interval is root-found between the two points where
# the walk first falls past the cutoff (.rgWalkEnd()), the grid's end where
# it never does: list(curve, interval), the curve relative to its largest.
# @noRd
.rgProfile <- function(grid, fit, refit, ...) {
    walk <- .rgProfileWalk(grid, fit, refit, ...)
    below <- rev(which(grid < fit$estimate))
    above <- which(grid > fit$estimate)
    list(
        curve = walk$value - max(walk$value),
        interval = c(
            .rgWalkEnd(below, grid[[1L]], walk, fit, refit, ...),
            .rgWalkEnd(above, grid[[length(grid)]], walk, fit, refit, ...)
        )
    )
}

# The profile walked out from the grid point nearest the fit either way, each
# point refitted from its neighbours' fits (.rgWalkStep()): the grid, and the
# value and parameters at each of its points.
# @noRd
#' @importFrom purrr accumulate
.rgProfileWalk <- function(grid, fit, refit, ...) {
    k0 <- which.min(abs(grid - fit$estimate))
    start <- list(par = fit$par)
    up <- accumulate(
        grid[k0:length(grid)],
        .rgWalkStep,
        .init = start,
        refit = refit,
        ...
    )
    down <- accumulate(
        rev(grid[seq_len(k0 - 1L)]),
        .rgWalkStep,
        .init = start,
        refit = refit,
        ...
    )
    points <- c(rev(down[-1L]), up[-1L])
    list(
        grid = grid,
        value = map_dbl(points, "value"),
        par = map(points, "par")
    )
}

# One point of .rgProfileWalk(): the refit at `value` started on the line
# through the last two points' fits, or at the last one's where no point
# comes before it (the fit starts the first point, too close to it to steer
# a line) or a refit failed; -Inf where this one fails, the last fit kept
# for the next to start from.
# @noRd
.rgWalkStep <- function(previous, value, refit, ...) {
    before <- previous$before
    start <- if (is.null(before)) {
        previous$par
    } else {
        reach <- (value - previous$at) / (previous$at - before$at)
        previous$par + reach * (previous$par - before$par)
    }
    fit <- refit(value, start, ...)
    if (is.null(fit)) {
        return(list(par = previous$par, value = -Inf))
    }
    list(
        at = value,
        par = fit$par,
        value = fit$value,
        before = if (!is.null(previous$at)) previous[c("at", "par")]
    )
}

# One end of .rgProfile()'s interval, over the walk's points `toward` its end
# of the grid (`end`) in walking order, the fit first: the root of
# .h2LrExcess() between the first point past the cutoff and the one before
# it, each refit started on the line between their fits (.rgWalkProfile());
# `end` where no point is past it.
# @noRd
#' @importFrom stats qchisq uniroot
.rgWalkEnd <- function(toward, end, walk, fit, refit, ...) {
    cutoff <- qchisq(0.95, 1) / 2
    path <- list(
        at = c(fit$estimate, walk$grid[toward]),
        par = c(list(fit$par), walk$par[toward]),
        excess = c(-cutoff, pmin(fit$value - walk$value[toward] - cutoff, 1e10))
    )
    past <- match(TRUE, path$excess > 0)
    if (is.na(past)) {
        return(end)
    }
    pair <- c(past - 1L, past)[order(path$at[c(past - 1L, past)])]
    uniroot(
        .h2LrExcess,
        path$at[pair],
        f.lower = path$excess[[pair[[1L]]]],
        f.upper = path$excess[[pair[[2L]]]],
        top = fit$value,
        profile = .rgWalkProfile,
        ends = list(at = path$at[pair], par = path$par[pair]),
        refit = refit,
        ...,
        tol = 1e-8
    )$root
}

# The profile at `value` between two points of a walk (`ends`), refitted from
# the line between their fits.
# @noRd
.rgWalkProfile <- function(value, ends, refit, ...) {
    share <- (value - ends$at[[1L]]) / (ends$at[[2L]] - ends$at[[1L]])
    start <- ends$par[[1L]] + share * (ends$par[[2L]] - ends$par[[1L]])
    refit(value, start, ...)$value %||% -Inf
}

# ---- the likelihood of a group at once ---------------------------------------

# A group's fit by the likelihood of all its conditions' rotated statistics at
# once (H2LikelihoodParam(jointFit = TRUE)): each condition's heritability and
# intercept from its own likelihood (step 1, as for a pair), then every pair's
# covariance and, where the window fits it, shared noise together (step 2), the
# genetic correlation matrix built from canonical partial correlations so that
# it is one. At two conditions this is a pair's own step 2. Each pair's p-value
# is the likelihood ratio of its covariance at zero with the rest refitted, and
# its interval is where the likelihood ratio along its correlation, the rest
# held, stays within the chi-square's cutoff.
# @noRd
.rgJointLikelihood <- function(param, ud) {
    P <- length(ud$z)
    pairs <- .rgGroupPairs(P)
    d <- .rgEigenData(param, ud)
    if (is.null(d)) {
        return(pmap(pairs, .rgPairResult, nSnps = length(ud$z[[1L]])))
    }
    proc <- param@procedures
    members <- map(
        seq_len(P),
        .rgLikelihoodMember,
        d = d,
        floor = proc@varianceFloor
    )
    steps <- map2(
        members,
        ud$heldIntercept,
        .rgLikelihoodStep1,
        param = param
    )
    h2 <- map_dbl(steps, "h2")
    intercepts <- map_dbl(steps, "intercept")
    base <- list(nSnps = d$M, nComponents = length(d$lam))
    fitted <- !is.na(h2) & !is.na(intercepts)
    if (!all(fitted & h2 > 0)) {
        return(.rgJointLive(param, ud, h2, fitted, base))
    }
    joint <- .rgJointData(
        d,
        members,
        h2,
        intercepts,
        param,
        .rgTieOf(param, ud, length(h2))
    )
    fit <- .rgJointFit(param, ud, d, members, joint)
    if (is.null(fit)) {
        return(pmap(pairs, .rgJointNoEstimate, h2 = h2, base = base))
    }
    pmap(
        pairs,
        .rgJointPairResult,
        fit = fit,
        joint = joint,
        ud = ud,
        base = base,
        param = param
    )
}

# A pair of the group without an estimate: its covariance at 0 where both
# sides have a step-1 fit and one's heritability sits at its bound of 0, as
# a pair's own fit has it (`fitted`: each condition's step 1 fitted).
# @noRd
.rgJointNoEstimate <- function(i, j, h2, base, fitted = !is.na(h2)) {
    at <- c(i, j)
    bounded <- all(fitted[at]) && any(h2[at] <= 0)
    exec(
        .rgPairResult,
        i = i,
        j = j,
        !!!base,
        h2 = h2[at],
        covariance = if (bounded) 0 else NA_real_
    )
}

# A group some of whose conditions step 1 leaves without a heritability (at
# its bound of 0, or without a fit): the others fitted without them, at once
# where three or more remain (.rgFit()); a pair with a side at 0 has its
# covariance at 0 too, and one with a side unfitted none
# (.rgJointNoEstimate()).
# @noRd
.rgJointLive <- function(param, ud, h2, fitted, base) {
    pairs <- .rgGroupPairs(length(h2))
    at <- which(fitted & h2 > 0)
    refit <- if (length(at) >= 2L) .rgFit(param, .rgSubUnit(ud, at))
    pmap(
        pairs,
        .rgJointLivePair,
        at = at,
        refit = refit,
        h2 = h2,
        fitted = fitted,
        base = base
    )
}

# Pair (i, j) of .rgJointLive(): the refit's where both sides are in it.
# @noRd
.rgJointLivePair <- function(i, j, at, refit, h2, fitted, base) {
    sides <- match(c(i, j), at)
    if (anyNA(sides)) {
        return(.rgJointNoEstimate(i, j, h2, base, fitted))
    }
    k <- detect_index(refit, .rgSamePair, sides = sides)
    list_assign(refit[[k]], i = i, j = j)
}

# @noRd
.rgSamePair <- function(pair, sides) {
    pair$i == sides[[1L]] && pair$j == sides[[2L]]
}

# What step 2 reads: every condition's rotated statistics (a column each) and
# own modelled variances at step 1's fits, the scale of the noise each pair
# shares, the noise terms' roots, the settings' floor, and where tied noise
# ties a pair's (.rgTieOf()).
# @noRd
.rgJointData <- function(d, members, h2, intercepts, param, tie) {
    variances <- pmap(list(h2, members, intercepts), .h2WindowVariance)
    rotated <- d$u
    P <- length(h2)
    list(
        P = P,
        nPairs = P * (P - 1L) / 2L,
        u = exec(cbind, !!!rotated),
        v = exec(cbind, !!!variances),
        lam = d$lam,
        lam2M = d$lam^2 / d$M,
        crossScale = 1 / sqrt(tcrossprod(d$n)),
        h2 = h2,
        intercepts = intercepts,
        root = sqrt(tcrossprod(intercepts)),
        floor = param@procedures@varianceFloor,
        tied = tie
    )
}

# The log-likelihood of every condition's statistics at canonical partial
# correlations theta[1:nPairs] (row by row of the lower triangle) and, where
# the window fits it, the pairs' shared noise after them (or `held`): the
# "joint" likelihood (.h2Scoring()), each condition's density given those
# before it through the LDL' decomposition of each direction's covariance
# (block by block on several references), each conditional variance floored
# at the procedures' floor. At two conditions the second's term is the
# pair's own, the "cross" likelihood's.
# @noRd
.rgJointLl <- function(theta, joint, held) {
    .h2Likelihood("joint", list(joint = joint, held = held), theta)$value
}

# ---- canonical partial correlations ------------------------------------------

# The Cholesky factor of the correlation matrix of canonical partial
# correlations `z`, row by row of the lower triangle (z_21, z_31, z_32, ...):
# W_i1 = z_i1, W_ij = z_ij times the root of the product of 1 - z_im^2 over m
# < j, and W_ii that root over m < i, so that every row has length 1. Any `z`
# within [-1, 1] gives a correlation matrix, and z_21 is the correlation of
# the first two conditions.
# @noRd
.rgCpcFactor <- function(z, P) {
    rows <- map(seq_len(P), .rgCpcRow, z = z, P = P)
    exec(rbind, !!!rows)
}

# @noRd
.rgCpcRow <- function(i, z, P) {
    zi <- z[(i - 1L) * (i - 2L) / 2L + seq_len(i - 1L)]
    lengths <- c(1, cumprod(sqrt(1 - zi^2)))
    c(zi * lengths[seq_len(i - 1L)], lengths[[i]], rep(0, P - i))
}

# The canonical partial correlations of a correlation matrix `R`: its
# Cholesky factor's entries over what each row has left. A matrix on the
# boundary is moved a hair inside it so that it has a factor.
# @noRd
.rgCorrelationCpc <- function(R) {
    P <- nrow(R)
    W <- t(chol((1 - 1e-10) * R + 1e-10 * diag(P)))
    z <- list_c(map(seq_len(P), .rgCpcOfRow, W = W))
    pmin(pmax(z, -1), 1)
}

# @noRd
.rgCpcOfRow <- function(i, W) {
    w <- W[i, seq_len(i - 1L)]
    used <- c(0, cumsum(w^2))[seq_len(i - 1L)]
    w / sqrt(pmax(1 - used, .Machine$double.eps))
}

# The symmetric matrix of the pairs' `values` (row by row of the lower
# triangle), 0 on the diagonal.
# @noRd
.rgPairMatrix <- function(values, P) {
    at <- outer(seq_len(P), seq_len(P), .rgPairPosition)
    index <- as.vector(replace(at, is.na(at), length(values) + 1L))
    matrix(c(values, 0)[index], P, P)
}

# The position of pair (k, l) row by row of the lower triangle; NA on the
# diagonal.
# @noRd
.rgPairPosition <- function(k, l) {
    hi <- pmax(k, l)
    lo <- pmin(k, l)
    if_else(hi == lo, NA_real_, (hi - 1) * (hi - 2) / 2 + lo)
}

# A symmetric matrix's lower triangle, row by row.
# @noRd
.rgPairValues <- function(M) {
    P <- nrow(M)
    M[cbind(rep(seq_len(P)[-1L], seq_len(P - 1L)), sequence(seq_len(P - 1L)))]
}

# ---- the joint fit -----------------------------------------------------------

# Step 2 at once under the settings' shared noise, or, where a window fit has
# none that converges and beats no covariance, under the supplied noise
# correlation; NULL when there is none.
# @noRd
.rgJointFit <- function(param, ud, d, members, joint) {
    mode <- param@crossNoise
    held <- switch(
        mode,
        supplied = ud$noiseCorrelation * joint$root,
        zero = matrix(0, joint$P, joint$P),
        tied = replace(
            ud$phenotypicCorrelation,
            is.na(ud$phenotypicCorrelation),
            0
        ),
        NULL
    )
    fit <- .rgJointBest(param, ud, d, members, joint, mode, held)
    if (
        !is.null(fit) ||
            !identical(mode, "window") ||
            is.null(ud$noiseCorrelation)
    ) {
        return(fit)
    }
    fallback <- ud$noiseCorrelation * joint$root
    .rgJointBest(param, ud, d, members, joint, "fallback", fallback)
}

# The best fit by Fisher scoring (.h2Scoring(), the "joint" likelihood) from the
# pairs' own fits put together, from no covariance and from the null's own fit
# (which nests in the full model, so a run from it ends no lower), that
# converges and, where the procedures say, beats no covariance; NULL when none
# does, or when the null's refit does not converge.
# @noRd
.rgJointBest <- function(param, ud, d, members, joint, mode, held) {
    proc <- param@procedures
    window <- identical(mode, "window")
    start <- .rgJointStart(param, ud, d, members, joint, held)
    sharedStart <- if (window) start$shared
    starts <- list(
        c(start$z, sharedStart),
        c(rep(0, joint$nPairs), sharedStart)
    )
    null <- .rgJointCoordinateNull(
        starts[[2L]],
        seq_len(joint$nPairs),
        joint,
        proc,
        held
    )
    if (is.null(null)) {
        return(NULL)
    }
    bounds <- .rgJointBounds(joint, proc, held)
    best <- .h2ScoringBest(
        .h2Scoring(
            "joint",
            .rgJointScoringData(joint, held, bounds),
            c(starts, list(null$par)),
            bounds$lower,
            bounds$upper,
            proc@maxit
        ),
        if (proc@mustBeatNull) null$value else -Inf
    )
    if (is.null(best)) {
        return(NULL)
    }
    list_assign(best, mode = mode, held = held, window = window)
}

# What the "joint" likelihood reads: the group's data, the held shared noise
# (NULL where fitted) and the parameters' bounds, within which its
# covariances' derivatives are differenced.
# @noRd
.rgJointScoringData <- function(joint, held, bounds) {
    list(
        joint = joint,
        held = held,
        lower = bounds$lower,
        upper = bounds$upper
    )
}

# The joint fit's start: every pair's own fit (its step 2 on the two
# conditions), their correlations put together and moved to the nearest
# correlation matrix shrunk a little inside the boundary, as canonical partial
# correlations; and their shared noise, the statistics' correlation over |z| < 8
# for a pair whose own fit has none.
# @noRd
.rgJointStart <- function(param, ud, d, members, joint, held) {
    if (isTRUE(joint$mixed)) {
        return(.rgMixedStart(param, ud, joint))
    }
    pairs <- .rgGroupPairs(joint$P)
    own <- pmap(
        pairs,
        .rgJointPairStart,
        param = param,
        ud = ud,
        d = d,
        members = members,
        joint = joint
    )
    R <- .rgPairsMatrix(map_dbl(own, "r"), pairs, joint$P) + diag(joint$P)
    nearest <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    list(
        z = .rgCorrelationCpc(0.98 * nearest + 0.02 * diag(joint$P)),
        shared = .rgPairValues(
            .rgPairsMatrix(map_dbl(own, "shared"), pairs, joint$P)
        )
    )
}

# The symmetric matrix of values given for `pairs` (.rgGroupPairs()'s
# order), 0 on the diagonal.
# @noRd
.rgPairsMatrix <- function(values, pairs, P) {
    .rgPairMatrix(values[order(.rgPairPosition(pairs$i, pairs$j))], P)
}

# One pair's own step-2 fit, as a start: its correlation and shared noise.
# @noRd
.rgJointPairStart <- function(i, j, param, ud, d, members, joint) {
    at <- c(i, j)
    pairD <- list_modify(d, u = d$u[at], n = d$n[at], r = d$r[at])
    pairUd <- .rgPairUnitData(ud, at)
    cross <- .rgCrossData(
        pairD,
        members[at],
        joint$h2[at],
        joint$intercepts[at],
        joint$floor,
        joint$tied[[i, j]] > 0
    )
    fit <- .rgLikelihoodCrossFit(param, pairUd, pairD, cross)
    fallback <- .rgStatisticCorrelation(pairUd$z)
    if (is.null(fit)) {
        return(list(r = 0, shared = fallback))
    }
    window <- identical(fit$mode, "window")
    list(
        r = fit$par[[1L]] / cross$bound,
        shared = if (window) fit$par[[2L]] else fallback
    )
}

# A unit's data for two of its conditions, `at`, with the LD they read, as a
# pair fitted on its own takes it (given no other condition).
# @noRd
.rgPairUnitData <- function(ud, at) {
    ld <- .rgSubLd(ud, at)
    c(
        list(z = ud$z[at], nVariant = ud$nVariant[at]),
        ld,
        list(
            overlap = ud$overlap[at],
            heldIntercept = ud$heldIntercept[at],
            noiseCorrelation = .rgSubMatrix(ud$noiseCorrelation, at),
            phenotypicCorrelation = .rgSubMatrix(ud$phenotypicCorrelation, at)
        )
    )
}

# A unit's data for its conditions `at`, fitted as a group of their own:
# .rgPairUnitData()'s, given the same conditions as the unit's.
# @noRd
.rgSubUnit <- function(ud, at) {
    c(.rgPairUnitData(ud, at), list(conditional = ud$conditional))
}

# The LD a unit's conditions `at` read: the unit's one source, the one
# reference they all read, or the references they read with each one's.
# @noRd
.rgSubLd <- function(ud, at) {
    if (is.null(ud$lds)) {
        return(list(ld = ud$ld))
    }
    of <- ud$ldOf[at]
    used <- unique(of)
    if (length(used) == 1L) {
        return(list(ld = ud$lds[[used]]))
    }
    list(ld = NULL, lds = ud$lds[used], ldOf = match(of, used))
}

# @noRd
.rgSubMatrix <- function(x, at) {
    if (is.null(x)) NULL else x[at, at, drop = FALSE]
}

# The joint fit's bounds: [-1, 1] for each partial correlation, and the cross
# intercept's bounds for each shared noise the window fits.
# @noRd
.rgJointBounds <- function(joint, proc, held) {
    window <- is.null(held)
    bounds <- proc@crossInterceptBounds
    shared <- if (window) joint$nPairs else 0L
    list(
        lower = c(rep(-1, joint$nPairs), rep(bounds[[1L]], shared)),
        upper = c(rep(1, joint$nPairs), rep(bounds[[2L]], shared))
    )
}

# The best fit with the parameters `at` held at `values` (0: the null) and
# the rest refitted from `theta`'s by Fisher scoring, its value alone read
# (its parameters only start a neighbour's refit); the held fit itself where
# nothing is left; NULL where the refit does not converge.
# @noRd
.rgJointCoordinateNull <- function(theta, at, joint, proc, held, values = 0) {
    fixed <- replace(theta, at, values)
    if (length(at) == length(theta)) {
        return(list(value = .rgJointLl(fixed, joint, held), par = fixed))
    }
    free <- setdiff(seq_along(theta), at)
    bounds <- .rgJointBounds(joint, proc, held)
    best <- .h2ScoringBest(.h2Scoring(
        "joint",
        .rgJointScoringData(joint, held, bounds),
        list(fixed[free]),
        bounds$lower[free],
        bounds$upper[free],
        proc@maxit,
        valueOnly = TRUE,
        full = fixed,
        free = free
    ))
    if (is.null(best)) {
        return(NULL)
    }
    list(value = best$value, par = replace(fixed, free, best$par))
}

# The joint fit's profile at partial correlation `at` = `value`, the rest
# refitted from `start`'s (.rgJointCoordinateNull()): a point of
# .rgProfile().
# @noRd
.rgJointProfile <- function(value, start, at, joint, held, proc) {
    .rgJointCoordinateNull(start, at, joint, proc, held, value)
}

# ---- a pair from the joint fit -----------------------------------------------

# Pair (i, j) of the joint fit: its correlation, covariance and the two
# heritabilities, or, given the group's other conditions (conditional =
# TRUE), its partial correlation, partial covariance and the two variances
# left, with each condition's r2; its test and interval (.rgJointPairTest());
# and the noise it shares.
# @noRd
.rgJointPairResult <- function(i, j, fit, joint, ud, base, param) {
    conditional <- isTRUE(ud$conditional)
    rest <- setdiff(seq_len(joint$P), c(i, j))
    R <- tcrossprod(.rgCpcFactor(fit$par[seq_len(joint$nPairs)], joint$P))
    G <- R * sqrt(tcrossprod(joint$h2))
    shared <- .rgJointShared(fit, R, joint)
    partial <- .rgPartial(G, i, j, if (conditional) rest else integer(0))
    test <- .rgJointPairTest(i, j, fit, joint, R, shared, param, conditional)
    source <- .rgPairNoiseSource(fit$mode, joint$tied, i, j)
    exec(
        .rgPairResult,
        i = i,
        j = j,
        !!!base,
        rg = .rgCorrelation(partial$cov, partial$var),
        rgLower = test$interval[[1L]],
        rgUpper = test$interval[[2L]],
        covariance = partial$cov,
        h2 = partial$var,
        noiseCorrelation = .rgReportedNoise(
            source,
            shared[[i, j]],
            joint$root[[i, j]],
            list(noiseCorrelation = .rgSubMatrix(ud$noiseCorrelation, c(i, j)))
        ),
        noiseSource = source,
        p = test$p,
        curve = test$curve,
        r2 = if (conditional) {
            c(.rgExplained(G, i), .rgExplained(G, j))
        } else {
            c(NA_real_, NA_real_)
        }
    )
}

# A pair's likelihood-ratio test and 95% interval on the coordinate that is
# its correlation: with the pair first, where its correlation is the first
# partial correlation, or, given the rest, last, where it is the last (the
# partial correlation given all before it). The test refits every other
# parameter with that coordinate at 0, and the interval and curve with it at
# each value they read (.rgJointProfile(), .rgProfile()).
# @noRd
#' @importFrom stats pchisq
.rgJointPairTest <- function(i, j, fit, joint, R, shared, param, conditional) {
    rest <- setdiff(seq_len(joint$P), c(i, j))
    order <- if (conditional) c(rest, i, j) else c(i, j, rest)
    at <- if (conditional) joint$nPairs else 1L
    perm <- .rgJointPermute(joint, order)
    held <- .rgSubMatrix(fit$held, order)
    theta <- c(
        .rgCorrelationCpc(R[order, order]),
        if (fit$window) .rgPairValues(shared[order, order])
    )
    proc <- param@procedures
    ll <- .rgJointLl(theta, perm, held)
    null <- .rgJointCoordinateNull(theta, at, perm, proc, held)
    profile <- .rgProfile(
        .RG_CURVE_GRID,
        list(estimate = theta[[at]], par = theta, value = ll),
        .rgJointProfile,
        at = at,
        joint = perm,
        held = held,
        proc = proc
    )
    list(
        p = if (identical(param@inference, "lrt") && !is.null(null)) {
            pchisq(max(2 * (ll - null$value), 0), 1, lower.tail = FALSE)
        } else {
            NA_real_
        },
        interval = profile$interval,
        curve = profile$curve
    )
}

# The noise the pairs share at the joint fit: fitted in the window, the held
# value, or, tied, the phenotypes' correlation less the covariance.
# @noRd
.rgJointShared <- function(fit, R, joint) {
    if (fit$window) {
        return(.rgPairMatrix(fit$par[-seq_len(joint$nPairs)], joint$P))
    }
    fit$held - joint$tied * R * sqrt(tcrossprod(joint$h2))
}

# The step-2 data with the conditions in `order`.
# @noRd
.rgJointPermute <- function(joint, order) {
    if (isTRUE(joint$mixed)) {
        return(list_assign(
            joint,
            u = joint$u[order],
            v = joint$v[order],
            basisOf = joint$basisOf[order],
            crossScale = joint$crossScale[order, order],
            h2 = joint$h2[order],
            intercepts = joint$intercepts[order],
            root = joint$root[order, order],
            tied = joint$tied[order, order]
        ))
    }
    list_modify(
        joint,
        u = joint$u[, order, drop = FALSE],
        v = joint$v[, order, drop = FALSE],
        crossScale = joint$crossScale[order, order],
        h2 = joint$h2[order],
        intercepts = joint$intercepts[order],
        root = joint$root[order, order],
        tied = joint$tied[order, order]
    )
}

# ---- the regressions: the cross-products on the model's columns -------------

# A pair's fit by the eigen regression: each condition's whitened squares u^2 /
# lambda on (lambda / M, 1 / n) as estimateLocalH2()'s eigen regression fits
# them, then the whitened cross-products u1 u2 / lambda on (lambda / M, 1 /
# sqrt(n1 n2)). No cross term for the reference's noise, as the likelihood has
# none.
# @noRd
.rgEigenRegression <- function(param, ud) {
    d <- .rgEigenData(param, ud)
    if (is.null(d) || length(d$lam) < 3L) {
        return(list(.rgPairResult(nSnps = length(ud$z[[1L]]))))
    }
    taper <- .h2Taper(param@procedures@eigenvalueTaper, d$lam)
    rows <- map(seq_len(2L), .rgEigenRows, d = d, taper = taper)
    cross <- list(
        y = d$u[[1L]] * d$u[[2L]] / d$lam,
        X = cbind(d$lam / d$M, 1 / sqrt(d$n[[1L]] * d$n[[2L]])),
        taper = taper,
        redundancy = rep(1, length(d$lam))
    )
    fit <- .rgRegressionPair(param, ud, rows, cross)
    list(.rgRegressionResult(param, fit, d$M, length(d$lam)))
}

# Condition `k`'s rows for the eigen regression, as .h2EigenRegression() forms
# them.
# @noRd
.rgEigenRows <- function(k, d, taper) {
    .rgEigenSquares(d$u[[k]], d$lam, d$M, d$nRef, d$n[[k]], taper)
}

# The whitened squares u^2 / lambda of rotated correlations `u` on the
# directions `lam`, on (lambda / M - 1 / nRef, 1 / n).
# @noRd
.rgEigenSquares <- function(u, lam, M, nRef, n, taper) {
    y <- u^2 / lam
    list(
        y = y,
        X = cbind(lam / M - 1 / nRef, 1 / n),
        chi2 = n * y,
        n = n,
        taper = taper,
        redundancy = rep(1, length(y)),
        lam = lam,
        weights = "variance"
    )
}

# A pair's fit by the variant regression: each condition's n r^2 on (n l / M,
# 1) as estimateLocalH2()'s variant regression fits them, then the
# cross-products sqrt(n1 n2) r1 r2 on (sqrt(n1 n2) l / M, 1), cross-trait LD
# score regression in the window, with the LD-redundancy weight where the
# procedures take it.
# @noRd
.rgVariantRegression <- function(param, ud) {
    M <- length(ud$z[[1L]])
    if (M < 3L) {
        return(list(.rgPairResult(nSnps = M)))
    }
    n <- map_dbl(ud$nVariant, .h2SampleSize, how = param@sampleSize)
    r <- map2(ud$z, ud$nVariant, .zToCor)
    l <- .h2WindowLdScores(param, ud$ld, M)
    redundancy <- .rgRedundancy(l, ud$ld, param)
    rows <- map(
        seq_len(2L),
        .rgVariantRows,
        r = r,
        n = n,
        l = l,
        w = redundancy
    )
    nn <- sqrt(n[[1L]] * n[[2L]])
    cross <- list(
        y = nn * r[[1L]] * r[[2L]],
        X = cbind(nn * l / M, 1),
        taper = rep(1, M),
        redundancy = redundancy,
        ldModel = .rgVariantModel(param, list(ud$ld, ud$ld), n)
    )
    fit <- .rgRegressionPair(param, ud, rows, cross)
    list(.rgRegressionResult(param, fit, M, NA_integer_))
}

# What a variant regression pair's sandwich needs to form its statistics'
# covariances under the model (.rgModelS()): each condition's LD source (the
# same one twice on one reference), the frequency weights, each condition's
# sample size and the window's size.
# @noRd
.rgVariantModel <- function(param, lds, n) {
    list(
        lds = lds,
        weights = map(lds, .h2SourceWeights, param = param),
        shared = .rgCrossWeights(param, lds),
        n = n,
        M = .ldSourceSize(lds[[1L]])
    )
}

# The LD-redundancy weight 1 / max(l, 1) where the procedures take it, l the
# window's scores with every variant counted once (.h2RedundancyScores()).
# @noRd
.rgRedundancy <- function(l, ld, param) {
    1 / pmax(.h2RedundancyScores(param, ld, l), 1)
}

# Condition `k`'s rows for the variant regression, as .h2VariantRegression()
# forms them.
# @noRd
.rgVariantRows <- function(k, r, n, l, w) {
    y <- n[[k]] * r[[k]]^2
    list(
        y = y,
        X = cbind(n[[k]] * l / length(l), 1),
        chi2 = y,
        n = n[[k]],
        taper = rep(1, length(y)),
        redundancy = w,
        ldTotal = l,
        wLd = 1 / w,
        M = length(l),
        hsqScale = 1,
        weights = "ldScore",
        index = seq_along(y)
    )
}

# The two steps of a regression pair on the rows both conditions' chi-square
# filters keep: each condition's fit (its intercept held where a GWAS side's
# is supplied), then the cross-products', the shared noise free, held or zero;
# a window fit that cannot be solved falls back to the supplied noise
# correlation where one is given. The fit carries what the standard errors
# read; for the jackknife, the rows themselves.
# @noRd
.rgRegressionPair <- function(param, ud, rows, cross) {
    proc <- param@procedures
    keep <- .h2Chi2Kept(proc, rows[[1L]]$chi2, rows[[1L]]$n) &
        .h2Chi2Kept(proc, rows[[2L]]$chi2, rows[[2L]]$n)
    kept <- map(rows, .h2SubsetRows, keep = keep)
    keptCross <- .rgSubsetCross(cross, keep)
    steps <- map2(kept, ud$heldIntercept, .rgRegressionStep1, param = param)
    fit <- .rgCrossRegressionFit(param, ud, kept, keptCross, steps)
    list_assign(fit, rows = kept, cross = keptCross, steps = steps, ud = ud)
}

# Rows `keep` of the cross-products' inputs.
# @noRd
.rgSubsetCross <- function(cross, keep) {
    list_modify(
        cross,
        y = cross$y[keep],
        X = cross$X[keep, , drop = FALSE],
        taper = cross$taper[keep],
        redundancy = cross$redundancy[keep]
    )
}

# Step 1 for one condition: the settings' own regression fit, or with the
# intercept held at a GWAS side's genome-wide value.
# @noRd
.rgRegressionStep1 <- function(rows, heldIntercept, param) {
    if (is.na(heldIntercept)) {
        return(.h2RegressionFit(param, rows))
    }
    .h2SolveRegression(param, rows, heldIntercept)
}

# Step 2: the cross-products' weighted fit, the shared noise as the settings
# hold it, or its fallback where a window fit cannot be solved.
# @noRd
.rgCrossRegressionFit <- function(param, ud, rows, cross, steps) {
    mode <- param@crossNoise
    intercepts <- map_dbl(steps, .rgStepNoise)
    held <- .rgHeldCross(param, ud, .rgRootProduct(intercepts))
    fit <- try_fetch(
        .rgCrossIrwls(param, rows, cross, steps, held),
        error = .rgSolveFailed
    )
    if (!is.null(fit) || !identical(mode, "window")) {
        return(list_assign(fit %||% list(), mode = mode))
    }
    if (is.null(ud$noiseCorrelation)) {
        return(list(mode = mode))
    }
    fallback <- ud$noiseCorrelation[[1L, 2L]] * .rgRootProduct(intercepts)
    list_assign(
        .rgCrossIrwls(param, rows, cross, steps, fallback),
        mode = "fallback"
    )
}

# A condition's noise term at its step-1 fit: the intercept, or 1 - h2 for
# tied noise.
# @noRd
.rgStepNoise <- function(step) {
    step$coef[[2L]]
}

# @noRd
.rgRootProduct <- function(x) {
    sqrt(max(x[[1L]], 0) * max(x[[2L]], 0))
}

# A cross-product fit that cannot be solved has no estimate.
# @noRd
.rgSolveFailed <- function(cnd) {
    NULL
}

# The cross-products' weighted fit by the regressions' engine
# (.h2SolveRegression()), each weight the inverse of a cross-product's
# variance at the current fit (.h2CrossWeights()), the first weights flat or
# at no covariance. The shared-noise column is fitted unless `held` holds it
# (NA: fitted), tied noise taking the covariance off it.
# @noRd
.rgCrossIrwls <- function(param, rows, cross, steps, held) {
    proc <- param@procedures
    mu <- map2(rows, steps, .rgStepMean, proc = proc)
    crossRows <- list_assign(
        cross,
        n = rows[[1L]]$n,
        weights = "cross",
        mu11 = mu[[1L]],
        mu22 = mu[[2L]]
    )
    window <- is.na(held)
    fit <- .h2SolveRegression(
        param,
        crossRows,
        if (window) NULL else held,
        tied = !window && identical(param@crossNoise, "tied")
    )
    list_assign(fit, mu = mu)
}

# A condition's modelled mean at its step-1 fit, as its weights read it: h2
# within [0, 1] and the noise term no less than 0 (or 1 where the weights hold
# it), floored above 0.
# @noRd
.rgStepMean <- function(rows, step, proc) {
    h2 <- min(max(step$coef[[1L]], 0), 1)
    noise <- if (proc@fitWeightIntercept) max(step$coef[[2L]], 0) else 1
    pmax(h2 * rows$X[, 1L] + noise * rows$X[, 2L], 1e-12)
}

# A regression pair's values: the covariance and the shared noise it fits,
# each condition's heritability from step 1, and the standard errors the
# settings' inference gives (.rgRegressionSe()), with rg's 95% Wald interval
# and the Wald test of no covariance.
# @noRd
.rgRegressionResult <- function(param, fit, nSnps, nComponents) {
    h2 <- map_dbl(fit$steps, .rgStepH2)
    base <- list(nSnps = nSnps, nComponents = nComponents, h2 = h2)
    if (is.null(fit$coef)) {
        return(exec(.rgPairResult, !!!base))
    }
    covariance <- fit$coef[[1L]]
    rg <- .rgCorrelation(covariance, h2)
    se <- .rgRegressionSe(param, fit, rg)
    root <- .rgRootProduct(map_dbl(fit$steps, .rgStepNoise))
    exec(
        .rgPairResult,
        !!!base,
        rg = rg,
        rgSe = se$rg,
        rgLower = rg - qnorm(0.975) * se$rg,
        rgUpper = rg + qnorm(0.975) * se$rg,
        covariance = covariance,
        covarianceSe = se$covariance,
        noiseCorrelation = .rgReportedNoise(
            fit$mode,
            fit$coef[[2L]],
            root,
            fit$ud
        ),
        noiseSource = fit$mode,
        p = .rgWaldP(covariance, se$covariance)
    )
}

# A condition's heritability at its step-1 fit.
# @noRd
.rgStepH2 <- function(step) {
    step$coef[[1L]]
}

# The two-sided Wald test of no covariance; NA without a standard error.
# @noRd
#' @importFrom stats pnorm qnorm
.rgWaldP <- function(estimate, se) {
    if (is.na(se) || se <= 0) {
        return(NA_real_)
    }
    2 * pnorm(-abs(estimate / se))
}

# The covariance's and rg's standard errors under the settings' inference:
# the model's sandwich over the three regressions jointly, the jackknife over
# contiguous blocks of the window, or none.
# @noRd
.rgRegressionSe <- function(param, fit, rg) {
    switch(
        param@inference,
        information = .rgSandwichSe(fit, rg),
        jackknife = .rgRegressionJackknife(param, fit),
        list(covariance = NA_real_, rg = NA_real_)
    )
}

# The covariance's and rg's standard errors from the regressions' joint
# sandwich (.rgSandwichCovariance()); rg's is the delta method's.
# @noRd
.rgSandwichSe <- function(fit, rg) {
    .rgDeltaSe(
        .rgSandwichCovariance(fit),
        map_dbl(fit$steps, .rgStepH2),
        rg
    )
}

# The covariance's standard error from a sandwich (`cov`, with h1, h2 and the
# covariance at `at`), and rg's by the delta method.
# @noRd
.rgDeltaSe <- function(sandwich, h2, rg) {
    at <- sandwich$at
    covarianceSe <- sqrt(sandwich$cov[[at[[3L]], at[[3L]]]])
    if (is.na(rg)) {
        return(list(covariance = covarianceSe, rg = NA_real_))
    }
    grad <- c(-rg / (2 * h2), 1 / sqrt(prod(h2)))
    list(
        covariance = covarianceSe,
        rg = sqrt(as.vector(crossprod(grad, sandwich$cov[at, at] %*% grad)))
    )
}

# The sandwich covariance of the three regressions' coefficients (the two
# conditions' and the cross-products'), their estimating equations stacked:
# each regression's own bread, and a meat that takes the covariances of a
# row's squares and cross-product under the model, Gaussian moments of the
# modelled means (Var(y11) = 2 mu11^2, Var(y12) = mu11 mu22 + mu12^2,
# Cov(y11, y22) = 2 mu12^2, Cov(y11, y12) = 2 mu11 mu12) where the rows are
# independent directions, and of the statistics' covariances across every
# pair of rows where they are variants, which LD couples (.rgLdMeat()). `at`
# gives where h1, h2 and the covariance sit; each condition's block is its
# own fit's sandwich (.h2RegressionSe()).
# @noRd
.rgSandwichCovariance <- function(fit) {
    parts <- c(fit$steps, list(fit))
    X <- map(parts, c("design", "X"))
    w <- map(parts, "w")
    bread <- .rgBlockDiagonal(map2(X, w, .rgBread))
    meat <- if (is.null(fit$cross$ldModel)) {
        .rgMeat(X, w, .rgModelMeans(fit))
    } else {
        .rgLdMeat(map2(X, w, `*`), .rgModelS(fit))
    }
    list(
        cov = solve(bread, t(solve(bread, meat))),
        at = cumsum(c(1L, map_int(X[seq_len(2L)], ncol)))
    )
}

# The modelled means of a row's two squares and cross-product at the fits.
# @noRd
.rgModelMeans <- function(fit) {
    list(
        mu11 = pmax(
            as.vector(fit$rows[[1L]]$X %*% fit$steps[[1L]]$coef),
            1e-12
        ),
        mu22 = pmax(
            as.vector(fit$rows[[2L]]$X %*% fit$steps[[2L]]$coef),
            1e-12
        ),
        mu12 = as.vector(fit$cross$X %*% fit$coef)
    )
}

# @noRd
.rgBread <- function(X, w) {
    crossprod(X, X * w)
}

# The stacked meat: block (a, b) is X_a' diag(w_a w_b C_ab) X_b, C_ab the
# covariance of rows' responses a and b under the model.
# @noRd
.rgMeat <- function(X, w, mu) {
    C <- list(
        c(2 * mu$mu11^2, 2 * mu$mu12^2, 2 * mu$mu11 * mu$mu12),
        c(2 * mu$mu12^2, 2 * mu$mu22^2, 2 * mu$mu22 * mu$mu12),
        c(
            2 * mu$mu11 * mu$mu12,
            2 * mu$mu22 * mu$mu12,
            mu$mu11 * mu$mu22 + mu$mu12^2
        )
    )
    n <- length(mu$mu11)
    rowsOf <- map(seq_len(3L), .rgMeatRow, X = X, w = w, C = C, n = n)
    exec(rbind, !!!rowsOf)
}

# Row block `a` of the stacked meat.
# @noRd
.rgMeatRow <- function(a, X, w, C, n) {
    blocks <- map(
        seq_len(3L),
        .rgMeatBlock,
        a = a,
        X = X,
        w = w,
        C = C,
        n = n
    )
    exec(cbind, !!!blocks)
}

# @noRd
.rgMeatBlock <- function(b, a, X, w, C, n) {
    cab <- C[[a]][(b - 1L) * n + seq_len(n)]
    crossprod(X[[a]], X[[b]] * (w[[a]] * w[[b]] * cab))
}

# The stacked meat over a variant regression pair's kept variants: block
# (a, b) is (X_a w_a)' C_ab (X_b w_b), C_ab the covariances of the rows'
# responses a and b across every pair of variants, the Gaussian moments of
# the statistics' covariances S at the fits (.rgRowCovariance()). `Xw` holds
# each regression's weighted design.
# @noRd
.rgLdMeat <- function(Xw, S) {
    rowsOf <- map(seq_len(3L), .rgLdMeatRow, Xw = Xw, S = S)
    exec(rbind, !!!rowsOf)
}

# Row block `a` of .rgLdMeat().
# @noRd
.rgLdMeatRow <- function(a, Xw, S) {
    blocks <- map(seq_len(3L), .rgLdMeatBlock, a = a, Xw = Xw, S = S)
    exec(cbind, !!!blocks)
}

# Block (a, b) of .rgLdMeat(); below the diagonal, the transpose of (b, a).
# @noRd
.rgLdMeatBlock <- function(b, a, Xw, S) {
    if (a > b) {
        return(t(.rgLdMeatBlock(a, b, Xw, S)))
    }
    crossprod(Xw[[a]], .rgRowCovariance(a, b, S) %*% Xw[[b]])
}

# The covariance of responses a and b (1: condition 1's squares, 2: condition
# 2's, 3: the cross-products) between every pair of variants, from S11, S22
# and S12 = Cov(z1, z2'): Cov(z1j^2, z1k^2) = 2 S11_jk^2, Cov(z1j^2, z2k^2) =
# 2 S12_jk^2, Cov(z1j^2, z1k z2k) = 2 S11_jk S12_jk, Cov(z2j^2, z1k z2k) = 2
# S12_kj S22_jk and Cov(z1j z2j, z1k z2k) = S11_jk S22_jk + S12_jk S12_kj.
# @noRd
.rgRowCovariance <- function(a, b, S) {
    switch(
        str_c(a, b),
        "11" = 2 * S$s11^2,
        "12" = 2 * S$s12^2,
        "13" = 2 * S$s11 * S$s12,
        "22" = 2 * S$s22^2,
        "23" = 2 * t(S$s12) * S$s22,
        "33" = S$s11 * S$s22 + S$s12 * t(S$s12)
    )
}

# A variant regression pair's statistics' covariances under the model at its
# fits, over its kept variants: each condition's (n_k h_k / M) R_k W_k R_k +
# c_k R_k from its step-1 fit, and S12 = (sqrt(n1 n2) g / M) R1 W12 R2 + c12 N
# from the cross-products' (covariance, shared noise), N the reference's R,
# or the two references' mean. The heritabilities and noise terms are held
# at no less than 0, the covariance and the shared noise within the roots of
# their products, so that the three are one Gaussian's on one reference.
# @noRd
.rgModelS <- function(fit) {
    model <- fit$cross$ldModel
    steps <- map(fit$steps, .rgStepModelCoef)
    root <- sqrt(steps[[1L]] * steps[[2L]])
    cross <- pmin(pmax(fit$coef, -root), root)
    idx <- fit$rows[[1L]]$index
    S <- list(
        s11 = .rgOwnModelS(1L, model, steps[[1L]]),
        s22 = .rgOwnModelS(2L, model, steps[[2L]]),
        s12 = .rgCrossModelS(model, cross)
    )
    map(S, .rgSubsetSquare, idx = idx)
}

# A step-1 fit's (h2, noise), each held at no less than 0.
# @noRd
.rgStepModelCoef <- function(step) {
    pmax(step$coef[seq_len(2L)], 0)
}

# Condition `k`'s statistics' covariance (n_k h / M) R W R + c R.
# @noRd
.rgOwnModelS <- function(k, model, coef) {
    .h2GlsS(
        model$lds[[k]],
        model$weights[[k]] %||% rep(1, model$M),
        model$n[[k]] * coef[[1L]] / model$M,
        coef[[2L]]
    )
}

# The statistics' cross-covariance (sqrt(n1 n2) g / M) R1 W12 R2 + c12 N, `coef`
# the cross-products' (g, c12): through .h2GlsS() on one reference, formed
# from the two references' matrices on two.
# @noRd
.rgCrossModelS <- function(model, coef) {
    scale <- sqrt(model$n[[1L]] * model$n[[2L]]) * coef[[1L]] / model$M
    shared <- rep_len(model$shared, model$M)
    if (identical(model$lds[[1L]], model$lds[[2L]])) {
        return(.h2GlsS(model$lds[[1L]], shared, scale, coef[[2L]]))
    }
    R1 <- .ldSourceMatrix(model$lds[[1L]])
    R2 <- .ldSourceMatrix(model$lds[[2L]])
    scale * (R1 %*% (shared * R2)) + coef[[2L]] * (R1 + R2) / 2
}

# Rows and columns `idx` of a square matrix.
# @noRd
.rgSubsetSquare <- function(S, idx) {
    S[idx, idx, drop = FALSE]
}

# A block-diagonal matrix of `blocks`.
# @noRd
.rgBlockDiagonal <- function(blocks) {
    as.matrix(Matrix::bdiag(blocks))
}

# The jackknife over contiguous blocks of the window (the procedures'
# jackknifeBlocks): both steps refitted without each block in turn.
# @noRd
.rgRegressionJackknife <- function(param, fit) {
    nBlocks <- param@procedures@jackknifeBlocks
    if (is.null(nBlocks)) {
        abort(str_c(
            "In one window, a jackknife needs `jackknifeBlocks` ",
            "(H2RegressionParam()): there are no LD blocks to leave out."
        ))
    }
    nRows <- length(fit$cross$y)
    blockOf <- cut(seq_len(nRows), min(nBlocks, nRows), labels = FALSE)
    loo <- map(
        unique(blockOf),
        .rgRegressionLoo,
        blockOf = blockOf,
        param = param,
        fit = fit
    )
    rgLoo <- map_dbl(loo, "rg")
    list(
        covariance = .jackknifeSeFromLoo(map_dbl(loo, "covariance")),
        rg = if (anyNA(rgLoo)) NA_real_ else .jackknifeSeFromLoo(rgLoo)
    )
}

# The pair refitted without block `b`.
# @noRd
.rgRegressionLoo <- function(b, blockOf, param, fit) {
    keep <- blockOf != b
    rows <- map(fit$rows, .h2SubsetRows, keep = keep)
    steps <- map2(
        rows,
        fit$ud$heldIntercept,
        .rgRegressionStep1,
        param = param
    )
    cross <- .rgCrossRegressionFit(
        param,
        fit$ud,
        rows,
        .rgSubsetCross(fit$cross, keep),
        steps
    )
    covariance <- if (is.null(cross$coef)) NA_real_ else cross$coef[[1L]]
    list(
        covariance = covariance,
        rg = .rgCorrelation(covariance, map_dbl(steps, .rgStepH2))
    )
}

# ---- the variant GLS with the cross-trait weight -----------------------------

# A pair's fit by the variant GLS: each condition's n r^2 - 1 on (1, l) as
# estimateLocalH2()'s variant GLS takes them, weighted by S o S with S = t R^2
# + R and t its first pass, then the cross-products sqrt(n1 n2) r1 r2 on (1,
# l), weighted by their covariance under the model, S11 o S22 + S12 o S12,
# with S12 = t12 R^2 + c12 R at the cross-products' own first pass. The
# standard errors come from the Gaussian covariances of the three responses
# through each fit's GLS map.
# @noRd
.rgVariantGls <- function(param, ud) {
    M <- length(ud$z[[1L]])
    if (M < 3L) {
        return(list(.rgPairResult(nSnps = M)))
    }
    g <- .rgGlsData(param, ud)
    steps <- map(seq_len(2L), .rgGlsStep1, g = g, held = ud$heldIntercept)
    base <- list(nSnps = M, h2 = map_dbl(steps, "h2"))
    cross <- .rgGlsCrossFit(param, ud, g, steps)
    if (is.null(cross$coef)) {
        return(list(exec(.rgPairResult, !!!base)))
    }
    list(.rgGlsResult(param, ud, g, steps, cross, base))
}

# Step 2 under the settings' shared noise, or under the supplied noise
# correlation where a window fit cannot be solved; without a coefficient
# when neither can.
# @noRd
.rgGlsCrossFit <- function(param, ud, g, steps) {
    mode <- param@crossNoise
    root <- .rgRootProduct(map_dbl(steps, "intercept"))
    held <- .rgHeldCross(param, ud, root)
    fit <- try_fetch(.rgGlsCross(param, g, steps, held), error = .rgSolveFailed)
    if (!is.null(fit) || !identical(mode, "window")) {
        return(list_assign(fit %||% list(), mode = mode))
    }
    if (is.null(ud$noiseCorrelation)) {
        return(list(mode = mode))
    }
    fallback <- ud$noiseCorrelation[[1L, 2L]] * root
    list_assign(.rgGlsCross(param, g, steps, fallback), mode = "fallback")
}

# What the pair's GLS fits read: each condition's response and the
# cross-products', the window's LD scores, R and R^2 (once, for the weights'
# S matrices), and the sample sizes.
# @noRd
.rgGlsData <- function(param, ud) {
    M <- length(ud$z[[1L]])
    n <- map_dbl(ud$nVariant, .h2SampleSize, how = param@sampleSize)
    r <- map2(ud$z, ud$nVariant, .zToCor)
    weights <- .h2SourceWeights(param, ud$ld) %||% rep(1, M)
    R <- .h2GlsS(ud$ld, rep(1, M), 0)
    list(
        y = list(n[[1L]] * r[[1L]]^2 - 1, n[[2L]] * r[[2L]]^2 - 1),
        y12 = sqrt(n[[1L]] * n[[2L]]) * r[[1L]] * r[[2L]],
        l = .h2WindowLdScores(param, ud$ld, M),
        R = R,
        R2 = .h2GlsS(ud$ld, weights, 1) - R,
        weights = weights,
        n = n,
        M = M,
        ld = ud$ld,
        proc = param@procedures
    )
}

# Step 1 for one condition, as .h2VariantGls() fits it: the weight's first
# pass t = sum(y) / sum(l), Omega = S o S cut at the procedures' share of its
# eigenvalue mass, then the GLS on (1, l); with a GWAS side's intercept held,
# on l alone. The S the cross-products' weight and the standard errors read
# takes the first pass at no less than 0, so that it is a covariance.
# @noRd
.rgGlsStep1 <- function(k, g, held) {
    free <- is.na(held[[k]])
    y <- g$y[[k]] - (if (free) 0 else held[[k]] - 1)
    t <- sum(y) / sum(g$l)
    om <- .h2GlsOmegaEigen(
        g$ld,
        matrix(g$weights, g$M, 1L),
        t,
        cutpoint = g$proc@omegaFraction
    )
    fit <- .h2GlsFit(
        y,
        if (free) cbind(1, g$l) else cbind(g$l),
        om,
        g$proc@ridge
    )
    slope <- fit$coef[[if (free) 2L else 1L]]
    list(
        h2 = slope * g$M / g$n[[k]],
        intercept = if (free) fit$coef[[1L]] + 1 else held[[k]],
        slope = slope,
        slopeAt = if (free) 2L else 1L,
        A = fit$A,
        t = max(t, 0),
        S = max(t, 0) * g$R2 + g$R
    )
}

# Step 2: the cross-products' GLS, the shared noise fitted on (1, l) or held
# (`held`; NA: fitted), weighted by S11 o S22 + S12 o S12' at their first pass
# (S12 is symmetric on one reference; on two, `g` carries the cross-population
# LD scores as l, R1 R2 as R2 and (R1 + R2) / 2 as R)
# (ordinary least squares for a fitted shared noise). Tied noise holds Ph12
# less the covariance, slope M / sqrt(n1 n2), so the slope's column is l less
# that. The first pass is held where the two conditions' S matrices make a
# joint covariance with S12, its genetic term within the root of the product
# of theirs and its noise term within [-1, 1] (theirs being 1), so that the
# weight is one: in one window the first pass of the shared noise can lie far
# outside it.
# @noRd
.rgGlsCross <- function(param, g, steps, held) {
    toNoise <- identical(param@crossNoise, "tied") * g$M / sqrt(prod(g$n))
    slopeColumn <- g$l - toNoise
    first <- if (is.na(held)) .rgOlsCoef(g$y12, g$l) else c(held, NA_real_)
    t12 <- if (is.na(held)) {
        first[[2L]]
    } else {
        sum(g$y12 - held) / sum(slopeColumn)
    }
    bound <- sqrt(steps[[1L]]$t * steps[[2L]]$t)
    S12 <- min(max(t12, -bound), bound) *
        g$R2 +
        min(max(first[[1L]] - toNoise * t12, -1), 1) * g$R
    om <- .h2TruncatedEigen(
        steps[[1L]]$S * steps[[2L]]$S + S12 * t(S12),
        g$proc@omegaFraction
    )
    y <- if (is.na(held)) g$y12 else g$y12 - held
    X <- if (is.na(held)) cbind(1, g$l) else cbind(slopeColumn)
    fit <- .h2GlsFit(y, X, om, g$proc@ridge)
    list(
        coef = if (is.na(held)) {
            fit$coef
        } else {
            c(held - toNoise * fit$coef, fit$coef)
        },
        A = fit$A,
        slopeAt = if (is.na(held)) 2L else 1L,
        S12 = S12
    )
}

# The ordinary least-squares coefficients of `y` on (1, l).
# @noRd
.rgOlsCoef <- function(y, l) {
    as.vector(coef(lm(y ~ l)))
}

# A GLS pair's values: the covariance on the heritability's scale, the shared
# noise as a correlation, and under the information the standard errors from
# the three fits' joint covariance (.rgGlsCovariance()).
# @noRd
.rgGlsResult <- function(param, ud, g, steps, cross, base) {
    nn <- sqrt(g$n[[1L]] * g$n[[2L]])
    covariance <- cross$coef[[2L]] * g$M / nn
    rg <- .rgCorrelation(covariance, base$h2)
    se <- if (identical(param@inference, "information")) {
        .rgGlsSe(g, steps, cross, rg, g$M / nn)
    } else {
        list(covariance = NA_real_, rg = NA_real_)
    }
    root <- .rgRootProduct(map_dbl(steps, "intercept"))
    exec(
        .rgPairResult,
        !!!base,
        rg = rg,
        rgSe = se$rg,
        rgLower = rg - qnorm(0.975) * se$rg,
        rgUpper = rg + qnorm(0.975) * se$rg,
        covariance = covariance,
        covarianceSe = se$covariance,
        noiseCorrelation = .rgReportedNoise(
            cross$mode,
            cross$coef[[1L]],
            root,
            ud
        ),
        noiseSource = cross$mode,
        p = .rgWaldP(covariance, se$covariance)
    )
}

# The covariance's and rg's standard errors from the three GLS fits' joint
# covariance, A_a Cov(y_a, y_b) A_b', with the responses' Gaussian
# covariances under the weights' S matrices: Var(y11) = 2 S11 o S11,
# Cov(y11, y22) = 2 S12 o S12, Cov(y11, y12) = 2 S11 o S12, Cov(y22, y12) = 2
# S22 o S12', Var(y12) = S11 o S22 + S12 o S12', S12 = Cov(z1, z2) (symmetric
# on one reference, not on two). rg = s12 / sqrt(s1 s2) of the three slopes.
# @noRd
.rgGlsSe <- function(g, steps, cross, rg, scale) {
    S <- list(steps[[1L]]$S, steps[[2L]]$S, cross$S12)
    A <- list(steps[[1L]]$A, steps[[2L]]$A, cross$A)
    at <- c(steps[[1L]]$slopeAt, steps[[2L]]$slopeAt, cross$slopeAt)
    cov <- .rgGlsSlopeCovariance(A, S, at)
    covarianceSe <- sqrt(cov[[3L, 3L]]) * scale
    if (is.na(rg)) {
        return(list(covariance = covarianceSe, rg = NA_real_))
    }
    s <- c(steps[[1L]]$slope, steps[[2L]]$slope)
    grad <- c(-rg / (2 * s), 1 / sqrt(prod(s)))
    list(
        covariance = covarianceSe,
        rg = sqrt(as.vector(crossprod(grad, cov %*% grad)))
    )
}

# The 3 x 3 covariance of the three fits' slopes.
# @noRd
#' @importFrom purrr pmap_dbl
.rgGlsSlopeCovariance <- function(A, S, at) {
    C <- list(
        list(
            2 * S[[1L]] * S[[1L]],
            2 * S[[3L]] * S[[3L]],
            2 * S[[1L]] * S[[3L]]
        ),
        list(NULL, 2 * S[[2L]] * S[[2L]], 2 * S[[2L]] * t(S[[3L]])),
        list(NULL, NULL, S[[1L]] * S[[2L]] + S[[3L]] * t(S[[3L]]))
    )
    pairs <- expand.grid(a = seq_len(3L), b = seq_len(3L))
    upper <- pairs[pairs$a <= pairs$b, ]
    values <- pmap_dbl(upper, .rgGlsSlopeCov, A = A, C = C, at = at)
    entry <- match(
        pmin(pairs$a, pairs$b) * 10L + pmax(pairs$a, pairs$b),
        upper$a * 10L + upper$b
    )
    matrix(values[entry], 3L, 3L)
}

# @noRd
.rgGlsSlopeCov <- function(a, b, A, C, at) {
    as.vector(A[[a]][at[[a]], ] %*% C[[a]][[b]] %*% A[[b]][at[[b]], ])
}

# ---- the moments: LAVA's Omega ----------------------------------------------

# A group's fit by the moments, as LAVA makes it: every condition's whitened
# components delta_k = Lambda^-1/2 Q'r_k on one K, and Omega = Delta'Delta / K
# - Sigma. Sigma's diagonal is each condition's noise per component tied to
# its heritability, (1 - h2 + h2 phi) / (n - 1) with the overlap factor phi
# (1 without it), so K Omega_kk is the heritability estimateLocalH2() gives;
# its off-diagonal is the supplied noise correlation times the root of the
# two diagonal terms, or 0. K Omega_kl is the covariance. Each pair's p-value
# is LAVA's simulation (.rgMomentsP()), and its interval the 2.5% and 97.5%
# quantiles of rg over noncentral Wishart draws at the estimate
# (.rgMomentsInterval()). That is the fixed-effects reading; under random
# effects, .rgRandomMoments().
# @noRd
.rgEigenMoments <- function(param, ud) {
    if (identical(param@effects, "random")) {
        return(.rgRandomMoments(param, ud))
    }
    d <- .rgEigenData(param, ud)
    pairs <- .rgGroupPairs(length(ud$z))
    if (is.null(d) || any(d$n <= length(d$lam) + 1L)) {
        return(pmap(pairs, .rgPairResult, nSnps = length(ud$z[[1L]])))
    }
    K <- length(d$lam)
    whitened <- map(d$u, .rgWhiten, lam = d$lam)
    delta <- exec(cbind, !!!whitened)
    Sigma <- .rgMomentsSigma(param, ud, d, delta)
    Omega <- crossprod(delta) / K - Sigma
    if (isTRUE(ud$conditional)) {
        return(pmap(
            pairs,
            .rgMomentsConditionalPair,
            Omega = Omega,
            Sigma = Sigma,
            K = K,
            M = d$M,
            param = param,
            draws = .rgMomentsDraws(Omega, Sigma, K)
        ))
    }
    pmap(
        pairs,
        .rgMomentsPair,
        Omega = Omega,
        Sigma = Sigma,
        K = K,
        M = d$M,
        param = param
    )
}

# Every pair (i, j), i < j, of a group's `P` conditions.
# @noRd
#' @importFrom utils combn
.rgGroupPairs <- function(P) {
    ij <- combn(P, 2L)
    list(i = ij[1L, ], j = ij[2L, ])
}

# Whether a group's conditions are fitted at once: under the moments, and
# under the likelihood with jointFit; otherwise pair by pair.
# @noRd
.rgFitsAtOnce <- function(param) {
    identical(param@criterion, "moments") ||
        (identical(param@criterion, "likelihood") &&
            isTRUE(param@procedures@jointFit))
}

# The rotated correlations whitened by the eigenvalues: LAVA's delta.
# @noRd
.rgWhiten <- function(u, lam) {
    u / sqrt(lam)
}

# Sigma: each condition's noise per component on the diagonal, tied to the
# heritability its own moments give (.h2EigenMoments()'s, with the overlap
# factor where the settings take it), and the noise the conditions share off
# it: the supplied correlation times the root of the two diagonal terms,
# none, or tied (.rgMomentsTiedNoise()).
# The moments on one reference under random effects: the group on the
# reference's eigenpairs at the smallest sample size (one K, as LAVA keeps),
# read as on several references (.rgMomentsGroupOf()), so that the estimates
# are LAVA's; each pair's p-value and interval, or each pair's given the
# group's others, from the components' Gaussian distribution with their
# genetic variance.
# @noRd
.rgRandomMoments <- function(param, ud) {
    P <- length(ud$z)
    M <- length(ud$z[[1L]])
    pairs <- .rgGroupPairs(P)
    one <- list_assign(ud, lds = list(ud$ld), ldOf = rep(1L, P))
    basis <- .rgBasisEigen(1L, param, one)
    if (is.null(basis)) {
        return(pmap(pairs, .rgPairResult, nSnps = M))
    }
    sides <- map(
        seq_len(P),
        .rgMemberSide,
        bases = list(basis),
        ud = one,
        param = param
    )
    if (any(map_dbl(sides, "n") <= length(basis$lam) + 1L)) {
        return(pmap(pairs, .rgPairResult, nSnps = M))
    }
    group <- .rgMomentsGroupOf(param, one, sides, list(basis))
    pairFit <- if (isTRUE(ud$conditional)) {
        .rgMixedConditionalPair
    } else {
        .rgGroupPair
    }
    pmap(pairs, pairFit, group = group, M = M)
}

# @noRd
.rgMomentsSigma <- function(param, ud, d, delta) {
    K <- length(d$lam)
    phi <- if (identical(param@referenceNoise, "overlapFactor")) {
        .h2OverlapFactor(d$n, ud$ld$nSamples, ud$overlap)
    } else {
        rep(1, length(d$n))
    }
    noise <- .rgMomentsOwn(
        colSums(delta^2),
        K,
        d$n,
        phi,
        .rgMomentsIntercepts(ud)
    )$noise
    shared <- switch(
        param@crossNoise,
        supplied = ud$noiseCorrelation * sqrt(tcrossprod(noise)),
        tied = .rgMomentsTiedNoise(ud, crossprod(delta), d$n, phi, K),
        diag(noise, length(noise))
    )
    shared - diag(diag(shared), length(noise)) + diag(noise, length(noise))
}

# Each condition's heritability K Omega_kk from its own components' sum of
# squares `dtd` on K components, and its noise per component tied to it,
# (c + h2 (phi - 1)) / (n - 1), LAVA's (1 - h2 + h2 phi) / (n - 1) at c = 1:
# the null intercept c (`intercept`: 1, or a GWAS side's genome-wide
# intercept where it is held) and the overlap factor phi (1 without it).
# @noRd
.rgMomentsOwn <- function(dtd, K, n, phi, intercept) {
    m <- K / (n - 1)
    h2 <- (dtd - m * intercept) / (1 + m * (phi - 1))
    list(h2 = h2, noise = pmax((intercept + h2 * (phi - 1)) / (n - 1), 0))
}

# Each condition's null intercept for the moments: a GWAS side's genome-wide
# intercept where `gwasH2Estimate` holds it, 1 otherwise.
# @noRd
.rgMomentsIntercepts <- function(ud) {
    held <- ud$heldIntercept %||% rep(NA_real_, length(ud$z))
    replace(held, is.na(held), 1)
}

# Tied shared noise: per component, the phenotypes' covariance less the part
# of the covariance the components carry, (Ph_kl - (1 - phi) G_kl) / (n - 1),
# the off-diagonal counterpart of the diagonal's (1 - (1 - phi) h2) / (n -
# 1). G_kl = Delta'Delta_kl - K Sigma_kl makes it linear in Sigma_kl, so
# Sigma_kl = (Ph_kl - f Delta'Delta_kl) / (a - f K), with a and f the pair's
# n - 1 and 1 - phi (alike on shared donors), taken as geometric means.
# @noRd
.rgMomentsTiedNoise <- function(ud, dtd, n, phi, K) {
    a <- sqrt(tcrossprod(n - 1))
    f <- 1 - sqrt(tcrossprod(phi))
    (ud$phenotypicCorrelation - f * dtd) / (a - f * K)
}

# One pair's values from the group's Omega: the covariance and the two
# heritabilities on the heritability's scale, rg, LAVA's simulation p-value
# and interval, and the noise correlation Sigma holds.
# @noRd
.rgMomentsPair <- function(i, j, Omega, Sigma, K, M, param) {
    at <- c(i, j)
    omega <- Omega[at, at]
    sigma <- Sigma[at, at]
    covariance <- K * omega[[1L, 2L]]
    h2 <- K * diag(omega)
    interval <- .rgMomentsInterval(omega, sigma, K)
    .rgPairResult(
        i = i,
        j = j,
        nSnps = M,
        nComponents = K,
        rg = .rgCorrelation(covariance, h2),
        rgLower = interval[[1L]],
        rgUpper = interval[[2L]],
        covariance = covariance,
        h2 = h2,
        noiseCorrelation = if (prod(diag(sigma)) > 0) {
            sigma[[1L, 2L]] / sqrt(prod(diag(sigma)))
        } else {
            NA_real_
        },
        noiseSource = param@crossNoise,
        p = .rgMomentsP(omega, sigma, K)
    )
}

# LAVA's simulation p-value for no covariance (Werme et al. 2022, Methods):
# under the null Omega_0, the observed variances and no covariance, the
# components' cross-products are noncentral Wishart, K degrees of freedom,
# scale Sigma, noncentrality K Omega_0. Partial integration: given one
# condition's components, the covariance statistic omega_12 is normal, so each
# draw of that condition contributes its conditional p-value, and the p-value
# is their mean, half the draws conditioning on each condition. 10,000 draws,
# then 100,000 below 1e-4 and 1,000,000 below 1e-6. NA unless both variances
# are positive and the noise is not perfectly correlated.
# @noRd
.rgMomentsP <- function(omega, sigma, K) {
    if (any(diag(omega) <= 0)) {
        return(NA_real_)
    }
    if (prod(diag(sigma)) - sigma[[1L, 2L]]^2 <= 0) {
        return(NA_real_)
    }
    p <- .rgMomentsSimulatedP(omega, sigma, K, 1e4)
    if (p >= 1e-4) {
        return(p)
    }
    p100k <- .rgMomentsSimulatedP(omega, sigma, K, 1e5)
    if (p100k >= 1e-6) {
        return(p100k)
    }
    .rgMomentsSimulatedP(omega, sigma, K, 1e6)
}

# The mean conditional p-value over `nDraws` draws, half conditioning on each
# condition.
# @noRd
.rgMomentsSimulatedP <- function(omega, sigma, K, nDraws) {
    half <- nDraws / 2
    mean(c(
        .rgConditionalP(1L, omega, sigma, K, half),
        .rgConditionalP(2L, omega, sigma, K, half)
    ))
}

# The conditional p-values given condition `given`'s components. With the null
# means of the two conditions' components on two orthogonal axes, lengths
# sqrt(K omega_kk), a draw of the given condition's components needs its two
# coordinates on those axes and the squared length of the rest; the other
# condition's components given them are normal with mean m + beta (e - m) and
# variance tau^2 per component, beta = sigma_12 / sigma_gg and tau^2 =
# sigma_oo - sigma_12^2 / sigma_gg, so the statistic delta_o' delta_g / K -
# sigma_12 is normal, and its two-sided tail beyond the observed omega_12 is
# the draw's p-value.
# @noRd
#' @importFrom stats rnorm rchisq pnorm
.rgConditionalP <- function(given, omega, sigma, K, nDraws) {
    other <- 3L - given
    sg <- sigma[[given, given]]
    s12 <- sigma[[1L, 2L]]
    own <- rnorm(nDraws, sd = sqrt(sg))
    cross <- rnorm(nDraws, sd = sqrt(sg))
    rest <- sg * rchisq(nDraws, df = K - 2L)
    meanGiven <- sqrt(K * omega[[given, given]])
    meanOther <- sqrt(K * omega[[other, other]])
    beta <- s12 / sg
    tau2 <- sigma[[other, other]] - s12^2 / sg
    noiseDot <- meanGiven * own + own^2 + cross^2 + rest
    mu <- (meanOther * cross + beta * noiseDot) / K - s12
    sd <- sqrt(tau2 * ((meanGiven + own)^2 + cross^2 + rest)) / K
    observed <- abs(omega[[1L, 2L]])
    pnorm(-observed, mu, sd) + pnorm(observed, mu, sd, lower.tail = FALSE)
}

# LAVA's 95% interval for rg: 10,000 draws of the components' cross-products,
# noncentral Wishart with K degrees of freedom, scale Sigma and noncentrality
# K Omega at the estimate (matrixsampling::rwishart(), whose noncentrality is
# M'M), each turned into Omega and rg; the 2.5% and 97.5% quantiles of the
# draws' rg. NA unless the estimate is a covariance matrix, or without
# matrixsampling.
# @noRd
#' @importFrom stats quantile
.rgMomentsInterval <- function(omega, sigma, K) {
    none <- c(NA_real_, NA_real_)
    if (any(diag(omega) <= 0) || det(omega) < 0) {
        return(none)
    }
    if (!requireNamespace("matrixsampling", quietly = TRUE)) {
        return(none)
    }
    draws <- matrixsampling::rwishart(1e4, K, sigma, K * omega)
    rg <- .rgDrawCorrelation(draws, sigma, K)
    unname(quantile(rg, c(0.025, 0.975), na.rm = TRUE))
}

# Each Wishart draw's rg, NA where a draw's variance is not positive.
# @noRd
.rgDrawCorrelation <- function(draws, sigma, K) {
    o11 <- draws[1L, 1L, ] / K - sigma[[1L, 1L]]
    o22 <- draws[2L, 2L, ] / K - sigma[[2L, 2L]]
    o12 <- draws[1L, 2L, ] / K - sigma[[1L, 2L]]
    if_else(o11 > 0 & o22 > 0, o12 / sqrt(pmax(o11 * o22, 0)), NA_real_)
}

# ---- the moments given a unit's other conditions -----------------------------

# One pair of a group given its other conditions Z (conditional = TRUE), from
# the group's Omega, as LAVA's partial correlation and multiple regression
# read it (Werme et al. 2022): the partial covariance Omega_ij|Z and the two
# partial variances, on the heritability's scale (times K), and rg their
# correlation; each condition's r2, the share of its genetic variance that
# the group's other conditions explain together; LAVA's simulation p-value of
# no partial covariance (.rgPartialP()); the 2.5% and 97.5% quantiles of rg
# and the r2s over the noncentral Wishart draws at the estimate. With no
# other condition, the pair's own values, r2 its rg squared.
# @noRd
.rgMomentsConditionalPair <- function(i, j, Omega, Sigma, K, M, param, draws) {
    rest <- setdiff(seq_len(nrow(Omega)), c(i, j))
    partial <- .rgPartial(Omega, i, j, rest)
    covariance <- K * partial$cov
    h2 <- K * partial$var
    drawn <- .rgConditionalDrawQuantiles(draws, Sigma, K, i, j, rest)
    sigma <- Sigma[c(i, j), c(i, j)]
    .rgPairResult(
        i = i,
        j = j,
        nSnps = M,
        nComponents = K,
        rg = .rgCorrelation(covariance, h2),
        rgLower = drawn$rg[[1L]],
        rgUpper = drawn$rg[[2L]],
        covariance = covariance,
        h2 = h2,
        noiseCorrelation = if (prod(diag(sigma)) > 0) {
            sigma[[1L, 2L]] / sqrt(prod(diag(sigma)))
        } else {
            NA_real_
        },
        noiseSource = param@crossNoise,
        p = .rgPartialP(Omega, Sigma, K, i, j, rest),
        r2 = c(.rgExplained(Omega, i), .rgExplained(Omega, j)),
        r2Lower = drawn$r2Lower,
        r2Upper = drawn$r2Upper
    )
}

# The covariance of conditions i and j and their two variances given the
# conditions `rest`, NA where rest's covariance is not positive definite (its
# inverse is then undefined, as LAVA refuses it).
# @noRd
.rgPartial <- function(Omega, i, j, rest) {
    ij <- c(i, j)
    if (length(rest) == 0L) {
        return(list(cov = Omega[[i, j]], var = diag(Omega)[ij]))
    }
    if (!.rgPositiveDefinite(Omega[rest, rest, drop = FALSE])) {
        return(list(cov = NA_real_, var = c(NA_real_, NA_real_)))
    }
    across <- Omega[rest, ij, drop = FALSE]
    left <- Omega[ij, ij] -
        crossprod(across, solve(Omega[rest, rest, drop = FALSE], across))
    list(cov = left[[1L, 2L]], var = diag(left))
}

# The share of condition k's genetic variance the others of its group explain
# together, the r2 of LAVA's multiple regression; NA where k's variance is not
# positive or the others' covariance is not positive definite.
# @noRd
.rgExplained <- function(Omega, k) {
    rest <- setdiff(seq_len(nrow(Omega)), k)
    if (
        Omega[[k, k]] <= 0 ||
            !.rgPositiveDefinite(Omega[rest, rest, drop = FALSE])
    ) {
        return(NA_real_)
    }
    across <- Omega[rest, k]
    sum(across * solve(Omega[rest, rest, drop = FALSE], across)) / Omega[[k, k]]
}

# Whether a symmetric matrix is positive definite, to its scale.
# @noRd
.rgPositiveDefinite <- function(x) {
    values <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    min(values) > sqrt(.Machine$double.eps) * max(abs(values), 1e-300)
}

# Whether a symmetric matrix is positive semi-definite, to its scale.
# @noRd
.rgPositiveSemidefinite <- function(x) {
    values <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    min(values) >= -sqrt(.Machine$double.eps) * max(abs(values), 1e-300)
}

# A factor L with L L' = x of a positive semi-definite matrix.
# @noRd
.rgPsdFactor <- function(x) {
    ei <- eigen(x, symmetric = TRUE)
    ei$vectors %*% diag(sqrt(pmax(ei$values, 0)), nrow(x))
}

# The group's 10,000 noncentral Wishart draws at the estimate (as
# .rgMomentsInterval() draws a pair's), or NULL where Omega is not a
# covariance matrix or matrixsampling is not installed.
# @noRd
.rgMomentsDraws <- function(Omega, Sigma, K) {
    if (!.rgPositiveSemidefinite(Omega) || any(diag(Omega) <= 0)) {
        return(NULL)
    }
    if (!requireNamespace("matrixsampling", quietly = TRUE)) {
        return(NULL)
    }
    matrixsampling::rwishart(1e4, K, Sigma, K * Omega)
}

# The 2.5% and 97.5% quantiles over the draws of the pair's rg given `rest`
# and of the two conditions' r2s; NA without draws.
# @noRd
#' @importFrom stats quantile
.rgConditionalDrawQuantiles <- function(draws, Sigma, K, i, j, rest) {
    none <- c(NA_real_, NA_real_)
    if (is.null(draws)) {
        return(list(rg = none, r2Lower = none, r2Upper = none))
    }
    stats <- map(
        seq_len(dim(draws)[[3L]]),
        .rgConditionalDraw,
        draws = draws,
        Sigma = Sigma,
        K = K,
        i = i,
        j = j,
        rest = rest
    )
    rg <- quantile(map_dbl(stats, 1L), c(0.025, 0.975), na.rm = TRUE)
    r2i <- quantile(map_dbl(stats, 2L), c(0.025, 0.975), na.rm = TRUE)
    r2j <- quantile(map_dbl(stats, 3L), c(0.025, 0.975), na.rm = TRUE)
    list(
        rg = unname(rg),
        r2Lower = unname(c(r2i[[1L]], r2j[[1L]])),
        r2Upper = unname(c(r2i[[2L]], r2j[[2L]]))
    )
}

# One draw's rg given `rest` and the two conditions' r2s.
# @noRd
.rgConditionalDraw <- function(d, draws, Sigma, K, i, j, rest) {
    omega <- draws[,, d] / K - Sigma
    partial <- .rgPartial(omega, i, j, rest)
    c(
        .rgCorrelation(partial$cov, partial$var),
        .rgExplained(omega, i),
        .rgExplained(omega, j)
    )
}

# LAVA's simulation p-value of no partial covariance of i and j given `rest`
# (Werme et al. 2022, Methods), written from the paper as .rgMomentsP() is:
# under the null Omega_0, the estimate with Omega_ij replaced by the value
# rest implies, the components' cross-products are noncentral Wishart; given
# the components of i and rest (or of j and rest), the partial covariance is
# normal, as it is linear in the other's components, so each draw gives its
# conditional p-value, and the p-value is their mean, half the draws
# conditioning each way, adaptively 10,000, 100,000 and 1,000,000 draws. NA
# where the partial covariance is undefined or Omega_0 is not a covariance
# matrix. With no other condition, the pair's own (.rgMomentsP()).
# @noRd
.rgPartialP <- function(Omega, Sigma, K, i, j, rest) {
    if (length(rest) == 0L) {
        return(.rgMomentsP(Omega[c(i, j), c(i, j)], Sigma[c(i, j), c(i, j)], K))
    }
    partial <- .rgPartial(Omega, i, j, rest)
    P <- nrow(Omega)
    at <- outer(seq_len(P) == i, seq_len(P) == j) |
        outer(seq_len(P) == j, seq_len(P) == i)
    null <- Omega - partial$cov * at
    if (is.na(partial$cov) || !.rgPositiveSemidefinite(null)) {
        return(NA_real_)
    }
    sim <- list(
        top = sqrt(K) * t(.rgPsdFactor(null)),
        Sigma = Sigma,
        K = K,
        observed = partial$cov
    )
    p <- .rgPartialSimulatedP(sim, i, j, rest, 1e4)
    if (p >= 1e-4) {
        return(p)
    }
    p100k <- .rgPartialSimulatedP(sim, i, j, rest, 1e5)
    if (p100k >= 1e-6) {
        return(p100k)
    }
    .rgPartialSimulatedP(sim, i, j, rest, 1e6)
}

# The mean conditional p-value over `nDraws` draws, half conditioning on each
# of the pair's two conditions (with rest).
# @noRd
.rgPartialSimulatedP <- function(sim, i, j, rest, nDraws) {
    mean(c(
        .rgPartialConditionalP(i, j, rest, sim, nDraws / 2),
        .rgPartialConditionalP(j, i, rest, sim, nDraws / 2)
    ))
}

# The conditional p-values given the components of condition `g` and `rest`
# (C), condition `f` free. The null's means lie in a span of as many
# dimensions as the group has conditions, the components' first rows here
# (`top`, M with M'M = K Omega_0), so a draw needs C's noise there and the
# cross-product of its noise on the remaining K - P rows, a central Wishart.
# f's components given C's are normal, mean M_f + (Delta_C - M_C) B and
# variance tau^2 per component, B = Sigma_CC^-1 Sigma_Cf; the partial
# covariance is w'Delta_f / K - (Sigma_gf - a'Sigma_rest,f), w = Delta_g -
# Delta_rest a with a = Omega_rest^-1 Omega_rest,g, all of C's.
# @noRd
#' @importFrom stats pnorm
.rgPartialConditionalP <- function(g, f, rest, sim, nDraws) {
    C <- c(g, rest)
    P <- ncol(sim$top)
    sigmaCC <- sim$Sigma[C, C, drop = FALSE]
    B <- solve(sigmaCC, sim$Sigma[C, f])
    tau2 <- sim$Sigma[[f, f]] - sum(sim$Sigma[f, C] * B)
    draw <- list(
        noise = .rgCorrelatedNormals(nDraws * P, sigmaCC),
        gram = .rgCentralGram(nDraws, sim$K - P, sigmaCC),
        C = C,
        f = f,
        g = g,
        rest = rest,
        B = B,
        tau2 = tau2
    )
    mean(
        map_dbl(seq_len(nDraws), .rgPartialDrawP, draw = draw, sim = sim),
        na.rm = TRUE
    )
}

# One draw's conditional p-value; NA where its rest's covariance is singular.
# @noRd
.rgPartialDrawP <- function(k, draw, sim, P = ncol(sim$top)) {
    rows <- (k - 1L) * P + seq_len(P)
    DC <- sim$top[, draw$C, drop = FALSE] + draw$noise[rows, , drop = FALSE]
    S <- crossprod(DC) + draw$gram[,, k]
    DM <- crossprod(DC, sim$top)
    omegaC <- S / sim$K - sim$Sigma[draw$C, draw$C, drop = FALSE]
    a <- try_fetch(
        solve(omegaC[-1L, -1L, drop = FALSE], omegaC[-1L, 1L]),
        error = .rgSolveFailed
    )
    if (is.null(a)) {
        return(NA_real_)
    }
    w <- c(1, -a)
    shared <- sim$Sigma[[draw$g, draw$f]] -
        sum(a * sim$Sigma[draw$rest, draw$f])
    mu <- (sum(w * DM[, draw$f]) +
        sum(w * ((S - DM[, draw$C, drop = FALSE]) %*% draw$B))) /
        sim$K -
        shared
    sd <- sqrt(draw$tau2 * sum(w * (S %*% w))) / sim$K
    observed <- abs(sim$observed)
    pnorm(-observed, mu, sd) + pnorm(observed, mu, sd, lower.tail = FALSE)
}

# `n` draws of rows with covariance `Sigma`.
# @noRd
#' @importFrom stats rnorm
.rgCorrelatedNormals <- function(n, Sigma) {
    matrix(rnorm(n * nrow(Sigma)), n) %*% chol(Sigma)
}

# `nDraws` central Wishart cross-products with `df` degrees of freedom and
# scale `Sigma`, as an array; from the rows themselves where df is below the
# dimension (stats::rWishart() needs it at least that), none where it is 0.
# @noRd
#' @importFrom stats rWishart
.rgCentralGram <- function(nDraws, df, Sigma) {
    q <- nrow(Sigma)
    if (df >= q) {
        return(rWishart(nDraws, df, Sigma))
    }
    if (df <= 0) {
        return(array(0, c(q, q, nDraws)))
    }
    rows <- .rgCorrelatedNormals(nDraws * df, Sigma)
    grams <- map(seq_len(nDraws), .rgBlockGram, rows = rows, df = df)
    array(unlist(grams), c(q, q, nDraws))
}

# @noRd
.rgBlockGram <- function(k, rows, df) {
    crossprod(rows[(k - 1L) * df + seq_len(df), , drop = FALSE])
}

# =============================================================================
# estimateLocalRg on QTL summary statistics
# =============================================================================

#' @rdname estimateLocalRg
#' @export
#' @importFrom BiocParallel SerialParam
setMethod(
    "estimateLocalRg",
    "QtlSumStats",
    function(
        data,
        gwas = NULL,
        param = NULL,
        jointSpecification = NULL,
        contexts = NULL,
        traitId = NULL,
        cisWindow = NULL,
        BPPARAM = SerialParam(),
        referenceOverlap = 0,
        noiseCorrelation = NULL,
        gwasH2Estimate = NULL,
        conditional = FALSE
    ) {
        settings <- .rgParam(
            param,
            gwas,
            noiseCorrelation,
            gwasH2Estimate,
            conditional
        )
        selected <- .rgSelectedRows(data, contexts, traitId, cisWindow)
        run <- .rgRun(gwas, noiseCorrelation, gwasH2Estimate, conditional)
        units <- .rgUnits(
            selected$data,
            selected$rows,
            jointSpecification,
            settings,
            run
        )
        fits <- .bpMap(
            units,
            .rgUnitFit,
            data = selected$data,
            param = settings,
            overlaps = .cisH2RowOverlaps(
                referenceOverlap,
                selected$rows,
                length(selected$data),
                "estimateLocalRg"
            ),
            run = run,
            BPPARAM = BPPARAM
        )
        .rgEstimate(fits, data, settings, run)
    }
)

# The QC'd collection's rows the selection keeps, on the collection as the
# selection subsets it.
# @noRd
.rgSelectedRows <- function(data, contexts, traitId, cisWindow) {
    .cisH2AssertQc(data, "estimateLocalRg")
    selection <- .qssSelection(
        data,
        contexts,
        traitId,
        cisWindow,
        "estimateLocalRg"
    )
    list(
        data = selection$data,
        rows = .cisH2SelectedRows(selection, "estimateLocalRg")
    )
}

# What every unit's fit reads besides its own data: the GWAS, the supplied
# noise correlations, each GWAS study's genome-wide intercept, and whether the
# pairs are reported given their unit's other conditions.
# @noRd
.rgRun <- function(gwas, noiseCorrelation, gwasH2Estimate, conditional) {
    list(
        gwas = gwas,
        noiseCorrelation = noiseCorrelation,
        gwasIntercepts = .rgGwasIntercepts(gwas, gwasH2Estimate),
        gwasLiability = .rgGwasLiability(gwas),
        conditional = conditional
    )
}

# Each GWAS study's factor to the liability scale (.h2LiabilityFactor()), NA
# for one that records no prevalence; NULL without a GWAS.
# @noRd
.rgGwasLiability <- function(gwas) {
    if (is.null(gwas)) {
        return(NULL)
    }
    studies <- unique(as.character(gwas$study))
    set_names(map_dbl(studies, .h2LiabilityFactor, data = gwas), studies)
}

# The settings the pairs are fitted under: `param`, or the provisional default
# (the eigen likelihood with the shared noise fitted in the window), refused
# where a pair cannot be fitted under them; `individual` for a QtlDataset's
# pairs.
# @noRd
.rgParam <- function(
    param,
    gwas,
    noiseCorrelation,
    gwasH2Estimate,
    conditional,
    individual = FALSE
) {
    .assertMethodParam(param %||% list(), "H2EstimateParam", "param")
    # An empty list carries no settings, as everywhere a Param is taken.
    settings <- if (length(param) == 0L) H2EstimateParam() else param
    .h2RefuseLocalSetting(settings, individual, "estimateLocalRg")
    if (!is.null(gwas)) {
        assertClass(gwas, "GwasSumStats", .var.name = "gwas")
    }
    if (!is.null(noiseCorrelation)) {
        .rgAssertNoiseCorrelation(noiseCorrelation)
    }
    if (!is.null(gwasH2Estimate)) {
        assertClass(gwasH2Estimate, "H2Estimate", .var.name = "gwasH2Estimate")
    }
    assertFlag(conditional, .var.name = "conditional")
    if (conditional && !.rgFitsAtOnce(settings)) {
        abort(str_c(
            "estimateLocalRg: conditional = TRUE reads each group's ",
            "conditions fitted at once: the moments, or the likelihood with ",
            "H2LikelihoodParam(jointFit = TRUE)."
        ))
    }
    problem <- .rgNoiseProblem(
        settings,
        gwas,
        noiseCorrelation,
        individual
    ) %||%
        .rgSettingProblem(settings, gwas, gwasH2Estimate)
    if (!is.null(problem)) {
        abort(glue("estimateLocalRg: {problem}."))
    }
    settings
}

# Where the shared noise a pair's settings hold cannot be had, or NULL.
# @noRd
.rgNoiseProblem <- function(settings, gwas, noiseCorrelation, individual) {
    tied <- identical(settings@crossNoise, "tied")
    if (tied && !individual) {
        return(str_c(
            "tied shared noise reads each pair's phenotypic correlation and ",
            "sample overlap, which summary statistics do not carry"
        ))
    }
    if (tied && !is.null(gwas)) {
        return(str_c(
            "tied shared noise reads the two phenotypes' correlation, which ",
            "a GWAS side does not have"
        ))
    }
    if (
        identical(settings@crossNoise, "supplied") && is.null(noiseCorrelation)
    ) {
        return("crossNoise = \"supplied\" needs `noiseCorrelation`")
    }
    NULL
}

# What else a pair cannot be fitted under, or NULL.
# @noRd
.rgSettingProblem <- function(settings, gwas, gwasH2Estimate) {
    if (identical(settings@criterion, "reml")) {
        return(str_c(
            "REML has no bivariate fit here (qgg::greml fits one trait); ",
            "tied shared noise takes what individual-level data add"
        ))
    }
    if (
        identical(settings@criterion, "likelihood") &&
            !identical(settings@noise, "free")
    ) {
        return(str_c(
            "the likelihood's covariance step reads each condition's own ",
            "intercept, so it fits free noise"
        ))
    }
    if (!is.null(gwasH2Estimate) && is.null(gwas)) {
        return("`gwasH2Estimate` holds a GWAS's intercept, so it needs `gwas`")
    }
    NULL
}

# A noise-correlation matrix: square, symmetric, named alike on both sides,
# 1 on the diagonal and within [-1, 1].
# @noRd
.rgAssertNoiseCorrelation <- function(x) {
    ok <- is.matrix(x) &&
        is.numeric(x) &&
        nrow(x) == ncol(x) &&
        !anyNA(x) &&
        !is.null(rownames(x)) &&
        identical(rownames(x), colnames(x)) &&
        isSymmetric(unname(x)) &&
        all(abs(diag(x) - 1) < 1e-8) &&
        all(abs(x) <= 1 + 1e-8)
    if (ok) {
        return(invisible(TRUE))
    }
    abort(str_c(
        "estimateLocalRg: `noiseCorrelation` must be a symmetric correlation ",
        "matrix, 1 on the diagonal, named by the conditions on both sides ",
        "(as mashResidualCorrelation() returns one)."
    ))
}

# Each GWAS study's genome-wide intercept from `gwasH2Estimate`, NULL without
# one.
# @noRd
.rgGwasIntercepts <- function(gwas, gwasH2Estimate) {
    if (is.null(gwasH2Estimate)) {
        return(NULL)
    }
    studies <- unique(as.character(gwas$study))
    set_names(
        map_dbl(studies, .rgGwasIntercept, h2 = gwasH2Estimate),
        studies
    )
}

# One GWAS study's genome-wide intercept: its row's in `h2`.
# @noRd
.rgGwasIntercept <- function(study, h2) {
    at <- which(
        as.character(h2$study) == study & is.na(h2$context) & is.na(h2$trait)
    )
    intercept <- if (length(at) == 1L) h2$intercept[[at]] else NA_real_
    if (!is.na(intercept)) {
        return(intercept)
    }
    abort(glue(
        "estimateLocalRg: `gwasH2Estimate` holds no genome-wide intercept ",
        "for GWAS study '{study}'."
    ))
}

# ---- the units: pairs, or groups fitted at once ------------------------------

# The units the pairs are fitted in. Without `gwas`, the selected rows are
# grouped along each jointSpecification's axes (contexts within each (study,
# trait) by default), and each group gives its pairs, or, where the settings
# fit a group at once (.rgFitsAtOnce()), one unit of all its rows; with
# `gwas`, each selected row pairs with each GWAS study.
# @noRd
.rgUnits <- function(data, rows, jointSpecification, settings, run) {
    if (!is.null(run$gwas) && !run$conditional) {
        .rgAssertGwasPairing(jointSpecification)
        return(.rgGwasUnits(data, rows, run$gwas))
    }
    specs <- parseJointSpecification(jointSpecification %||% "context", data)
    units <- list_c(map(
        specs,
        .rgSpecUnits,
        data = data,
        rows = rows,
        atOnce = .rgFitsAtOnce(settings),
        gwas = run$gwas
    ))
    if (length(units) == 0L) {
        abort(str_c(
            "estimateLocalRg: no two selected rows share a group along the ",
            "joint axes, so there is no pair to estimate."
        ))
    }
    units
}

# One joint spec's units: its selected rows in scope, grouped by the identity
# columns its axes leave out; with `gwas`, a group of one row too.
# @noRd
.rgSpecUnits <- function(spec, data, rows, atOnce, gwas) {
    inScope <- rows[.rgRowsInScope(data, rows, spec$scope)]
    complement <- setdiff(c("study", "context", "trait"), spec$axes)
    keys <- .rgJoinedColumns(data, inScope, complement, "\r")
    groups <- unname(split(inScope, keys))
    list_c(map(
        groups[lengths(groups) >= (if (is.null(gwas)) 2L else 1L)],
        .rgGroupUnits,
        data = data,
        axes = intersect(c("study", "context", "trait"), spec$axes),
        atOnce = atOnce,
        gwas = gwas
    ))
}

# Whether each of `rows` lies in the spec's `scope` (every row without one).
# @noRd
.rgRowsInScope <- function(data, rows, scope) {
    inKeys <- map(
        names(scope),
        .rgRowsInScopeKey,
        data = data,
        rows = rows,
        scope = scope
    )
    reduce(inKeys, `&`, .init = rep(TRUE, length(rows)))
}

# @noRd
.rgRowsInScopeKey <- function(key, data, rows, scope) {
    is_in(.rgColumnAt(key, data, rows), scope[[key]])
}

# The identity columns `cols` of `rows`, joined by `sep`; one key for all
# when there are none.
# @noRd
.rgJoinedColumns <- function(data, rows, cols, sep) {
    if (length(cols) == 0L) {
        return(rep("all", length(rows)))
    }
    parts <- map(cols, .rgColumnAt, data = data, rows = rows)
    exec(str_c, !!!parts, sep = sep)
}

# Identity column `col` of `rows`: an mcols column, as `[[` addresses a
# collection's elements.
# @noRd
.rgColumnAt <- function(col, data, rows) {
    as.character(mcols(data, use.names = FALSE)[[col]][rows])
}

# A group's units: one of all its rows where the settings fit a group at
# once, its pairs otherwise. Rows on different LD references are fitted each
# on its own (.rgUnitLd()); each row is labelled by its values on the joint
# axes, as `noiseCorrelation` names them.
# @noRd
.rgGroupUnits <- function(group, data, axes, atOnce, gwas) {
    labels <- .rgJoinedColumns(data, group, axes, "|")
    if (!is.null(gwas)) {
        return(map(
            unique(as.character(gwas$study)),
            .rgGroupGwasUnit,
            group = group,
            labels = labels
        ))
    }
    if (atOnce) {
        return(list(.rgUnit(group, labels)))
    }
    pairs <- .rgGroupPairs(length(group))
    map2(pairs$i, pairs$j, .rgPairUnit, group = group, labels = labels)
}

# @noRd
.rgPairUnit <- function(i, j, group, labels) {
    .rgUnit(group[c(i, j)], labels[c(i, j)])
}

# A group with a GWAS study, for its pairs given the group's other rows: the
# rows labelled by their values on the joint axes, the GWAS by its study.
# @noRd
.rgGroupGwasUnit <- function(study, group, labels) {
    .rgUnit(group, c(labels, study), study)
}

# A unit: its QTL rows, their labels, and the GWAS study it pairs them with.
# @noRd
.rgUnit <- function(rows, labels, gwasStudy = NA_character_) {
    list(rows = rows, labels = labels, gwasStudy = gwasStudy)
}

# Each selected row with each GWAS study, labelled by the two studies, as the
# noise they share is the two studies' sample overlap.
# @noRd
.rgGwasUnits <- function(data, rows, gwas) {
    grid <- expand.grid(
        row = rows,
        study = unique(as.character(gwas$study)),
        stringsAsFactors = FALSE
    )
    pmap(grid, .rgGwasUnit, data = data)
}

# @noRd
.rgGwasUnit <- function(row, study, data) {
    .rgUnit(row, c(as.character(data$study[[row]]), study), study)
}

# ---- one unit's fit ---------------------------------------------------------

# One unit's pairs: its conditions' statistics on the variants they all carry,
# matched by allele to the first condition's coding (and, where they read
# different LD references, held by every one, a strand-ambiguous variant
# kept only where it needs no strand decision), each condition's LD reference
# for those variants in the first QTL row's window, its reference overlap and
# held intercept, and the noise correlations its labels name; then the
# settings' fit, of its QTL-GWAS pairs alone where it has a GWAS.
# @noRd
.rgUnitFit <- function(unit, data, param, overlaps, run) {
    window <- .rowRegionAt(data, unit$rows[[1L]])
    members <- .rgUnitMembers(unit, data, run$gwas, window)
    label <- .rgUnitLabel(members)
    refs <- .rgDistinctRefs(.rgUnitRefs(unit, data, run$gwas, window))
    oneReference <- length(refs$refs) == 1L
    common <- .rgCommonHeldBy(
        .rgCommonVariants(members, is.na(unit$gwasStudy) && oneReference),
        if (oneReference) list() else refs$refs,
        label
    )
    out <- list(window = window, members = members)
    if (length(common$ids) < 2L) {
        empty <- .rgNoPairs(length(members), length(common$ids))
        return(list_assign(out, pairs = .rgUnitPairs(empty, unit, members)))
    }
    ud <- c(
        list(z = common$z, nVariant = common$N),
        .rgUnitLd(common$ids, refs, param, label),
        .rgUnitTerms(unit, overlaps, run, label)
    )
    .rgAssertUnitOverlap(param, ud, label)
    pairs <- .rgFit(param, ud)
    list_assign(out, pairs = .rgUnitPairs(pairs, unit, members))
}

# A QtlSumStats unit's per-condition terms: each condition's reference
# overlap (a GWAS side's none) and held intercept, the noise correlations its
# labels name, and whether its pairs are given its other conditions.
# @noRd
.rgUnitTerms <- function(unit, overlaps, run, label) {
    list(
        overlap = c(overlaps[unit$rows], if (!is.na(unit$gwasStudy)) 0),
        heldIntercept = .rgHeldIntercepts(
            length(unit$rows),
            unit$gwasStudy,
            run$gwasIntercepts
        ),
        noiseCorrelation = .rgUnitNoise(
            run$noiseCorrelation,
            unit$labels,
            label
        ),
        conditional = run$conditional
    )
}

# A unit's pairs without an estimate, `nSnps` the variants it holds.
# @noRd
.rgNoPairs <- function(P, nSnps) {
    pmap(.rgGroupPairs(P), .rgPairResult, nSnps = nSnps)
}

# The pairs a unit reports: all of them, or, where it has a GWAS (last among
# its members), each QTL condition's with it.
# @noRd
.rgUnitPairs <- function(pairs, unit, members) {
    if (is.na(unit$gwasStudy)) {
        return(pairs)
    }
    keep(pairs, .rgIsGwasPair, P = length(members))
}

# @noRd
.rgIsGwasPair <- function(pair, P) {
    pair$j == P
}

# Each condition's held intercept: none for the `nQtl` QTL conditions, and for
# a GWAS side its genome-wide intercept where `gwasH2Estimate` gives one.
# @noRd
.rgHeldIntercepts <- function(nQtl, gwasStudy, gwasIntercepts) {
    qtl <- rep(NA_real_, nQtl)
    if (is.na(gwasStudy)) {
        return(qtl)
    }
    c(qtl, gwasIntercepts[[gwasStudy]] %||% NA_real_)
}

# The unit's conditions: each QTL row's statistics, and the GWAS study's in
# the first row's window.
# @noRd
.rgUnitMembers <- function(unit, data, gwas, window) {
    qtl <- map(unit$rows, .rgQtlMember, data = data)
    if (is.na(unit$gwasStudy)) {
        return(qtl)
    }
    c(qtl, list(.rgGwasMember(gwas, unit$gwasStudy, window)))
}

# @noRd
.rgQtlMember <- function(row, data) {
    st <- as.character(data$study[[row]])
    ctx <- as.character(data$context[[row]])
    tr <- as.character(data$trait[[row]])
    list(
        study = st,
        context = ctx,
        trait = tr,
        df = .sumStatsDf(
            data,
            studyName = st,
            context = ctx,
            trait = tr,
            require = c("Z", "N"),
            derive = "zFromBetaSe"
        )
    )
}

# A GWAS study's statistics in `window`.
# @noRd
#' @importFrom IRanges subsetByOverlaps
.rgGwasMember <- function(gwas, study, window) {
    entry <- sumStats(gwas, studyName = study)
    inWindow <- if (is.null(window)) {
        entry[0L]
    } else {
        subsetByOverlaps(entry, window)
    }
    list(
        study = study,
        context = NA_character_,
        trait = NA_character_,
        df = .entryToSumstatDf(
            inWindow,
            require = c("Z", "N"),
            derive = "zFromBetaSe",
            label = glue("GwasSumStats[{study}]")
        )
    )
}

# The variants every condition carries, matched by allele to the first
# condition's coding: their ids, and each condition's z (sign-flipped where
# its alleles are swapped) and sample sizes. QTL rows harmonized to one LD
# reference keep their strand-ambiguous variants; a GWAS side keeps one only
# where its alleles agree as given, as matchVariants() judges each pair (a
# swapped A/T could be a strand flip).
# @noRd
.rgCommonVariants <- function(members, keepAmbiguous) {
    dfs <- map(members, "df")
    ids <- dfs[[1L]]$variant_id
    matches <- map(
        dfs[-1L],
        .rgMatchVariants,
        ids = ids,
        removeStrandAmbiguous = !keepAmbiguous
    )
    keep <- reduce(map(matches, "idxA"), intersect, .init = seq_along(ids))
    list(
        ids = ids[keep],
        z = c(
            list(dfs[[1L]]$z[keep]),
            map2(matches, dfs[-1L], .rgAligned, keep = keep, column = "z")
        ),
        N = c(
            list(dfs[[1L]]$N[keep]),
            map2(matches, dfs[-1L], .rgAligned, keep = keep, column = "N")
        )
    )
}

# @noRd
.rgMatchVariants <- function(df, ids, removeStrandAmbiguous) {
    matchVariants(
        ids,
        df$variant_id,
        removeStrandAmbiguous = removeStrandAmbiguous
    )
}

# A condition's `column` on the kept variants, in the first condition's
# order; z sign-flipped where the alleles are swapped.
# @noRd
.rgAligned <- function(m, df, keep, column) {
    at <- match(keep, m$idxA)
    values <- df[[column]][m$idxB[at]]
    if (identical(column, "z")) values * m$sign[at] else values
}

# The unit's conditions, named for messages.
# @noRd
.rgUnitLabel <- function(members) {
    tuples <- map_chr(members, .rgMemberTuple)
    glue("estimateLocalRg: {str_flatten(tuples, ' and ')}")
}

# @noRd
.rgMemberTuple <- function(member) {
    parts <- c(member$study, member$context, member$trait)
    str_c("(", str_flatten(parts[!is.na(parts)], ", "), ")")
}

# The unit's noise correlations, by its conditions' labels; NULL without a
# matrix.
# @noRd
.rgUnitNoise <- function(noiseCorrelation, labels, label) {
    if (is.null(noiseCorrelation)) {
        return(NULL)
    }
    missing <- setdiff(labels, rownames(noiseCorrelation))
    if (length(missing) > 0L) {
        abort(glue(
            "{label}: `noiseCorrelation` names no condition ",
            "{str_flatten(missing, ', ')}; it is named by the conditions' ",
            "values on the joint axes, or by the two studies for a GWAS pair."
        ))
    }
    noiseCorrelation[labels, labels, drop = FALSE]
}

# A reference cannot hold more of a QTL study's samples than it has, where
# the overlap factor reads the overlap.
# @noRd
.rgAssertUnitOverlap <- function(param, ud, label) {
    if (!identical(param@referenceNoise, "overlapFactor")) {
        return(invisible(NULL))
    }
    pwalk(
        list(
            n = map_dbl(ud$nVariant, mean),
            overlap = ud$overlap,
            nRef = map_dbl(seq_along(ud$z), .rgMemberLdSize, ud = ud)
        ),
        .cisH2AssertOverlap,
        label = label
    )
}

# ---- conditions on different LD references -----------------------------------

# Each condition's LD reference: a QTL row's own, and a GWAS study's own where
# its rows in the window name one (.rgGwasRef()), the first QTL row's
# otherwise.
# @noRd
.rgUnitRefs <- function(unit, data, gwas, window) {
    qtl <- map(unit$rows, .rowLdRef, x = data)
    if (is.na(unit$gwasStudy)) {
        return(qtl)
    }
    c(qtl, list(.rgGwasRef(gwas, unit$gwasStudy, window) %||% qtl[[1L]]))
}

# A GWAS study's own LD reference in `window`: the one its rows overlapping
# the window name; NULL where they name none, refused where they name two.
# @noRd
.rgGwasRef <- function(gwas, study, window) {
    if (is.null(window)) {
        return(NULL)
    }
    rows <- which(as.character(gwas$study) == study)
    inWindow <- rows[map_lgl(rows, .rgRowOverlaps, x = gwas, window = window)]
    refs <- .rgDistinctRefs(compact(map(inWindow, .rowLdRef, x = gwas)))$refs
    if (length(refs) > 1L) {
        abort(glue(
            "estimateLocalRg: GWAS study {study}'s rows in the window name ",
            "different LD references; give them one."
        ))
    }
    if (length(refs) == 0L) NULL else refs[[1L]]
}

# Whether row `i` of `x` overlaps `window`.
# @noRd
.rgRowOverlaps <- function(i, x, window) {
    region <- .rowRegionAt(x, i)
    !is.null(region) && overlapsAny(region, window)
}

# The distinct LD references among `refs`, identical ones counted once, and
# which one each condition reads.
# @noRd
.rgDistinctRefs <- function(refs) {
    first <- map_int(refs, .rgFirstIdentical, refs = refs)
    keys <- unique(first)
    list(refs = refs[keys], of = match(first, keys))
}

# @noRd
.rgFirstIdentical <- function(ref, refs) {
    detect_index(refs, identical, ref)
}

# The common variants every reference in `refs` holds, so that each gives the
# LD of every variant the pairs read; `common` as it is for none.
# @noRd
.rgCommonHeldBy <- function(common, refs, label) {
    if (length(refs) == 0L || length(common$ids) == 0L) {
        return(common)
    }
    keep <- which(reduce(
        map(refs, .rgRefHolds, ids = common$ids, label = label),
        `&`
    ))
    list(
        ids = common$ids[keep],
        z = map(common$z, `[`, keep),
        N = map(common$N, `[`, keep)
    )
}

# Which of `ids` reference `ref` holds, matched by allele.
# @noRd
.rgRefHolds <- function(ref, ids, label) {
    matched <- .ldFromRefMatch(ref, ids, label, "drop")
    is_in(ids, matched$keptIds %||% character(0))
}

# The unit's LD: one source where its conditions read one reference, else
# one per distinct reference (`lds`) and the one each condition reads
# (`ldOf`).
# @noRd
.rgUnitLd <- function(ids, refs, param, label) {
    sources <- map(
        refs$refs,
        .cisH2LdSource,
        ids = ids,
        param = param,
        label = label
    )
    if (length(sources) == 1L) {
        return(list(ld = sources[[1L]]))
    }
    list(ld = NULL, lds = sources, ldOf = refs$of)
}

# Condition `k`'s LD source.
# @noRd
.rgMemberLd <- function(k, ud) {
    ud$ld %||% ud$lds[[ud$ldOf[[k]]]]
}

# @noRd
.rgMemberLdSize <- function(k, ud) {
    .rgMemberLd(k, ud)$nSamples %||% NA_real_
}

# Tied noise reads the conditions' phenotypes on the donors they share, which
# a summary row does not carry.
# @noRd
.rgAssertTiedDonors <- function(param, ud, label) {
    if (
        !identical(param@crossNoise, "tied") ||
            !is.null(ud$phenotypicCorrelation)
    ) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: tied noise reads the conditions' phenotypes on the donors ",
        "they share, which a summary row does not carry."
    ))
}

# The settings' fit of one unit: its pairs' results. A pair the likelihood
# fits on its own, given no other condition, has each r2 its rg squared.
# @noRd
.rgFit <- function(param, ud) {
    if (!is.null(ud$lds)) {
        return(.rgTwoLdFit(param, ud))
    }
    fit <- switch(
        str_c(param@statistic, "/", param@criterion),
        "eigen/likelihood" = if (length(ud$z) > 2L) {
            .rgJointLikelihood
        } else {
            .rgEigenLikelihood
        },
        "eigen/regression" = .rgEigenRegression,
        "eigen/moments" = .rgEigenMoments,
        "variant/regression" = .rgVariantRegression,
        "variant/gls" = .rgVariantGls
    )
    pairs <- fit(param, ud)
    ownPair <- length(ud$z) == 2L && identical(param@criterion, "likelihood")
    if (isTRUE(ud$conditional) && ownPair) map(pairs, .rgPairR2) else pairs
}

# @noRd
.rgPairR2 <- function(pair) {
    rg2 <- pair$values$rg^2
    list_assign(pair, values = list_assign(pair$values, r2_1 = rg2, r2_2 = rg2))
}

# =============================================================================
# Conditions on different LD references
# -----------------------------------------------------------------------------
# Two conditions whose statistics come from populations with different LD
# (across ancestries, or a GWAS on its own panel) share the genetic term
# R_1 D_12 R_2 and, with no shared samples, no noise (section 4c of the
# design notes): r_k = R_k beta_k + e_k, Cov(r_1, r_2) = G_12 R_1 R_2 / M on
# the per-SD scale. Each condition's heritability is its own fit on its own
# reference, as estimateLocalH2() makes it; the covariance takes each
# setting's two-reference form. A pair whose two sides read one reference is
# fitted on it as before.
# =============================================================================

# A unit whose conditions read different LD references: each pair on the LD
# its two sides read (.rgMixedPair()), or every condition of a group at once
# where the settings fit one (.rgMixedGroup()).
# @noRd
.rgTwoLdFit <- function(param, ud) {
    P <- length(ud$z)
    if (P == 2L) {
        return(list(.rgTwoLdPair(param, ud)))
    }
    if (.rgFitsAtOnce(param)) {
        return(.rgMixedGroup(param, ud))
    }
    pairs <- .rgGroupPairs(P)
    map2(pairs$i, pairs$j, .rgMixedPair, param = param, ud = ud)
}

# Pair (i, j) of a unit on several references, on the LD its two sides read.
# @noRd
.rgMixedPair <- function(i, j, param, ud) {
    pair <- .rgFit(param, .rgPairUnitData(ud, c(i, j)))[[1L]]
    list_assign(pair, i = i, j = j)
}

# A pair whose two sides read different references, by the settings'
# criterion's two-reference form.
# @noRd
.rgTwoLdPair <- function(param, ud) {
    fit <- switch(
        str_c(param@statistic, "/", param@criterion),
        "eigen/likelihood" = .rgTwoLdLikelihood,
        "eigen/regression" = .rgTwoLdEigenRegression,
        "eigen/moments" = .rgTwoLdMoments,
        "variant/regression" = .rgTwoLdVariantRegression,
        "variant/gls" = .rgTwoLdGls
    )
    fit(.rgPairNoiseSettings(param, ud), ud)
}

# A pair's settings as its two sides allow: under tied noise, two conditions
# that share no donors (their phenotypes' correlation NA, as for two
# studies') share no noise, the overlap it scales with being none.
# @noRd
.rgPairNoiseSettings <- function(param, ud) {
    unshared <- identical(param@crossNoise, "tied") &&
        is.na(ud$phenotypicCorrelation[[1L, 2L]])
    if (unshared) methods::initialize(param, crossNoise = "zero") else param
}

# Where tied noise ties a pair's shared noise to its phenotypes: every pair
# whose two conditions share donors (their phenotypes' correlation known);
# none under other settings.
# @noRd
.rgTieOf <- function(param, ud, P) {
    if (!identical(param@crossNoise, "tied")) {
        return(matrix(0, P, P))
    }
    1 * !is.na(ud$phenotypicCorrelation)
}

# The noise a pair reports it shares: the fit's, or none for tied noise
# between two conditions that share no donors.
# @noRd
.rgPairNoiseSource <- function(mode, tie, i, j) {
    if (identical(mode, "tied") && tie[[i, j]] == 0) "zero" else mode
}

# Condition `k`'s eigen statistics on its own reference: the eigenpairs kept
# under the settings' truncation at its own sample size, and its correlations
# rotated onto them; NULL when fewer than the truncation's minimum are kept.
# @noRd
.rgSideEigen <- function(k, param, ud) {
    ld <- .rgMemberLd(k, ud)
    n <- .h2SampleSize(param@sampleSize, ud$nVariant[[k]])
    M <- length(ud$z[[k]])
    dec <- .h2Basis(param, ld, n, M)
    if (length(dec$values) < param@truncation@minimum) {
        return(NULL)
    }
    r <- .zToCor(ud$z[[k]], ud$nVariant[[k]])
    list(
        lam = dec$values,
        vectors = dec$vectors,
        plain = dec$plain,
        rotation = dec$rotation,
        weights = dec$weights,
        M = M,
        n = n,
        r = r,
        u = .rgRotate(r, dec$vectors),
        nRef = if (identical(param@referenceNoise, "eigenTerm")) {
            ld$nSamples
        } else {
            Inf
        },
        ld = ld
    )
}

# A side's data for its own likelihood (.h2LikelihoodFree()).
# @noRd
.rgSideMember <- function(side, floor) {
    list(
        lam = side$lam,
        bstar = side$u,
        n = side$n,
        nRef = side$nRef,
        M = side$M,
        floor = floor
    )
}

# ---- the likelihood on two references ----------------------------------------

# A pair's likelihood on two references: each side's heritability and intercept
# from its own likelihood on its own reference (step 1), then the covariance.
# With no shared noise, from the reduced likelihood (.rgReduce()) fitted as a
# pair's step 2 on one reference; with shared noise, which on two references
# takes the shape (R1 + R2) / 2, from the Gaussian of both sides' statistics at
# once, fitted as a group of two (.rgMixedLikelihood()).
# @noRd
.rgTwoLdLikelihood <- function(param, ud) {
    if (!identical(param@crossNoise, "zero")) {
        return(.rgMixedLikelihood(param, ud)[[1L]])
    }
    M <- length(ud$z[[1L]])
    sides <- map(seq_len(2L), .rgSideEigen, param = param, ud = ud)
    if (any(map_lgl(sides, is.null))) {
        return(.rgPairResult(nSnps = M))
    }
    floor <- param@procedures@varianceFloor
    members <- map(sides, .rgSideMember, floor = floor)
    steps <- map2(
        members,
        ud$heldIntercept,
        .rgLikelihoodStep1,
        param = param
    )
    h2 <- map_dbl(steps, "h2")
    intercepts <- map_dbl(steps, "intercept")
    base <- list(
        nSnps = M,
        nComponents = min(lengths(map(sides, "lam"))),
        h2 = h2
    )
    if (anyNA(h2) || anyNA(intercepts)) {
        return(exec(.rgPairResult, !!!base))
    }
    if (any(h2 <= 0)) {
        return(exec(.rgPairResult, !!!base, covariance = 0))
    }
    v <- pmap(
        list(h2 = h2, data = members, intercept = intercepts),
        .h2WindowVariance
    )
    reduced <- .rgReduce(sides, v, floor, sqrt(prod(h2)))
    .rgReducedResult(param, reduced, base)
}

# The two sides' statistics as independent pairs (section 4c): each whitened
# by its step-1 variances, D_k^-1/2 u_k, and rotated by the SVD of W =
# D_1^-1/2 B D_2^-1/2, B the genetic covariance's shape (Lambda_1 C Lambda_2
# / M with C = Q_1'Q_2 the two references' eigenvectors' overlap on the
# shared variants; .rgCrossPieces()), so that a_i and b_i have correlation g
# s_i at covariance g and the rest are independent. The shared noise's shape
# N = (Lambda_1 C + C Lambda_2) / 2, the LD of the samples the two share
# taken as the average of the two references' (Lambda on one), is kept with
# its reduced form Nred, which the rotation does not make diagonal.
# @noRd
.rgReduce <- function(sides, v, floor, bound) {
    pieces <- .rgCrossPieces(sides[[1L]], sides[[2L]])
    roots <- tcrossprod(sqrt(sides[[1L]]$lam), sqrt(sides[[2L]]$lam))
    B <- pieces$genetic * roots / sides[[1L]]$M
    N <- pieces$noise * roots
    whitened <- sweep(B / sqrt(v[[1L]]), 2L, sqrt(v[[2L]]), `/`)
    dec <- svd(whitened)
    noise <- sweep(N / sqrt(v[[1L]]), 2L, sqrt(v[[2L]]), `/`)
    list(
        a = as.vector(crossprod(dec$u, sides[[1L]]$u / sqrt(v[[1L]]))),
        b = as.vector(crossprod(dec$v, sides[[2L]]$u / sqrt(v[[2L]]))),
        s = dec$d,
        U = dec$u,
        V = dec$v,
        C = pieces$C,
        B = B,
        N = N,
        Nred = crossprod(dec$u, noise %*% dec$v),
        floor = floor,
        bound = bound
    )
}

# Two kept bases' overlaps on their components, delta = lambda^-1/2 u: C =
# Q_a'Q_b of their plain eigenvectors; the genetic cross term Lambda_a^1/2
# Q_a' W_ab Q_b Lambda_b^1/2, W_ab = diag(sqrt(w_a w_b)) of the two
# references' frequency weights (C o sqrt(lambda_a lambda_b') without them);
# and the shared noise's shape (Lambda_a^1/2 C Lambda_b^-1/2 + Lambda_a^-1/2
# C Lambda_b^1/2) / 2. Where weights made pseudo-eigenpairs (.h2Weighted()),
# each is turned onto them, U_a' X U_b.
# @noRd
.rgCrossPieces <- function(a, b) {
    Qa <- (a$plain %||% a)$vectors
    Qb <- (b$plain %||% b)$vectors
    ra <- sqrt(if (is.null(a$plain)) a$lam else a$plain$values)
    rb <- sqrt(if (is.null(b$plain)) b$lam else b$plain$values)
    C <- crossprod(Qa, Qb)
    overlap <- if (is.null(a$weights)) {
        C
    } else {
        crossprod(Qa * sqrt(a$weights * b$weights), Qb)
    }
    pieces <- list(
        C = C,
        genetic = overlap * tcrossprod(ra, rb),
        noise = C * (tcrossprod(ra, 1 / rb) + tcrossprod(1 / ra, rb)) / 2
    )
    map(pieces, .rgRotatePair, Ua = a$rotation, Ub = b$rotation)
}

# @noRd
.rgRotatePair <- function(X, Ua, Ub) {
    if (is.null(Ua)) X else crossprod(Ua, X %*% Ub)
}

# The reduced pairs' log-likelihood at covariance g: each b_i given a_i, mean
# g s_i a_i and variance 1 - (g s_i)^2, floored; a_i's own terms, and the
# unpaired coordinates', do not move with g.
# @noRd
.rgReducedLl <- function(g, red) {
    .h2Likelihood("reduced", red, g)$value
}

# The reduced pairs' log-likelihood at covariance g as a point of its profile
# (.rgProfile()), nothing else to refit.
# @noRd
.rgReducedRefit <- function(g, start, red) {
    list(value = .rgReducedLl(g, red), par = numeric(0))
}

# ---- a group on several references at once -----------------------------------

# A group of three or more on several references, fitted at once as the
# settings fit one: the likelihood's joint fit (.rgMixedLikelihood()), or the
# moments (.rgMixedMoments()).
# @noRd
.rgMixedGroup <- function(param, ud) {
    if (identical(param@criterion, "moments")) {
        return(.rgMixedMoments(param, ud))
    }
    .rgMixedLikelihood(param, ud)
}

# The joint likelihood on several references: each reference's eigenpairs
# kept at the smallest sample size among the conditions that read it, each
# condition's statistics rotated onto its reference's and its step-1 fit
# there (.rgMixedStep1()), then the correlation matrix fitted at once as on
# one reference (.rgJointFit()), the log-likelihood the Gaussian of every
# condition's rotated statistics together (.rgJointLl()), the pairs' shared
# noise as the settings take it.
# @noRd
.rgMixedLikelihood <- function(param, ud) {
    P <- length(ud$z)
    pairs <- .rgGroupPairs(P)
    M <- length(ud$z[[1L]])
    bases <- map(seq_along(ud$lds), .rgBasisEigen, param = param, ud = ud)
    if (any(map_lgl(bases, is.null))) {
        return(pmap(pairs, .rgPairResult, nSnps = M))
    }
    step1 <- .rgMixedStep1(param, ud, bases)
    h2 <- step1$h2
    base <- list(nSnps = M, nComponents = min(lengths(map(bases, "lam"))))
    fitted <- !is.na(h2) & !is.na(step1$intercepts)
    if (!all(fitted & h2 > 0)) {
        return(.rgJointLive(param, ud, h2, fitted, base))
    }
    joint <- .rgMixedData(
        bases,
        step1,
        ud$ldOf,
        param@procedures@varianceFloor,
        .rgTieOf(param, ud, P)
    )
    fit <- .rgJointFit(param, ud, NULL, step1$members, joint)
    if (is.null(fit)) {
        return(pmap(pairs, .rgJointNoEstimate, h2 = h2, base = base))
    }
    pmap(
        pairs,
        .rgJointPairResult,
        fit = fit,
        joint = joint,
        ud = ud,
        base = base,
        param = param
    )
}

# Step 1 on several references: each condition's statistics rotated onto
# its reference's eigenpairs (.rgMemberSide()) and its heritability and
# intercept from its own likelihood there.
# @noRd
.rgMixedStep1 <- function(param, ud, bases) {
    sides <- map(
        seq_along(ud$z),
        .rgMemberSide,
        bases = bases,
        ud = ud,
        param = param
    )
    members <- map(sides, .rgSideMember, floor = param@procedures@varianceFloor)
    steps <- map2(
        members,
        ud$heldIntercept,
        .rgLikelihoodStep1,
        param = param
    )
    list(
        sides = sides,
        members = members,
        h2 = map_dbl(steps, "h2"),
        intercepts = map_dbl(steps, "intercept")
    )
}

# Reference `a`'s eigenpairs, kept under the settings' truncation at the
# smallest sample size among the conditions that read it; NULL when fewer
# than the minimum are kept.
# @noRd
.rgBasisEigen <- function(a, param, ud) {
    readers <- which(ud$ldOf == a)
    n <- min(map_dbl(
        ud$nVariant[readers],
        .h2SampleSize,
        how = param@sampleSize
    ))
    M <- length(ud$z[[1L]])
    ld <- ud$lds[[a]]
    dec <- .h2Basis(param, ld, n, M)
    if (length(dec$values) < param@truncation@minimum) {
        return(NULL)
    }
    list(
        lam = dec$values,
        vectors = dec$vectors,
        plain = dec$plain,
        rotation = dec$rotation,
        weights = dec$weights,
        ld = ld,
        nRef = if (identical(param@referenceNoise, "eigenTerm")) {
            ld$nSamples
        } else {
            Inf
        }
    )
}

# Condition `k`'s statistics on its reference's kept eigenpairs.
# @noRd
.rgMemberSide <- function(k, bases, ud, param) {
    basis <- bases[[ud$ldOf[[k]]]]
    r <- .zToCor(ud$z[[k]], ud$nVariant[[k]])
    list(
        lam = basis$lam,
        vectors = basis$vectors,
        plain = basis$plain,
        rotation = basis$rotation,
        weights = basis$weights,
        M = length(r),
        n = .h2SampleSize(param@sampleSize, ud$nVariant[[k]]),
        r = r,
        u = .rgRotate(r, basis$vectors),
        nRef = basis$nRef,
        ld = basis$ld
    )
}

# What the joint fit on several references reads: each condition's rotated
# statistics and own modelled variances at step 1 (.rgMixedStep1()), the
# reference it reads, and for each two references a and b the genetic
# covariance's shape B_ab = Lambda_a C_ab Lambda_b / M and the shared noise's
# N_ab = (Lambda_a C_ab + C_ab Lambda_b) / 2 (C_ab = Q_a'Q_b; on one
# reference Lambda^2 / M and Lambda), so that conditions k and l covary G_kl
# B_ab + c_kl N_ab, the noise c_kl on the intercepts' scale over the root of
# the two sample sizes as on one reference; the procedures' variance floor;
# and where tied noise ties a pair's (.rgTieOf()).
# @noRd
.rgMixedData <- function(bases, step1, basisOf, floor, tie) {
    P <- length(step1$h2)
    grid <- expand.grid(a = seq_along(bases), b = seq_along(bases))
    list(
        mixed = TRUE,
        P = P,
        nPairs = P * (P - 1L) / 2L,
        u = map(step1$sides, "u"),
        v = pmap(
            list(
                h2 = step1$h2,
                data = step1$members,
                intercept = step1$intercepts
            ),
            .h2WindowVariance
        ),
        basisOf = basisOf,
        B = matrix(
            pmap(grid, .rgMixedShape, bases = bases, M = step1$sides[[1L]]$M),
            length(bases)
        ),
        N = matrix(
            pmap(grid, .rgMixedNoiseShape, bases = bases),
            length(bases)
        ),
        crossScale = 1 / sqrt(tcrossprod(map_dbl(step1$sides, "n"))),
        h2 = step1$h2,
        intercepts = step1$intercepts,
        root = sqrt(tcrossprod(pmax(step1$intercepts, 0))),
        floor = floor,
        tied = tie
    )
}

# @noRd
.rgMixedNoiseShape <- function(a, b, bases) {
    pieces <- .rgCrossPieces(bases[[a]], bases[[b]])
    pieces$noise * tcrossprod(sqrt(bases[[a]]$lam), sqrt(bases[[b]]$lam))
}

# @noRd
.rgMixedShape <- function(a, b, bases, M) {
    pieces <- .rgCrossPieces(bases[[a]], bases[[b]])
    pieces$genetic *
        tcrossprod(sqrt(bases[[a]]$lam), sqrt(bases[[b]]$lam)) /
        M
}

# The joint fit's start on several references: every pair's own fit on the
# LD its two sides read (.rgMixedPair()), its correlations put together and
# moved to the nearest correlation matrix shrunk a little inside the
# boundary, as canonical partial correlations; and its shared noise, on the
# intercepts' scale, where the settings fit it: a pair on one reference
# (one study's conditions) its own window fit's, a pair on two none (their
# statistics' correlation carries their genetic covariance, and two studies
# seldom share many samples).
# @noRd
.rgMixedStart <- function(param, ud, joint) {
    pairs <- .rgGroupPairs(joint$P)
    own <- map2(
        pairs$i,
        pairs$j,
        .rgMixedPairStart,
        param = methods::initialize(
            param,
            procedures = methods::initialize(param@procedures, jointFit = FALSE)
        ),
        ud = ud,
        joint = joint
    )
    R <- .rgPairsMatrix(map_dbl(own, "r"), pairs, joint$P) + diag(joint$P)
    nearest <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    list(
        z = .rgCorrelationCpc(0.98 * nearest + 0.02 * diag(joint$P)),
        shared = .rgPairValues(
            .rgPairsMatrix(map_dbl(own, "shared"), pairs, joint$P)
        )
    )
}

# One pair's own fit, as a start: its rg (0 where it has none), with no
# shared noise on two references; and its shared noise on the intercepts'
# scale, its own window fit's on one reference, else none.
# @noRd
.rgMixedPairStart <- function(i, j, param, ud, joint) {
    oneLd <- ud$ldOf[[i]] == ud$ldOf[[j]]
    window <- identical(param@crossNoise, "window")
    own <- if (oneLd) param else methods::initialize(param, crossNoise = "zero")
    values <- .rgMixedPair(i, j, own, ud)$values
    noise <- if (oneLd && window) values$noise_correlation else NA_real_
    list(
        r = if (is.na(values$rg)) 0 else values$rg,
        shared = if (is.na(noise)) 0 else noise * joint$root[[i, j]]
    )
}

# The moments of a group on several references at once. Given the group's
# other conditions (conditional = TRUE), each pair's partial covariance and
# r2s from the group's matrix H (.rgMomentsGroupOf()), with the p-value and
# intervals from draws of every condition's components together
# (.rgMixedDrawn()), as the settings' effects read them. Otherwise each pair
# on its own (.rgMixedPair()): a pair's statistic reads its two sides alone.
# @noRd
.rgMixedMoments <- function(param, ud) {
    P <- length(ud$z)
    pairs <- .rgGroupPairs(P)
    if (!isTRUE(ud$conditional)) {
        return(map2(pairs$i, pairs$j, .rgMixedPair, param = param, ud = ud))
    }
    M <- length(ud$z[[1L]])
    bases <- map(seq_along(ud$lds), .rgBasisEigen, param = param, ud = ud)
    if (any(map_lgl(bases, is.null))) {
        return(pmap(pairs, .rgPairResult, nSnps = M))
    }
    sides <- map(
        seq_len(P),
        .rgMemberSide,
        bases = bases,
        ud = ud,
        param = param
    )
    if (any(map_dbl(sides, "n") <= lengths(map(sides, "lam")) + 1L)) {
        return(pmap(pairs, .rgPairResult, nSnps = M))
    }
    group <- .rgMomentsGroupOf(param, ud, sides, bases)
    pmap(pairs, .rgMixedConditionalPair, group = group, M = M)
}

# The shared noise per component a group's settings hold: the supplied
# correlation times the root of the two conditions' own; tied, each pair that
# shares donors as on one reference (.rgMomentsTiedNoise()), whose two
# conditions read one reference, and none between conditions that share no
# donors; or none.
# @noRd
.rgMixedSigma <- function(param, ud, own, shapes) {
    P <- length(own)
    switch(
        param@crossNoise,
        supplied = ud$noiseCorrelation *
            sqrt(tcrossprod(map_dbl(own, "noise"))),
        tied = matrix(
            pmap_dbl(
                expand.grid(k = seq_len(P), l = seq_len(P)),
                .rgMixedTiedCell,
                ud = ud,
                own = own,
                shapes = shapes
            ),
            P
        ),
        matrix(0, P, P)
    )
}

# Tied noise per component between conditions k and l, (Ph_kl - f
# Delta'Delta_kl) / (a - f K) with a and f the pair's n - 1 and 1 - phi as
# geometric means; none on the diagonal and between conditions that share no
# donors.
# @noRd
.rgMixedTiedCell <- function(k, l, ud, own, shapes) {
    ph <- ud$phenotypicCorrelation[[k, l]]
    if (k == l || is.na(ph)) {
        return(0)
    }
    C <- shapes[[ud$ldOf[[k]], ud$ldOf[[l]]]]$C
    raw <- sum(own[[k]]$delta * (C %*% own[[l]]$delta))
    a <- sqrt((own[[k]]$n - 1) * (own[[l]]$n - 1))
    f <- 1 - sqrt(own[[k]]$phi * own[[l]]$phi)
    (ph - f * raw) / (a - f * own[[k]]$K)
}

# The group's H: each condition's heritability K Omega_kk on the diagonal,
# and between conditions k and l on references a and b the cross-projection
# delta_k' C_ab delta_l less the held noise's part (sigma, per component)
# times the shapes' scale (.rgCrossMoments()); on one reference C_aa = I and
# the scale is 1, LAVA's K Omega_kl.
# @noRd
.rgMixedH <- function(own, basisOf, shapes, sigma) {
    P <- length(own)
    cells <- expand.grid(k = seq_len(P), l = seq_len(P))
    H <- matrix(
        pmap_dbl(
            cells,
            .rgMixedHEntry,
            own = own,
            shapes = shapes,
            basisOf = basisOf,
            sigma = sigma
        ),
        P
    )
    list(H = H, own = own, shapes = shapes, basisOf = basisOf, sigma = sigma)
}

# Every two references' moment shapes, as a matrix list, scaled as the
# settings' effects read the covariance.
# @noRd
.rgMixedMomentShapes <- function(bases, effects) {
    grid <- expand.grid(a = seq_along(bases), b = seq_along(bases))
    matrix(
        pmap(grid, .rgMixedMomentShape, bases = bases, effects = effects),
        length(bases)
    )
}

# References a and b's overlap on the components (.rgCrossPieces()): C, the
# genetic cross term (sqrt(lambda_a lambda_b') o C without frequency
# weights), the shared noise's shape and its trace against C; and the scale
# that puts a cross-projection on the covariance's scale: under fixed effects
# 1, the components' cross-products themselves; under random effects, the
# kept components' share of the per-SD genetic covariance, sqrt(tr Lambda_a
# tr Lambda_b) / sum(C o genetic) (1 on one reference).
# @noRd
.rgMixedMomentShape <- function(a, b, bases, effects) {
    pieces <- .rgCrossPieces(bases[[a]], bases[[b]])
    traces <- sqrt(sum(bases[[a]]$lam) * sum(bases[[b]]$lam))
    list(
        C = pieces$C,
        genetic = pieces$genetic,
        noiseShape = pieces$noise,
        noiseTrace = sum(pieces$C * pieces$noise),
        perM = 1 / traces,
        scale = if (identical(effects, "fixed")) {
            1
        } else {
            traces / sum(pieces$C * pieces$genetic)
        }
    )
}

# @noRd
.rgMixedHEntry <- function(k, l, own, shapes, basisOf, sigma) {
    if (k == l) {
        return(own[[k]]$h2)
    }
    shape <- shapes[[basisOf[[k]], basisOf[[l]]]]
    raw <- sum(own[[k]]$delta * (shape$C %*% own[[l]]$delta))
    (raw - sigma[[k, l]] * shape$noiseTrace) * shape$scale
}

# One pair of a mixed group given the group's other conditions: the partial
# covariance and variances from H, rg their correlation, each condition's r2,
# and the p-value and 95% intervals from draws of the components.
# @noRd
.rgMixedConditionalPair <- function(i, j, group, M) {
    H <- group$H
    rest <- setdiff(seq_len(nrow(H)), c(i, j))
    partial <- .rgPartial(H, i, j, rest)
    drawn <- .rgMixedDrawn(group, i, j, rest, partial$cov)
    .rgPairResult(
        i = i,
        j = j,
        nSnps = M,
        nComponents = min(map_dbl(group$own, "K")),
        rg = .rgCorrelation(partial$cov, partial$var),
        rgLower = drawn$rg[[1L]],
        rgUpper = drawn$rg[[2L]],
        covariance = partial$cov,
        h2 = partial$var,
        noiseCorrelation = .rgSigmaCorrelation(group, i, j),
        noiseSource = .rgPairNoiseSource(group$mode, group$tie, i, j),
        p = drawn$p,
        r2 = c(.rgExplained(H, i), .rgExplained(H, j)),
        r2Lower = drawn$r2Lower,
        r2Upper = drawn$r2Upper
    )
}

# A pair's shared noise per component as a correlation of the two
# conditions' own, NA where either has none.
# @noRd
.rgSigmaCorrelation <- function(group, i, j) {
    own <- group$own[[i]]$noise * group$own[[j]]$noise
    if (own > 0) group$sigma[[i, j]] / sqrt(own) else NA_real_
}

# The pair's p-value of no partial covariance and its 95% intervals, from
# Gaussian draws of every condition's components: under the null, H with the
# pair's entry H_iZ H_ZZ^-1 H_Zj (no partial covariance, LAVA's Omega_0),
# the share of draws whose partial covariance is at least the observed in
# size (10,000 draws, 100,000 below 1e-3); and at the estimate, the 2.5% and
# 97.5% quantiles of the draws' partial rg and r2s. NA where H does not give
# the components a covariance matrix.
# @noRd
.rgMixedDrawn <- function(group, i, j, rest, observed) {
    none <- c(NA_real_, NA_real_)
    empty <- list(p = NA_real_, rg = none, r2Lower = none, r2Upper = none)
    if (is.na(observed)) {
        return(empty)
    }
    null <- .rgMixedNullH(group$H, i, j, rest)
    p <- .rgMixedDrawP(group, null, i, j, rest, observed, 1e4)
    stats <- .rgMixedDrawStats(group, group$H, i, j, rest, 1e4)
    if (is.null(stats) || is.na(p)) {
        return(list_assign(empty, p = p))
    }
    r2i <- quantile(stats$r2i, c(0.025, 0.975), na.rm = TRUE)
    r2j <- quantile(stats$r2j, c(0.025, 0.975), na.rm = TRUE)
    list(
        p = if (p < 1e-3) {
            .rgMixedDrawP(group, null, i, j, rest, observed, 1e5)
        } else {
            p
        },
        rg = unname(quantile(stats$rg, c(0.025, 0.975), na.rm = TRUE)),
        r2Lower = unname(c(r2i[[1L]], r2j[[1L]])),
        r2Upper = unname(c(r2i[[2L]], r2j[[2L]]))
    )
}

# H with the pair's partial covariance given `rest` removed.
# @noRd
.rgMixedNullH <- function(H, i, j, rest) {
    implied <- if (length(rest) == 0L) {
        0
    } else {
        as.vector(
            H[i, rest, drop = FALSE] %*%
                solve(H[rest, rest, drop = FALSE], H[rest, j])
        )
    }
    replace(H, rbind(c(i, j), c(j, i)), implied)
}

# The share of `nDraws` draws under `H` whose partial covariance is at least
# `observed` in size, counted with the observed draw itself; NA where `H`
# gives the components no covariance matrix.
# @noRd
.rgMixedDrawP <- function(group, H, i, j, rest, observed, nDraws) {
    stats <- .rgMixedDrawStats(group, H, i, j, rest, nDraws)
    if (is.null(stats)) {
        return(NA_real_)
    }
    (1 + sum(abs(stats$cov) >= abs(observed), na.rm = TRUE)) / (1 + nDraws)
}

# `nDraws` draws of every condition's components at `H` (.rgMixedDrawSources()),
# each turned into the group's H as the estimate was, then the pair's partial
# covariance and rg and the two r2s; NULL where `H` gives no covariance
# matrix.
# @noRd
.rgMixedDrawStats <- function(group, H, i, j, rest, nDraws) {
    sources <- .rgMixedDrawSources(group, H)
    if (is.null(sources)) {
        return(NULL)
    }
    parts <- list_c(map(
        sources,
        .rgMixedSourceDraws,
        nDraws = nDraws / length(sources),
        group = group,
        i = i,
        j = j,
        rest = rest
    ))
    list(
        cov = list_c(map(parts, "cov")),
        rg = list_c(map(parts, "rg")),
        r2i = list_c(map(parts, "r2i")),
        r2j = list_c(map(parts, "r2j"))
    )
}

# Where a group's components are drawn from at `H`: under random effects,
# about no mean with the covariance H implies (.rgMixedStacked()); under
# fixed effects, about genetic means with H's cross-products
# (.rgFixedMeans()) with the noise's covariance, half the draws building the
# means in the conditions' order and half in reverse. NULL where the
# covariance is not one.
# @noRd
.rgMixedDrawSources <- function(group, H) {
    L <- try_fetch(chol(.rgMixedStacked(group, H)), error = .rgSolveFailed)
    if (is.null(L)) {
        return(NULL)
    }
    if (!identical(group$effects, "fixed")) {
        return(list(list(L = L, mean = 0)))
    }
    P <- nrow(H)
    map(
        list(seq_len(P), rev(seq_len(P))),
        .rgFixedSource,
        L = L,
        H = H,
        group = group
    )
}

# @noRd
.rgFixedSource <- function(order, L, H, group) {
    means <- .rgFixedMeans(H, group$own, group$shapes, group$basisOf, order)
    list(L = L, mean = list_c(means))
}

# Genetic means for a group's components under fixed effects whose
# cross-products are H's: built in `order`, each condition's in its own
# basis (.rgFixedMeanStep()); a list of each condition's.
# @noRd
.rgFixedMeans <- function(H, own, shapes, basisOf, order) {
    reduce(
        order,
        .rgFixedMeanStep,
        .init = vector("list", length(own)),
        H = H,
        own = own,
        shapes = shapes,
        basisOf = basisOf,
        order = order
    )
}

# Condition k's mean given those built before it in `order`: within the span
# of their images in its basis (C_ab mu_l), the shortest combination whose
# cross-products with them are H's; outside it, its observed components'
# direction, at the length H_kk leaves. Where the combination alone is
# longer than H_kk allows, it is shortened to it, and the cross-products
# fall short.
# @noRd
.rgFixedMeanStep <- function(mus, k, H, own, shapes, basisOf, order) {
    earlier <- order[seq_len(match(k, order) - 1L)]
    images <- map(
        earlier,
        .rgFixedImage,
        k = k,
        mus = mus,
        shapes = shapes,
        basisOf = basisOf
    )
    span <- .rgSpanFit(exec(cbind, !!!images), H[earlier, k], own[[k]]$delta)
    total <- max(H[[k, k]], 0)
    inSpan <- sum(span$inSpan^2)
    mu <- if (inSpan >= total) {
        span$inSpan * sqrt(total / max(inSpan, .Machine$double.xmin))
    } else {
        span$inSpan + sqrt(total - inSpan) * span$direction
    }
    replace(mus, k, list(as.vector(mu)))
}

# Condition j's mean seen from condition k's basis.
# @noRd
.rgFixedImage <- function(j, k, mus, shapes, basisOf) {
    as.vector(shapes[[basisOf[[k]], basisOf[[j]]]]$C %*% mus[[j]])
}

# The shortest combination of W's columns with cross-products `t` with them
# (none where W has none), and the direction of `delta`'s part outside W's
# span (none where it has no such part).
# @noRd
.rgSpanFit <- function(W, t, delta) {
    if (is.null(W)) {
        return(list(inSpan = 0 * delta, direction = .rgUnitVector(delta)))
    }
    dec <- svd(W)
    keep <- dec$d > max(dec$d) * 1e-10
    A <- dec$u[, keep, drop = FALSE]
    inSpan <- A %*% (crossprod(dec$v[, keep, drop = FALSE], t) / dec$d[keep])
    rest <- delta - A %*% crossprod(A, delta)
    list(inSpan = as.vector(inSpan), direction = .rgUnitVector(as.vector(rest)))
}

# `x` scaled to length 1, or none where it has no length.
# @noRd
.rgUnitVector <- function(x) {
    size <- sqrt(sum(x^2))
    if (size > 0) x / size else 0 * x
}

# One source's draws, 10,000 at a time.
# @noRd
.rgMixedSourceDraws <- function(source, nDraws, group, i, j, rest) {
    sizes <- diff(unique(c(seq(0, nDraws, by = 1e4), nDraws)))
    map(
        sizes,
        .rgMixedDrawChunk,
        source = source,
        group = group,
        i = i,
        j = j,
        rest = rest
    )
}

# The components' joint covariance at `H`: under random effects, each
# condition's variances (.rgSideComponentVariance()) on its diagonal block,
# and between k and l on references a and b, (H_kl / sqrt(tr Lambda_a tr
# Lambda_b)) times sqrt(lambda_a lambda_b') o C_ab; under fixed effects the
# noise's alone, each condition's noise per component and none between.
# Either way the held shared noise between, sigma_kl times its shape.
# @noRd
.rgMixedStacked <- function(group, H) {
    P <- nrow(H)
    rowsOf <- map(seq_len(P), .rgMixedStackedRow, group = group, H = H)
    exec(rbind, !!!rowsOf)
}

# @noRd
.rgMixedStackedRow <- function(k, group, H) {
    blocks <- map(
        seq_len(nrow(H)),
        .rgMixedStackedBlock,
        k = k,
        group = group,
        H = H
    )
    exec(cbind, !!!blocks)
}

# @noRd
.rgMixedStackedBlock <- function(l, k, group, H) {
    fixed <- identical(group$effects, "fixed")
    if (k == l) {
        own <- group$own[[k]]
        variance <- if (fixed) {
            rep(own$noise, own$K)
        } else {
            .rgSideComponentVariance(own)
        }
        return(diag(variance, length(variance)))
    }
    shape <- group$shapes[[group$basisOf[[k]], group$basisOf[[l]]]]
    genetic <- if (fixed) 0 else H[[k, l]] * shape$perM * shape$genetic
    genetic + group$sigma[[k, l]] * shape$noiseShape
}

# One chunk of draws' statistics: every draw's H cell by cell over the chunk
# (.rgMixedDrawCells()), then each draw's pair.
# @noRd
#' @importFrom stats rnorm
.rgMixedDrawChunk <- function(size, source, group, i, j, rest) {
    L <- source$L
    draws <- source$mean +
        crossprod(L, matrix(rnorm(nrow(L) * size), nrow(L)))
    K <- map_dbl(group$own, "K")
    ends <- cumsum(K)
    deltas <- map2(ends - K + 1, ends, .rgDrawRows, draws = draws)
    cells <- .rgMixedDrawCells(deltas, group)
    stats <- map(
        seq_len(size),
        .rgMixedDrawPairAt,
        cells = cells,
        P = length(K),
        i = i,
        j = j,
        rest = rest
    )
    list(
        cov = map_dbl(stats, "cov"),
        rg = map_dbl(stats, "rg"),
        r2i = map_dbl(stats, "r2i"),
        r2j = map_dbl(stats, "r2j")
    )
}

# The draws' H, one row per cell (column-major) and one column per draw: each
# condition's heritability from its components' squares with its noise per
# component held at the estimate's, and each pair's cross-projection less the
# held noise's part, times the shapes' scale (.rgMixedHEntry()'s).
# @noRd
.rgMixedDrawCells <- function(deltas, group) {
    P <- length(deltas)
    values <- pmap(
        expand.grid(k = seq_len(P), l = seq_len(P)),
        .rgMixedDrawCell,
        deltas = deltas,
        group = group
    )
    exec(rbind, !!!values)
}

# @noRd
.rgMixedDrawCell <- function(k, l, deltas, group) {
    if (k == l) {
        own <- group$own[[k]]
        return(colSums(deltas[[k]]^2) - own$K * own$noise)
    }
    shape <- group$shapes[[group$basisOf[[k]], group$basisOf[[l]]]]
    raw <- colSums(deltas[[k]] * (shape$C %*% deltas[[l]]))
    (raw - group$sigma[[k, l]] * shape$noiseTrace) * shape$scale
}

# Draw d's pair from the cells.
# @noRd
.rgMixedDrawPairAt <- function(d, cells, P, i, j, rest) {
    .rgMixedDrawPair(matrix(cells[, d], P), i, j, rest)
}

# @noRd
.rgDrawRows <- function(from, to, draws) {
    draws[from:to, , drop = FALSE]
}

# One draw's partial covariance and rg for the pair and the two r2s.
# @noRd
.rgMixedDrawPair <- function(H, i, j, rest) {
    partial <- .rgPartial(H, i, j, rest)
    list(
        cov = partial$cov,
        rg = .rgCorrelation(partial$cov, partial$var),
        r2i = .rgExplained(H, i),
        r2j = .rgExplained(H, j)
    )
}

# ---- the regressions and the GLS on two references --------------------------

# A pair's variant regression on two references: each condition's n r^2 on
# its own LD scores, as on one reference, then the cross-products sqrt(n1 n2)
# r1 r2 on the cross-population LD scores sum_k r1_jk r2_jk and, for the
# shared noise, an intercept (the diagonal of (R1 + R2) / 2 being 1), as the
# settings take it; a cross-product's LD-redundancy weight is the root of the
# two conditions'.
# @noRd
.rgTwoLdVariantRegression <- function(param, ud) {
    M <- length(ud$z[[1L]])
    if (M < 3L) {
        return(.rgPairResult(nSnps = M))
    }
    n <- map_dbl(ud$nVariant, .h2SampleSize, how = param@sampleSize)
    r <- map2(ud$z, ud$nVariant, .zToCor)
    lds <- map(seq_len(2L), .rgMemberLd, ud = ud)
    l <- map(lds, .h2WindowLdScores, param = param, M = M)
    redundancy <- map2(l, lds, .rgRedundancy, param = param)
    rows <- pmap(
        list(k = seq_len(2L), l = l, w = redundancy),
        .rgVariantRows,
        r = r,
        n = n
    )
    nn <- sqrt(n[[1L]] * n[[2L]])
    cross <- list(
        y = nn * r[[1L]] * r[[2L]],
        X = cbind(
            nn * .rgCrossLdScores(lds, M, .rgCrossWeights(param, lds)) / M,
            1
        ),
        taper = rep(1, M),
        redundancy = sqrt(redundancy[[1L]] * redundancy[[2L]]),
        ldModel = .rgVariantModel(param, lds, n)
    )
    fit <- .rgRegressionPair(param, ud, rows, cross)
    .rgRegressionResult(param, fit, M, NA_integer_)
}

# The root of two references' frequency weights, sqrt(w1 w2) (each 1 at alpha
# = -1).
# @noRd
.rgCrossWeights <- function(param, lds) {
    weights <- map(lds, .h2SourceWeights, param = param)
    if (is.null(weights[[1L]])) 1 else sqrt(weights[[1L]] * weights[[2L]])
}

# The cross-population LD scores of two conditions' references, sum_k r1_jk
# r2_jk w_k, `weights` the root of the two references' frequency weights
# (each 1 without them).
# @noRd
.rgCrossLdScores <- function(lds, M, weights) {
    as.vector(.ldSquaredProduct(lds[[1L]], matrix(weights, M, 1L), lds[[2L]]))
}

# A pair's eigen regression on two references: each condition's whitened
# squares on its own reference, as on one, then the reduced pairs'
# cross-products a_i b_i (.rgReduce() at the step-1 fits) on s_i and, for the
# shared noise, the reduced noise shape's diagonal over the root of the two
# sample sizes, weighted 1 / (1 + mu_i^2) as the procedures update them; the
# standard errors from the three regressions' joint sandwich
# (.rgTwoLdSandwich()).
# @noRd
.rgTwoLdEigenRegression <- function(param, ud) {
    M <- length(ud$z[[1L]])
    sides <- map(seq_len(2L), .rgSideEigen, param = param, ud = ud)
    if (any(map_lgl(sides, is.null)) || min(lengths(map(sides, "lam"))) < 3L) {
        return(.rgPairResult(nSnps = M))
    }
    proc <- param@procedures
    rows <- map(sides, .rgSideRows, proc = proc)
    kept <- map(rows, .rgChi2KeptRows, proc = proc)
    steps <- map2(kept, ud$heldIntercept, .rgRegressionStep1, param = param)
    v <- map2(
        sides,
        map2(rows, steps, .rgStepMean, proc = proc),
        .rgSideVariance
    )
    red <- .rgReduce(sides, v, 1e-12, NA_real_)
    crossScale <- 1 / sqrt(prod(map_dbl(sides, "n")))
    root <- .rgRootProduct(map_dbl(steps, .rgStepNoise))
    fit <- list(
        steps = steps,
        rows = kept,
        cross = .rgReducedIrwls(
            param,
            red,
            crossScale,
            .rgHeldCross(param, ud, root)
        ),
        red = red,
        mu = map2(kept, steps, .rgStepMean, proc = proc),
        lam = map(sides, "lam"),
        keep = map(rows, .rgChi2Keep, proc = proc),
        M = M,
        crossScale = crossScale,
        root = root,
        mode = param@crossNoise,
        ud = ud
    )
    .rgTwoLdRegressionResult(param, fit, M)
}

# A side's eigen-regression rows.
# @noRd
.rgSideRows <- function(side, proc) {
    .rgEigenSquares(
        side$u,
        side$lam,
        side$M,
        side$nRef,
        side$n,
        .h2Taper(proc@eigenvalueTaper, side$lam)
    )
}

# The rows a side's chi-square filter keeps, and which.
# @noRd
.rgChi2KeptRows <- function(rows, proc) {
    .h2SubsetRows(rows, .rgChi2Keep(rows, proc))
}

# @noRd
.rgChi2Keep <- function(rows, proc) {
    .h2Chi2Kept(proc, rows$chi2, rows$n)
}

# A side's modelled variance per direction at its step-1 fit, lambda mu.
# @noRd
.rgSideVariance <- function(side, mu) {
    side$lam * mu
}

# The reduced pairs' cross-products a_i b_i regressed on (s_i, crossScale
# Nred_ii) with no intercept by the regressions' engine
# (.h2SolveRegression()), the shared noise's column fitted unless `held`
# holds it (NA: fitted): the pairs are whitened, so each weight is the
# inverse of a product's variance, 1 + mu_i^2, at the current fit, from no
# covariance, as the procedures update them. They carry no eigenvalue for
# the eigenvalue taper to read, s_i being the cross-shape's singular values,
# and are of one variance; each side's own regression reads the taper.
# @noRd
.rgReducedIrwls <- function(param, red, crossScale, held) {
    n <- length(red$a)
    rows <- list(
        y = red$a * red$b,
        X = cbind(red$s, crossScale * diag(red$Nred)),
        taper = rep(1, n),
        redundancy = rep(1, n),
        weights = "cross",
        mu11 = rep(1, n),
        mu22 = rep(1, n)
    )
    .h2SolveRegression(
        param,
        rows,
        if (is.na(held)) NULL else held,
        tied = FALSE
    )
}

# A two-reference eigen regression's values: the covariance and the shared
# noise it fits, each condition's heritability from step 1, the sandwich's
# standard errors where the settings ask (rg's interval Wald's, its test of
# no covariance too).
# @noRd
.rgTwoLdRegressionResult <- function(param, fit, M) {
    h2 <- map_dbl(fit$steps, .rgStepH2)
    covariance <- fit$cross$coef[[1L]]
    rg <- .rgCorrelation(covariance, h2)
    se <- if (identical(param@inference, "information")) {
        .rgDeltaSe(.rgTwoLdSandwich(fit), h2, rg)
    } else {
        list(covariance = NA_real_, rg = NA_real_)
    }
    .rgPairResult(
        nSnps = M,
        nComponents = length(fit$red$s),
        h2 = h2,
        rg = rg,
        rgSe = se$rg,
        rgLower = rg - qnorm(0.975) * se$rg,
        rgUpper = rg + qnorm(0.975) * se$rg,
        covariance = covariance,
        covarianceSe = se$covariance,
        noiseCorrelation = .rgReportedNoise(
            fit$mode,
            fit$cross$coef[[2L]],
            fit$root,
            fit$ud
        ),
        noiseSource = fit$mode,
        p = .rgWaldP(covariance, se$covariance)
    )
}

# The three regressions' joint sandwich on two references: each one's bread,
# and a meat of their responses' Gaussian covariances under the model
# (.rgTwoLdCovariances()). The cross fit's dependence on step 1 through the
# whitening has mean zero at the model, so the bread is block-diagonal.
# @noRd
.rgTwoLdSandwich <- function(fit) {
    X <- list(
        fit$rows[[1L]]$X,
        fit$rows[[2L]]$X,
        fit$cross$design$X
    )
    w <- list(fit$steps[[1L]]$w, fit$steps[[2L]]$w, fit$cross$w)
    bread <- .rgBlockDiagonal(map2(X, w, .rgBread))
    meat <- .rgTwoLdMeat(X, w, .rgTwoLdCovariances(fit))
    list(
        cov = solve(bread, t(solve(bread, meat))),
        at = cumsum(c(1L, map_int(X[seq_len(2L)], ncol)))
    )
}

# The responses' covariance blocks under the model at the fits, on the rows
# each step kept. With X = g S + c Nred the reduced pairs' cross-covariance
# and G12 = g B + c N the rotated statistics' (c the shared noise on the
# statistics' scale): Cov(y1_i, y2_j) = 2 G12_ij^2 / (lambda1_i lambda2_j),
# Cov(y1_i, y12_k) = 2 mu1_i U_ik (U X)_ik, Cov(y2_j, y12_k) = 2 mu2_j V_jk (V
# X')_jk, Cov(y12_k, y12_l) = [k = l] + X_kl X_lk, each side's squares 2 mu^2
# apart.
# @noRd
.rgTwoLdCovariances <- function(fit) {
    g <- fit$cross$coef[[1L]]
    c12 <- fit$cross$coef[[2L]] * fit$crossScale
    red <- fit$red
    X <- g * diag(red$s, length(red$s)) + c12 * red$Nred
    keep1 <- fit$keep[[1L]]
    keep2 <- fit$keep[[2L]]
    U <- red$U[keep1, , drop = FALSE]
    V <- red$V[keep2, , drop = FALSE]
    G12 <- (g * red$B + c12 * red$N)[keep1, keep2, drop = FALSE]
    lam <- map2(fit$lam, fit$keep, `[`)
    list(
        c11 = 2 * fit$mu[[1L]]^2,
        c22 = 2 * fit$mu[[2L]]^2,
        c33 = diag(length(red$s)) + X * t(X),
        c12 = 2 * G12^2 / tcrossprod(lam[[1L]], lam[[2L]]),
        c13 = 2 * fit$mu[[1L]] * U * (U %*% X),
        c23 = 2 * fit$mu[[2L]] * V * (V %*% t(X))
    )
}

# The stacked meat of the three regressions on two references: each side's
# diagonal block X' diag(w^2 c) X, and every other X_a' (w_a o C_ab o w_b')
# X_b.
# @noRd
.rgTwoLdMeat <- function(X, w, cv) {
    m11 <- .rgMeatDiagonal(X[[1L]], w[[1L]], cv$c11)
    m22 <- .rgMeatDiagonal(X[[2L]], w[[2L]], cv$c22)
    m33 <- .rgMeatCross(X[[3L]], X[[3L]], w[[3L]], w[[3L]], cv$c33)
    m12 <- .rgMeatCross(X[[1L]], X[[2L]], w[[1L]], w[[2L]], cv$c12)
    m13 <- .rgMeatCross(X[[1L]], X[[3L]], w[[1L]], w[[3L]], cv$c13)
    m23 <- .rgMeatCross(X[[2L]], X[[3L]], w[[2L]], w[[3L]], cv$c23)
    rbind(
        cbind(m11, m12, m13),
        cbind(t(m12), m22, m23),
        cbind(t(m13), t(m23), m33)
    )
}

# @noRd
.rgMeatDiagonal <- function(X, w, c) {
    crossprod(X, X * (w^2 * c))
}

# @noRd
.rgMeatCross <- function(Xa, Xb, wa, wb, cab) {
    crossprod(Xa, (wa * cab * rep(wb, each = length(wa))) %*% Xb)
}

# A pair's variant GLS on two references: each condition's fit on its own LD
# (.rgGlsStep1()), then the cross-products sqrt(n1 n2) r1 r2 as on one
# reference (.rgGlsCrossFit()), on the cross-population LD scores, S12's
# genetic term t12 R1 R2 and its shared noise's shape (R1 + R2) / 2: the LD
# of the samples the two share, taken as the average of the two references'.
# @noRd
.rgTwoLdGls <- function(param, ud) {
    M <- length(ud$z[[1L]])
    if (M < 3L) {
        return(.rgPairResult(nSnps = M))
    }
    lds <- map(seq_len(2L), .rgMemberLd, ud = ud)
    gs <- map(lds, .rgGlsSideData, param = param, ud = ud)
    steps <- map2(seq_len(2L), gs, .rgGlsStep1, held = ud$heldIntercept)
    base <- list(nSnps = M, h2 = map_dbl(steps, "h2"))
    roots <- sqrt(gs[[1L]]$weights * gs[[2L]]$weights)
    crossG <- list_assign(
        gs[[1L]],
        l = .rgCrossLdScores(lds, M, roots),
        R2 = gs[[1L]]$R %*% (roots * gs[[2L]]$R),
        R = (gs[[1L]]$R + gs[[2L]]$R) / 2
    )
    cross <- .rgGlsCrossFit(param, ud, crossG, steps)
    if (is.null(cross$coef)) {
        return(exec(.rgPairResult, !!!base))
    }
    .rgGlsResult(param, ud, crossG, steps, cross, base)
}

# One side's GLS data, read on its own reference.
# @noRd
.rgGlsSideData <- function(ld, param, ud) {
    .rgGlsData(param, list_assign(ud, ld = ld))
}

# ---- the moments on two references -------------------------------------------

# A pair's moments on two references: each condition's whitened components
# delta_k = Lambda_k^-1/2 Q_k'r_k on its own reference's K_k and its
# heritability K_k Omega_kk as on one, and the covariance from the
# components' cross-projection delta_1' C delta_2, C = Q_1'Q_2, less the
# supplied shared noise's part, as the settings' effects read it: the pair as
# a group of two (.rgMomentsGroupOf()), its inference .rgGroupPair()'s.
# @noRd
.rgTwoLdMoments <- function(param, ud) {
    M <- length(ud$z[[1L]])
    sides <- map(seq_len(2L), .rgSideEigen, param = param, ud = ud)
    if (
        any(map_lgl(sides, is.null)) ||
            any(map_dbl(sides, "n") <= lengths(map(sides, "lam")) + 1L)
    ) {
        return(.rgPairResult(nSnps = M))
    }
    group <- .rgMomentsGroupOf(
        param,
        list_assign(ud, ldOf = c(1L, 2L)),
        sides,
        sides
    )
    .rgGroupPair(1L, 2L, group, M)
}

# A group's moments as its fits and draws read them: each condition's
# components, heritability and noise per component (.rgSideMoments()), every
# two references' shapes as the settings' effects scale them, the shared
# noise the settings hold (.rgMixedSigma()), and H (.rgMixedH()); with the
# settings' shared noise, its tie and their effects.
# @noRd
.rgMomentsGroupOf <- function(param, ud, sides, bases) {
    own <- pmap(
        list(sides, ud$overlap, .rgMomentsIntercepts(ud)),
        .rgSideMoments,
        param = param
    )
    shapes <- .rgMixedMomentShapes(bases, param@effects)
    sigma <- .rgMixedSigma(param, ud, own, shapes)
    list_assign(
        .rgMixedH(own, ud$ldOf, shapes, sigma),
        mode = param@crossNoise,
        tie = .rgTieOf(param, ud, length(own)),
        effects = param@effects
    )
}

# One pair of a group on its own, from the group's H: the covariance and the
# two heritabilities, rg, and the p-value and interval as the group's effects
# read them (.rgGroupPairInference()).
# @noRd
.rgGroupPair <- function(i, j, group, M) {
    H <- group$H
    h2 <- c(H[[i, i]], H[[j, j]])
    inference <- .rgGroupPairInference(.rgSubGroup(group, c(i, j)))
    .rgPairResult(
        i = i,
        j = j,
        nSnps = M,
        nComponents = min(group$own[[i]]$K, group$own[[j]]$K),
        rg = .rgCorrelation(H[[i, j]], h2),
        rgLower = inference$interval[[1L]],
        rgUpper = inference$interval[[2L]],
        covariance = H[[i, j]],
        h2 = h2,
        noiseCorrelation = .rgSigmaCorrelation(group, i, j),
        noiseSource = .rgPairNoiseSource(group$mode, group$tie, i, j),
        p = inference$p
    )
}

# A group's conditions `at` as a group of their own.
# @noRd
.rgSubGroup <- function(group, at) {
    list_assign(
        group,
        H = group$H[at, at],
        own = group$own[at],
        basisOf = group$basisOf[at],
        sigma = group$sigma[at, at],
        tie = group$tie[at, at]
    )
}

# A pair's p-value of no covariance and its 95% interval for rg. Under random
# effects, from the components' Gaussian distribution with their genetic
# variance (.rgTwoLdMomentsP(), .rgTwoLdMomentsInterval()); under fixed
# effects, with the genetic components held at means (.rgFixedPairP(), and
# the draws' quantiles at the estimate).
# @noRd
.rgGroupPairInference <- function(pair) {
    if (identical(pair$effects, "random")) {
        shape <- pair$shapes[[pair$basisOf[[1L]], pair$basisOf[[2L]]]]
        cross <- .rgCrossMoments(shape, pair$own, pair$sigma[[1L, 2L]])
        return(list(
            p = .rgTwoLdMomentsP(cross, pair$own),
            interval = .rgTwoLdMomentsInterval(cross, pair$own)
        ))
    }
    stats <- if (all(diag(pair$H) > 0)) {
        .rgMixedDrawStats(pair, pair$H, 1L, 2L, integer(0), 1e4)
    }
    list(
        p = .rgFixedPairP(pair),
        interval = if (is.null(stats)) {
            c(NA_real_, NA_real_)
        } else {
            unname(quantile(stats$rg, c(0.025, 0.975), na.rm = TRUE))
        }
    )
}

# A side's components, heritability and noise per component, at its null
# intercept (.rgMomentsIntercepts()).
# @noRd
.rgSideMoments <- function(side, overlap, intercept, param) {
    K <- length(side$lam)
    delta <- side$u / sqrt(side$lam)
    phi <- if (identical(param@referenceNoise, "overlapFactor")) {
        .h2OverlapFactor(side$n, side$ld$nSamples, overlap)
    } else {
        1
    }
    own <- .rgMomentsOwn(sum(delta^2), K, side$n, phi, intercept)
    list(
        delta = delta,
        lam = side$lam,
        K = K,
        n = side$n,
        phi = phi,
        h2 = own$h2,
        noise = own$noise
    )
}

# The cross-projection of a pair's components on their references' shape
# (.rgMixedMomentShape()) and the covariance it reads: delta_1' C delta_2
# less the shared noise's part, sigma12 sum(C o Nd), with Nd the noise's
# shape on the components (the identity on one reference) and sigma12 the
# noise per component, times the shape's scale.
# @noRd
.rgCrossMoments <- function(shape, own, sigma12) {
    raw <- as.vector(crossprod(own[[1L]]$delta, shape$C %*% own[[2L]]$delta))
    offset <- sigma12 * shape$noiseTrace
    list(
        C = shape$C,
        genetic = shape$genetic,
        noiseShape = shape$noiseShape,
        sigma12 = sigma12,
        offset = offset,
        statistic = raw - offset,
        scale = shape$scale,
        covariance = (raw - offset) * shape$scale
    )
}

# A side's components' variances under the per-SD random-effects model at its
# heritability: (h2 / tr Lambda) lambda + its noise, per component.
# @noRd
.rgSideComponentVariance <- function(side) {
    max(side$h2, 0) / sum(side$lam) * side$lam + side$noise
}

# The cross-projection's p-value for no covariance. With no shared noise the
# two sides' components are independent, the statistic is sum_i rho_i x_i y_i
# over the singular values rho of V_1^1/2 C V_2^1/2 and standard normal x, y,
# and given y it is normal with variance sum rho_i^2 y_i^2, so each draw of y
# contributes its two-sided normal tail (partial integration, as LAVA's):
# 10,000 draws, then 100,000 below 1e-4 and 1,000,000 below 1e-6. With shared
# noise the components are drawn together under no covariance
# (.rgTwoLdDrawP()).
# @noRd
.rgTwoLdMomentsP <- function(cross, own) {
    v <- map(own, .rgSideComponentVariance)
    if (cross$sigma12 != 0) {
        return(.rgTwoLdDrawP(cross, own, v))
    }
    rho <- svd(
        sqrt(v[[1L]]) * cross$C * rep(sqrt(v[[2L]]), each = length(v[[1L]])),
        nu = 0L,
        nv = 0L
    )$d
    observed <- abs(cross$statistic)
    p <- .rgBilinearTail(rho, observed, 1e4, 0, 0)
    if (p >= 1e-4) {
        return(p)
    }
    p100k <- .rgBilinearTail(rho, observed, 1e5, 0, 0)
    if (p100k >= 1e-6) {
        return(p100k)
    }
    .rgBilinearTail(rho, observed, 1e6, 0, 0)
}

# The fixed-effects p-value of a pair's covariance. Under the null each
# condition's genetic components are held at means with its observed
# direction and heritability and no cross-projection with the other's
# (.rgFixedMeans()), and only the noise is drawn. With no shared noise, by
# partial integration: in the singular vectors of C = U S V', the statistic
# is sum_i s_i xi_i eta_i, normal given eta, the components' noise isotropic;
# half the draws hold the first condition's direction and draw the second,
# half the reverse. 10,000 draws, then 100,000 below 1e-4 and 1,000,000
# below 1e-6. With shared noise, the share of draws of both
# (.rgMixedDrawP()), 100,000 below 1e-3. NA unless both heritabilities are
# positive.
# @noRd
.rgFixedPairP <- function(pair) {
    if (any(diag(pair$H) <= 0)) {
        return(NA_real_)
    }
    null <- replace(pair$H, rbind(c(1L, 2L), c(2L, 1L)), 0)
    observed <- pair$H[[1L, 2L]]
    if (pair$sigma[[1L, 2L]] != 0) {
        drawn <- .rgMixedDrawP(pair, null, 1L, 2L, integer(0), observed, 1e4)
        if (is.na(drawn) || drawn >= 1e-3) {
            return(drawn)
        }
        return(.rgMixedDrawP(pair, null, 1L, 2L, integer(0), observed, 1e5))
    }
    tails <- .rgFixedTails(pair, null)
    p <- .rgFixedTail(tails, observed, 1e4)
    if (p >= 1e-4) {
        return(p)
    }
    p100k <- .rgFixedTail(tails, observed, 1e5)
    if (p100k >= 1e-6) {
        return(p100k)
    }
    .rgFixedTail(tails, observed, 1e6)
}

# The two partial integrations' terms: rho = s sqrt(sigma_1 sigma_2), and for
# each order the integrated side's means and the drawn side's, on the
# singular vectors, over their noise's root.
# @noRd
.rgFixedTails <- function(pair, null) {
    shape <- pair$shapes[[pair$basisOf[[1L]], pair$basisOf[[2L]]]]
    dec <- svd(shape$C)
    roots <- sqrt(map_dbl(pair$own, "noise"))
    orders <- list(forward = c(1L, 2L), reverse = c(2L, 1L))
    means <- map(
        orders,
        .rgFixedMeans,
        H = null,
        own = pair$own,
        shapes = pair$shapes,
        basisOf = pair$basisOf
    )
    c(
        list(rho = dec$d * prod(roots)),
        map(means, .rgFixedOnSingular, dec = dec, roots = roots)
    )
}

# A pair's means on C's singular vectors, over their noise's root.
# @noRd
.rgFixedOnSingular <- function(means, dec, roots) {
    list(
        first = as.vector(crossprod(dec$u, means[[1L]])) / roots[[1L]],
        second = as.vector(crossprod(dec$v, means[[2L]])) / roots[[2L]]
    )
}

# The mean conditional tail over `nDraws` draws, half each way: holding the
# first condition's direction, the second drawn; and the reverse.
# @noRd
.rgFixedTail <- function(tails, observed, nDraws) {
    half <- nDraws / 2
    mean(c(
        .rgBilinearTail(
            tails$rho,
            observed,
            half,
            tails$forward$first,
            tails$forward$second
        ),
        .rgBilinearTail(
            tails$rho,
            observed,
            half,
            tails$reverse$second,
            tails$reverse$first
        )
    ))
}

# With shared noise: the share of draws of both sides' components under no
# covariance (the noise's cross-covariance kept) whose statistic is at least
# the observed in size, counted with the observed draw itself; 10,000 draws,
# 100,000 below 1e-3.
# @noRd
.rgTwoLdDrawP <- function(cross, own, v) {
    L <- try_fetch(
        chol(.rgTwoLdJoint(v, cross$sigma12 * cross$noiseShape)),
        error = .rgSolveFailed
    )
    if (is.null(L)) {
        return(NA_real_)
    }
    p <- .rgTwoLdDrawShare(1e4, L, cross, own[[1L]]$K)
    if (p >= 1e-3) p else .rgTwoLdDrawShare(1e5, L, cross, own[[1L]]$K)
}

# @noRd
.rgTwoLdDrawShare <- function(nDraws, L, cross, K1) {
    sizes <- diff(unique(c(seq(0, nDraws, by = 1e4), nDraws)))
    stats <- list_c(map(
        sizes,
        .rgTwoLdDrawStatistic,
        L = L,
        cross = cross,
        K1 = K1
    ))
    (1 + sum(abs(stats) >= abs(cross$statistic))) / (1 + nDraws)
}

# `size` draws' statistics, the noise's part taken off as the estimate's is.
# @noRd
#' @importFrom stats rnorm
.rgTwoLdDrawStatistic <- function(size, L, cross, K1) {
    draws <- crossprod(L, matrix(rnorm(nrow(L) * size), nrow(L)))
    d1 <- draws[seq_len(K1), , drop = FALSE]
    d2 <- draws[-seq_len(K1), , drop = FALSE]
    colSums(d1 * (cross$C %*% d2)) - cross$offset
}

# The two sides' components' joint covariance: each one's variances, and the
# cross-covariance `cross12` between them.
# @noRd
.rgTwoLdJoint <- function(v, cross12) {
    rbind(
        cbind(diag(v[[1L]], length(v[[1L]])), cross12),
        cbind(t(cross12), diag(v[[2L]], length(v[[2L]])))
    )
}

# The mean conditional tail over `nDraws` draws, 100,000 at a time, of the
# statistic sum_i rho_i (a_i + x_i)(b_i + y_i) with standard normal x and y:
# given y it is normal, mean sum_i rho_i a_i (b_i + y_i) and variance sum_i
# rho_i^2 (b_i + y_i)^2.
# @noRd
.rgBilinearTail <- function(rho, observed, nDraws, a, b) {
    sizes <- diff(unique(c(seq(0, nDraws, by = 1e5), nDraws)))
    tails <- map_dbl(
        sizes,
        .rgBilinearTailChunk,
        rho = rho,
        observed = observed,
        a = a,
        b = b
    )
    sum(tails * sizes) / nDraws
}

# @noRd
#' @importFrom stats rnorm pnorm
.rgBilinearTailChunk <- function(size, rho, observed, a, b) {
    y <- matrix(rnorm(size * length(rho)), size) + rep(b, each = size)
    centre <- as.vector(y %*% (rho * a))
    sd <- sqrt(as.vector(y^2 %*% rho^2))
    mean(
        pnorm(-abs(observed), centre, sd) +
            pnorm(abs(observed), centre, sd, lower.tail = FALSE)
    )
}

# The 95% interval for rg: 10,000 draws of the two sides' components at the
# estimate (their variances, and the cross-covariance (G12 / M) Lambda_1^1/2
# C Lambda_2^1/2 the covariance implies plus the shared noise's), each turned
# into the two heritabilities and the covariance as the estimate was, and the
# 2.5% and 97.5% quantiles of the draws' rg. NA where the implied joint
# covariance is not positive definite.
# @noRd
#' @importFrom stats quantile
.rgTwoLdMomentsInterval <- function(cross, own) {
    v <- map(own, .rgSideComponentVariance)
    perM <- cross$covariance / sqrt(sum(own[[1L]]$lam) * sum(own[[2L]]$lam))
    cross12 <- perM * cross$genetic + cross$sigma12 * cross$noiseShape
    L <- try_fetch(chol(.rgTwoLdJoint(v, cross12)), error = .rgSolveFailed)
    if (is.null(L)) {
        return(c(NA_real_, NA_real_))
    }
    rg <- list_c(map(
        rep(1e3, 10L),
        .rgTwoLdDrawRg,
        L = L,
        cross = cross,
        own = own
    ))
    unname(quantile(rg, c(0.025, 0.975), na.rm = TRUE))
}

# `size` draws' rg at the estimate.
# @noRd
.rgTwoLdDrawRg <- function(size, L, cross, own) {
    K1 <- own[[1L]]$K
    draws <- crossprod(L, matrix(rnorm(nrow(L) * size), nrow(L)))
    d1 <- draws[seq_len(K1), , drop = FALSE]
    d2 <- draws[-seq_len(K1), , drop = FALSE]
    h1 <- map_dbl(colSums(d1^2), .rgDrawH2, side = own[[1L]])
    h2 <- map_dbl(colSums(d2^2), .rgDrawH2, side = own[[2L]])
    covariance <- (colSums(d1 * (cross$C %*% d2)) - cross$offset) * cross$scale
    if_else(h1 > 0 & h2 > 0, covariance / sqrt(pmax(h1 * h2, 0)), NA_real_)
}

# A draw's heritability from its components' sum of squares, the side's noise
# per component held at the estimate's.
# @noRd
.rgDrawH2 <- function(dtd, side) {
    dtd - side$K * side$noise
}

# The pair's values from the reduced likelihood: the covariance at its
# maximum by Fisher scoring (.h2Scoring()) from 0 and half its bound either
# way, rg its ratio to that bound, rg's 95% interval from inverting the
# likelihood ratio and the curve the RgEstimate keeps (.rgProfile()), and
# the ratio's p-value against no covariance; no shared noise. NA where no
# fit converges.
# @noRd
#' @importFrom stats pchisq
.rgReducedResult <- function(param, red, base) {
    bound <- red$bound
    fit <- .h2ScoringBest(.h2Scoring(
        "reduced",
        red,
        as.list(c(0, -bound / 2, bound / 2)),
        -bound,
        bound,
        param@procedures@maxit
    ))
    if (is.null(fit)) {
        return(exec(.rgPairResult, !!!base))
    }
    profile <- .rgProfile(
        .RG_CURVE_GRID * bound,
        list(estimate = fit$par, par = numeric(0), value = fit$value),
        .rgReducedRefit,
        red = red
    )
    lr <- 2 * (fit$value - .rgReducedLl(0, red))
    exec(
        .rgPairResult,
        !!!base,
        curve = profile$curve,
        rg = fit$par / bound,
        rgLower = profile$interval[[1L]] / bound,
        rgUpper = profile$interval[[2L]] / bound,
        covariance = fit$par,
        noiseCorrelation = 0,
        noiseSource = "zero",
        p = if (identical(param@inference, "lrt")) {
            pchisq(max(lr, 0), 1, lower.tail = FALSE)
        } else {
            NA_real_
        }
    )
}

# =============================================================================
# estimateLocalRg on a QtlDataset
# =============================================================================

#' @rdname estimateLocalRg
#' @param region,genotypeFilterParam,residualizationParam,phenotypeFilterParam
#'   For a \code{QtlDataset} (each study of a \code{MultiStudyQtlDataset}),
#'   as in \code{\link{estimateLocalH2}}: one region
#'   every selected trait is analysed over (otherwise each trait's window),
#'   and the genotype filters, covariate residualization and phenotype
#'   filters, as \code{\link{twasWeightsPipeline}} takes them.
#' @export
setMethod(
    "estimateLocalRg",
    "QtlDataset",
    function(
        data,
        gwas = NULL,
        param = NULL,
        jointSpecification = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        genotypeFilterParam = GenotypeFilterParam(),
        residualizationParam = ResidualizationParam(),
        phenotypeFilterParam = PhenotypeFilterParam(),
        noiseCorrelation = NULL,
        gwasH2Estimate = NULL,
        conditional = FALSE
    ) {
        settings <- .rgParam(
            param,
            gwas,
            noiseCorrelation,
            gwasH2Estimate,
            conditional,
            TRUE
        )
        .h2AssertOneRegion(region, cisWindow, "estimateLocalRg(QtlDataset)")
        filter <- .qtlResolveFilter(genotypeFilterParam)
        qd <- .qtlRestrictSamples(data, filter$keepSamples)
        run <- .rgRun(gwas, noiseCorrelation, gwasH2Estimate, conditional)
        args <- .rgDatasetArgs(
            region,
            cisWindow,
            residualizationParam,
            phenotypeFilterParam,
            jointSpecification
        )
        units <- .rgDatasetUnits(qd, settings, contexts, traitId, args, run)
        work <- .rgDatasetWork(units, settings, run)
        fits <- .qtlRunWork(
            work,
            qd,
            filter,
            residualizationParam,
            withAf = settings@alpha != -1
        )
        .rgEstimate(fits, .qtlGenotypeSe(qd), settings, run)
    }
)

# What a QtlDataset's units are built from: the one region or each trait's
# window, the residualization and the resolved phenotype filter, as
# estimateLocalH2() takes them, and the pairing.
# @noRd
.rgDatasetArgs <- function(
    region,
    cisWindow,
    residualizationParam,
    phenotypeFilterParam,
    jointSpecification
) {
    list(
        region = .makeXRegions(region, FALSE)[[1L]],
        cisWindow = cisWindow,
        residualizationParam = residualizationParam,
        phenotypeFilter = .qtlResolvePhenotypeFilter(
            phenotypeFilterParam,
            naDefault = "drop",
            naAllowed = c("drop", "impute"),
            caller = "estimateLocalRg"
        ),
        jointSpecification = jointSpecification
    )
}

# The units' fits as a plan: each unit's genotypes read and residualized as
# its spec says (one read per window), then its pairs fitted on them.
# @noRd
.rgDatasetWork <- function(units, param, run) {
    .qtlWork(
        units,
        map(units, .rgDatasetSpec),
        .rgDatasetUnitFit,
        list(param = param, run = run)
    )
}

# ---- a QtlDataset's units ----------------------------------------------------

# The units a QtlDataset's pairs are fitted in, as on summary statistics: the
# selected (context, trait) conditions grouped along each jointSpecification's
# axes (each trait's contexts by default), each group giving its pairs or,
# where the settings fit a group at once, one unit of all its conditions;
# with `gwas`, each condition with each GWAS study.
# @noRd
.rgDatasetUnits <- function(qd, param, contexts, traitId, args, run) {
    conditions <- .rgDatasetConditions(qd, contexts, traitId, args)
    if (!is.null(run$gwas) && !run$conditional) {
        .rgAssertGwasPairing(args$jointSpecification)
        return(list_c(map(
            conditions,
            .rgDatasetGwasUnits,
            qd = qd,
            args = args,
            gwas = run$gwas
        )))
    }
    specs <- parseJointSpecification(args$jointSpecification %||% "context", qd)
    units <- list_c(map(
        specs,
        .rgDatasetSpecUnits,
        qd = qd,
        conditions = conditions,
        atOnce = .rgFitsAtOnce(param),
        args = args,
        gwas = run$gwas
    ))
    if (length(units) == 0L) {
        abort(str_c(
            "estimateLocalRg(QtlDataset): no two selected conditions share a ",
            "group along the joint axes, so there is no pair to estimate."
        ))
    }
    units
}

# Every selected (context, trait) condition, as estimateLocalH2() selects
# them, with its residualized phenotype.
# @noRd
.rgDatasetConditions <- function(qd, contexts, traitId, args) {
    selected <- parseTraitIds(traitId, qd)
    useCtx <- .twasQdsResolveContexts(qd, parseContexts(contexts, qd))
    scope <- .twasQdsScope(
        studyName(qd),
        useCtx,
        .twasQdsResolveTraits(
            qd,
            useCtx,
            selected,
            args$region,
            "estimateLocalRg(QtlDataset)"
        ),
        if (.hasNestedTraitSelection(selected)) selected
    )
    list_flatten(map(
        scope$contexts[[studyName(qd)]],
        .rgDatasetContextConditions,
        qd = qd,
        scope = scope,
        args = args
    ))
}

# With `gwas`, each QTL condition pairs with each GWAS study, so a
# jointSpecification, which pairs QTL conditions, is refused.
# @noRd
.rgAssertGwasPairing <- function(jointSpecification) {
    if (is.null(jointSpecification)) {
        return(invisible(NULL))
    }
    abort(str_c(
        "estimateLocalRg: with `gwas`, each QTL row pairs with each GWAS ",
        "study; `jointSpecification` pairs QTL rows, so give one or the other."
    ))
}

# One context's selected conditions: each trait the scope keeps, with its
# phenotype residualized as estimateLocalH2() residualizes it.
# @noRd
.rgDatasetContextConditions <- function(cx, qd, scope, args) {
    map(
        .enumContextTraits(cx, qd, scope, studyName(qd)),
        .rgDatasetCondition,
        cx = cx,
        qd = qd,
        args = args
    )
}

# One condition: its identity and its residualized phenotype, named by the
# samples that have one.
# @noRd
.rgDatasetCondition <- function(tid, cx, qd, args) {
    Y <- .fmResidPheno(
        qd,
        contexts = cx,
        traitId = tid,
        phenotypeFilter = args$phenotypeFilter,
        residualizationParam = args$residualizationParam
    )
    list(
        study = studyName(qd),
        context = cx,
        trait = tid,
        y = set_names(as.vector(Y[, 1L]), rownames(Y))
    )
}

# One joint spec's units: its conditions in scope, grouped by the identity
# axes it leaves out; with `gwas`, a group of one condition too. A QtlDataset
# holds one study, so the study axis pairs nothing and is refused.
# @noRd
.rgDatasetSpecUnits <- function(spec, qd, conditions, atOnce, args, gwas) {
    if (is_in("study", spec$axes)) {
        abort(str_c(
            "estimateLocalRg(QtlDataset): a QtlDataset holds one study, so ",
            "the study axis pairs nothing; pair across studies on their ",
            "summary statistics (a QtlSumStats)."
        ))
    }
    inScope <- keep(conditions, .rgConditionInScope, scope = spec$scope)
    complement <- setdiff(c("context", "trait"), spec$axes)
    keys <- map_chr(inScope, .rgConditionLabel, cols = complement)
    groups <- unname(split(inScope, keys))
    list_c(map(
        groups[lengths(groups) >= (if (is.null(gwas)) 2L else 1L)],
        .rgDatasetGroupUnits,
        qd = qd,
        axes = intersect(c("context", "trait"), spec$axes),
        atOnce = atOnce,
        args = args,
        gwas = gwas
    ))
}

# Whether a condition lies in a spec's `scope` (every one without one).
# @noRd
.rgConditionInScope <- function(condition, scope) {
    all(map_lgl(
        names(scope),
        .rgConditionInScopeKey,
        condition = condition,
        scope = scope
    ))
}

# @noRd
.rgConditionInScopeKey <- function(key, condition, scope) {
    is_in(condition[[key]], scope[[key]])
}

# A condition's values on the identity columns `cols`, joined by "|".
# @noRd
.rgConditionLabel <- function(condition, cols) {
    str_flatten(as.character(unlist(condition[cols])), "|")
}

# A group's units: one of all its conditions where the settings fit a group
# at once, its pairs otherwise, each condition labelled by its values on the
# joint axes; with `gwas`, the group with each GWAS study, labelled by its
# study.
# @noRd
.rgDatasetGroupUnits <- function(group, qd, axes, atOnce, args, gwas) {
    labels <- map_chr(group, .rgConditionLabel, cols = axes)
    if (!is.null(gwas)) {
        return(map(
            unique(as.character(gwas$study)),
            .rgDatasetGroupGwasUnit,
            group = group,
            labels = labels,
            qd = qd,
            args = args
        ))
    }
    if (atOnce) {
        return(list(.rgDatasetUnit(group, labels, qd, args)))
    }
    pairs <- .rgGroupPairs(length(group))
    map2(
        pairs$i,
        pairs$j,
        .rgDatasetPairUnit,
        group = group,
        labels = labels,
        qd = qd,
        args = args
    )
}

# @noRd
.rgDatasetPairUnit <- function(i, j, group, labels, qd, args) {
    .rgDatasetUnit(group[c(i, j)], labels[c(i, j)], qd, args)
}

# @noRd
.rgDatasetGroupGwasUnit <- function(study, group, labels, qd, args) {
    .rgDatasetUnit(group, c(labels, study), qd, args, gwasStudy = study)
}

# Each GWAS study's unit with one condition, labelled by the two studies.
# @noRd
.rgDatasetGwasUnits <- function(condition, qd, args, gwas) {
    map(
        unique(as.character(gwas$study)),
        .rgDatasetGwasUnit,
        condition = condition,
        qd = qd,
        args = args
    )
}

# @noRd
.rgDatasetGwasUnit <- function(study, condition, qd, args) {
    .rgDatasetUnit(
        list(condition),
        c(condition$study, study),
        qd,
        args,
        gwasStudy = study
    )
}

# A unit: its QTL conditions and their labels, the GWAS study it pairs them
# with, the window both or all of its traits' windows hold, the contexts whose
# covariates its genotypes are residualized on, and the donors all its
# conditions share.
# @noRd
.rgDatasetUnit <- function(conditions, labels, qd, args, gwasStudy = NA) {
    traits <- unique(map_chr(conditions, "trait"))
    windows <- map(traits, .rgTraitWindow, qd = qd, args = args)
    list(
        conditions = conditions,
        labels = labels,
        gwasStudy = gwasStudy,
        window = reduce(windows, GenomicRanges::intersect),
        contexts = unique(map_chr(conditions, "context")),
        samples = reduce(map(conditions, .rgConditionSamples), intersect)
    )
}

# @noRd
.rgTraitWindow <- function(tid, qd, args) {
    .fmGenoWindow(qd, tid, args$region, args$cisWindow)
}

# @noRd
.rgConditionSamples <- function(condition) {
    names(condition$y)
}

# The genotypes a unit's fit is on: its window, residualized on its contexts'
# covariates, filtered and scaled on its shared donors; none for a window its
# traits do not share or fewer than three shared donors.
# @noRd
.rgDatasetSpec <- function(unit) {
    if (length(unit$window) == 0L || length(unit$samples) < 3L) {
        return(NULL)
    }
    .qtlFitSpec(unit$window, unit$contexts, unit$samples)
}

# ---- one QtlDataset unit's fit -----------------------------------------------

# One unit's pairs on its block: each condition's marginal correlations with
# its genotypes on the unit's shared donors, those genotypes as the LD (a GWAS
# side's too, unless it names its own reference in the window, which it then
# reads), every condition's reference overlap 1 (a GWAS side's 0), the
# phenotypes' correlation for tied noise, and the variants a GWAS side
# carries matched by allele; then the settings' fit.
# @noRd
.rgDatasetUnitFit <- function(unit, param, run, block = NULL) {
    qtl <- .rgDatasetMembers(unit, block$X, block$af)
    members <- c(
        qtl$members,
        if (!is.na(unit$gwasStudy)) {
            list(.rgGwasMember(run$gwas, unit$gwasStudy, unit$window))
        }
    )
    out <- list(window = unit$window, members = members)
    label <- .rgUnitLabel(members)
    gwasRef <- .rgDatasetGwasRef(unit, run$gwas)
    common <- if (is.null(qtl$Xs)) {
        list(ids = character(0))
    } else {
        .rgCommonHeldBy(
            .rgCommonVariants(members, is.na(unit$gwasStudy)),
            compact(list(gwasRef)),
            label
        )
    }
    if (length(common$ids) < 2L) {
        empty <- .rgNoPairs(length(members), length(common$ids))
        return(list_assign(out, pairs = .rgUnitPairs(empty, unit, members)))
    }
    ud <- c(
        list(z = common$z, nVariant = common$N),
        .rgDatasetLd(qtl, common$ids, gwasRef, param, label),
        .rgDatasetTerms(unit, qtl, run, label)
    )
    pairs <- .rgFit(param, ud)
    list_assign(out, pairs = .rgUnitPairs(pairs, unit, members))
}

# A GWAS side's own reference in a QtlDataset unit's window, if it names one.
# @noRd
.rgDatasetGwasRef <- function(unit, gwas) {
    if (is.na(unit$gwasStudy)) {
        return(NULL)
    }
    .rgGwasRef(gwas, unit$gwasStudy, unit$window)
}

# A QtlDataset unit's LD: the genotypes on its shared donors, for a GWAS side
# too unless the GWAS names its own reference in the window, which that side
# then reads.
# @noRd
.rgDatasetLd <- function(qtl, ids, gwasRef, param, label) {
    genotypes <- list(
        ldFactor = qtl$Xs[, ids, drop = FALSE],
        nSamples = qtl$n,
        het = unname(qtl$het[ids])
    )
    if (is.null(gwasRef)) {
        return(list(ld = genotypes))
    }
    list(
        ld = NULL,
        lds = list(genotypes, .cisH2LdSource(ids, gwasRef, param, label)),
        ldOf = c(rep(1L, length(qtl$members)), 2L)
    )
}

# A QtlDataset unit's per-condition terms: each condition's overlap with its
# own genotypes and none for a GWAS side, the held intercepts, the noise
# correlations its labels name, the phenotypes' correlation on the shared
# donors (NA with a GWAS side), and whether its pairs are given its other
# conditions.
# @noRd
.rgDatasetTerms <- function(unit, qtl, run, label) {
    list(
        overlap = c(rep(1, length(qtl$members)), if (!is.na(unit$gwasStudy)) 0),
        heldIntercept = .rgHeldIntercepts(
            length(qtl$members),
            unit$gwasStudy,
            run$gwasIntercepts
        ),
        noiseCorrelation = .rgUnitNoise(
            run$noiseCorrelation,
            unit$labels,
            label
        ),
        phenotypicCorrelation = .rgWithGwasSide(qtl$Ph, unit$gwasStudy),
        conditional = run$conditional
    )
}

# The conditions' phenotypes' correlation with a GWAS side's row and column
# NA: a GWAS shares no phenotypes, or donors, with them.
# @noRd
.rgWithGwasSide <- function(Ph, gwasStudy) {
    if (is.null(Ph) || is.na(gwasStudy)) {
        return(Ph)
    }
    rbind(cbind(Ph, NA_real_), NA_real_)
}

# The unit's QTL conditions on its shared donors that `X` (the residualized
# block) holds: each one's identity and statistics, the genotypes standardized
# on those donors over the variants that vary there, their count and the
# phenotypes' correlation. No statistics with fewer than three such donors.
# @noRd
#' @importFrom matrixStats colSds
#' @importFrom stats cor
.rgDatasetMembers <- function(unit, X, af = NULL) {
    ids <- map(unit$conditions, .rgDatasetMemberId)
    donors <- if (is.null(X)) {
        character(0)
    } else {
        intersect(unit$samples, rownames(X))
    }
    if (length(donors) < 3L) {
        return(list(members = ids))
    }
    Xd <- .qtlRowsOf(X, donors)
    sds <- colSds(Xd)
    Xs <- scale(Xd[, !is.na(sds) & sds > 0, drop = FALSE])
    phenotypes <- map(unit$conditions, .rgConditionY, donors = donors)
    Y <- exec(cbind, !!!phenotypes)
    n <- length(donors)
    r <- cor(Xs, Y)
    list(
        members = map2(
            ids,
            seq_len(ncol(Y)),
            .rgDatasetMember,
            r = r,
            variants = colnames(Xs),
            n = n
        ),
        Xs = Xs,
        n = n,
        Ph = cor(Y),
        het = .rgBlockHet(af, X, colnames(Xs))
    )
}

# Each of `ids`' heterozygosity 2p(1 - p) from the block's allele
# frequencies (one per column of `X`), named by variant; NULL where the
# block carries none.
# @noRd
.rgBlockHet <- function(af, X, ids) {
    if (is.null(af)) {
        return(NULL)
    }
    p <- af[match(ids, colnames(X))]
    stats::setNames(2 * p * (1 - p), ids)
}

# @noRd
.rgDatasetMemberId <- function(condition) {
    condition[c("study", "context", "trait")]
}

# @noRd
.rgConditionY <- function(condition, donors) {
    condition$y[donors]
}

# One condition's statistics: its marginal t statistics on the n donors.
# @noRd
.rgDatasetMember <- function(id, k, r, variants, n) {
    c(
        id,
        list(df = tibble(variant_id = variants, z = .corToZ(r[, k], n), N = n))
    )
}

# =============================================================================
# estimateLocalRg on a MultiStudyQtlDataset
# =============================================================================

#' @rdname estimateLocalRg
#' @export
setMethod(
    "estimateLocalRg",
    "MultiStudyQtlDataset",
    function(
        data,
        gwas = NULL,
        param = NULL,
        jointSpecification = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        genotypeFilterParam = GenotypeFilterParam(),
        residualizationParam = ResidualizationParam(),
        phenotypeFilterParam = PhenotypeFilterParam(),
        referenceOverlap = 0,
        noiseCorrelation = NULL,
        gwasH2Estimate = NULL,
        conditional = FALSE
    ) {
        settings <- .rgParam(
            param,
            gwas,
            noiseCorrelation,
            gwasH2Estimate,
            conditional,
            TRUE
        )
        label <- "estimateLocalRg(MultiStudyQtlDataset)"
        .h2AssertOneRegion(region, cisWindow, label)
        filter <- .qtlResolveFilter(genotypeFilterParam)
        run <- .rgRun(gwas, noiseCorrelation, gwasH2Estimate, conditional)
        args <- .rgDatasetArgs(
            region,
            cisWindow,
            residualizationParam,
            phenotypeFilterParam,
            jointSpecification
        )
        selection <- list(
            contexts = contexts,
            traitId = traitId,
            args = args,
            referenceOverlap = referenceOverlap
        )
        plan <- .rgMsPlan(data, filter, settings, run, selection)
        fits <- .rgMsFits(plan$units, plan$studies, plan$ctx)
        .rgEstimate(fits, plan$seqinfoOf, plan$ctx$param, run)
    }
)

# A MultiStudyQtlDataset run's studies (their samples filtered), units, what
# every unit's fit reads (`ctx`: the settings, the run, the genotype filter
# and the residualization), and what the RgEstimate takes its genome from:
# the first study's genotypes, as every collection holds one.
# @noRd
.rgMsPlan <- function(data, filter, settings, run, selection) {
    studies <- map(as.list(data), .qtlRestrictSamples, filter$keepSamples)
    conditions <- .rgMsConditions(studies, sumStats(data), selection)
    list(
        studies = studies,
        units = .rgMsUnits(conditions, data, settings, selection$args, run),
        ctx = list(
            param = settings,
            run = run,
            filter = filter,
            residualizationParam = selection$args$residualizationParam
        ),
        seqinfoOf = .qtlGenotypeSe(studies[[1L]])
    )
}

# ---- a MultiStudyQtlDataset's conditions and units ---------------------------

# Every selected condition of the collection: each study's (context, trait)
# conditions as the QtlDataset method selects them, and each selected row of
# its summary statistics, each with the window it sits in. A selection of
# contexts or traits applies to the studies that hold them.
# @noRd
.rgMsConditions <- function(studies, ss, selection) {
    individual <- list_c(imap(
        studies,
        .rgMsStudyConditions,
        selection = selection
    ))
    c(individual %||% list(), .rgMsSummaryConditions(ss, selection))
}

# @noRd
.rgMsStudyConditions <- function(qd, name, selection) {
    useContexts <- .rgMsHeld(selection$contexts, contexts(qd))
    useTraits <- .rgMsHeld(selection$traitId, .spListTraits(qd))
    if (
        identical(useContexts, character(0)) ||
            identical(useTraits, character(0))
    ) {
        return(list())
    }
    conditions <- .rgDatasetConditions(
        qd,
        useContexts,
        useTraits,
        selection$args
    )
    map(conditions, .rgMsIndividualCondition, qd = qd, args = selection$args)
}

# The selected values a component holds; NULL (all) without a selection, and
# a selection that is not a plain vector of names left to the component.
# @noRd
.rgMsHeld <- function(selected, held) {
    if (is.null(selected) || !is.character(selected)) {
        return(selected)
    }
    intersect(selected, held)
}

# @noRd
.rgMsIndividualCondition <- function(condition, qd, args) {
    list_assign(
        condition,
        kind = "individual",
        window = .rgTraitWindow(condition$trait, qd, args)
    )
}

# The summary statistics' selected rows as conditions, each with its row,
# window and reference overlap.
# @noRd
.rgMsSummaryConditions <- function(ss, selection) {
    if (is.null(ss)) {
        return(list())
    }
    useContexts <- .rgMsHeld(
        selection$contexts,
        unique(as.character(ss$context))
    )
    useTraits <- .rgMsHeld(selection$traitId, unique(as.character(ss$trait)))
    if (
        identical(useContexts, character(0)) ||
            identical(useTraits, character(0))
    ) {
        return(list())
    }
    selected <- .rgSelectedRows(
        ss,
        useContexts,
        useTraits,
        selection$args$cisWindow
    )
    overlaps <- .cisH2RowOverlaps(
        selection$referenceOverlap,
        selected$rows,
        length(selected$data),
        "estimateLocalRg"
    )
    map(
        selected$rows,
        .rgMsSummaryCondition,
        data = selected$data,
        overlaps = overlaps
    )
}

# @noRd
.rgMsSummaryCondition <- function(row, data, overlaps) {
    list(
        kind = "summary",
        study = as.character(data$study[[row]]),
        context = as.character(data$context[[row]]),
        trait = as.character(data$trait[[row]]),
        row = row,
        data = data,
        window = .rowRegionAt(data, row),
        overlap = overlaps[[row]]
    )
}

# The units the collection's pairs are fitted in: its conditions grouped
# along each jointSpecification's axes (each trait's contexts by default; the
# study axis pairs across studies), each group giving its pairs or, where the
# settings fit a group at once, one unit; with `gwas`, each condition with
# each GWAS study.
# @noRd
.rgMsUnits <- function(conditions, data, settings, args, run) {
    if (!is.null(run$gwas) && !run$conditional) {
        .rgAssertGwasPairing(args$jointSpecification)
        return(list_c(map(conditions, .rgMsGwasUnits, gwas = run$gwas)))
    }
    specs <- parseJointSpecification(
        args$jointSpecification %||% "context",
        data
    )
    units <- list_c(map(
        specs,
        .rgMsSpecUnits,
        conditions = conditions,
        atOnce = .rgFitsAtOnce(settings),
        gwas = run$gwas
    ))
    if (length(units) == 0L) {
        abort(str_c(
            "estimateLocalRg: no two selected conditions share a group along ",
            "the joint axes, so there is no pair to estimate."
        ))
    }
    units
}

# One joint spec's units: its conditions in scope, grouped by the identity
# columns its axes leave out.
# @noRd
.rgMsSpecUnits <- function(spec, conditions, atOnce, gwas) {
    inScope <- keep(conditions, .rgConditionInScope, scope = spec$scope)
    complement <- setdiff(c("study", "context", "trait"), spec$axes)
    keys <- map_chr(inScope, .rgConditionLabel, cols = complement)
    groups <- unname(split(inScope, keys))
    list_c(map(
        groups[lengths(groups) >= (if (is.null(gwas)) 2L else 1L)],
        .rgMsGroupUnits,
        axes = intersect(c("study", "context", "trait"), spec$axes),
        atOnce = atOnce,
        gwas = gwas
    ))
}

# A group's units, each condition labelled by its values on the joint axes.
# @noRd
.rgMsGroupUnits <- function(group, axes, atOnce, gwas) {
    labels <- map_chr(group, .rgConditionLabel, cols = axes)
    if (!is.null(gwas)) {
        return(map(
            unique(as.character(gwas$study)),
            .rgMsGroupGwasUnit,
            group = group,
            labels = labels
        ))
    }
    if (atOnce) {
        return(list(.rgMsUnit(group, labels)))
    }
    pairs <- .rgGroupPairs(length(group))
    map2(pairs$i, pairs$j, .rgMsPairUnit, group = group, labels = labels)
}

# @noRd
.rgMsPairUnit <- function(i, j, group, labels) {
    .rgMsUnit(group[c(i, j)], labels[c(i, j)])
}

# @noRd
.rgMsGroupGwasUnit <- function(study, group, labels) {
    .rgMsUnit(group, c(labels, study), study)
}

# Each GWAS study's unit with one condition, labelled by the two studies.
# @noRd
.rgMsGwasUnits <- function(condition, gwas) {
    map(unique(as.character(gwas$study)), .rgMsGwasUnit, condition = condition)
}

# @noRd
.rgMsGwasUnit <- function(study, condition) {
    .rgMsUnit(list(condition), c(condition$study, study), study)
}

# A unit: its conditions and labels, the GWAS study it pairs them with, and
# the window all its conditions' windows hold.
# @noRd
.rgMsUnit <- function(conditions, labels, gwasStudy = NA_character_) {
    windows <- compact(map(conditions, "window"))
    list(
        conditions = conditions,
        labels = labels,
        gwasStudy = gwasStudy,
        window = if (length(windows) > 0L) {
            reduce(windows, GenomicRanges::intersect)
        }
    )
}

# The LD sources a unit's conditions read, identical ones counted once
# (.rgDistinctRefs()): a study's genotypes for its conditions, a summary
# row's reference, and a GWAS study's own reference in the window (the first
# condition's source where it names none).
# @noRd
.rgMsSources <- function(unit, gwas) {
    keys <- map(unit$conditions, .rgMsSourceKey)
    if (is.na(unit$gwasStudy)) {
        return(.rgDistinctRefs(keys))
    }
    gwasRef <- .rgGwasRef(gwas, unit$gwasStudy, unit$window)
    .rgDistinctRefs(c(keys, list(gwasRef %||% keys[[1L]])))
}

# @noRd
.rgMsSourceKey <- function(condition) {
    if (identical(condition$kind, "individual")) {
        return(list(genotypes = condition$study))
    }
    .rowLdRef(condition$data, condition$row)
}

# ---- a MultiStudyQtlDataset's fits -------------------------------------------

# The units' fits a window at a time: each study's block of the window read
# once, for every unit of the window, then each unit's pairs; the fits in the
# units' order.
# @noRd
.rgMsFits <- function(units, studies, ctx) {
    keys <- map_chr(units, .rgMsWindowKey)
    groups <- unname(split(seq_along(units), factor(keys, unique(keys))))
    fits <- list_c(map(
        groups,
        .rgMsWindowFits,
        units = units,
        studies = studies,
        ctx = ctx
    ))
    fits[order(list_c(groups))]
}

# @noRd
.rgMsWindowKey <- function(unit) {
    if (length(unit$window) == 0L) "none" else as.character(unit$window)[[1L]]
}

# One window's units: each study's parts of them, then each unit's fit.
# @noRd
.rgMsWindowFits <- function(at, units, studies, ctx) {
    parts <- imap(studies, .rgMsStudyParts, units = units[at], ctx = ctx)
    map(
        seq_along(at),
        .rgMsUnitFitAt,
        units = units[at],
        parts = parts,
        ctx = ctx
    )
}

# One study's parts of `units`: each unit's conditions in the study on their
# shared donors off the study's block (.rgDatasetMembers()); NULL for a unit
# with none there.
# @noRd
.rgMsStudyParts <- function(qd, name, units, ctx) {
    studyUnits <- map(units, .rgMsStudyUnit, study = name)
    present <- which(!map_lgl(studyUnits, is.null))
    none <- rep(list(NULL), length(units))
    if (length(present) == 0L) {
        return(none)
    }
    work <- .qtlWork(
        studyUnits[present],
        map(studyUnits[present], .rgDatasetSpec),
        .rgMsStudyPart
    )
    fitted <- .qtlRunWork(
        work,
        qd,
        ctx$filter,
        ctx$residualizationParam,
        withAf = ctx$param@alpha != -1
    )
    replace(none, present, fitted)
}

# A unit's conditions in one study, on the unit's window and their shared
# donors.
# @noRd
.rgMsStudyUnit <- function(unit, study) {
    mine <- keep(unit$conditions, .rgMsInStudy, study = study)
    if (length(mine) == 0L) {
        return(NULL)
    }
    list(
        conditions = mine,
        window = unit$window,
        contexts = unique(map_chr(mine, "context")),
        samples = reduce(map(mine, .rgConditionSamples), intersect)
    )
}

# @noRd
.rgMsInStudy <- function(condition, study) {
    identical(condition$kind, "individual") && identical(condition$study, study)
}

# @noRd
.rgMsStudyPart <- function(unit, block = NULL) {
    .rgDatasetMembers(unit, block$X, block$af)
}

# @noRd
.rgMsUnitFitAt <- function(k, units, parts, ctx) {
    .rgMsUnitFit(units[[k]], map(parts, k), ctx)
}

# One unit's pairs: each condition's statistics (a study's on its shared
# donors, a summary row's as given, a GWAS study's in the window) on the
# variants they all carry, matched by allele to the first condition's coding
# and held by every reference, each condition reading its own LD source;
# then the settings' fit.
# @noRd
.rgMsUnitFit <- function(unit, parts, ctx) {
    members <- .rgMsMembers(unit, parts, ctx$run$gwas)
    out <- list(window = unit$window, members = members)
    label <- .rgUnitLabel(members)
    sources <- .rgMsSources(unit, ctx$run$gwas)
    common <- .rgMsCommon(members, unit, sources, label)
    if (length(common$ids) < 2L) {
        empty <- .rgNoPairs(length(members), length(common$ids))
        return(list_assign(out, pairs = .rgUnitPairs(empty, unit, members)))
    }
    ud <- c(
        list(z = common$z, nVariant = common$N),
        .rgMsLd(common$ids, sources, parts, ctx$param, label),
        .rgMsTerms(unit, parts, ctx$run, label)
    )
    .rgAssertTiedDonors(ctx$param, ud, label)
    pairs <- .rgFit(ctx$param, ud)
    list_assign(out, pairs = .rgUnitPairs(pairs, unit, members))
}

# The unit's members in its conditions' order, then its GWAS side.
# @noRd
.rgMsMembers <- function(unit, parts, gwas) {
    qtl <- map(unit$conditions, .rgMsMember, parts = parts)
    if (is.na(unit$gwasStudy)) {
        return(qtl)
    }
    c(qtl, list(.rgGwasMember(gwas, unit$gwasStudy, unit$window)))
}

# A condition's statistics: a summary row's, or its study part's (none where
# the study had too few shared donors).
# @noRd
.rgMsMember <- function(condition, parts) {
    if (identical(condition$kind, "summary")) {
        return(.rgQtlMember(condition$row, condition$data))
    }
    found <- detect(
        parts[[condition$study]]$members,
        .rgMsSameCondition,
        condition = condition
    )
    found %||% condition[c("study", "context", "trait")]
}

# @noRd
.rgMsSameCondition <- function(member, condition) {
    identical(member$context, condition$context) &&
        identical(member$trait, condition$trait)
}

# The variants every condition carries, matched by allele (strand-ambiguous
# ones kept only on one LD source) and held by every reference; none where a
# condition has no statistics.
# @noRd
.rgMsCommon <- function(members, unit, sources, label) {
    if (any(map_lgl(members, .rgMemberLacksStats))) {
        return(list(ids = character(0)))
    }
    oneSource <- length(sources$refs) == 1L && is.na(unit$gwasStudy)
    .rgCommonHeldBy(
        .rgCommonVariants(members, oneSource),
        discard(sources$refs, .rgIsGenotypesSource),
        label
    )
}

# @noRd
.rgMemberLacksStats <- function(member) {
    is.null(member$df)
}

# @noRd
.rgIsGenotypesSource <- function(source) {
    is.list(source) && identical(names(source), "genotypes")
}

# The unit's LD: one source where its conditions read one, else one per
# source (`lds`) and the one each condition reads (`ldOf`).
# @noRd
.rgMsLd <- function(ids, sources, parts, param, label) {
    lds <- map(
        sources$refs,
        .rgMsSourceLd,
        ids = ids,
        parts = parts,
        param = param,
        label = label
    )
    if (length(lds) == 1L) {
        return(list(ld = lds[[1L]]))
    }
    list(ld = NULL, lds = lds, ldOf = sources$of)
}

# One source's LD over `ids`: a study's standardized genotypes on its shared
# donors, matched by allele and sign-flipped into the ids' coding, or a
# reference's.
# @noRd
.rgMsSourceLd <- function(source, ids, parts, param, label) {
    if (!.rgIsGenotypesSource(source)) {
        return(.cisH2LdSource(ids, source, param, label))
    }
    part <- parts[[source$genotypes]]
    factor <- .rgAlignedFactor(part$Xs, ids)
    list(
        ldFactor = factor,
        nSamples = part$n,
        het = .rgAlignedHet(part$het, colnames(part$Xs), ids)
    )
}

# A study's heterozygosities in the ids' order, matched by allele as its
# genotypes are (.rgAlignedFactor()); a swap leaves each as it is. NULL
# where the study's block carried no frequencies.
# @noRd
.rgAlignedHet <- function(het, studyIds, ids) {
    if (is.null(het)) {
        return(NULL)
    }
    m <- matchVariants(ids, studyIds, removeStrandAmbiguous = FALSE)
    at <- order(m$idxA)
    unname(het[m$idxB[at]])
}

# @noRd
.rgAlignedFactor <- function(Xs, ids) {
    m <- matchVariants(ids, colnames(Xs), removeStrandAmbiguous = FALSE)
    at <- order(m$idxA)
    named <- sweep(Xs[, m$idxB[at], drop = FALSE], 2L, m$sign[at], `*`)
    `colnames<-`(named, ids[m$idxA[at]])
}

# A multi-study unit's per-condition terms: overlap 1 for a study's conditions
# (their own genotypes), a summary row's own and none for a GWAS side, the
# held intercepts, the noise correlations its labels name, the phenotypes'
# correlation where every condition is one study's, and whether its pairs are
# given its other conditions.
# @noRd
.rgMsTerms <- function(unit, parts, run, label) {
    list(
        overlap = c(
            map_dbl(unit$conditions, .rgMsOverlap),
            if (!is.na(unit$gwasStudy)) 0
        ),
        heldIntercept = .rgHeldIntercepts(
            length(unit$conditions),
            unit$gwasStudy,
            run$gwasIntercepts
        ),
        noiseCorrelation = .rgUnitNoise(
            run$noiseCorrelation,
            unit$labels,
            label
        ),
        phenotypicCorrelation = .rgWithGwasSide(
            .rgMsPhenotypic(unit, parts),
            unit$gwasStudy
        ),
        conditional = run$conditional
    )
}

# @noRd
.rgMsOverlap <- function(condition) {
    if (identical(condition$kind, "summary")) condition$overlap else 1
}

# The phenotypes' correlation on each study's shared donors, for tied noise:
# each study's block among its own conditions, and NA between conditions of
# different studies, which share no donors; NULL where a condition is a
# summary row's.
# @noRd
.rgMsPhenotypic <- function(unit, parts) {
    kinds <- unique(map_chr(unit$conditions, "kind"))
    if (!identical(kinds, "individual")) {
        return(NULL)
    }
    studies <- map_chr(unit$conditions, "study")
    within <- map_int(seq_along(studies), .rgMsRankInStudy, studies = studies)
    P <- length(studies)
    cells <- expand.grid(k = seq_len(P), l = seq_len(P))
    matrix(
        pmap_dbl(
            cells,
            .rgMsPhenotypicCell,
            studies = studies,
            within = within,
            parts = parts
        ),
        P
    )
}

# Condition k's place among its study's conditions in the unit.
# @noRd
.rgMsRankInStudy <- function(k, studies) {
    sum(studies[seq_len(k)] == studies[[k]])
}

# @noRd
.rgMsPhenotypicCell <- function(k, l, studies, within, parts) {
    if (!identical(studies[[k]], studies[[l]])) {
        return(NA_real_)
    }
    parts[[studies[[k]]]]$Ph[[within[[k]], within[[l]]]]
}

# ---- the RgEstimate ----------------------------------------------------------

# The units' pairs as an RgEstimate on the data's genome, recording whether
# the run (.rgRun()) gives them their groups' other conditions, with a GWAS
# side on the liability scale too where the run gives its study a factor,
# and one warning listing the pairs that have no rg.
# @noRd
.rgEstimate <- function(fits, data, settings, run) {
    rowsOf <- list_c(map(fits, .rgFitRows, liability = run$gwasLiability))
    values <- list_rbind(map(rowsOf, "values"))
    out <- .rgEstimateNew(
        map(rowsOf, "window"),
        as.list(values),
        data,
        settings,
        run$conditional
    )
    .rgWarnNoEstimate(out)
    .rgWarnNoInterval(settings, out)
    out
}

# The moments' interval on one reference under fixed effects draws from the
# noncentral Wishart with matrixsampling; one warning where it is not
# installed and a pair with an rg has no interval (random effects' draws, and
# two references', need no package).
# @noRd
.rgWarnNoInterval <- function(settings, out) {
    if (
        !identical(settings@criterion, "moments") ||
            identical(settings@effects, "random") ||
            requireNamespace("matrixsampling", quietly = TRUE)
    ) {
        return(invisible(NULL))
    }
    md <- mcols(out, use.names = FALSE)
    if (!any(is.finite(md$rg) & is.na(md$rg_lower))) {
        return(invisible(NULL))
    }
    warn(str_c(
        "estimateLocalRg: the moments' interval draws from the noncentral ",
        "Wishart with the matrixsampling package, which is not installed; ",
        "rg_lower and rg_upper are NA."
    ))
}

# One unit's rows: each pair's window, its two conditions' identities and its
# values, on the liability scale too where a side converts (.rgPairRow()).
# @noRd
.rgFitRows <- function(fit, liability) {
    map(
        fit$pairs,
        .rgPairRow,
        members = fit$members,
        window = fit$window %||% GenomicRanges::GRanges(),
        liability = liability
    )
}

# @noRd
.rgPairRow <- function(pair, members, window, liability) {
    one <- members[[pair$i]]
    two <- members[[pair$j]]
    identity <- list(
        study_1 = one$study,
        context_1 = one$context,
        trait_1 = one$trait,
        study_2 = two$study,
        context_2 = two$context,
        trait_2 = two$trait
    )
    factors <- map_dbl(list(one, two), .rgSideLiability, liability = liability)
    values <- list_assign(
        pair$values,
        !!!.rgLiabilityValues(pair$values, factors)
    )
    list(window = window, values = as_tibble(c(identity, values)))
}

# A side's factor to the liability scale: its GWAS study's in `liability`
# (NA where the study records no prevalence); NA for a QTL condition, which
# has a context and a trait.
# @noRd
.rgSideLiability <- function(member, liability) {
    gwas <- is.na(member$context) && is.na(member$trait)
    if (!gwas || is.null(liability)) {
        return(NA_real_)
    }
    liability[[member$study]]
}

# A pair's values with each side that has a factor (`factors`, NA for one
# without) on the liability scale: its heritability times its factor, and
# the covariance and its standard error times the root of each side's
# factor, a side without one taking 1 as ldsc's gencov_obs_to_liab() does
# for a quantitative trait; NA where neither side converts. rg, its
# interval and the r2s, ratios of these, do not change.
# @noRd
.rgLiabilityValues <- function(values, factors) {
    scale <- if (all(is.na(factors))) {
        NA_real_
    } else {
        prod(sqrt(replace(factors, is.na(factors), 1)))
    }
    list(
        covariance_liability = scale * values$covariance,
        covariance_liability_se = scale * values$covariance_se,
        h2_1_liability = factors[[1L]] * values$h2_1,
        h2_2_liability = factors[[2L]] * values$h2_2
    )
}

# The pairs without an rg, listed in one warning: a window holding too few
# variants both conditions carry, a heritability that is not positive, or a
# fit that did not converge.
# @noRd
.rgWarnNoEstimate <- function(out) {
    none <- is.na(out$rg)
    if (!any(none)) {
        return(invisible(NULL))
    }
    keys <- str_c(
        .rgJoinedColumns(out, which(none), .RG_ID_COLS[seq_len(3L)], "|"),
        " ~ ",
        .rgJoinedColumns(out, which(none), .RG_ID_COLS[-seq_len(3L)], "|")
    )
    warn(glue(
        "estimateLocalRg: {sum(none)} of {length(none)} pair(s) have no rg ",
        "(too few shared variants, a heritability that is not positive, or ",
        "no converged fit): {str_flatten(keys, '; ')}."
    ))
}

# =============================================================================
# summarizeLocalRg: rg across genes
# =============================================================================

#' @rdname summarizeLocalRg
#' @export
#' @importFrom checkmate assertNumber
setMethod(
    "summarizeLocalRg",
    "RgEstimate",
    function(
        x,
        summary = c("pooled", "mean", "shared"),
        mashModel = NULL,
        minH2 = 0
    ) {
        summaries <- arg_match(
            summary,
            c("pooled", "mean", "shared"),
            multiple = TRUE
        )
        assertNumber(minH2, finite = TRUE, .var.name = "minH2")
        implied <- .rgMashImplied(mashModel, x)
        .rgWarnBoundedPooling(x, summaries)
        md <- mcols(x, use.names = FALSE)
        keys <- .rgGroupKeys(md)
        groups <- unname(split(seq_along(x), factor(keys, unique(keys))))
        rows <- list_c(map(
            groups,
            .rgSummaryGroup,
            md = md,
            summaries = summaries,
            minH2 = minH2,
            implied = implied
        ))
        .rgSummaryFrame(rows, x)
    }
)

# The identity columns a summary's groups share: everything but the traits.
# @noRd
.RG_GROUP_COLS <- c("study_1", "context_1", "study_2", "context_2")

# Each pair's group, NA kept apart from every value.
# @noRd
.rgGroupKeys <- function(md) {
    parts <- map(.RG_GROUP_COLS, .rgGroupKeyPart, md = md)
    exec(str_c, !!!parts, sep = "\r")
}

# @noRd
#' @importFrom stringr str_replace_na
.rgGroupKeyPart <- function(col, md) {
    str_replace_na(as.character(md[[col]]), "\n")
}

# One group's rows, one per summary, each led by the group's identity and
# ending in mash's correlation of its two conditions.
# @noRd
.rgSummaryGroup <- function(at, md, summaries, minH2, implied) {
    first <- at[[1L]]
    identity <- map(
        set_names(.RG_GROUP_COLS),
        .rgGroupValue,
        md = md,
        k = first
    )
    mash <- .rgMashCorrelation(implied, identity)
    map(
        summaries,
        .rgSummaryRow,
        at = at,
        md = md,
        minH2 = minH2,
        identity = identity,
        mash = mash
    )
}

# @noRd
.rgGroupValue <- function(col, md, k) {
    as.character(md[[col]][[k]])
}

# @noRd
.rgSummaryRow <- function(summary, at, md, minH2, identity, mash) {
    values <- switch(
        summary,
        pooled = .rgPooledSummary(md, at),
        mean = .rgMeanSummary(md, at, minH2),
        shared = .rgSharedSummary(md, at)
    )
    c(identity, list(summary = summary), values, list(mash_rg = mash))
}

# A summary's values: rg, its standard error, its 95% interval (Wald's unless
# given) and p-value (Wald's unless given), and the genes it used.
# @noRd
#' @importFrom stats qnorm
.rgSummaryValues <- function(
    rg,
    se,
    n,
    interval = rg + c(-1, 1) * qnorm(0.975) * se,
    p = .rgWaldP(rg, se)
) {
    list(
        rg = rg,
        rg_se = se,
        rg_lower = interval[[1L]],
        rg_upper = interval[[2L]],
        p_value = p,
        n_genes = as.integer(n)
    )
}

# The pooled ratio of sums over the group's genes with all three estimates,
# and its leave-one-gene-out jackknife.
# @noRd
.rgPooledSummary <- function(md, at) {
    keep <- at[
        is.finite(md$covariance[at]) &
            is.finite(md$h2_1[at]) &
            is.finite(md$h2_2[at])
    ]
    parts <- list(
        covariance = md$covariance[keep],
        h1 = md$h2_1[keep],
        h2 = md$h2_2[keep]
    )
    rg <- .rgPooledRatio(seq_along(keep), parts)
    loo <- map_dbl(seq_along(keep), .rgPooledLoo, parts = parts)
    se <- if (length(keep) >= 2L && !anyNA(loo)) {
        .jackknifeSeFromLoo(loo)
    } else {
        NA_real_
    }
    .rgSummaryValues(rg, se, length(keep))
}

# The ratio of sums over genes `at` of `parts`; NA unless both sums of
# heritabilities are positive.
# @noRd
.rgPooledRatio <- function(at, parts) {
    s1 <- sum(parts$h1[at])
    s2 <- sum(parts$h2[at])
    if (length(at) == 0L || s1 <= 0 || s2 <= 0) {
        return(NA_real_)
    }
    sum(parts$covariance[at]) / sqrt(s1 * s2)
}

# @noRd
.rgPooledLoo <- function(g, parts) {
    .rgPooledRatio(setdiff(seq_along(parts$h1), g), parts)
}

# The mean of the group's genes' rg, over those where it is defined and both
# heritabilities exceed `minH2`.
# @noRd
#' @importFrom stats sd
.rgMeanSummary <- function(md, at, minH2) {
    passes <- is.finite(md$rg[at]) &
        is_in(md$h2_1[at] > minH2, TRUE) &
        is_in(md$h2_2[at] > minH2, TRUE)
    rg <- md$rg[at[passes]]
    n <- length(rg)
    .rgSummaryValues(
        if (n > 0L) mean(rg) else NA_real_,
        if (n >= 2L) sd(rg) / sqrt(n) else NA_real_,
        n
    )
}

# One rg for the group: the maximum of its genes' summed log-likelihoods, the
# interval where their likelihood ratio stays within the cutoff and its
# p-value at 0, and the sandwich standard error from each gene's slope at the
# maximum over the sum's curvature there. NA where no gene keeps a
# log-likelihood.
# @noRd
#' @importFrom stats pchisq
.rgSharedSummary <- function(md, at) {
    curves <- compact(md$rg_loglik[at])
    if (length(curves) == 0L) {
        return(.rgSummaryValues(NA_real_, NA_real_, 0L))
    }
    total <- reduce(curves, `+`)
    peak <- .rgCurvePeak(total)
    slopes <- map_dbl(curves, .rgCurveSlope, peak = peak)
    zero <- which.min(abs(.RG_CURVE_GRID))
    .rgSummaryValues(
        peak$at,
        sqrt(sum(slopes^2)) / abs(peak$curvature),
        length(curves),
        interval = .rgCurveInterval(total, peak),
        p = pchisq(
            max(2 * (peak$value - total[[zero]]), 0),
            1,
            lower.tail = FALSE
        )
    )
}

# A curve's maximum on the rg grid, refined by the parabola through it and
# its two neighbours (the grid's end where it lies there): where it lies, its
# value and the curvature there, and the grid point it is nearest.
# @noRd
.rgCurvePeak <- function(curve) {
    k <- which.max(curve)
    h <- diff(.RG_CURVE_GRID)[[1L]]
    mid <- min(max(k, 2L), length(curve) - 1L)
    f <- curve[mid + c(-1L, 0L, 1L)]
    second <- (f[[1L]] - 2 * f[[2L]] + f[[3L]]) / h^2
    first <- (f[[3L]] - f[[1L]]) / (2 * h)
    interior <- k == mid && second < 0
    offset <- if (interior) -first / second else 0
    list(
        at = .RG_CURVE_GRID[[k]] + offset,
        value = if (interior) f[[2L]] - first^2 / (2 * second) else curve[[k]],
        curvature = second,
        mid = mid
    )
}

# A gene's slope at the shared maximum: its parabola's through the three grid
# points around it.
# @noRd
.rgCurveSlope <- function(curve, peak) {
    h <- diff(.RG_CURVE_GRID)[[1L]]
    f <- curve[peak$mid + c(-1L, 0L, 1L)]
    first <- (f[[3L]] - f[[1L]]) / (2 * h)
    second <- (f[[1L]] - 2 * f[[2L]] + f[[3L]]) / h^2
    first + second * (peak$at - .RG_CURVE_GRID[[peak$mid]])
}

# Where the summed likelihood ratio from the maximum reaches the chi-square's
# 95% cutoff on each side, between grid points by linear interpolation; the
# grid's end where it never does.
# @noRd
#' @importFrom stats qchisq
.rgCurveInterval <- function(total, peak) {
    excess <- 2 * (peak$value - total) - qchisq(0.95, 1)
    below <- which(.RG_CURVE_GRID < peak$at)
    above <- which(.RG_CURVE_GRID > peak$at)
    c(
        .rgCurveCrossing(rev(below), excess, -1),
        .rgCurveCrossing(above, excess, 1)
    )
}

# The first point along `path` (grid indices from the maximum outward) where
# `excess` turns positive, interpolated from the point before; `end` when it
# never does.
# @noRd
.rgCurveCrossing <- function(path, excess, end) {
    out <- path[excess[path] > 0]
    if (length(out) == 0L) {
        return(end)
    }
    k <- out[[1L]]
    inner <- k - sign(end)
    if (!is_in(inner, path) && excess[[inner]] > 0) {
        return(.RG_CURVE_GRID[[k]])
    }
    w <- excess[[inner]] / (excess[[inner]] - excess[[k]])
    .RG_CURVE_GRID[[inner]] +
        w * (.RG_CURVE_GRID[[k]] - .RG_CURVE_GRID[[inner]])
}

# The likelihood bounds each pair (a covariance within the root of the
# product of the heritabilities, each at no less than 0), so a pooled ratio
# of its estimates reads low; one warning says so.
# @noRd
.rgWarnBoundedPooling <- function(x, summaries) {
    if (
        !is_in("pooled", summaries) ||
            !identical(h2EstimateParam(x)@criterion, "likelihood")
    ) {
        return(invisible(NULL))
    }
    warn(str_c(
        "summarizeLocalRg: the likelihood bounds each pair's covariance ",
        "within the root of the product of its heritabilities, so their ",
        "pooled ratio reads low; the moments' or a regression's unbounded ",
        "estimates pool without that."
    ))
}

# ---- mash's cross-check ------------------------------------------------------

# The correlation matrix of a fitted mash model's prior, sum_k (sum_l pi_kl
# grid_l^2) U_k (fitted_g$pi in mashr's expand_cov() order: the null first
# where it has one, then each grid value's components); NULL without a model.
# A conditional RgEstimate takes none: mash's correlation is not given a
# group's other conditions.
# @noRd
#' @importFrom stats cov2cor
.rgMashImplied <- function(mashModel, x) {
    if (is.null(mashModel)) {
        return(NULL)
    }
    if (isConditional(x)) {
        abort(str_c(
            "summarizeLocalRg: mash's prior correlation is not given a ",
            "group's other conditions, so it does not cross-check a ",
            "conditional RgEstimate."
        ))
    }
    g <- mashModel$fitted_g
    names <- rownames(g$Ulist[[1L]])
    if (is.null(g$pi) || is.null(g$grid) || is.null(names)) {
        abort(str_c(
            "summarizeLocalRg: `mashModel` must be a fitted mash model ",
            "(mashModelFit() or mashr::mash()), its covariances named by ",
            "condition."
        ))
    }
    pi <- if (isTRUE(g$usepointmass)) g$pi[-1L] else g$pi
    weights <- as.vector(matrix(pi, nrow = length(g$Ulist)) %*% g$grid^2)
    cov2cor(reduce(map2(g$Ulist, weights, `*`), `+`))
}

# mash's correlation of the group's two conditions, each named by its values
# on the study and context columns the two sides differ on; NA without a
# model or where either is not one of its conditions.
# @noRd
.rgMashCorrelation <- function(implied, identity) {
    if (is.null(implied)) {
        return(NA_real_)
    }
    one <- c(identity$study_1, identity$context_1)
    two <- c(identity$study_2, identity$context_2)
    differ <- is.na(one) | is.na(two) | one != two
    labels <- c(
        str_flatten(str_replace_na(one[differ]), "|"),
        str_flatten(str_replace_na(two[differ]), "|")
    )
    if (!any(differ) || !all(is_in(labels, rownames(implied)))) {
        return(NA_real_)
    }
    implied[[labels[[1L]], labels[[2L]]]]
}

# The summaries' rows as a DataFrame, its metadata the settings and the
# reading of the pairs summarized.
# @noRd
.rgSummaryFrame <- function(rows, x) {
    columns <- as.list(list_rbind(map(rows, as_tibble)))
    frame <- exec(S4Vectors::DataFrame, !!!columns, check.names = FALSE)
    S4Vectors::`metadata<-`(
        frame,
        value = list(
            param = h2EstimateParam(x),
            conditional = isConditional(x)
        )
    )
}
