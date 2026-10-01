#' @title P-value Combination Methods
#' @description Functions for combining p-values across multiple tests or
#'   methods: Cauchy combination (ACAT), harmonic mean p-value (HMP), poolr,
#'   GBJ, and aSPU methods. Also includes the null correlation matrix for TWAS
#'   z-scores used by correlation-adjusted methods.
#' @name pecotmr-pval-combine
#' @keywords internal
#' @importFrom stats pt pcauchy
NULL

# Internal: two-tailed normal p-value from a signed Z. The single shared helper
# across the sumstats-QC, TWAS, h2, mash and S-LDSC paths. Uses the numerically
# stable 2 * pnorm(-|z|) form (1 - pnorm(|z|) underflows earlier at large |z|).
# Returns NA where z is NA; |z| > ~37 underflows to 0 (R's pnorm limit).
# @noRd
.zToPvalue <- function(z) 2 * pnorm(-abs(z))

# Derive effect size + standard error from a signed Z, allele frequency (maf),
# and sample size (n) for a standardized phenotype, using the model-exact
# relationship  se = 1 / sqrt(2*maf*(1-maf)*(n + z^2)),  beta = z * se. The
# (n + z^2) term accounts for the residual-variance reduction from the SNP's own
# effect (residual = n/(n + z^2)); dropping z^2 is the weak-effect limit. Shared
# by the RAISS imputation (summaryStatsQc) and the CIP MR helpers.
# @noRd
.zToBetaSe <- function(z, maf, n) {
    se <- 1 / sqrt(2 * maf * (1 - maf) * (n + z * z))
    list(beta = z * se, se = se)
}

#' Wald-test two-sided p-value
#'
#' Two-sided p-value for a Wald test of \code{beta = 0} using the t-distribution
#' with \code{n - 2} degrees of freedom.
#'
#' @param beta Numeric vector of effect-size estimates.
#' @param se Numeric vector of standard errors (same length as \code{beta}).
#' @param n Integer sample size (used for the degrees of freedom).
#' @return Numeric vector of two-sided p-values.
#' @examples
#' waldTestPval(beta = 0.3, se = 0.1, n = 1000)
#' @importFrom checkmate assertNumeric
#' @export
waldTestPval <- function(beta, se, n) {
    assertNumeric(beta)
    assertNumeric(se)
    assertNumeric(n)
    # Calculate the t statistic
    tValue <- beta / se
    # Degrees of freedom
    df <- n - 2
    # Calculate two-tailed p-value
    pValue <- 2 * pt(-abs(tValue), df = df, lower.tail = TRUE)

    return(pValue)
}

pvalAcat <- function(pvals, naRm = TRUE) {
    # ACAT (Aggregated Cauchy Association Test) -- Liu & Xie (2020).
    # T = mean(tan(pi * (0.5 - p_i)))
    #
    # Robustness handling merged from the former pvalCauchy implementation:
    #   - naRm: optionally drop NAs
    #   - clip p > 0.99 to 0.99 to bound the contribution of near-1 p-values
    #   - small-p asymptotic: tan(pi*(0.5 - p)) ~ 1/(pi*p) for p < 1e-15 to
    #     avoid Inf from floating-point precision loss in pi*0.5
    #   - large-stat asymptotic: when the mean Cauchy variate is > 1e15 the
    #     CDF tail collapses to (1/T) / pi (Cauchy survival expansion)
    present <- if (naRm) pvals[!is.na(pvals)] else pvals
    if (length(present) == 0L) {
        return(NA_real_)
    }
    if (length(present) == 1L) {
        return(present[[1]])
    }
    # Capped below 1: tan(pi * (0.5 - 1)) is -Inf, which would sink the
    # combined statistic regardless of the other p-values.
    capped <- pmin(present, 0.99)
    cauchyVals <- if_else(
        capped < 1e-15,
        1 / (capped * pi),
        tan(pi * (0.5 - capped))
    )
    stat <- mean(cauchyVals)
    if (!is.finite(stat)) {
        return(NA_real_)
    }
    if (stat > 1e15) {
        return((1 / stat) / pi)
    }
    pcauchy(stat, lower.tail = FALSE)
}

pvalHmp <- function(pvals) {
    # Make sure harmonicmeanp is installed
    if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
        abort("Package 'harmonicmeanp' is required for this function.")
    }
    # https://search.r-project.org/CRAN/refmans/harmonicmeanp/html/pLandau.html
    L <- length(pvals)
    HMP <- L / sum(pvals^-1)

    LOC_L1 <- 0.874367040387922
    SCALE <- 1.5707963267949

    return(harmonicmeanp::pLandau(
        1 / HMP,
        mu = log(L) + LOC_L1,
        sigma = SCALE,
        lower.tail = FALSE
    ))
}

# Raise the "unknown method" error from a switch() default slot (which cannot
# host the two-statement glue + abort form).
# @noRd
.abortUnknownMethod <- function(kind, method) {
    msg <- glue("Unknown {kind} method: '{method}'")
    abort(msg)
}

pvalPoolr <- function(pvals, method, R, methodArgs = poolrConfig()) {
    if (!requireNamespace("poolr", quietly = TRUE)) {
        abort("Package 'poolr' is required for this method.")
    }
    fn <- switch(
        method,
        fisher = poolr::fisher,
        stouffer = poolr::stouffer,
        invchisq = poolr::invchisq,
        .abortUnknownMethod("poolr", method)
    )
    exec(fn, pvals, R = R, !!!as.list(methodArgs))$p
}

pvalGbj <- function(zScores, R, method, methodArgs = gbjConfig()) {
    if (!requireNamespace("GBJ", quietly = TRUE)) {
        abort("Package 'GBJ' is required for this method.")
    }
    fn <- switch(
        method,
        gbj = GBJ::GBJ,
        bj = GBJ::BJ,
        hc = GBJ::HC,
        ghc = GBJ::GHC,
        minp = GBJ::minP,
        gbj_omni = GBJ::OMNI_ss,
        .abortUnknownMethod("GBJ", method)
    )
    result <- exec(
        fn,
        test_stats = zScores,
        cor_mat = R,
        !!!as.list(methodArgs)
    )
    pvalName <- switch(
        method,
        gbj = "GBJ_pvalue",
        bj = "BJ_pvalue",
        hc = "HC_pvalue",
        ghc = "GHC_pvalue",
        minp = "minP_pvalue",
        gbj_omni = "OMNI_pvalue"
    )
    result[[pvalName]]
}

pvalAspu <- function(
    zScores = NULL,
    pvals = NULL,
    R,
    method,
    methodArgs = aspuConfig()
) {
    if (!requireNamespace("aSPU", quietly = TRUE)) {
        abort(
            "Package 'aSPU' is required for this method."
        )
    }
    switch(
        method,
        aspu = exec(
            aSPU::aSPUs,
            Zs = zScores,
            corSNP = R,
            !!!as.list(methodArgs)
        )$pvs["aSPUs"],
        # GATES2's only formals are the two pecotmr supplies, so there is
        # nothing for a caller to configure here.
        gates = aSPU::GATES2(ldmatrix = R, p = pvals)[["Pg"]],
        .abortUnknownMethod("aSPU", method)
    )
}

# =============================================================================
# combinePValues -- unified dispatcher
# =============================================================================

# Methods that require an R correlation matrix.
.combinePvalMethodsNeedingR <- c(
    "fisher",
    "stouffer",
    "invchisq",
    "gbj",
    "bj",
    "hc",
    "ghc",
    "minp",
    "gbj_omni",
    "aspu",
    "gates"
)

# Methods that require signed zScores (and therefore cannot be derived from
# pvals alone). All other methods can work from pvals.
.combinePvalMethodsNeedingZ <- c(
    "gbj",
    "bj",
    "hc",
    "ghc",
    "minp",
    "gbj_omni",
    "aspu"
)

.combinePvalKnownMethods <- c(
    "acat",
    "hmp",
    "bonferroni",
    "fisher",
    "stouffer",
    "invchisq",
    "gbj",
    "bj",
    "hc",
    "ghc",
    "minp",
    "gbj_omni",
    "aspu",
    "gates"
)

# Internal: align an R correlation matrix to a target order. If R has
# rownames/colnames, reorder to match `targetNames`; require every target
# name to be present. If R is unnamed, only length check.
#' @importFrom checkmate assertMatrix
.combinePvalAlignR <- function(R, targetNames) {
    if (is.null(R)) {
        return(NULL)
    }
    assertMatrix(R)
    if (nrow(R) != ncol(R)) {
        abort("`R` must be square.")
    }
    rNames <- rownames(R)
    cNames <- colnames(R)
    hasNames <- !is.null(rNames) && !is.null(cNames)
    if (hasNames) {
        if (!identical(rNames, cNames)) {
            abort("`R` rownames and colnames must be identical.")
        }
        missing <- setdiff(targetNames, rNames)
        if (length(missing) > 0L) {
            shown <- str_flatten(utils::head(missing, 5), ", ")
            moreSuffix <- if (length(missing) > 5L) {
                nMore <- length(missing) - 5L
                glue(" (and {nMore} more)")
            } else {
                ""
            }
            msg <- glue("`R` is missing entries for: {shown}{moreSuffix}")
            abort(msg)
        }
        R <- R[targetNames, targetNames, drop = FALSE]
    } else {
        if (nrow(R) != length(targetNames)) {
            nR <- nrow(R)
            nTarget <- length(targetNames)
            msg <- glue(
                "Unnamed `R` must have nrow = length(pvals); got ",
                "nrow(R) = {nR}, length(pvals) = {nTarget}."
            )
            abort(msg)
        }
    }
    R
}

# Internal: compute one method's combined p-value. Assumes inputs have
# already been validity-filtered (positive, finite, < 1) and aligned with R.
.combinePvalSingle <- function(method, pvals, zScores, R, methodArgs = list()) {
    args <- .combinePvalMethodArgs(methodArgs, method)
    switch(
        method,
        acat = pvalAcat(pvals),
        hmp = pvalHmp(pvals),
        bonferroni = min(length(pvals) * min(pvals), 1.0),
        fisher = ,
        stouffer = ,
        invchisq = pvalPoolr(pvals, method = method, R = R, methodArgs = args),
        gbj = ,
        bj = ,
        hc = ,
        ghc = ,
        minp = ,
        gbj_omni = pvalGbj(
            zScores,
            R = R,
            method = method,
            methodArgs = args
        ),
        aspu = pvalAspu(
            zScores = zScores,
            R = R,
            method = "aspu",
            methodArgs = args
        ),
        gates = pvalAspu(pvals = pvals, R = R, method = "gates"),
        .abortUnknownMethod("combination", method)
    )
}

# One method's entry out of the pvalMethodsConfig() record, or that method's
# own empty record when the caller configured nothing for it. Returning the
# constructor's result rather than list() keeps each engine's pecotmr
# defaults -- poolr's `adjust = "generalized"` -- in force either way.
# @noRd
.combinePvalMethodArgs <- function(methodArgs, method) {
    i <- match(method, names(methodArgs))
    if (!is.na(i)) {
        return(methodArgs[[i]])
    }
    ctors <- .pvalMethodCtors()
    if (is_in(method, names(ctors))) ctors[[method]]() else list()
}

#' Combine P-values via Any of a Menu of Methods
#'
#' Unified dispatcher for combining a vector of p-values (and/or z-scores) into
#' a single combined p-value. Supports independent-test methods (ACAT, HMP,
#' Bonferroni) and correlation-adjusted methods (Fisher / Stouffer /
#' inverse-chi-square via \code{poolr}; GBJ / BJ / HC / GHC / minP / GBJ-omnibus
#' via \code{GBJ}; aSPU / GATES via \code{aSPU}). Multiple methods may be
#' requested in a single call; the function returns a per-method result list
#' keyed by method name.
#'
#' Either \code{pvals} or \code{zScores} may be supplied. If only \code{zScores}
#' is given, two-sided p-values are derived as \code{p = 2 * (1 - pnorm(|z|))}.
#' Methods that require signed \code{zScores} (\code{gbj}, \code{bj}, \code{hc},
#' \code{ghc}, \code{minp}, \code{gbj_omni}, \code{aspu}) cannot be derived from
#' p-values alone and error if \code{zScores} is missing.
#'
#' Methods that require a correlation matrix \code{R} (\code{fisher},
#' \code{stouffer}, \code{invchisq}, \code{gbj}, \code{bj}, \code{hc},
#' \code{ghc}, \code{minp}, \code{gbj_omni}, \code{aspu}, \code{gates}) error if
#' \code{R} is missing. Methods that do not use \code{R} silently ignore it.
#' When \code{R} is named, it is realigned to match the order of \code{pvals} /
#' \code{zScores}, with a hard error if any entry is missing from \code{R}'s
#' names.
#'
#' An internal validity filter drops entries where \code{!is.finite(pvals) |
#' pvals <= 0 | pvals >= 1}; a warning is emitted when any are dropped.
#'
#' @param pvals Optional numeric vector of p-values. Required for methods that
#'   work on p-values; derivable from \code{zScores}.
#' @param zScores Optional numeric vector of signed z-scores.
#' @param methods Character vector of combination method names; see above for
#'   the menu (lowercase).
#' @param R Optional correlation matrix aligned to \code{pvals} /
#'   \code{zScores}. Required for the correlation-adjusted methods.
#' @param naRm Logical; if \code{TRUE} (default), drop NA p-values before
#'   combination.
#' @param methodArgs Per-method engine options, built with
#'   \code{\link{pvalMethodsConfig}} -- for example
#'   \code{pvalMethodsConfig(aspu = aspuConfig(n.perm = 5000))}. Configuring
#'   a method that \code{methods} does not request is an error rather than a
#'   silent no-op. \code{acat}, \code{hmp}, \code{bonferroni} and
#'   \code{gates} take no options.
#' @return A list with two elements:
#'   \describe{
#'     \item{input}{Summary of the call: \code{nPvalsIn},
#'       \code{nZScoresIn}, \code{nValid}, and the aligned
#'       \code{Raligned} matrix (or \code{NULL}).}
#'     \item{results}{Named list keyed by method, each element a list
#'       with \code{method} and \code{pval}.}
#'   }
#' @examples
#' combinePValues(pvals = c(0.01, 0.2, 0.5), methods = "fisher", R = diag(3))
#' @export
#' @importFrom checkmate assertFlag
combinePValues <- function(
    pvals = NULL,
    zScores = NULL,
    methods,
    R = NULL,
    naRm = TRUE,
    methodArgs = pvalMethodsConfig()
) {
    assertFlag(naRm)
    .assertMethodConfig(methodArgs, "pvalMethodsConfig", "methodArgs")
    methods <- .combinePvalCheckMethods(methods)
    .combinePvalCheckMethodArgs(methodArgs, methods)
    nPvalsIn <- if (is.null(pvals)) 0L else length(pvals)
    nZScoresIn <- if (is.null(zScores)) 0L else length(zScores)
    .combinePvalCheckPrereqs(methods, zScores, R)
    pvals <- .combinePvalDerivePvals(pvals, zScores)
    .combinePvalCheckLengths(pvals, zScores)
    filt <- .combinePvalFilter(pvals, zScores, naRm)
    targetNames <- .combinePvalTargetNames(pvals, zScores, R, filt$keep)
    Raligned <- .combinePvalAlignR(R, targetNames)
    input <- list(
        nPvalsIn = nPvalsIn,
        nZScoresIn = nZScoresIn,
        nValid = length(filt$pvalsK),
        Raligned = Raligned
    )
    if (length(filt$pvalsK) < 1L) {
        return(.combinePvalResult(.combinePvalNaMethods(methods), input))
    }
    perMethod <- .combinePvalRunMethods(
        methods,
        filt$pvalsK,
        filt$zScoresK,
        Raligned,
        methodArgs
    )
    .combinePvalResult(perMethod, input)
}

# Coerce + validate the requested methods against the known set.
# @noRd
.combinePvalCheckMethods <- function(methods) {
    if (missing(methods) || length(methods) == 0L) {
        known <- str_flatten(.combinePvalKnownMethods, ", ")
        msg <- glue("`methods` is required (one or more of: {known}).")
        abort(msg)
    }
    methods <- as.character(methods)
    unknown <- setdiff(methods, .combinePvalKnownMethods)
    if (length(unknown) > 0L) {
        unknownStr <- str_flatten(unknown, ", ")
        known <- str_flatten(.combinePvalKnownMethods, ", ")
        msg <- glue("Unknown method(s): {unknownStr}. Known: {known}")
        abort(msg)
    }
    methods
}

# Options for a method that was not requested would be dropped in silence --
# exactly the failure the constructors exist to prevent -- so say so instead.
# @noRd
.combinePvalCheckMethodArgs <- function(methodArgs, methods) {
    configured <- names(methodArgs)
    unused <- setdiff(configured, methods)
    if (length(unused) == 0L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "combinePValues: `methodArgs` configures ",
        "{str_flatten(unused, ', ')}, which `methods` does not request ",
        "({str_flatten(methods, ', ')}), so those options would be ignored."
    )
    abort(msg)
}

# Method-level prerequisites: z-score-requiring + R-requiring methods.
# @noRd
.combinePvalCheckPrereqs <- function(methods, zScores, R) {
    needZ <- intersect(methods, .combinePvalMethodsNeedingZ)
    if (length(needZ) > 0L && is.null(zScores)) {
        needZStr <- str_flatten(needZ, ", ")
        msg <- glue(
            "Method(s) {needZStr} require `zScores`; supplied input only ",
            "has pvals. Signed z-scores cannot be recovered from p-values ",
            "alone."
        )
        abort(msg)
    }
    needR <- intersect(methods, .combinePvalMethodsNeedingR)
    if (length(needR) > 0L && is.null(R)) {
        needRStr <- str_flatten(needR, ", ")
        msg <- glue(
            "Method(s) {needRStr} require an `R` correlation matrix; ",
            "got NULL."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Derive pvals from zScores (two-sided) when only zScores were supplied.
# @noRd
.combinePvalDerivePvals <- function(pvals, zScores) {
    if (is.null(pvals) && !is.null(zScores)) {
        pvals <- 2 * stats::pnorm(-abs(as.numeric(zScores)))
    }
    if (is.null(pvals)) {
        abort("Either `pvals` or `zScores` must be supplied.")
    }
    pvals
}

# @noRd
.combinePvalCheckLengths <- function(pvals, zScores) {
    if (!is.null(zScores) && length(zScores) != length(pvals)) {
        msg <- glue(
            "`pvals` and `zScores` must have the same length when both ",
            "supplied."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Internal validity filter + optional NA drop. Returns list(keep, pvalsK,
# zScoresK).
# @noRd
.combinePvalFilter <- function(pvals, zScores, naRm) {
    naMask <- is.na(pvals) | if (is.null(zScores)) FALSE else is.na(zScores)
    invalidMask <- !naMask & (!is.finite(pvals) | pvals <= 0 | pvals >= 1)
    dropMask <- (naRm & naMask) | invalidMask
    if (any(dropMask)) {
        nDrop <- sum(dropMask)
        nNa <- sum(naMask & dropMask)
        nInvalid <- sum(invalidMask)
        msg <- glue(
            "combinePValues: dropped {nDrop} entry/entries ",
            "({nNa} NA, {nInvalid} invalid)."
        )
        warn(msg)
    }
    keep <- !dropMask
    list(
        keep = keep,
        pvalsK = pvals[keep],
        zScoresK = if (is.null(zScores)) NULL else zScores[keep]
    )
}

# Stable target-name vector: pvals names, else zScores names, else R rownames,
# else positional labels.
# @noRd
.combinePvalTargetNames <- function(pvals, zScores, R, keep) {
    if (!is.null(names(pvals))) {
        return(names(pvals)[keep])
    }
    if (!is.null(zScores) && !is.null(names(zScores))) {
        return(names(zScores)[keep])
    }
    if (!is.null(R) && !is.null(rownames(R)) && nrow(R) == length(pvals)) {
        return(rownames(R)[keep])
    }
    as.character(seq_len(sum(keep)))
}

# All-NA per-method result (used when no valid entries survive filtering).
# @noRd
.combinePvalNaMethods <- function(methods) {
    set_names(map(methods, .combinePvalNaEntry), methods)
}

# @noRd
.combinePvalNaEntry <- function(m) {
    list(method = m, pval = NA_real_)
}

# Run each method, warning + NA on per-method failure.
# @noRd
.combinePvalRunMethods <- function(
    methods,
    pvalsK,
    zScoresK,
    Raligned,
    methodArgs
) {
    set_names(
        map(
            methods,
            .combinePvalRunOne,
            pvalsK = pvalsK,
            zScoresK = zScoresK,
            Raligned = Raligned,
            methodArgs = methodArgs
        ),
        methods
    )
}

# @noRd
#' @importFrom rlang try_fetch
.combinePvalRunOne <- function(m, pvalsK, zScoresK, Raligned, methodArgs) {
    p <- try_fetch(
        .combinePvalSingle(
            m,
            pvals = pvalsK,
            zScores = zScoresK,
            R = Raligned,
            methodArgs = methodArgs
        ),
        error = function(cnd) {
            msg <- glue("combinePValues: method '{m}' failed")
            warn(msg, parent = cnd)
            NA_real_
        }
    )
    list(method = m, pval = as.numeric(p))
}

# @noRd
.combinePvalResult <- function(perMethod, input) {
    list(input = input, results = perMethod)
}

# =============================================================================
# Engine option records
# =============================================================================

#' @title Options for the GBJ Family of Combination Tests
#' @description Build a checked record of extra arguments for the \pkg{GBJ}
#'   set-based tests behind \code{\link{combinePValues}} methods
#'   \code{"gbj"}, \code{"bj"}, \code{"hc"}, \code{"ghc"} and \code{"minp"}.
#'   All five share one signature, so one constructor covers them.
#' @param ... Arguments for \code{GBJ::GBJ()} and its siblings -- in practice
#'   \code{pairwise_cors}, the vector of pairwise correlations the test would
#'   otherwise derive from \code{cor_mat}. Names are checked against the
#'   engine's formals, so a misspelling is an error rather than a silent drop.
#'   \code{test_stats} and \code{cor_mat} are supplied by pecotmr and refused.
#' @return A \code{MethodConfig} record for
#'   \code{\link{pvalMethodsConfig}}.
#' @seealso \code{\link{combinePValues}}, \code{\link{gbjOmniConfig}}
#' @examples
#' gbjConfig(pairwise_cors = rep(0.1, 3))
#' @export
gbjConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            test_stats = "the caller's `zScores`",
            cor_mat = "the caller's `R`"
        ),
        "gbjConfig"
    )
    .newMethodConfig(
        c("GBJ::GBJ", "GBJ::BJ", "GBJ::HC", "GBJ::GHC", "GBJ::minP"),
        defaults = list(),
        extra = extra,
        label = "gbjConfig",
        engine = "gbj"
    )
}

#' @title Options for the GBJ Omnibus Combination Test
#' @description Build a checked record of extra arguments for
#'   \code{GBJ::OMNI_ss()}, the engine behind \code{\link{combinePValues}}
#'   method \code{"gbj_omni"}. It takes a different option from its siblings
#'   (\code{num_boots} rather than \code{pairwise_cors}), so it has its own
#'   constructor.
#' @param ... Arguments for \code{GBJ::OMNI_ss()} -- in practice
#'   \code{num_boots}, the bootstrap count the omnibus null is built from.
#'   \code{test_stats} and \code{cor_mat} are supplied by pecotmr and refused.
#' @return A \code{MethodConfig} record for
#'   \code{\link{pvalMethodsConfig}}.
#' @seealso \code{\link{combinePValues}}, \code{\link{gbjConfig}}
#' @examples
#' gbjOmniConfig(num_boots = 200)
#' @export
gbjOmniConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            test_stats = "the caller's `zScores`",
            cor_mat = "the caller's `R`"
        ),
        "gbjOmniConfig"
    )
    .newMethodConfig(
        "GBJ::OMNI_ss",
        defaults = list(),
        extra = extra,
        label = "gbjOmniConfig",
        engine = "gbjOmni"
    )
}

#' @title Options for the aSPU Combination Test
#' @description Build a checked record of extra arguments for
#'   \code{aSPU::aSPUs()}, the engine behind \code{\link{combinePValues}}
#'   method \code{"aspu"}.
#' @param ... Arguments for \code{aSPU::aSPUs()}: \code{pow} (the SPU powers
#'   the omnibus minimises over), \code{n.perm} (permutation count -- aSPU's
#'   p-value resolution is bounded by it), \code{Ps} and \code{prune}.
#'   \code{Zs} and \code{corSNP} are supplied by pecotmr and refused.
#' @return A \code{MethodConfig} record for
#'   \code{\link{pvalMethodsConfig}}.
#' @seealso \code{\link{combinePValues}}
#' @examples
#' aspuConfig(n.perm = 5000)
#' @export
aspuConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            Zs = "the caller's `zScores`",
            corSNP = "the caller's `R`"
        ),
        "aspuConfig"
    )
    .newMethodConfig(
        "aSPU::aSPUs",
        defaults = list(),
        extra = extra,
        label = "aspuConfig",
        engine = "aspu"
    )
}

#' @title Options for the poolr Combination Tests
#' @description Build a record of extra arguments for \code{poolr::fisher()},
#'   \code{poolr::stouffer()} and \code{poolr::invchisq()} -- the engines
#'   behind \code{\link{combinePValues}} methods \code{"fisher"},
#'   \code{"stouffer"} and \code{"invchisq"}. All three share one signature.
#' @param ... Arguments for the chosen \pkg{poolr} test: \code{adjust} (the
#'   correlation adjustment; pecotmr defaults it to \code{"generalized"}),
#'   \code{m}, \code{size}, \code{threshold}, \code{side}, \code{batchsize}
#'   and \code{nearpd}. The \pkg{poolr} signatures end in \code{...}, so
#'   names cannot be checked here and a misspelling is passed through.
#'   \code{p} and \code{R} are supplied by pecotmr and refused.
#' @return A \code{MethodConfig} record for
#'   \code{\link{pvalMethodsConfig}}.
#' @seealso \code{\link{combinePValues}}
#' @examples
#' poolrConfig(adjust = "empirical", side = 1)
#' @export
poolrConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            p = "the caller's `pvals`",
            R = "the caller's `R`"
        ),
        "poolrConfig"
    )
    .newMethodConfig(
        c("poolr::fisher", "poolr::stouffer", "poolr::invchisq"),
        defaults = list(adjust = "generalized"),
        extra = extra,
        label = "poolrConfig",
        engine = "poolr"
    )
}

# Which constructor configures each combinePValues() method. Keys that share
# an engine share a constructor: the five GBJ set tests have one signature,
# and so do the three poolr tests. A method absent from this map takes no
# options -- acat, hmp and bonferroni are pecotmr's own arithmetic, and
# gates's aSPU::GATES2() has no formal beyond the two pecotmr supplies.
# @noRd
.pvalMethodCtors <- function() {
    list(
        fisher = poolrConfig,
        stouffer = poolrConfig,
        invchisq = poolrConfig,
        gbj = gbjConfig,
        bj = gbjConfig,
        hc = gbjConfig,
        ghc = gbjConfig,
        minp = gbjConfig,
        gbj_omni = gbjOmniConfig,
        aspu = aspuConfig
    )
}

#' @title Per-Method Options for combinePValues
#' @description Bundle per-method engine options for
#'   \code{\link{combinePValues}}. Each entry is named for the method it
#'   configures and built with that method's constructor, so an option meant
#'   for one engine cannot silently reach another.
#' @param ... Named entries, one per method: \code{fisher}, \code{stouffer},
#'   \code{invchisq} (\code{\link{poolrConfig}}); \code{gbj}, \code{bj},
#'   \code{hc}, \code{ghc}, \code{minp} (\code{\link{gbjConfig}});
#'   \code{gbj_omni} (\code{\link{gbjOmniConfig}}); \code{aspu}
#'   (\code{\link{aspuConfig}}). A plain named list is accepted and passed
#'   through that method's constructor. \code{acat}, \code{hmp},
#'   \code{bonferroni} and \code{gates} take no options and are not valid
#'   names here.
#' @return A \code{MethodConfig} record for \code{combinePValues(methodArgs =)}.
#' @seealso \code{\link{combinePValues}}
#' @examples
#' pvalMethodsConfig(aspu = aspuConfig(n.perm = 5000))
#' @export
pvalMethodsConfig <- function(...) {
    .newNestedConfig(list(...), .pvalMethodCtors(), "pvalMethodsConfig")
}
