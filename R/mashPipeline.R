#' @title Run mashr Across Multi-Context QTL or GWAS Summary Statistics
#' @description End-to-end driver: from `(strong, random, null)` sumstats
#'   collections, builds the variant x context Bhat / Shat matrices, estimates
#'   the residual correlation (\code{Vhat}), and fits the mash model with
#'   canonical + PCA + flash + ED covariance components, returning the fitted
#'   covariance list and the estimated mixture weights.
#' @param sumStatsList Named list (or \code{SimpleList}) of
#'   \code{\link{QtlSumStats}} or \code{\link{GwasSumStats}} objects. Required
#'   names: \code{"strong"} (discovery variants), \code{"random"} (random
#'   background). Optional: \code{"null"} (for residual correlation estimation).
#' @param alpha Numeric (length 1). Variance-stabilising-transform exponent
#'   forwarded to \code{mashr::mash_set_data()}. Use \code{alpha = 0} on the
#'   BETA scale, \code{alpha = 1} on the Z scale.
#' @param residualCorrelation Optional pre-computed residual correlation matrix
#'   (\code{Vhat}). When supplied, replaces the inline
#'   \code{mashr::estimate_null_correlation_simple()} call entirely and the
#'   \code{"random"} slot of \code{sumStatsList} becomes optional (the function
#'   does not need it for anything else). Useful when \code{Vhat} was estimated
#'   previously on a larger reference and shipped as a static artefact (the
#'   legacy MWE pattern).
#' @param priorCovariances Optional named list of square covariance matrices
#'   (the \code{Ulist} \code{mashr::mash()} consumes). When supplied, replaces
#'   the canonical + PCA + flash + ED chain (\code{cov_canonical} /
#'   \code{cov_pca} / \code{cov_flash} / \code{cov_ed}) entirely; mash sees only
#'   the supplied matrices. Every entry must be a \code{ncol(Bhat) x ncol(Bhat)}
#'   matrix. Useful when \code{U} was learnt on a larger reference and shipped
#'   as a static artefact (the legacy MWE pattern).
#' @param nPcs Optional integer; number of principal components seeded into
#'   \code{mashr::cov_pca()}. Defaults to \code{ncol(Bhat) - 1}. Ignored when
#'   \code{priorCovariances} is supplied.
#' @param inputScale One of \code{"auto"} (default), \code{"beta"},
#'   \code{"z"}. Controls which (Bhat, Shat) pair is extracted from each
#'   sumstats entry:
#'   \describe{
#'     \item{\code{"beta"}}{Bhat = BETA, Shat = SE -- the standard
#'       effect-size scale mashr was designed around. Requires every
#'       entry to carry BETA + SE mcols.}
#'     \item{\code{"z"}}{Bhat = Z, Shat = 1 -- z-score scale. Requires Z.}
#'     \item{\code{"auto"}}{Use BETA + SE when every entry carries both;
#'       otherwise fall back to (Z, 1) when every entry carries Z. Mixed
#'       inputs (some entries missing BETA, others missing Z) are a
#'       hard error.}
#'   }
#'   \code{alpha} should be chosen consistently with the resolved scale:
#'   typically \code{alpha = 0} for beta, \code{alpha = 1} for z.
#' @param setSeed Integer. RNG seed for reproducibility of
#'   \code{mashr::cov_flash} and \code{mashr::cov_ed}. Default 999.
#' @return A list with elements \code{U} (the combined covariance list:
#'   canonical + PCA + flash + ED) and \code{w} (the estimated mixture weights).
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' ss <- qtlSumStatsMulticontextExample
#' sumStatsList <- list(strong = ss, random = ss)
#' mashPipeline(sumStatsList, alpha = 0, nPcs = 2L)
#' @export
mashPipeline <- function(
    sumStatsList,
    alpha,
    residualCorrelation = NULL,
    priorCovariances = NULL,
    nPcs = NULL,
    inputScale = c("auto", "beta", "z"),
    setSeed = 999
) {
    inputScale <- arg_match(inputScale)
    .mashRequirePriorPackages()
    # Accept either a base list or a S4Vectors::SimpleList.
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    .mashValidateSumStatsList(sumStatsList, residualCorrelation)
    withr::local_seed(setSeed)
    vhat <- .mashResolveVhat(
        sumStatsList,
        alpha,
        inputScale,
        residualCorrelation
    )
    # mashPriorCovariances() owns the cov_* chain, the supplied-prior bypass,
    # and the mash() weight fit; mashPipeline just forwards its arguments.
    prior <- mashPriorCovariances(
        sumStatsList,
        alpha,
        vhat = vhat,
        priorCovariances = priorCovariances,
        nPcs = nPcs,
        inputScale = inputScale,
        setSeed = NULL
    )
    list(U = prior$U, w = prior$w)
}

# `sumStatsList` must be a named list of strong[/random/null] partitions;
# `random` is required only when vhat must be derived from data (no supplied
# residualCorrelation).
# @noRd
.mashValidateSumStatsList <- function(sumStatsList, residualCorrelation) {
    if (!is.list(sumStatsList) || is.null(names(sumStatsList))) {
        msg <- glue(
            "mashPipeline: `sumStatsList` must be a named list (or ",
            "SimpleList) of QtlSumStats / GwasSumStats objects, named with at ",
            "least 'strong' and 'random' (optionally 'null')."
        )
        abort(msg)
    }
    required <- if (is.null(residualCorrelation)) {
        c("strong", "random")
    } else {
        "strong"
    }
    missingNames <- setdiff(required, names(sumStatsList))
    if (length(missingNames) > 0L) {
        entrWord <- if (length(missingNames) == 1L) "entry" else "entries"
        msg <- glue(
            "mashPipeline: `sumStatsList` is missing required {entrWord}: ",
            "{str_flatten(shQuote(missingNames), ', ')}."
        )
        abort(msg)
    }
    extraNames <- setdiff(names(sumStatsList), c("strong", "random", "null"))
    if (length(extraNames) > 0L) {
        msg <- glue(
            "mashPipeline: `sumStatsList` has unrecognised entries: ",
            "{str_flatten(shQuote(extraNames), ', ')}. ",
            "Only 'strong', 'random', and 'null' are accepted."
        )
        abort(msg)
    }
}

# Vhat: the supplied residualCorrelation, else estimate it via
# mashResidualCorrelation ('simple' when a null set exists, else 'identity').
# setSeed = NULL leaves the RNG stream (seeded by the caller) untouched, so
# delegated calls consume it in the original order.
# @noRd
.mashResolveVhat <- function(
    sumStatsList,
    alpha,
    inputScale,
    residualCorrelation
) {
    if (!is.null(residualCorrelation)) {
        return(residualCorrelation)
    }
    hasNull <- is_in("null", names(sumStatsList)) && !is.null(sumStatsList$null)
    mashResidualCorrelation(
        sumStatsList,
        alpha,
        method = if (hasNull) "simple" else "identity",
        inputScale = inputScale,
        setSeed = NULL
    )
}

#' @title Estimate the mash Residual Correlation Matrix (Vhat)
#' @description Estimate the residual (null) correlation matrix \eqn{\hat V}
#'   that \code{mashr::mash_set_data()} consumes as \code{V}. Split out of
#'   \code{\link{mashPipeline}} so the estimation logic lives in one place and
#'   is reusable by the mixture-prior workflows.
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) of
#'   \code{\link{QtlSumStats}} / \code{\link{GwasSumStats}}. \code{method =
#'   "simple"} needs a \code{"null"} entry; \code{method = "mle"} needs
#'   \code{"random"}; \code{method = "identity"} reads only the condition count
#'   off \code{"strong"}.
#' @param alpha mash \code{alpha} (forwarded to \code{mashr::mash_set_data()}).
#' @param method Estimator, all on the \code{"null"} partition unless noted.
#'   \code{"simple"} = mashr \code{estimate_null_correlation_simple()};
#'   \code{"identity"} = \code{diag(nConditions)} (reads only \code{"strong"});
#'   \code{"simple_specific"} = \code{Matrix::nearPD(cov(nullZ), corr = TRUE)};
#'   \code{"corshrink"} = \code{CorShrink::CorShrinkData()} adaptive-shrinkage
#'   correlation (needs the \pkg{CorShrink} package); \code{"mle"} =
#'   \code{mashr::mash_estimate_corr_em()} on a random subset (needs
#'   \code{"random"} + \code{priorCovariances}).
#' @param priorCovariances Prior \code{U} list (required by \code{method =
#'   "mle"}).
#' @param nSubset,maxIter \code{method = "mle"} controls (random-subset size and
#'   EM iterations).
#' @param inputScale SumStats -> matrix conversion scale (\code{"auto"} /
#'   \code{"beta"} / \code{"z"}).
#' @param setSeed Integer seed, or \code{NULL} to leave the ambient RNG stream
#'   untouched (how \code{mashPipeline} keeps one continuous stream across its
#'   delegated calls).
#' @return A conditions x conditions residual correlation matrix.
#' @seealso \code{\link{mashPipeline}}, \code{\link{mashPriorCovariances}}
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' ss <- qtlSumStatsMulticontextExample
#' mashResidualCorrelation(list(strong = ss, random = ss), alpha = 0,
#'   method = "identity")
#' @export
mashResidualCorrelation <- function(
    sumStatsList,
    alpha,
    method = c("simple", "identity", "mle", "corshrink", "simple_specific"),
    priorCovariances = NULL,
    nSubset = 6000L,
    maxIter = 6L,
    inputScale = c("auto", "beta", "z"),
    setSeed = 999
) {
    method <- arg_match(method)
    inputScale <- arg_match(inputScale)
    if (!requireNamespace("mashr", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "To use this function, please install mashr: ",
            "https://cran.r-project.org/web/packages/mashr/index.html"
        )
        abort(msg)
        # nocov end
    }
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    if (!is.null(setSeed)) {
        withr::local_seed(setSeed)
    }
    if (method == "identity") {
        strongMats <- .mashSumStatsToMatrices(
            sumStatsList$strong,
            "strong",
            inputScale = inputScale
        )
        return(diag(rep(1, ncol(strongMats$b))))
    }
    if (method == "simple") {
        return(.mashResidCorSimple(sumStatsList, alpha, inputScale))
    }
    if (method == "mle") {
        return(.mashResidCorMle(
            sumStatsList,
            alpha,
            inputScale,
            priorCovariances,
            nSubset,
            maxIter
        ))
    }
    .mashResidCorNullBased(sumStatsList, alpha, inputScale, method)
}

# method 'simple': mashr's estimate_null_correlation_simple on the null set.
# @noRd
.mashResidCorSimple <- function(sumStatsList, alpha, inputScale) {
    if (is.null(sumStatsList$null)) {
        msg <- glue(
            "mashResidualCorrelation: method 'simple' requires a 'null' entry ",
            "in `sumStatsList`."
        )
        abort(msg)
    }
    nullMats <- .mashSumStatsToMatrices(
        sumStatsList$null,
        "null",
        inputScale = inputScale
    )
    mashr::estimate_null_correlation_simple(
        mashr::mash_set_data(
            nullMats$b,
            Shat = nullMats$s,
            alpha,
            zero_Bhat_Shat_reset = 1000
        )
    )
}

# method 'mle': EM refinement of V against the prior U over a random subset.
# @noRd
.mashResidCorMle <- function(
    sumStatsList,
    alpha,
    inputScale,
    priorCovariances,
    nSubset,
    maxIter
) {
    if (is.null(sumStatsList$random)) {
        msg <- glue(
            "mashResidualCorrelation: method 'mle' requires a 'random' entry ",
            "in `sumStatsList`."
        )
        abort(msg)
    }
    if (is.null(priorCovariances)) {
        msg <- glue(
            "mashResidualCorrelation: method 'mle' requires ",
            "`priorCovariances` ",
            "(the prior U the EM refines V against)."
        )
        abort(msg)
    }
    randomMats <- .mashSumStatsToMatrices(
        sumStatsList$random,
        "random",
        inputScale = inputScale
    )
    n <- nrow(randomMats$b)
    idx <- sample(seq_len(n), min(nSubset, n))
    dsub <- mashr::mash_set_data(
        randomMats$b[idx, , drop = FALSE],
        Shat = randomMats$s[idx, , drop = FALSE],
        alpha,
        zero_Bhat_Shat_reset = 1000
    )
    fit <- mashr::mash_estimate_corr_em(
        dsub,
        priorCovariances,
        max_iter = maxIter,
        details = TRUE
    )
    fit$V
}

# methods 'corshrink' / 'simple_specific': estimate V on the null z-matrix. The
# `null` partition is already the null variants (max|z| < 2), so no
# re-thresholding is needed.
# @noRd
.mashResidCorNullBased <- function(sumStatsList, alpha, inputScale, method) {
    if (is.null(sumStatsList$null)) {
        msg <- glue(
            "mashResidualCorrelation: method '{method}' requires a 'null' ",
            "entry in `sumStatsList` (the null variants V is estimated on)."
        )
        abort(msg)
    }
    nullMats <- .mashSumStatsToMatrices(
        sumStatsList$null,
        "null",
        inputScale = inputScale
    )
    nullZ <- nullMats$b / nullMats$s
    if (method == "simple_specific") {
        return(as.matrix(
            Matrix::nearPD(
                stats::cov(nullZ),
                conv.tol = 1e-06,
                doSym = TRUE,
                corr = TRUE
            )$mat
        ))
    }
    if (!requireNamespace("CorShrink", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "mashResidualCorrelation: method 'corshrink' needs the CorShrink ",
            "package. Install it, or use 'simple' / 'simple_specific'."
        )
        abort(msg)
        # nocov end
    }
    as.matrix(
        CorShrink::CorShrinkData(
            nullZ,
            ash.control = list(mixcompdist = "halfuniform"),
            image = "null"
        )$cor
    )
}

# Internal: build the requested data-driven covariance components off a prepared
# mashr data object, in a fixed order (canonical, pca, flash, flash_nonneg) so a
# given `components` set reproduces the same ordering everywhere. RNG is
# consumed only by cov_flash / cov_flash(nonneg). Shared by mashPriorCovariances
# and the exported mashCovarianceComponents.
# @noRd
.mashBuildComponents <- function(mashData, components, nPcs = NULL) {
    comps <- list()
    if (is_in("canonical", components)) {
        comps <- c(comps, mashr::cov_canonical(mashData))
    }
    if (is_in("pca", components)) {
        if (is.null(nPcs)) {
            nPcs <- ncol(mashData$Bhat) - 1
        }
        comps <- c(comps, mashr::cov_pca(mashData, npc = nPcs))
    }
    if (is_in("flash", components)) {
        comps <- c(comps, mashr::cov_flash(mashData))
    }
    if (is_in("flash_nonneg", components)) {
        comps <- c(comps, mashr::cov_flash(mashData, factors = "nonneg"))
    }
    comps
}

#' @title Build mash Data-Driven Covariance Components
#' @description Build the raw per-method data-driven covariance components
#'   (\code{cov_canonical} / \code{cov_pca} / \code{cov_flash} /
#'   \code{cov_flash(factors = "nonneg")}) for a \code{"strong"} SumStats set --
#'   the reusable building block \code{\link{mashPriorCovariances}} then refines
#'   with an ED / udr engine. Exposed on its own so a workflow can build or
#'   inspect a single component (the mixture-prior notebook's per-method steps)
#'   without running the full prior.
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) with at
#'   least \code{"strong"} (the discovery set the covariances are learned on).
#' @param alpha mash \code{alpha}.
#' @param vhat Residual correlation matrix (\code{V}); \code{NULL} -> identity.
#' @param components Any of \code{"canonical"}, \code{"pca"}, \code{"flash"},
#'   \code{"flash_nonneg"}. Built in that fixed order.
#' @param nPcs PCs seeded into \code{cov_pca}. Default \code{ncol(Bhat) - 1}.
#' @param inputScale SumStats -> matrix conversion scale.
#' @param setSeed Integer seed (\code{cov_flash} is stochastic), or \code{NULL}
#'   to leave the ambient RNG untouched.
#' @return A named list of covariance matrices (the concatenated components).
#' @seealso \code{\link{mashPriorCovariances}}
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' ss <- qtlSumStatsMulticontextExample
#' sumStatsList <- list(strong = ss, random = ss)
#' mashCovarianceComponents(sumStatsList, alpha = 0)
#' @export
mashCovarianceComponents <- function(
    sumStatsList,
    alpha,
    vhat = NULL,
    components = c("canonical", "pca", "flash", "flash_nonneg"),
    nPcs = NULL,
    inputScale = c("auto", "beta", "z"),
    setSeed = 999
) {
    inputScale <- arg_match(inputScale)
    .mashValidateComponents(components, "mashCovarianceComponents")
    .mashRequirePriorPackages()
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    if (!is.null(setSeed)) {
        withr::local_seed(setSeed)
    }
    mashData <- .mashMakeMashData(
        sumStatsList$strong,
        "strong",
        vhat,
        alpha,
        inputScale
    )
    .mashBuildComponents(mashData, components = components, nPcs = nPcs)
}

#' @title Estimate mash Prior Covariances and Mixture Weights
#' @description Build the mash prior covariance list \code{U} and mixture
#'   weights \code{w} from a \code{{strong[, random, null]}} SumStats
#'   collection. Split out of \code{\link{mashPipeline}} so the
#'   covariance-estimation chain lives in one place. The default builds every
#'   non-udr covariance component (\code{canonical + pca + flash +
#'   flash_nonneg}) and refines them with mashr extreme deconvolution
#'   (\code{cov_ed}).
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) with at
#'   least \code{"strong"} (the discovery set the covariances are learned on).
#' @param alpha mash \code{alpha}.
#' @param vhat Residual correlation matrix (\code{V}); typically from
#'   \code{\link{mashResidualCorrelation}}. \code{NULL} -> identity.
#' @param components Data-driven covariance components, any of
#'   \code{"canonical"}, \code{"pca"}, \code{"flash"} (default
#'   \code{cov_flash}), \code{"flash_nonneg"} (\code{cov_flash(factors =
#'   "nonneg")}). Built in that fixed order. Ignored by the \code{"ud"} /
#'   \code{"ud_ted"} engines.
#' @param engine Covariance-refinement engine. \code{"cov_ed"} (default; mashr's
#'   exported \code{cov_ed()} extreme deconvolution, whose default
#'   \code{algorithm = "bovy"} IS the Bovy et al. 2011 method -- weights from a
#'   final \code{mash()}); \code{"ud"} / \code{"ud_ted"} (\pkg{udr} ED / TED
#'   updates, returning weights directly -- OPT-IN, known numerical issues, so
#'   not the default; \code{"ud_ted"} additionally needs i.i.d. (z-scale) data).
#' @param nPcs PCs seeded into \code{cov_pca}. Default \code{ncol(Bhat) - 1}.
#' @param priorCovariances Optional caller-supplied prior \code{U}: a non-empty
#'   named list of \code{nCond x nCond} matrices. When supplied, the covariance
#'   chain is bypassed and only the mixture-weight fit runs.
#' @param priorComponents Optional caller-supplied raw covariance components (a
#'   non-empty named list, e.g. the concatenated
#'   \code{\link{mashCovarianceComponents}} outputs). When supplied, these are
#'   refined by \code{engine} instead of being rebuilt internally -- the
#'   mixture-prior pipeline where separate steps built the components.
#'   \code{components} / \code{nPcs} are then ignored. Distinct from
#'   \code{priorCovariances}, which bypasses the engine entirely.
#' @param udControl Named list overriding the \pkg{udr} controls for the
#'   \code{"ud"} / \code{"ud_ted"} engines: \code{n_unconstrained} (data-driven
#'   matrices to fit; default 50 -- the dominant cost, so reduce it for
#'   few-condition data), \code{maxiter} (default 1000), \code{tol},
#'   \code{tol.lik}. Ignored by \code{cov_ed}.
#' @param inputScale SumStats -> matrix conversion scale.
#' @param setSeed Integer seed, or \code{NULL} to leave the ambient RNG
#'   untouched.
#' @return \code{list(U, w, loglik)}: the covariance list, the
#'   \code{mashr::get_estimated_pi()} mixture weights, and the fit
#'   log-likelihood (\code{NULL} for the \code{cov_ed} engine).
#' @seealso \code{\link{mashPipeline}}, \code{\link{mashResidualCorrelation}}
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' ss <- qtlSumStatsMulticontextExample
#' sumStatsList <- list(strong = ss, random = ss)
#' mashPriorCovariances(sumStatsList, alpha = 0)
#' @export
mashPriorCovariances <- function(
    sumStatsList,
    alpha,
    vhat = NULL,
    components = c("canonical", "pca", "flash", "flash_nonneg"),
    engine = c("cov_ed", "ud", "ud_ted"),
    nPcs = NULL,
    priorCovariances = NULL,
    priorComponents = NULL,
    udControl = list(),
    inputScale = c("auto", "beta", "z"),
    setSeed = 999
) {
    engine <- arg_match(engine)
    inputScale <- arg_match(inputScale)
    .mashValidateComponents(components)
    .mashRequirePriorPackages()
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    if (!is.null(setSeed)) {
        withr::local_seed(setSeed)
    }
    mashData <- .mashMakeMashData(
        sumStatsList$strong,
        "strong",
        vhat,
        alpha,
        inputScale
    )
    result <- if (!is.null(priorCovariances)) {
        .mashUserPriorCovariances(priorCovariances, mashData)
    } else {
        .mashDataDrivenCovariances(
            mashData,
            priorComponents,
            components,
            nPcs,
            udControl,
            engine
        )
    }
    if (is.null(result$w)) {
        m <- mashr::mash(mashData, Ulist = result$U, outputlevel = 1)
        result$w <- mashr::get_estimated_pi(m)
    }
    list(U = result$U, w = result$w, loglik = result$loglik)
}

# Reject unknown prior-covariance component names.
# @noRd
.mashValidateComponents <- function(
    components,
    caller = "mashPriorCovariances"
) {
    valid <- c("canonical", "pca", "flash", "flash_nonneg")
    bad <- setdiff(as.character(components), valid)
    if (length(bad) > 0L) {
        msg <- glue(
            "{caller}: unknown component(s): ",
            "{str_flatten(bad, ', ')}. ",
            "Valid: {str_flatten(valid, ', ')}."
        )
        abort(msg)
    }
}

# mashr + flashier are required for the prior-covariance chain.
# @noRd
.mashRequirePriorPackages <- function() {
    if (!requireNamespace("mashr", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "To use this function, please install mashr: ",
            "https://cran.r-project.org/web/packages/mashr/index.html"
        )
        abort(msg)
        # nocov end
    }
    if (!requireNamespace("flashier", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "To use this function, please install flashier: ",
            "https://github.com/willwerscheid/flashier"
        )
        abort(msg)
        # nocov end
    }
}

# mashr mash_set_data over the STRONG effects (V defaults to identity).
# @noRd
.mashMakeMashData <- function(partition, label, vhat, alpha, inputScale) {
    mats <- .mashSumStatsToMatrices(partition, label, inputScale = inputScale)
    if (is.null(vhat)) {
        vhat <- diag(rep(1, ncol(mats$b)))
    }
    mashr::mash_set_data(
        mats$b,
        Shat = mats$s,
        V = vhat,
        alpha,
        zero_Bhat_Shat_reset = 1000
    )
}

# Caller-supplied prior covariance matrices (bypasses the cov_* chain; mashr
# sees only these). Validates they are a non-empty named list of nCond square
# matrices. Returns list(U, w = NULL, loglik = NULL).
# @noRd
.mashUserPriorCovariances <- function(priorCovariances, mashData) {
    if (
        !is.list(priorCovariances) ||
            length(priorCovariances) == 0L ||
            is.null(names(priorCovariances)) ||
            any(names(priorCovariances) == "")
    ) {
        msg <- glue(
            "mashPriorCovariances: `priorCovariances` must be a non-empty ",
            "named list of square covariance matrices."
        )
        abort(msg)
    }
    nCond <- ncol(mashData$Bhat)
    bad <- !map_lgl(priorCovariances, .mashCovIsSquareN, nCond = nCond)
    if (any(bad)) {
        msg <- glue(
            "mashPriorCovariances: every `priorCovariances` entry must ",
            "be a {nCond} x {nCond} matrix; offenders: ",
            "{str_flatten(names(priorCovariances)[bad], ', ')}."
        )
        abort(msg)
    }
    list(U = priorCovariances, w = NULL, loglik = NULL)
}

# Data-driven covariance chain: resolve the raw components (caller-supplied or
# built here), then refine + weight them with the chosen engine.
# @noRd
.mashDataDrivenCovariances <- function(
    mashData,
    priorComponents,
    components,
    nPcs,
    udControl,
    engine
) {
    comps <- .mashResolveComponents(
        priorComponents,
        mashData,
        components,
        nPcs
    )
    if (engine == "cov_ed") {
        .mashEngineCovEd(mashData, comps)
    } else {
        .mashEngineUd(mashData, engine, udControl)
    }
}

# Raw covariance components: caller-supplied `priorComponents` (validated) or
# freshly built via .mashBuildComponents.
# @noRd
.mashResolveComponents <- function(
    priorComponents,
    mashData,
    components,
    nPcs
) {
    if (is.null(priorComponents)) {
        return(.mashBuildComponents(
            mashData,
            components = components,
            nPcs = nPcs
        ))
    }
    if (
        !is.list(priorComponents) ||
            length(priorComponents) == 0L ||
            is.null(names(priorComponents)) ||
            any(names(priorComponents) == "")
    ) {
        msg <- glue(
            "mashPriorCovariances: `priorComponents` must be a non-empty ",
            "named list of covariance matrices (e.g. ",
            "mashCovarianceComponents() output)."
        )
        abort(msg)
    }
    priorComponents
}

# cov_ed engine: mashr's exported extreme-deconvolution wrapper (default
# algorithm = bovy) refines the components; mash() (upstream) learns the
# mixture weights. Returns list(U, w = NULL, loglik = NULL).
# @noRd
.mashEngineCovEd <- function(mashData, comps) {
    U.ed <- mashr::cov_ed(mashData, Ulist_init = comps)
    list(U = c(comps, U.ed), w = NULL, loglik = NULL)
}

# udr control list for the ud / ud_ted engines.
# @noRd
.mashUdControl <- function(engine, udControl, nCond) {
    list(
        unconstrained.update = if (engine == "ud_ted") "ted" else "ed",
        scaled.update = "fa",
        resid.update = "none",
        lambda = nCond,
        penalty.type = "iw",
        maxiter = udControl$maxiter,
        tol = udControl$tol,
        tol.lik = udControl$tol.lik
    )
}

# ud_fit with a directed error for the ud_ted / non-i.i.d. (per-variant SE)
# incompatibility.
# @noRd
.mashUdFit <- function(fit0, mashData, engine, udControl) {
    control <- .mashUdControl(engine, udControl, ncol(mashData$Bhat))
    tryCatch(
        udr::ud_fit(fit0, control = control, verbose = FALSE),
        error = function(e) {
            if (
                engine == "ud_ted" && str_detect(conditionMessage(e), "i.i.d")
            ) {
                msg <- glue(
                    "mashPriorCovariances: engine 'ud_ted' (udr TED update) ",
                    "needs i.i.d. data (a single shared V), which the beta ",
                    "scale does not provide (per-variant SE). Use engine 'ud' ",
                    "(ED update), or a z-scale input."
                )
                abort(msg)
            }
            # Not the ud_ted i.i.d. case rewrapped above -- re-raise the
            # original condition unchanged so unrelated udr failures surface
            # (and aren't swallowed as a NULL fit).
            cnd_signal(e)
        }
    )
}

# ud / ud_ted engine (OPT-IN; udr with known numerical issues). Seeds canonical
# as the scaled prior and generates n_unconstrained data-driven matrices,
# returning U + weights + loglik directly. Returns list(U, w, loglik).
# @noRd
.mashEngineUd <- function(mashData, engine, udControl) {
    if (!requireNamespace("udr", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "mashPriorCovariances: engine '{engine}' needs the udr package. ",
            "Install it, or use the default 'cov_ed'."
        )
        abort(msg)
        # nocov end
    }
    udControl <- utils::modifyList(
        list(
            n_unconstrained = 50L,
            maxiter = 1000L,
            tol = 1e-2,
            tol.lik = 1e-2
        ),
        udControl
    )
    U.can <- mashr::cov_canonical(mashData)
    fit0 <- udr::ud_init(
        mashData,
        n_unconstrained = udControl$n_unconstrained,
        U_scaled = U.can
    )
    fit <- .mashUdFit(fit0, mashData, engine, udControl)
    list(U = map(fit$U, "mat"), w = fit$w, loglik = fit$loglik)
}


# =============================================================================
# Mash model fit + posterior (mash_fit / mash_posterior notebooks)
# =============================================================================

#' @title Fit a mash Model for Posterior Computation
#' @description Fit a \pkg{mashr} model on a chosen partition using a supplied
#'   prior covariance list and residual correlation, returning the fitted model.
#'   This is the fit step (\code{mash_fit}'s first step): following Urbut et al.
#'   2019 the mixture weights are learned on the representative \code{"random"}
#'   partition, and the resulting model is then applied to the strong / target
#'   set by \code{\link{mashPosterior}}.
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) of
#'   \code{\link{QtlSumStats}} / \code{\link{GwasSumStats}}; must contain the
#'   \code{fitOn} entry.
#' @param alpha mash \code{alpha} (forwarded to \code{mashr::mash_set_data()}).
#' @param priorCovariances The prior covariance list (\code{U}) to fit with --
#'   e.g. the \code{$U} of a \code{\link{mashPriorCovariances}} result.
#' @param vhat Residual correlation matrix (\code{V}); \code{NULL} -> identity.
#' @param fitOn Partition to learn the mixture weights on: \code{"random"}
#'   (default, the standard unbiased choice) or \code{"strong"}.
#' @param outputLevel \code{mashr::mash()} \code{outputlevel} (default 4 -- the
#'   full model \code{\link{mashPosterior}} consumes).
#' @param inputScale SumStats -> matrix conversion scale.
#' @param setSeed Integer seed, or \code{NULL} to leave the ambient RNG
#'   untouched.
#' @return The fitted \pkg{mashr} model (the \code{mashr::mash()} object).
#' @seealso \code{\link{mashPosterior}}, \code{\link{mashPriorCovariances}}
#' @examples
#' data(mashInputExample)
#' mi <- mashInputExample
#' mk <- function(b, s) {
#'   qtlSumStatsFromBetaMatrix(as.matrix(mi[[b]]), as.matrix(mi[[s]]),
#'     study = "mash")
#' }
#' ssl <- list(strong = mk("strong.b", "strong.s"),
#'   random = mk("random.b", "random.s"), null = mk("null.b", "null.s"))
#' conds <- colnames(mi$strong.b)
#' vhat <- diag(length(conds))
#' dimnames(vhat) <- list(conds, conds)
#' prior <- mashPriorCovariances(ssl, alpha = 0, vhat = vhat,
#'   components = "canonical")
#' model <- mashModelFit(ssl, alpha = 0, priorCovariances = prior$U,
#'   vhat = vhat)
#' @export
mashModelFit <- function(
    sumStatsList,
    alpha,
    priorCovariances,
    vhat = NULL,
    fitOn = c("random", "strong"),
    outputLevel = 4L,
    inputScale = c("auto", "beta", "z"),
    setSeed = 999
) {
    fitOn <- arg_match(fitOn)
    inputScale <- arg_match(inputScale)
    if (!requireNamespace("mashr", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "To use this function, please install mashr: ",
            "https://cran.r-project.org/web/packages/mashr/index.html"
        )
        abort(msg)
        # nocov end
    }
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    .mashValidatePriorCovList(priorCovariances)
    if (is.null(sumStatsList[[fitOn]])) {
        msg <- glue(
            "mashModelFit: `sumStatsList` has no '{fitOn}' entry to fit on."
        )
        abort(msg)
    }
    if (!is.null(setSeed)) {
        withr::local_seed(setSeed)
    }
    mashData <- .mashMakeMashData(
        sumStatsList[[fitOn]],
        fitOn,
        vhat,
        alpha,
        inputScale
    )
    mashr::mash(mashData, Ulist = priorCovariances, outputlevel = outputLevel)
}

# `priorCovariances` must be a non-empty named list of covariance matrices.
# @noRd
.mashValidatePriorCovList <- function(priorCovariances) {
    if (
        is.null(priorCovariances) ||
            !is.list(priorCovariances) ||
            length(priorCovariances) == 0L ||
            is.null(names(priorCovariances))
    ) {
        msg <- glue(
            "mashModelFit: `priorCovariances` must be a non-empty named list ",
            "of covariance matrices (e.g. mashPriorCovariances()$U)."
        )
        abort(msg)
    }
}

#' @title Compute mash Posterior Matrices for a Target Set
#' @description Apply a fitted \pkg{mashr} model (from
#'   \code{\link{mashModelFit}}) to a target SumStats set, returning the
#'   posterior matrices (\code{PosteriorMean}, \code{PosteriorSD}, \code{lfsr},
#'   ...). This is the posterior step of the mash workflow -- \code{mash_fit}'s
#'   second step and \code{mash_posterior}'s per-analysis-unit computation.
#' @param model A fitted \pkg{mashr} model (the \code{\link{mashModelFit}}
#'   output).
#' @param sumStats A single \code{\link{QtlSumStats}} /
#'   \code{\link{GwasSumStats}} -- the target (e.g. strong) set to compute
#'   posteriors on.
#' @param alpha mash \code{alpha} (match the fit).
#' @param vhat Residual correlation matrix (\code{V}); \code{NULL} -> identity.
#' @param excludeCondition Character vector of condition (column) names to drop
#'   from the target AND the model before computing posteriors (the model's
#'   covariances are resized via \code{\link{updateMashModelCov}}). Default
#'   none.
#' @param outputPosteriorCov Return the full posterior covariance array (needed
#'   by \code{\link{fitMashContrast}}). Default \code{TRUE}.
#' @param inputScale SumStats -> matrix conversion scale.
#' @return The \code{mashr::mash_compute_posterior_matrices()} result: a list of
#'   \code{PosteriorMean} / \code{PosteriorSD} / \code{lfsr} /
#'   \code{NegativeProb} (+ \code{PosteriorCov} when \code{outputPosteriorCov =
#'   TRUE}).
#' @seealso \code{\link{mashModelFit}}, \code{\link{fitMashContrast}}
#' @examples
#' data(mashInputExample)
#' mi <- mashInputExample
#' mk <- function(b, s) {
#'   qtlSumStatsFromBetaMatrix(as.matrix(mi[[b]]), as.matrix(mi[[s]]),
#'     study = "mash")
#' }
#' ssl <- list(strong = mk("strong.b", "strong.s"),
#'   random = mk("random.b", "random.s"), null = mk("null.b", "null.s"))
#' conds <- colnames(mi$strong.b)
#' vhat <- diag(length(conds))
#' dimnames(vhat) <- list(conds, conds)
#' prior <- mashPriorCovariances(ssl, alpha = 0, vhat = vhat,
#'   components = "canonical")
#' model <- mashModelFit(ssl, alpha = 0, priorCovariances = prior$U,
#'   vhat = vhat)
#' mashPosterior(model, mk("strong.b", "strong.s"), alpha = 0, vhat = vhat)
#' @export
mashPosterior <- function(
    model,
    sumStats,
    alpha,
    vhat = NULL,
    excludeCondition = character(0),
    outputPosteriorCov = TRUE,
    inputScale = c("auto", "beta", "z")
) {
    inputScale <- arg_match(inputScale)
    if (!requireNamespace("mashr", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "To use this function, please install mashr: ",
            "https://cran.r-project.org/web/packages/mashr/index.html"
        )
        abort(msg)
        # nocov end
    }
    mats <- .mashSumStatsToMatrices(sumStats, "target", inputScale = inputScale)
    ex <- .mashExcludeConditions(
        mats$b,
        mats$s,
        vhat,
        model,
        as.character(excludeCondition)
    )
    vhat <- if (is.null(ex$vhat)) diag(rep(1, ncol(ex$b))) else ex$vhat
    mashData <- mashr::mash_set_data(
        ex$b,
        Shat = ex$s,
        V = vhat,
        alpha,
        zero_Bhat_Shat_reset = 1000
    )
    mashr::mash_compute_posterior_matrices(
        ex$model,
        mashData,
        output_posterior_cov = outputPosteriorCov
    )
}

# Drop `excludeCondition` columns from the target b/s (positionally, since a
# supplied vhat may have no dimnames), subsetting vhat + the model's covariances
# to match. Returns list(b, s, vhat, model); a no-op when nothing is excluded.
# @noRd
.mashExcludeConditions <- function(b, s, vhat, model, excludeCondition) {
    allConditions <- colnames(b)
    if (length(excludeCondition) == 0L) {
        return(list(b = b, s = s, vhat = vhat, model = model))
    }
    bad <- setdiff(excludeCondition, allConditions)
    if (length(bad) > 0L) {
        msg <- glue(
            "mashPosterior: excludeCondition not found in the target: ",
            "{str_flatten(bad, ', ')}."
        )
        abort(msg)
    }
    keep <- setdiff(allConditions, excludeCondition)
    if (length(keep) == 0L) {
        abort("mashPosterior: excludeCondition drops every condition.")
    }
    keepIdx <- match(keep, allConditions)
    if (!is.null(vhat)) {
        vhat <- vhat[keepIdx, keepIdx, drop = FALSE]
    }
    list(
        b = b[, keepIdx, drop = FALSE],
        s = s[, keepIdx, drop = FALSE],
        vhat = vhat,
        model = updateMashModelCov(model, allConditions, keep)
    )
}

# =============================================================================
# Mash pairwise contrast functions
# =============================================================================

#' Create a pairwise contrast column
#'
#' Sets +1 for the first condition and -1 for the second in a zero vector. Used
#' as a building block for contrast design matrices.
#'
#' @param pair A length-2 character vector naming the two conditions to
#'   contrast.
#' @param template A named numeric vector of zeros with names matching all
#'   conditions.
#' @return The template vector with +1 at \code{pair[1]} and -1 at
#'   \code{pair[2]}.
#' @examples
#' makePairwiseContrastCol(c("a", "b"), "mean_contrast_")
#' @export
makePairwiseContrastCol <- function(pair, template) {
    template[pair[1]] <- 1
    template[pair[2]] <- -1
    template
}

#' Compute pairwise contrasts from mash posterior
#'
#' For a single variant (row index), computes deviation contrasts (each
#' condition vs grand mean) and all pairwise contrasts from the mash posterior
#' mean and covariance. Supports condition grouping for weighted contrasts.
#'
#' @param index Integer row index of the variant in the posterior matrices.
#' @param origMean Matrix of original effect sizes (variants x conditions). Used
#'   to determine which conditions are "tested" (non-zero).
#' @param posteriorMean Matrix of mash posterior means (variants x conditions).
#' @param posteriorVcov 3D array of posterior covariance matrices (conditions x
#'   conditions x variants).
#' @param grouping Named integer vector mapping condition names to group IDs.
#'   Conditions with the same positive group ID are treated as replicates (e.g.,
#'   multiple datasets for the same cell type). Use 0 for ungrouped. If NULL
#'   (default), all conditions are treated independently.
#' @return A single-row data.frame with columns \code{mean_contrast_*},
#'   \code{se_contrast_*}, \code{p_contrast_*} for both deviation and pairwise
#'   contrasts. Returns NULL if fewer than 2 tested conditions.
#' @importFrom stringr str_remove_all
#' @importFrom utils combn
#' @examples
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3, dimnames = list("v1", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("v1", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' fitMashContrast(1L, om, pm, pv)
#' @export
fitMashContrast <- function(
    index,
    origMean,
    posteriorMean,
    posteriorVcov,
    grouping = NULL
) {
    populationNames <- colnames(posteriorMean)
    if (!is.null(populationNames)) {
        populationNames <- str_remove_all(populationNames, "BETA_")
    }
    origMeanVector <- origMean[index, ]
    names(origMeanVector) <- populationNames
    tested <- names(origMeanVector[origMeanVector != 0])
    if (length(tested) < 2) {
        return(NULL)
    }
    nPop <- length(tested)
    grouping <- if (is.null(grouping)) {
        set_names(rep(0L, nPop), tested)
    } else {
        grouping[tested]
    }
    contrastDesign <- .mashContrastDesign(tested, nPop, grouping)
    pm <- posteriorMean[index, tested]
    pv <- posteriorVcov[tested, tested, index]
    contrastDiff <- drop(t(contrastDesign) %*% pm)
    contrastVcov <- t(contrastDesign) %*% pv %*% contrastDesign
    contrastSe <- sqrt(diag(contrastVcov))
    contrastP <- .zToPvalue(contrastDiff / contrastSe)
    .mashContrastDf(
        index,
        posteriorMean,
        colnames(contrastDesign),
        contrastDiff,
        contrastSe,
        contrastP
    )
}

# Contrast design matrix over the `tested` conditions: a single pairwise column
# for two conditions, else deviation + grouped pairwise contrasts.
# @noRd
.mashContrastDesign <- function(tested, nPop, grouping) {
    pairwiseVector <- set_names(rep(0, nPop), tested)
    if (nPop <= 2) {
        pairwiseVector[tested[1]] <- 1
        pairwiseVector[tested[2]] <- -1
        return(matrix(
            pairwiseVector,
            ncol = 1,
            dimnames = list(tested, str_c(tested[1], "_vs_", tested[2]))
        ))
    }
    dev <- .mashDeviationContrast(tested, nPop, grouping)
    pwAdj <- .mashPairwiseContrast(tested, grouping, pairwiseVector)
    cbind(dev / (nPop - 1), pwAdj)
}

# Deviation contrasts (each condition vs the mean of the rest), with grouped
# conditions sharing their deviation weight.
# @noRd
.mashDeviationContrast <- function(tested, nPop, grouping) {
    dev <- matrix(-1, nPop, nPop, dimnames = list(tested, tested))
    diag(dev) <- nPop - 1
    uniqueGroups <- unique(grouping)
    for (grp in uniqueGroups[uniqueGroups > 0]) {
        grpMask <- grouping == grp
        grpSize <- sum(grpMask)
        diag(dev)[grpMask] <- (nPop - 1) / grpSize
        dev[grpMask, grpMask] <- (nPop - 1) / grpSize
    }
    colnames(dev) <- str_c(tested, "_deviation")
    dev
}

# Pairwise (all-pairs) contrasts, with grouped conditions' contributions split
# evenly across the group.
# @noRd
.mashPairwiseContrast <- function(tested, grouping, pairwiseVector) {
    twoCombn <- combn(tested, 2)
    pwNames <- apply(twoCombn, 2, str_flatten, collapse = "_vs_")
    pw <- apply(twoCombn, 2, makePairwiseContrastCol, pairwiseVector)
    colnames(pw) <- pwNames
    pwAdj <- pw
    for (col in colnames(pw)) {
        groups <- str_split(col, "_vs_")[[1]]
        groupValues <- grouping[is_in(names(grouping), groups)]
        relevant <- names(groupValues[groupValues > 0])
        if (n_distinct(groupValues) > 1 && length(relevant) > 0) {
            for (dg in unique(groupValues[groupValues > 0])) {
                rowsInGroup <- names(grouping[grouping == dg])
                matchedRow <- rowsInGroup[is_in(rowsInGroup, groups)]
                if (length(matchedRow) > 0) {
                    pwAdj[rowsInGroup, col] <- pw[matchedRow, col] /
                        length(rowsInGroup)
                }
            }
        }
    }
    pwAdj
}

# One-row contrast table (mean / se / p per contrast column) for `index`.
# @noRd
.mashContrastDf <- function(
    index,
    posteriorMean,
    cnames,
    contrastDiff,
    contrastSe,
    contrastP
) {
    fid <- rownames(posteriorMean)[index]
    if (is.null(fid)) {
        fid <- as.character(index)
    }
    df <- tibble(feature_id = fid)
    # unname: contrast* carry contrastDesign colnames; tibble (unlike
    # data.frame) preserves a named scalar's name on the column.
    for (i in seq_along(cnames)) {
        df[[str_c("mean_contrast_", cnames[i])]] <- unname(contrastDiff[i])
        df[[str_c("se_contrast_", cnames[i])]] <- unname(contrastSe[i])
        df[[str_c("p_contrast_", cnames[i])]] <- unname(contrastP[i])
    }
    df
}

#' Posterior contrast table over an entire mash posterior
#'
#' Orchestrates \code{\link{fitMashContrast}} across every feature (variant) of
#' a mash posterior: aligns the original effect matrix to the posterior columns,
#' runs the per-feature deviation + pairwise contrast, row-binds the results
#' (aligning the union of contrast columns), and orders them
#' (\code{mean}/\code{se}/\code{p}, deviation before pairwise). Features with
#' fewer than two tested conditions are dropped.
#'
#' @param posteriorMean Numeric matrix (features x conditions) of posterior
#'   means (\code{PosteriorMean} from \code{\link{mashPosterior}}).
#' @param posteriorVcov Numeric array (conditions x conditions x features) of
#'   posterior covariances (\code{PosteriorCov}).
#' @param origMean Numeric matrix (features x conditions) of the original effect
#'   estimates (e.g. \code{bhat}); used to decide which conditions were tested
#'   per feature. Aligned to \code{posteriorMean}'s columns by name; \code{NaN}
#'   is treated as 0 (untested).
#' @param grouping Optional named integer vector assigning conditions to groups
#'   (0 = independent); forwarded to \code{\link{fitMashContrast}} so replicate
#'   populations of one cell type share weight. Names are condition labels.
#' @return A \code{tibble} (features x contrasts) with a \code{feature_id}
#'   column followed by \code{mean_contrast_*}, \code{se_contrast_*},
#'   \code{p_contrast_*} columns. Empty when nothing is testable.
#' @seealso \code{\link{fitMashContrast}},
#'   \code{\link{metaAnalysisPerCondition}}
#' @importFrom dplyr bind_rows select matches
#' @examples
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("v1", c("a", "b", "c")))
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3, dimnames = list("v1", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' mashPosteriorContrast(pm, pv, om)
#' @export
mashPosteriorContrast <- function(
    posteriorMean,
    posteriorVcov,
    origMean,
    grouping = NULL
) {
    origMean <- origMean[, colnames(posteriorMean), drop = FALSE]
    origMean[is.nan(origMean)] <- 0

    parts <- map(
        seq_len(nrow(posteriorMean)),
        fitMashContrast,
        origMean = origMean,
        posteriorMean = posteriorMean,
        posteriorVcov = posteriorVcov,
        grouping = grouping
    )
    parts <- compact(parts)
    if (length(parts) == 0L) {
        return(tibble())
    }

    # Each part carries a feature_id column; bind_rows aligns the union of
    # contrast columns across features and preserves feature_id.
    res <- bind_rows(parts)
    dplyr::select(
        res,
        "feature_id",
        dplyr::matches("mean_contrast.*deviation"),
        dplyr::matches("mean_contrast.*_vs_"),
        dplyr::matches("se_contrast.*deviation"),
        dplyr::matches("se_contrast.*_vs_"),
        dplyr::matches("p_contrast.*deviation"),
        dplyr::matches("p_contrast.*_vs_")
    )
}

# =============================================================================
# Mash model subsetting functions
# =============================================================================

#' Subset a fitted mash model to a subset of conditions
#'
#' Updates the prior covariance matrices (\code{Ulist}) and mixture weights
#' (\code{pi}) in a fitted \code{mashr} model to match a reduced set of
#' conditions. Handles condition-specific, identity, and data-driven covariance
#' components.
#'
#' @param mashModel A fitted mash model object (from \code{mashr::mash}).
#' @param allSamples Character vector of all original condition names.
#' @param samples Character vector of the conditions to retain.
#' @return The updated mash model with resized covariance matrices and pruned
#'   mixture weights.
#' @examples
#' data(mashInputExample)
#' mi <- mashInputExample
#' mk <- function(b, s) {
#'   qtlSumStatsFromBetaMatrix(as.matrix(mi[[b]]), as.matrix(mi[[s]]),
#'     study = "mash")
#' }
#' ssl <- list(strong = mk("strong.b", "strong.s"),
#'   random = mk("random.b", "random.s"), null = mk("null.b", "null.s"))
#' conds <- colnames(mi$strong.b)
#' vhat <- diag(length(conds))
#' dimnames(vhat) <- list(conds, conds)
#' prior <- mashPriorCovariances(ssl, alpha = 0, vhat = vhat,
#'   components = "canonical")
#' model <- mashModelFit(ssl, alpha = 0, priorCovariances = prior$U,
#'   vhat = vhat)
#' updateMashModelCov(model, allSamples = conds, samples = conds[1:3])
#' @export
updateMashModelCov <- function(mashModel, allSamples, samples) {
    cov <- mashModel$fitted_g$Ulist

    # Remove matrices for dropped conditions
    unwanted <- setdiff(allSamples, samples)
    for (d in names(cov)) {
        if (is_in(d, unwanted) || is_in(d, str_c("ED_", unwanted))) {
            cov[[d]] <- NULL
        }
    }

    # Resize remaining matrices to match retained conditions
    for (d in names(cov)) {
        if (is_in(d, samples)) {
            # Condition-specific: single 1 on diagonal
            m <- matrix(0, length(samples), length(samples))
            m[which(samples == d), which(samples == d)] <- 1
            cov[[d]] <- m
        } else if (d == "identity") {
            m <- matrix(0, length(samples), length(samples))
            m[1, 1] <- 1
            cov[[d]] <- m
        } else if (is.null(colnames(cov[[d]]))) {
            cov[[d]] <- cov[[d]][
                seq_along(samples),
                seq_along(samples)
            ]
        } else {
            cov[[d]] <- cov[[d]][samples, samples]
        }
        cov[[d]] <- as.matrix(cov[[d]])
    }

    mashModel$fitted_g$Ulist <- cov

    # Prune mixture weights for removed conditions
    for (s in unwanted) {
        dropIdx <- grep(s, names(mashModel$fitted_g$pi), fixed = TRUE)
        if (length(dropIdx) > 0) {
            mashModel$fitted_g$pi <- mashModel$fitted_g$pi[-dropIdx]
        }
    }

    mashModel
}

#' Subset mash data matrices to specific SNPs and conditions
#'
#' Slices the \code{bhat}, \code{sbhat}, and \code{Z} matrices by row (SNPs) and
#' column (samples/conditions), and correspondingly subsets the \code{vhat}
#' covariance matrix.
#'
#' @param data A mash data list with elements \code{bhat}, \code{sbhat},
#'   \code{Z} (matrices), and \code{snp} (character vector).
#' @param vhat A square covariance matrix (conditions x conditions).
#' @param snps Character vector of SNP IDs to retain (row names).
#' @param samples Character vector of condition names to retain (column names).
#' @return A list with \code{data} (sliced data list) and \code{vhat} (sliced
#'   covariance matrix).
#' @examples
#' cond <- c("brain", "blood", "muscle")
#' p <- 8
#' bhat <- matrix(rnorm(p * 3), p, 3,
#'   dimnames = list(paste0("v", 1:p), cond))
#' sbhat <- matrix(abs(rnorm(p * 3)) + 0.1, p, 3,
#'   dimnames = list(paste0("v", 1:p), cond))
#' dat <- list(bhat = bhat, sbhat = sbhat, Z = bhat / sbhat,
#'   snp = paste0("v", 1:p))
#' vhat <- diag(3)
#' dimnames(vhat) <- list(cond, cond)
#' sliceMashData(dat, vhat = vhat, snps = 1:4, samples = NULL)
#' @export
sliceMashData <- function(data, vhat, snps, samples) {
    data$bhat <- as.matrix(data$bhat[snps, samples])
    data$sbhat <- as.matrix(data$sbhat[snps, samples])
    data$Z <- as.matrix(data$Z[snps, samples])
    vhat <- as.matrix(vhat[samples, samples])
    data$snp <- data$snp[is_in(data$snp, snps)]
    colnames(data$bhat) <- colnames(data$sbhat) <- colnames(data$Z) <- colnames(
        vhat
    ) <- samples
    list(data = data, vhat = vhat)
}

#' Sanitize NaN/Inf values in mash data
#'
#' Replaces NaN in \code{bhat} with 0 and NaN/Inf in \code{sbhat} with 1e3
#' (indicating high uncertainty).
#'
#' @param data A mash data list with \code{bhat} and \code{sbhat} matrices.
#' @return The data list with sanitized values.
#' @examples
#' sanitizeMashData(list(strong = list(z = matrix(rnorm(9), 3, 3))))
#' @export
sanitizeMashData <- function(data) {
    data$bhat[is.nan(data$bhat)] <- 0
    data$sbhat[is.nan(data$sbhat) | is.infinite(data$sbhat)] <- 1e3
    data
}

#' Random-Effects Meta-Analysis of Mash Pairwise Contrasts, per Condition
#'
#' For each condition (context), gathers all pairwise contrast effect sizes and
#' standard errors involving that condition, then runs a random-effects
#' meta-analysis (via \code{metafor::rma}; pecotmr does not implement its own
#' meta-analysis). Intended to be run on the output of
#' \code{\link{fitMashContrast}} / \code{\link{mashPosteriorContrast}}.
#'
#' @param effectSizes Numeric matrix (features x contrasts) of contrast effect
#'   sizes. Column names must follow the pattern
#'   \code{mean_contrast_<conditionA>_vs_<conditionB>}.
#' @param seValues Numeric matrix (features x contrasts) of contrast standard
#'   errors. Must have the same dimensions and column names as
#'   \code{effectSizes}.
#' @param seCutoff Numeric; minimum SE below which a contrast is excluded from
#'   the meta-analysis for a given feature (default 0).
#' @param metaMethod Between-study variance estimator for the random-effects
#'   meta-analysis, forwarded to \code{metafor::rma(method = )}. Default
#'   \code{"DL"} (DerSimonian-Laird); other options include \code{"REML"},
#'   \code{"ML"}, \code{"EB"}.
#' @return A tibble with columns:
#'   \describe{
#'     \item{condition}{The condition (context) name.}
#'     \item{contrast}{Pairwise contrast name (without prefix), e.g.
#'       \code{conditionA_vs_conditionB}.}
#'     \item{meta_pvalue}{P-value from the random-effects meta-analysis.}
#'     \item{meta_effect}{Pooled absolute effect size estimate.}
#'     \item{meta_se}{Standard error of the pooled estimate.}
#'     \item{tau2}{Between-study variance estimate.}
#'     \item{I2}{Heterogeneity measure (proportion of variance due to
#'       between-study variance), in [0, 1].}
#'   }
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @examples
#' effectSizes <- matrix(rnorm(12), 4, 3)
#' seValues <- matrix(abs(rnorm(12)) + 0.1, 4, 3)
#' metaAnalysisPerCondition(effectSizes, seValues)
#' @export
metaAnalysisPerCondition <- function(
    effectSizes,
    seValues,
    seCutoff = 0,
    metaMethod = "DL"
) {
    stopifnot(identical(dim(effectSizes), dim(seValues)))
    stopifnot(identical(colnames(effectSizes), colnames(seValues)))
    contrasts <- str_remove(colnames(effectSizes), "^mean_contrast_")
    conditions <- unique(c(
        str_remove(contrasts, "_vs_.*"),
        str_remove(contrasts, ".*_vs_")
    ))
    rows <- list_flatten(map(
        conditions,
        .metaConditionRows,
        effectSizes = effectSizes,
        seValues = seValues,
        contrasts = contrasts,
        seCutoff = seCutoff,
        metaMethod = metaMethod
    ))
    bind_rows(rows)
}

# Meta-analysis tibbles for every contrast column involving `condition`.
# @noRd
.metaConditionRows <- function(
    condition,
    effectSizes,
    seValues,
    contrasts,
    seCutoff,
    metaMethod
) {
    idx <- which(str_detect(colnames(effectSizes), condition))
    if (length(idx) == 0) {
        return(list())
    }
    condEffects <- effectSizes[, idx, drop = FALSE]
    condSes <- seValues[, idx, drop = FALSE]
    condContrasts <- contrasts[idx]
    map(
        seq_along(condContrasts),
        .metaContrastRow,
        condition = condition,
        condContrasts = condContrasts,
        condEffects = condEffects,
        condSes = condSes,
        seCutoff = seCutoff,
        metaMethod = metaMethod
    )
}

# One contrast's random-effects meta-analysis row. Fewer than two finite,
# above-cutoff points -> a degenerate row (single point passed through, else
# all-NA).
# @noRd
.metaOneContrast <- function(
    condition,
    contrast,
    es,
    se,
    seCutoff,
    metaMethod
) {
    keep <- se > seCutoff & is.finite(es) & is.finite(se)
    es <- es[keep]
    se <- se[keep]
    if (length(es) < 2) {
        return(tibble(
            condition = condition,
            contrast = contrast,
            meta_pvalue = if (length(es) == 1) {
                .zToPvalue(es / se)
            } else {
                NA_real_
            },
            meta_effect = if (length(es) == 1) es else NA_real_,
            meta_se = if (length(es) == 1) se else NA_real_,
            tau2 = NA_real_,
            I2 = NA_real_
        ))
    }
    ma <- .rmaMeta(es, se, method = metaMethod)
    tibble(
        condition = condition,
        contrast = contrast,
        meta_pvalue = .zToPvalue(ma$mean / ma$se),
        meta_effect = ma$mean,
        meta_se = ma$se,
        tau2 = ma$tau2,
        I2 = ma$I2
    )
}

#' Feature score from deviation contrasts (random-effects meta per condition)
#'
#' For each condition's deviation contrast, meta-analyzes the per-variant
#' absolute effect sizes (random-effects, via \code{metafor::rma}) and returns
#' the pooled Z-score (\eqn{\hat\mu / \mathrm{se}}). One score per condition --
#' the "meta" feature score of \code{mash_posterior.ipynb}.
#'
#' @param contrastResult A contrast table from
#'   \code{\link{mashPosteriorContrast}} (variants x contrasts) carrying
#'   \code{mean_contrast_*_deviation} and \code{se_contrast_*_deviation}
#'   columns.
#' @param metaMethod Between-study variance estimator forwarded to
#'   \code{metafor::rma} (default \code{"REML"}).
#' @return A \code{data.frame} with \code{condition} and \code{zScore}.
#' @seealso \code{\link{nSignificantScore}}, \code{\link{scoreFromCs}}
#' @examples
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3, dimnames = list("v1", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("v1", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' cr <- fitMashContrast(1L, om, pm, pv)
#' calculateFeatureScores(cr, metaMethod = "mean")
#' @export
calculateFeatureScores <- function(contrastResult, metaMethod = "REML") {
    cr <- as_tibble(contrastResult)
    effCols <- names(cr)[str_detect(names(cr), "mean_contrast_.*deviation")]
    if (length(effCols) == 0L) {
        return(tibble(condition = character(0), zScore = numeric(0)))
    }
    scores <- map_dbl(
        effCols,
        .metaContrastZScore,
        cr = cr,
        metaMethod = metaMethod
    )
    tibble(
        condition = str_remove(
            str_remove(effCols, "^mean_contrast_"),
            "_deviation$"
        ),
        zScore = as.numeric(scores)
    )
}

#' Feature score from the fraction of significant deviation contrasts
#'
#' For each condition's deviation contrast, the proportion of variants whose
#' contrast p-value falls below \code{pCutoff} -- the "n-significant" feature
#' score of \code{mash_posterior.ipynb}. No meta-analysis is involved.
#'
#' @param contrastResult A contrast table from
#'   \code{\link{mashPosteriorContrast}} carrying \code{p_contrast_*_deviation}
#'   columns.
#' @param pCutoff Significance threshold (default 1e-5).
#' @return A \code{data.frame} with \code{condition} and \code{ratio}
#'   (\eqn{n_{sig} / n}).
#' @examples
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3, dimnames = list("v1", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("v1", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' cr <- fitMashContrast(1L, om, pm, pv)
#' nSignificantScore(cr, pCutoff = 0.05)
#' @export
nSignificantScore <- function(contrastResult, pCutoff = 1e-5) {
    cr <- as_tibble(contrastResult)
    pCols <- names(cr)[str_detect(names(cr), "p_contrast_.*deviation")]
    if (length(pCols) == 0L) {
        return(tibble(condition = character(0), ratio = numeric(0)))
    }
    ratios <- map_dbl(pCols, .metaContrastSigRatio, cr = cr, pCutoff = pCutoff)
    tibble(
        condition = str_remove(
            str_remove(pCols, "^p_contrast_"),
            "_deviation$"
        ),
        ratio = as.numeric(ratios)
    )
}

#' Feature score from fine-mapped credible sets
#'
#' Scores a condition using a fine-mapping table: takes the lead (max-PIP)
#' variant of each credible set, intersects with the contrast variants, picks
#' the least-significant lead by deviation p-value, and returns the maximum
#' pairwise \eqn{|\mathrm{mean}| / \mathrm{se}} at that variant -- the "finemap"
#' feature score of \code{mash_posterior.ipynb}.
#'
#' @param fineMapping A fine-mapping \code{data.frame} with \code{cs_order},
#'   \code{pip} and \code{variants} columns (credible-set index 0 = not in a
#'   CS).
#' @param contrastResults A contrast table from
#'   \code{\link{mashPosteriorContrast}} with a \code{feature_id} column.
#' @param condition Condition label whose deviation p-value column is used; if
#'   that column is absent and exactly one pairwise contrast exists, the
#'   pairwise p-value is used instead.
#' @return A single numeric score, or \code{NA} when nothing overlaps.
#' @examples
#' data(mashPosteriorExample)
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3, dimnames = list("v1", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("v1", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' cr <- fitMashContrast(1L, om, pm, pv)
#' scoreFromCs(fineMapping = mashPosteriorExample$fineMapping,
#'   contrastResults = cr, condition = "a")
#' @export
scoreFromCs <- function(fineMapping, contrastResults, condition) {
    css <- setdiff(unique(fineMapping$cs_order), 0)
    if (length(css) == 0L) {
        return(NA_real_)
    }
    leadRows <- bind_rows(map(css, .csLeadRow, fineMapping = fineMapping))
    cr <- as_tibble(contrastResults)
    cr <- filter(cr, is_in(.data$feature_id, leadRows$variants))
    if (nrow(cr) == 0L) {
        return(NA_real_)
    }

    pDevCol <- names(cr)[str_detect(
        names(cr),
        str_c("p_contrast_", condition, "_deviation")
    )]
    if (length(pDevCol) > 0L) {
        pCol <- pDevCol[1L]
    } else {
        pvCols <- names(cr)[str_detect(names(cr), "p_contrast.*_vs_")]
        if (length(pvCols) != 1L) {
            return(NA_real_)
        }
        pCol <- pvCols
    }
    maxRow <- filter(cr, .data[[pCol]] == max(.data[[pCol]], na.rm = TRUE))
    meanCols <- names(cr)[str_detect(names(cr), "mean_contrast.*_vs_")]
    seCols <- names(cr)[str_detect(names(cr), "se_contrast.*_vs_")]
    meanVs <- as.numeric(as.matrix(select(maxRow, all_of(meanCols))))
    seVs <- as.numeric(as.matrix(select(maxRow, all_of(seCols))))
    max(abs(meanVs / seVs), na.rm = TRUE)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# TRUE when a prior covariance is a square matrix of dimension nCond.
# @noRd
.mashCovIsSquareN <- function(M, nCond) {
    is.matrix(M) && nrow(M) == nCond && ncol(M) == nCond
}

# One meta-analysis tibble for contrast `i` of a condition's contrast columns.
# @noRd
.metaContrastRow <- function(
    i,
    condition,
    condContrasts,
    condEffects,
    condSes,
    seCutoff,
    metaMethod
) {
    .metaOneContrast(
        condition,
        condContrasts[i],
        abs(as.numeric(condEffects[, i])),
        as.numeric(condSes[, i]),
        seCutoff,
        metaMethod
    )
}

# Meta-analysis z-score (mean/se) for one mean-contrast column, NA when empty.
# @noRd
.metaContrastZScore <- function(ec, cr, metaMethod) {
    seCol <- str_replace(ec, "^mean_contrast", "se_contrast")
    if (!is_in(seCol, names(cr))) {
        return(NA_real_)
    }
    es <- abs(as.numeric(cr[[ec]]))
    se <- as.numeric(cr[[seCol]])
    keep <- is.finite(es) & is.finite(se) & se > 0
    es <- es[keep]
    se <- se[keep]
    if (length(es) < 1L) {
        return(NA_real_)
    }
    ma <- .rmaMeta(es, se, method = metaMethod)
    ma$mean / ma$se
}

# Fraction of a p-contrast column below `pCutoff` (NA when no observations).
# @noRd
.metaContrastSigRatio <- function(pc, cr, pCutoff) {
    p <- as.numeric(cr[[pc]])
    nTot <- sum(!is.na(p))
    if (nTot == 0L) NA_real_ else sum(p < pCutoff, na.rm = TRUE) / nTot
}

# The lead (max-PIP) variant row(s) of credible set `cs`.
# @noRd
.csLeadRow <- function(cs, fineMapping) {
    tmp <- filter(fineMapping, .data$cs_order == cs)
    filter(tmp, .data$pip == max(.data$pip))
}
