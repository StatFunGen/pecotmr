# =============================================================================
# mash / udr engine interface
# -----------------------------------------------------------------------------
# Every call into mashr, udr, CorShrink and flashier lives here;
# mashPipeline.R orchestrates the stages, holds the Param bundles and
# prepares the input. The split is mechanical, so a reader can check it: a
# block belongs here if it calls one of those packages, or is the Options bag
# that forwards arguments to such a call.
#
# `mashResidualCorrelation()`, `mashCovarianceComponents()` and
# `mashPriorCovariances()` deliberately stay in mashPipeline.R: they name
# pipeline stages and dispatch among the helpers below without calling an
# engine themselves.
# =============================================================================

# method 'simple': mashr's estimate_null_correlation_simple on the null set.
# @noRd
.mashResidCorSimple <- function(
    sumStatsList,
    alpha,
    inputScale,
    mashDataArgs,
    corArgs
) {
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
    exec(
        mashr::estimate_null_correlation_simple,
        .mashSetData(nullMats$b, nullMats$s, alpha, NULL, mashDataArgs),
        !!!as.list(corArgs)
    )
}

# mash_set_data with the caller's own settings spliced in. V is passed only
# when supplied: mash_set_data's own default differs from an explicit NULL.
# @noRd
.mashSetData <- function(b, s, alpha, vhat, mashDataArgs) {
    args <- c(
        list(Bhat = b, Shat = s, alpha = alpha),
        if (!is.null(vhat)) list(V = vhat),
        as.list(mashDataArgs)
    )
    exec(mashr::mash_set_data, !!!args)
}

# method 'mle': EM refinement of V against the prior U over a random subset.
# @noRd
.mashResidCorMle <- function(
    sumStatsList,
    alpha,
    inputScale,
    priorCovariances,
    nSubset,
    maxIter,
    mashDataArgs,
    corArgs
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
    priorCovariances <- .mashAsUlist(priorCovariances)
    randomMats <- .mashSumStatsToMatrices(
        sumStatsList$random,
        "random",
        inputScale = inputScale
    )
    n <- nrow(randomMats$b)
    idx <- sample(seq_len(n), min(nSubset, n))
    dsub <- .mashSetData(
        randomMats$b[idx, , drop = FALSE],
        randomMats$s[idx, , drop = FALSE],
        alpha,
        NULL,
        mashDataArgs
    )
    fit <- exec(
        mashr::mash_estimate_corr_em,
        dsub,
        priorCovariances,
        max_iter = maxIter,
        details = TRUE,
        !!!as.list(corArgs)
    )
    fit$V
}

# methods 'corshrink' / 'simpleSpecific': estimate V on the null z-matrix. The
# `null` partition is already the null variants (max|z| < 2), so no
# re-thresholding is needed.
# @noRd
.mashResidCorNullBased <- function(
    sumStatsList,
    inputScale,
    method,
    corArgs
) {
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
    if (method == "simpleSpecific") {
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
        msg <- glue(
            "mashResidualCorrelation: method 'corshrink' needs the CorShrink ",
            "package. Install it, or use 'simple' / 'simpleSpecific'."
        )
        abort(msg)
    }
    as.matrix(
        exec(CorShrink::CorShrinkData, nullZ, !!!as.list(corArgs))$cor
    )
}

# Internal: build the requested data-driven covariance components off a prepared
# mashr data object, in a fixed order (canonical, pca, flash, flashNonneg) so a
# given `components` set reproduces the same ordering everywhere. RNG is
# consumed only by cov_flash / cov_flash(nonneg). Shared by mashPriorCovariances
# and the exported mashCovarianceComponents.
# @noRd
# Build the requested components, keeping the two roles apart: `canonical`
# are fixed structural hypotheses that go to the mash fit unrefined, while
# `dataDriven` are the generator output the engine refines. mashr's own eQTL
# vignette draws exactly this line -- cov_ed sees only the data-driven set.
# @noRd
.mashBuildComponents <- function(
    mashData,
    components,
    nPcs = NULL,
    componentArgs = NULL
) {
    npc <- nPcs %||% (ncol(mashData$Bhat) - 1)
    argsFor <- function(key) as.list(componentArgs[[key]] %||% list())
    canonical <- if (is_in("canonical", components)) {
        exec(
            mashr::cov_canonical,
            mashData,
            !!!compact(argsFor("canonical"))
        )
    }
    pca <- if (is_in("pca", components)) {
        exec(mashr::cov_pca, mashData, npc = npc, !!!compact(argsFor("pca")))
    }
    flash <- if (is_in("flash", components)) {
        exec(mashr::cov_flash, mashData, !!!compact(argsFor("flash")))
    }
    flashNonneg <- if (is_in("flashNonneg", components)) {
        exec(
            mashr::cov_flash,
            mashData,
            factors = "nonneg",
            !!!compact(argsFor("flashNonneg"))
        )
    }
    list(
        canonical = canonical,
        dataDriven = c(pca, flash, flashNonneg)
    )
}

# The covariance list and its mixture weights. A supplied `priorCovariances`
# short-circuits the cov_* chain; either way the weights come from a
# weights-only mash() fit unless the covariance step already produced them
# (the ED/UDR engines return their own).
# @noRd
.mashPriorFit <- function(
    mashData,
    priorCovariances,
    priorComponents,
    components,
    nPcs,
    engine,
    engineArgs,
    componentArgs,
    mashArgs
) {
    result <- if (!is.null(priorCovariances)) {
        .mashUserPriorCovariances(priorCovariances, mashData)
    } else {
        .mashDataDrivenCovariances(
            mashData,
            priorComponents,
            components,
            nPcs,
            engine,
            engineArgs = engineArgs,
            componentArgs = componentArgs
        )
    }
    w <- result$w %||%
        mashr::get_estimated_pi(exec(
            mashr::mash,
            mashData,
            Ulist = result$U,
            outputlevel = 1,
            !!!as.list(mashArgs)
        ))
    list(U = result$U, w = w, loglik = result$loglik)
}

#' @title Arguments For mashr's Canonical Covariance Components
#' @description Options for \code{mashr::cov_canonical}, which builds the
#'   fixed structural hypotheses (identity, singletons, equal effects, simple
#'   heterogeneity) from the condition count. These are hypotheses about
#'   shape, so they go to the mash fit unrefined.
#' @param cov_methods Which canonical components to build. \code{NULL}
#'   (default) leaves mashr's own selection in place.
#' @param ... Any other \code{mashr::cov_canonical} argument.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CovCanonicalOptions(cov_methods = c("identity", "equal_effects"))
#' @export
CovCanonicalOptions <- function(cov_methods = NULL, ...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "supplied from the data by the pipeline"
        ),
        "CovCanonicalOptions"
    )
    .newMethodOptions(
        "mashr::cov_canonical",
        defaults = list(cov_methods = cov_methods),
        extra = list(...),
        label = "CovCanonicalOptions",
        engine = "canonical"
    )
}

#' @title Arguments For mashr's PCA Covariance Components
#' @description Options for \code{mashr::cov_pca}, a data-driven generator.
#'   The number of components is not settable here: it comes from
#'   \code{mashPriorCovariances(nPcs = )}, which the caller also reports on.
#' @param subset Rows of the data to use. \code{NULL} (default) uses all.
#' @param ... Any other \code{mashr::cov_pca} argument.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CovPcaOptions(subset = 1:50)
#' @export
CovPcaOptions <- function(subset = NULL, ...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "supplied from the data by the pipeline"
        ),
        "CovPcaOptions"
    )
    .newMethodOptions(
        "mashr::cov_pca",
        defaults = list(subset = subset),
        extra = list(...),
        label = "CovPcaOptions",
        engine = "pca"
    )
}

#' @title Arguments For mashr's FLASH Covariance Components
#' @description Options for \code{mashr::cov_flash}, a data-driven generator.
#'   \code{factors} is not settable here: pecotmr sets it to distinguish the
#'   \code{"flash"} and \code{"flashNonneg"} components.
#' @param subset Rows of the data to use. \code{NULL} (default) uses all.
#' @param remove_singleton Drop singleton effect vectors. Default
#'   \code{FALSE}.
#' @param tag,output_model Passed through to \code{mashr::cov_flash}.
#' @param greedy_args,backfit_args Option lists forwarded to flashier's greedy
#'   and backfit stages; mashr nests these itself.
#' @param ... Any other \code{mashr::cov_flash} argument.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CovFlashOptions(remove_singleton = TRUE)
#' @export
CovFlashOptions <- function(
    subset = NULL,
    remove_singleton = FALSE,
    tag = NULL,
    output_model = NULL,
    greedy_args = list(),
    backfit_args = list(),
    ...
) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(data = "supplied from the data by the pipeline"),
        "CovFlashOptions"
    )
    .newMethodOptions(
        "mashr::cov_flash",
        defaults = list(
            subset = subset,
            remove_singleton = remove_singleton,
            tag = tag,
            output_model = output_model,
            greedy_args = greedy_args,
            backfit_args = backfit_args
        ),
        extra = list(...),
        label = "CovFlashOptions",
        engine = "flash"
    )
}

#' @title Arguments For Extreme Deconvolution Refinement
#' @description Options for \code{mashr::cov_ed}, which refines the
#'   data-driven components by Extreme Deconvolution. The components to refine
#'   are not settable here -- they are whatever \code{components} produced.
#'
#'   \code{mashr::cov_ed} takes \code{...}, so argument names \strong{cannot
#'   be checked} against it; see \code{\link{MethodOptions}}.
#' @param subset Rows of the data to use. \code{NULL} (default) uses all.
#' @param algorithm Deconvolution algorithm, \code{"bovy"} (default) or
#'   \code{"teem"}.
#' @param ... Any other argument \code{mashr::cov_ed} forwards.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CovEdOptions(algorithm = "teem")
#' @export
CovEdOptions <- function(subset = NULL, algorithm = "bovy", ...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "supplied from the data by the pipeline",
            Ulist_init = "the pipeline's data-driven priors"
        ),
        "CovEdOptions"
    )
    .newMethodOptions(
        "mashr::cov_ed",
        defaults = list(subset = subset, algorithm = algorithm),
        extra = list(...),
        label = "CovEdOptions",
        engine = "covEd"
    )
}

# udr's tunables live in a `control` list rather than in ud_fit's formals, so
# the valid names come from udr::ud_fit_control_default(). Still read live, so
# the check cannot drift from the installed udr.
# @noRd
.udrControlNames <- function() {
    if (!requireNamespace("udr", quietly = TRUE)) {
        return(NULL)
    }
    names(udr::ud_fit_control_default())
}

#' @title Arguments For udr's Initialisation Step
#' @description Arguments for \code{udr::ud_init()}, which sets up the
#'   mixture before \code{udr::ud_fit()} refines it. The data and the
#'   component matrices are pecotmr's to supply and are refused here; the
#'   remaining names --- \code{n_unconstrained}, \code{n_rank1} --- are
#'   checked against that function's live formals.
#'
#'   udr's own spellings, not pecotmr's: this is somebody else's signature,
#'   so \code{n_unconstrained} is written the way udr writes it.
#' @param ... Arguments for \code{udr::ud_init()}.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' UdInitOptions(n_unconstrained = 20L)
#' @export
UdInitOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            dat = "the mash data being refined",
            V = "the mash data's residual covariance",
            U_scaled = "the canonical components",
            U_unconstrained = "the data-driven components",
            control = "UdFitOptions()"
        ),
        "UdInitOptions"
    )
    .newMethodOptions(
        "udr::ud_init",
        defaults = list(),
        extra = extra,
        label = "UdInitOptions",
        engine = "udInit"
    )
}

#' @title Arguments For udr's Fitting Step
#' @description Arguments for \code{udr::ud_fit()}'s \code{control} list.
#'   udr's tunables do not live in \code{ud_fit}'s formals --- it takes an
#'   opaque \code{control} --- so the valid names come from
#'   \code{udr::ud_fit_control_default()} instead. Still a live source, not
#'   a transcription.
#' @param ... Any \code{udr} control field: \code{unconstrained.update},
#'   \code{maxiter}, \code{tol}, \code{tol.lik}, \code{penalty.type}, and
#'   the rest of \code{udr::ud_fit_control_default()}. The one exception is
#'   \code{lambda}, which pecotmr derives from the number of conditions and
#'   so refuses here.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' UdFitOptions(unconstrained.update = "ted", maxiter = 500)
#' @export
UdFitOptions <- function(...) {
    extra <- list(...)
    # `lambda` is the one control field pecotmr derives rather than fixes:
    # .mashUdControl() sets it to the number of conditions, and the merge
    # there lets a user value replace it. The rest are plain tunables, so
    # this bundle has exactly one owned name.
    .configRefuseOwned(
        extra,
        c(lambda = "the number of conditions in the mash data"),
        "UdFitOptions"
    )
    .udAssertUnconstrainedUpdate(extra)
    .newMethodOptions(
        NULL,
        defaults = list(),
        extra = extra,
        label = "UdFitOptions",
        engine = "udFit",
        accepted = .udrControlNames()
    )
}

#' @title Arguments For Unconstrained Deconvolution Refinement
#' @description Refining the data-driven components with \code{udr}.
#'
#'   The engine is two functions, so each one's arguments travel in their
#'   own bundle and are checked against their own live source:
#'   \code{init} against \code{udr::ud_init()}'s formals, \code{fit}
#'   against \code{udr::ud_fit_control_default()}'s names. Splitting them
#'   also retired two pecotmr spellings --- \code{unconstrainedUpdate} and
#'   \code{nUnconstrained} --- that had to be translated back to
#'   \code{unconstrained.update} and \code{n_unconstrained} on the way out.
#' @param init Arguments for \code{udr::ud_init()}, built with
#'   \code{\link{UdInitOptions}}.
#' @param fit Arguments for \code{udr::ud_fit()}'s control list, built with
#'   \code{\link{UdFitOptions}}.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CovUdrOptions(
#'     fit = UdFitOptions(unconstrained.update = "ted", maxiter = 500)
#' )
#' @export
CovUdrOptions <- function(
    init = UdInitOptions(),
    fit = UdFitOptions()
) {
    .assertMethodOptions(init, "UdInitOptions", "init")
    .assertMethodOptions(fit, "UdFitOptions", "fit")
    .newMethodOptions(
        NULL,
        defaults = list(init = init, fit = fit),
        extra = list(),
        label = "CovUdrOptions",
        engine = "covUdr"
    )
}

#' @title Settings For The mash Mixture-Weight Fit
#' @description Options forwarded to \code{mashr::mash}, the fit that
#'   estimates the mixture weights over the prior covariances. Checked
#'   against that function's live formals.
#'
#'   \code{data}, \code{Ulist} and \code{outputlevel} are refused: they are
#'   the inputs pecotmr assembles and the output level it needs.
#'   \code{seed} is refused too --- use the caller's \code{setSeed}, so one
#'   setting governs the whole run.
#' @param ... Arguments for \code{mashr::mash}, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MashOptions(nullweight = 10, optmethod = "mixSQP")
#' @export
MashOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "the sumstats pecotmr assembles",
            Ulist = "the prior covariances",
            outputlevel = "fixed by the calling entry point",
            seed = "the caller's `setSeed`"
        ),
        "MashOptions"
    )
    .newMethodOptions(
        "mashr::mash",
        defaults = list(),
        extra = extra,
        label = "MashOptions",
        engine = "mash"
    )
}

#' @title Settings For mashr Data Assembly
#' @description Options forwarded to \code{mashr::mash_set_data}, which turns
#'   the effect-size and standard-error matrices into mashr's data object.
#'
#'   \code{Bhat}, \code{Shat}, \code{alpha} and \code{V} are refused: the
#'   first two are the matrices pecotmr builds, and the last two are the
#'   caller's own \code{alpha} / \code{vhat} arguments.
#'
#'   \code{zero_Bhat_Shat_reset} defaults to \code{1000} here because that
#'   is the value pecotmr has always passed; it is now a setting rather than
#'   a literal.
#' @param ... Arguments for \code{mashr::mash_set_data}, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MashDataOptions(zero_Bhat_Shat_reset = 500)
#' @export
MashDataOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            Bhat = "the effect-size matrix pecotmr builds",
            Shat = "the standard-error matrix pecotmr builds",
            alpha = "the caller's `alpha`",
            V = "the caller's `vhat`"
        ),
        "MashDataOptions"
    )
    .newMethodOptions(
        "mashr::mash_set_data",
        defaults = list(zero_Bhat_Shat_reset = 1000),
        extra = extra,
        label = "MashDataOptions",
        engine = "mashData"
    )
}

#' @title Settings For mashr Posterior Computation
#' @description Options forwarded to
#'   \code{mashr::mash_compute_posterior_matrices}.
#'
#'   \code{g}, \code{data} and \code{output_posterior_cov} are refused: the
#'   first two are the fitted model and the data pecotmr assembles, and the
#'   third is the caller's \code{outputPosteriorCov}. \code{A} is refused as
#'   well --- the contrast matrix follows from \code{excludeCondition}.
#'   \code{seed} is refused; use \code{setSeed}.
#' @param ... Arguments for the posterior computation, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MashPosteriorOptions(pi_thresh = 1e-8)
#' @export
MashPosteriorOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            g = "the fitted mash model",
            data = "the sumstats pecotmr assembles",
            A = "derived from `excludeCondition`",
            output_posterior_cov = "the caller's `outputPosteriorCov`",
            seed = "the caller's `setSeed`"
        ),
        "MashPosteriorOptions"
    )
    .newMethodOptions(
        "mashr::mash_compute_posterior_matrices",
        defaults = list(),
        extra = extra,
        label = "MashPosteriorOptions",
        engine = "mashPosterior"
    )
}

#' @title Settings For The Simple Null-Correlation Estimator
#' @description Options forwarded to
#'   \code{mashr::estimate_null_correlation_simple}. \code{data} is refused:
#'   it is the null-partition data object pecotmr assembles.
#' @param ... Arguments for the estimator, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MashCorSimpleOptions(z_thresh = 3)
#' @export
MashCorSimpleOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(data = "the null-partition data pecotmr assembles"),
        "MashCorSimpleOptions"
    )
    .newMethodOptions(
        "mashr::estimate_null_correlation_simple",
        defaults = list(),
        extra = extra,
        label = "MashCorSimpleOptions",
        engine = "simple"
    )
}

#' @title Settings For The EM Null-Correlation Estimator
#' @description Options forwarded to \code{mashr::mash_estimate_corr_em}, the
#'   \code{method = "mle"} estimator.
#'
#'   \code{data}, \code{Ulist} and \code{max_iter} are refused: the first
#'   two are what pecotmr assembles, and \code{max_iter} is the caller's
#'   \code{maxIter}.
#'
#'   Not checkable: \code{mash_estimate_corr_em} ends in \code{...}, so any
#'   name is legal there and nothing can be rejected.
#' @param ... Arguments for the estimator, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MashCorEmOptions(tol = 1e-5)
#' @export
MashCorEmOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "the random-subset data pecotmr assembles",
            Ulist = "the prior covariances",
            max_iter = "the caller's `maxIter`"
        ),
        "MashCorEmOptions"
    )
    .newMethodOptions(
        "mashr::mash_estimate_corr_em",
        defaults = list(),
        extra = extra,
        label = "MashCorEmOptions",
        engine = "mle"
    )
}

#' @title Settings For The CorShrink Null-Correlation Estimator
#' @description Options forwarded to \code{CorShrink::CorShrinkData}, the
#'   \code{method = "corshrink"} estimator. \code{data} is refused: it is
#'   the null z-matrix pecotmr assembles.
#'
#'   \code{ash.control} and \code{image} default to what pecotmr has always
#'   passed, so they are settings now rather than literals.
#' @param ... Arguments for \code{CorShrinkData}, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CorShrinkOptions(nboot = 100)
#' @export
CorShrinkOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(data = "the null z-matrix pecotmr assembles"),
        "CorShrinkOptions"
    )
    .newMethodOptions(
        "CorShrink::CorShrinkData",
        defaults = list(
            ash.control = list(mixcompdist = "halfuniform"),
            image = "null"
        ),
        extra = extra,
        label = "CorShrinkOptions",
        engine = "corshrink"
    )
}

# cov_ed engine: mashr's exported extreme-deconvolution wrapper (default
# algorithm = bovy) refines the components; mash() (upstream) learns the
# mixture weights. Returns list(U, w = NULL, loglik = NULL).
# @noRd
.mashEngineCovEd <- function(mashData, dataDriven, engineArgs = NULL) {
    U.ed <- exec(
        mashr::cov_ed,
        mashData,
        Ulist_init = dataDriven,
        !!!compact(as.list(engineArgs %||% list()))
    )
    list(U = U.ed, w = NULL, loglik = NULL)
}

# ud_fit with a directed error for the ud_ted / non-i.i.d. (per-variant SE)
# incompatibility.
# @noRd
#' @importFrom rlang try_fetch
.mashUdFit <- function(fit0, fitOpts, nCond) {
    # is.matrix(fit$V) is the same test udr itself applies, and it reads the
    # POST-init fit: ud_init() is what expands a shared V per variant.
    control <- .mashUdControl(fitOpts, nCond, iid = is.matrix(fit0$V))
    try_fetch(
        udr::ud_fit(fit0, control = control, verbose = FALSE),
        error = function(cnd) {
            if (
                identical(control$unconstrained.update, "ted") &&
                    str_detect(conditionMessage(cnd), "i.i.d")
            ) {
                msg <- glue(
                    "mashPriorCovariances: the udr TED update needs ",
                    "i.i.d. data (a single shared V), which the beta scale ",
                    "does not provide (per-variant SE). Use ",
                    "CovUdrOptions(fit = UdFitOptions(unconstrained.update ",
                    "= 'ed')), or a z-scale input."
                )
                abort(msg, parent = cnd)
            }
            # Not the TED i.i.d. case rewrapped above -- re-raise the
            # original condition unchanged so unrelated udr failures surface
            # (and aren't swallowed as a NULL fit).
            cnd_signal(cnd)
        }
    )
}

# ud / ud_ted engine (OPT-IN; udr with known numerical issues). Seeds canonical
# as the scaled prior and generates n_unconstrained data-driven matrices,
# returning U + weights + loglik directly. Returns list(U, w, loglik).
# @noRd
.mashEngineUd <- function(mashData, dataDriven, canonical, engineArgs) {
    if (!requireNamespace("udr", quietly = TRUE)) {
        msg <- glue(
            "mashPriorCovariances: engine 'covUdr' needs the udr package. ",
            "Install it, or use the default 'covEd'."
        )
        abort(msg)
    }
    initOpts <- .mashUdHalf(engineArgs, "init", UdInitOptions)
    fitOpts <- .mashUdHalf(engineArgs, "fit", UdFitOptions)
    # udr distinguishes scaled components from unconstrained ones, which is
    # the same line mashr draws: the fixed hypotheses initialise `U_scaled`,
    # the data-driven components `U_unconstrained`. Previously the
    # data-driven set was discarded here and udr always started from
    # canonical alone.
    scaled <- canonical %||% mashr::cov_canonical(mashData)
    # n_unconstrained only applies when udr has to invent the unconstrained
    # components, so it is dropped outright when data-driven ones exist --
    # including one the caller set, which would otherwise be refused by udr.
    initArgs <- compact(list(
        U_scaled = scaled,
        U_unconstrained = if (length(dataDriven) > 0L) dataDriven,
        n_unconstrained = if (length(dataDriven) == 0L) {
            initOpts$n_unconstrained %||% 50L
        }
    ))
    userInit <- as.list(initOpts)[setdiff(names(initOpts), "n_unconstrained")]
    fit0 <- exec(udr::ud_init, mashData, !!!initArgs, !!!userInit)
    fit <- .mashUdFit(fit0, fitOpts, ncol(mashData$Bhat))
    list(U = map(fit$U, "mat"), w = fit$w, loglik = fit$loglik)
}

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
#' @param priorCovariances The prior covariance list (\code{U}) to fit with.
#'   Either shape the producers return is accepted: a bare named list of
#'   covariance matrices (\code{\link{mashCovarianceComponents}}) or the
#'   \code{list(U, w, loglik)} a \code{\link{mashPriorCovariances}} result
#'   carries, which is unwrapped to its \code{U}.
#' @param vhat Residual correlation matrix (\code{V}); \code{NULL} -> identity.
#' @param fitOn Partition to learn the mixture weights on: \code{"random"}
#'   (default, the standard unbiased choice) or \code{"strong"}.
#' @param outputLevel \code{mashr::mash()} \code{outputlevel} (default 4 -- the
#'   full model \code{\link{mashPosterior}} consumes).
#' @param inputScale SumStats -> matrix conversion scale.
#' @param mashDataArgs Extra arguments for
#'   \code{mashr::mash_set_data()}, built with
#'   \code{\link{MashDataOptions}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
#' @param mashArgs Extra arguments for \code{mashr::mash()}, built with
#'   \code{\link{MashOptions}} -- for example \code{nullweight},
#'   \code{optmethod} or \code{verbose}. \code{data}, \code{Ulist},
#'   \code{outputlevel} and \code{seed} are owned by pecotmr and are
#'   refused by the constructor.
#' @param setSeed Integer seed, or \code{NULL} to leave the ambient RNG
#'   untouched.
#' @return The fitted \pkg{mashr} model (the \code{mashr::mash()} object).
#' @seealso \code{\link{mashPosterior}}, \code{\link{mashPriorCovariances}}
#' @examples
#' data(mashInputExample)
#' mi <- mashInputExample
#' mk <- function(b, s) {
#'   qtlSumStatsFromBetaMatrix(as.matrix(mi[[b]]), as.matrix(mi[[s]]),
#'     studyName = "mash")
#' }
#' ssl <- list(strong = mk("strong.b", "strong.s"),
#'   random = mk("random.b", "random.s"), null = mk("null.b", "null.s"))
#' conds <- colnames(mi$strong.b)
#' vhat <- diag(length(conds))
#' dimnames(vhat) <- list(conds, conds)
#' prior <- mashPriorCovariances(ssl, alpha = 0, vhat = vhat,
#'   components = "canonical")
#' model <- mashModelFit(ssl, alpha = 0, priorCovariances = prior,
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
    mashDataArgs = MashDataOptions(),
    mashArgs = MashOptions(),
    setSeed = 999
) {
    .assertMethodOptions(mashDataArgs, "MashDataOptions", "mashDataArgs")
    .assertMethodOptions(mashArgs, "MashOptions", "mashArgs")
    fitOn <- arg_match(fitOn)
    inputScale <- arg_match(inputScale)
    if (!requireNamespace("mashr", quietly = TRUE)) {
        abort("Package 'mashr' is required for this function.")
    }
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    priorCovariances <- .mashAsUlist(priorCovariances)
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
        inputScale,
        mashDataArgs
    )
    exec(
        mashr::mash,
        mashData,
        Ulist = priorCovariances,
        outputlevel = outputLevel,
        !!!as.list(mashArgs)
    )
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
#' @param mashDataArgs Extra arguments for
#'   \code{mashr::mash_set_data()}, built with
#'   \code{\link{MashDataOptions}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
#' @param posteriorArgs Extra arguments for
#'   \code{mashr::mash_compute_posterior_matrices()}, built with
#'   \code{\link{MashPosteriorOptions}} -- for example
#'   \code{algorithm.version}. \code{g}, \code{data}, \code{A},
#'   \code{output_posterior_cov} and \code{seed} are owned by pecotmr
#'   and are refused by the constructor.
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
#'     studyName = "mash")
#' }
#' ssl <- list(strong = mk("strong.b", "strong.s"),
#'   random = mk("random.b", "random.s"), null = mk("null.b", "null.s"))
#' conds <- colnames(mi$strong.b)
#' vhat <- diag(length(conds))
#' dimnames(vhat) <- list(conds, conds)
#' prior <- mashPriorCovariances(ssl, alpha = 0, vhat = vhat,
#'   components = "canonical")
#' model <- mashModelFit(ssl, alpha = 0, priorCovariances = prior,
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
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = MashDataOptions(),
    posteriorArgs = MashPosteriorOptions()
) {
    .assertMethodOptions(mashDataArgs, "MashDataOptions", "mashDataArgs")
    .assertMethodOptions(posteriorArgs, "MashPosteriorOptions", "posteriorArgs")
    inputScale <- arg_match(inputScale)
    if (!requireNamespace("mashr", quietly = TRUE)) {
        abort("Package 'mashr' is required for this function.")
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
    mashData <- .mashSetData(ex$b, ex$s, alpha, vhat, mashDataArgs)
    exec(
        mashr::mash_compute_posterior_matrices,
        ex$model,
        mashData,
        output_posterior_cov = outputPosteriorCov,
        !!!as.list(posteriorArgs)
    )
}
