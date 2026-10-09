#' @include MethodParam.R
NULL

# Declared here rather than beside its constructor (far below) because
# MashComponentSelection references it and MashPriorParam needs that union.
#' @rdname MashComponentParam
#' @aliases MashComponentParam-class
#' @exportClass MashComponentParam
setClass(
    "MashComponentParam",
    contains = "MethodParam",
    slots = c(components = "list_OR_NULL")
)

# Both take either a name or that component/engine's own record, which is why
# neither is arg_match()ed in the constructor: collapsing one to a string
# would discard the options travelling with it.
#
# `components` is now exact -- only a MashComponentParam, not any
# MethodOptions. `engine` stays open to MethodOptions and always will:
# CovEdOptions() and CovUdrOptions() are external-engine argument bags, not
# pecotmr settings, so they never become Params.
# No character arm: a name is translated to the record it stands for by the
# constructor, so each slot has ONE shape.

# A name is a user-facing convenience, so it is translated to the record it
# stands for at the boundary rather than carried as a second slot shape.
# Both entry points use these: MashPriorParam() and mashPriorCovariances(),
# which is exported with the same convenience arguments.
#
# The alternative -- storing the string and resolving at use -- is what this
# replaced, and it made every consumer handle two shapes for one setting.
# @noRd
.mashNormalizeComponents <- function(v) {
    if (is(v, "MashComponentParam")) {
        return(v)
    }
    if (!is.character(v)) {
        abort(glue(
            "`components` must be component names or a ",
            "MashComponentParam() record; got {class(v)[[1L]]}"
        ))
    }
    .mashValidateComponents(v, "MashPriorParam")
    exec(
        MashComponentParam,
        !!!set_names(rep(list(list()), length(v)), v)
    )
}

# "none" has no constructor, so it is the absent engine: NULL. The default
# resolves to CovEdOptions(), so NULL in the slot is unambiguous.
# @noRd
.mashNormalizeEngine <- function(v) {
    if (.isMethodOptions(v)) {
        return(.mashAssertEngineRecord(v))
    }
    if (!is.character(v)) {
        abort(glue(
            "`engine` must be 'covEd', 'covUdr', 'none', or that engine's ",
            "Options record; got {class(v)[[1L]]}"
        ))
    }
    # arg_match() needs a symbol and would then report on the local name;
    # arg_match0() takes the user-facing name explicitly, so the message
    # says `engine` rather than whatever this binding is called.
    nm <- arg_match0(
        v[[1L]],
        c(names(.mashEngineCtors()), "none"),
        arg_nm = "engine"
    )
    if (identical(nm, "none")) {
        return(NULL)
    }
    .mashEngineCtors()[[nm]]()
}

# @noRd
.mashAssertEngineRecord <- function(v) {
    known <- names(.mashEngineCtors())
    got <- metadata(v)$engine
    if (is.null(got) || !is_in(got, known)) {
        abort(glue(
            "`engine`: these options are for '{got %||% 'unknown'}', which ",
            "does not refine mash prior covariances. Use ",
            "{str_flatten(str_c(known, 'Options()'), ' or ')}."
        ))
    }
    v
}

#' @rdname MashPriorParam
#' @aliases MashPriorParam-class
#' @exportClass MashPriorParam
setClass(
    "MashPriorParam",
    contains = "MethodParam",
    slots = c(
        priorCovariances = "list_OR_NULL",
        components = "MashComponentParam",
        engine = "MethodOptions_OR_NULL",
        nPcs = "numeric_OR_NULL"
    )
)

#' @title mash Prior-Covariance Settings
#' @description How the prior covariance matrices (the \code{Ulist} mash
#'   consumes) are obtained: supplied outright, or built from components and
#'   refined by an engine. Exactly the arguments
#'   \code{\link{mashPriorCovariances}} takes for that job, so the bundle
#'   travels there whole.
#'
#'   Supplying \code{priorCovariances} short-circuits the rest: the
#'   components are not built and the engine does not run.
#' @param priorCovariances Optional named list of square covariance matrices,
#'   or a \code{\link{MashPrior}}. \code{NULL} (default) builds them.
#' @param components Which prior-covariance components to build: any of
#'   \code{"canonical"}, \code{"pca"}, \code{"flash"},
#'   \code{"flashNonneg"}, or a \code{\link{MashComponentParam}} record to
#'   configure them. Ignored when \code{priorCovariances} is supplied.
#' @param engine How the data-driven components are refined:
#'   \code{"covEd"} (default), \code{"covUdr"} or \code{"none"}, or the
#'   matching constructor --- \code{\link{CovEdOptions}} /
#'   \code{\link{CovUdrOptions}} --- to configure it at the same time.
#'   Ignored when \code{priorCovariances} is supplied.
#' @param nPcs Optional integer; principal components seeded into
#'   \code{mashr::cov_pca()}. Defaults to \code{ncol} of the data. Read
#'   only when \code{components} includes \code{"pca"}.
#' @return A \code{MashPriorParam} object, a \code{\link{MethodParam}}.
#' @examples
#' MashPriorParam(components = c("canonical", "pca"), nPcs = 3)
#' @importFrom rlang arg_match0
#' @export
MashPriorParam <- function(
    priorCovariances = NULL,
    components = c("canonical", "pca", "flash", "flashNonneg"),
    engine = c("covEd", "covUdr", "none"),
    nPcs = NULL
) {
    new(
        "MashPriorParam",
        priorCovariances = priorCovariances,
        components = .mashNormalizeComponents(components),
        engine = .mashNormalizeEngine(engine),
        nPcs = nPcs
    )
}

#' @title Run mashr Across Multi-Context QTL or GWAS Summary Statistics
#' @description End-to-end driver: from `(strong, random, null)` sumstats
#'   collections, builds the variant x context Bhat / Shat matrices, estimates
#'   the residual correlation (\code{Vhat}), and fits the mash model with
#'   canonical + PCA + flash + ED covariance components, returning the fitted
#'   covariance list and the estimated mixture weights.
#' @param sumStatsList Named list (or \code{SimpleList}) of
#'   \code{\link{QtlSumStats}} or \code{\link{GwasSumStats}} objects. Required
#'   names: \code{"strong"} (discovery variants), \code{"random"} (random
#'   background, the partition the mixture weights are fit on). Optional:
#'   \code{"null"} (null variants), needed only by the
#'   \code{residualCorrelationMethod} values that estimate \eqn{\hat V} from
#'   them.
#' @param alpha Numeric (length 1). Variance-stabilising-transform exponent
#'   forwarded to \code{mashr::mash_set_data()}. Use \code{alpha = 0} on the
#'   BETA scale, \code{alpha = 1} on the Z scale.
#' @param residualCorrelation Optional pre-computed residual correlation matrix
#'   (\code{Vhat}). When supplied, it is used as-is:
#'   \code{residualCorrelationMethod} is not consulted, and the
#'   \code{"random"} slot of \code{sumStatsList} becomes optional (the
#'   function does not need it for anything else). Useful when \code{Vhat} was
#'   estimated previously on a larger reference and shipped as a static
#'   artefact (the legacy MWE pattern).
#' @param residualCorrelationMethod How to estimate \eqn{\hat V} when
#'   \code{residualCorrelation} is not supplied; forwarded to
#'   \code{\link{mashResidualCorrelation}}, which is where each estimator is
#'   described.
#'
#'   \code{"identity"} (the default) takes the residual correlation to be the
#'   identity, needs no extra partition, and never depends on which partitions
#'   were supplied. It is the safe default rather than the best one: it
#'   assumes conditions share no residual correlation, which
#'   overlapping-sample designs -- the usual multi-context QTL case, where the
#'   same donors are measured in every context -- violate. Supplying a
#'   \code{"null"} partition while leaving this at \code{"identity"} reports
#'   that the partition is unused, since that combination is more often an
#'   oversight than an intent.
#'
#'   The estimators differ in what they require: \code{"simple"},
#'   \code{"simpleSpecific"} and \code{"corshrink"} need a \code{"null"}
#'   entry; \code{"mle"} needs \code{"random"} plus a supplied
#'   \code{priorCovariances} to refine against. A named method whose
#'   requirement is unmet is a hard error, not a silent fallback.
#' @param prior How the prior covariance matrices are obtained, built with
#'   \code{\link{MashPriorParam}}: \code{priorCovariances} supplies them
#'   outright (short-circuiting the rest), otherwise \code{components}
#'   chooses which to build, \code{engine} how the data-driven ones are
#'   refined, and \code{nPcs} parametrises the \code{pca} component.
#'   Forwarded whole to \code{\link{mashPriorCovariances}}.
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
#' @param setSeed Integer. RNG seed for reproducibility of
#'   \code{mashr::cov_flash} and \code{mashr::cov_ed}. Default 999.
#' @return A list with elements \code{U} (the combined covariance list:
#'   canonical + PCA + flash + ED) and \code{w} (the estimated mixture weights).
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' ss <- qtlSumStatsMulticontextExample
#' sumStatsList <- list(strong = ss, random = ss)
#' mashPipeline(sumStatsList, alpha = 0, prior = MashPriorParam(nPcs = 2L))
#' @export
mashPipeline <- function(
    sumStatsList,
    alpha,
    residualCorrelation = NULL,
    residualCorrelationMethod = c(
        "identity",
        "simple",
        "simpleSpecific",
        "corshrink",
        "mle"
    ),
    prior = MashPriorParam(),
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = MashDataOptions(),
    mashArgs = MashOptions(),
    setSeed = 999
) {
    inputScale <- arg_match(inputScale)
    residualCorrelationMethod <- arg_match(residualCorrelationMethod)
    .mashPipelineAssert(prior, mashDataArgs, mashArgs)
    # Needed here as well as in the prior fit: the 'mle' Vhat estimator
    # refines against the supplied prior.
    priorCovariances <- prior$priorCovariances
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
        residualCorrelation,
        residualCorrelationMethod,
        priorCovariances,
        mashDataArgs
    )
    fitted <- .mashPipelinePrior(
        sumStatsList,
        alpha,
        vhat = vhat,
        prior = prior,
        inputScale = inputScale,
        mashDataArgs = mashDataArgs,
        mashArgs = mashArgs
    )
    list(U = fitted$U, w = fitted$w)
}

# Every bundle mashPipeline() takes, checked together.
# @noRd
.mashPipelineAssert <- function(prior, mashDataArgs, mashArgs) {
    .assertMethodParam(prior, "MashPriorParam", "prior")
    .assertMethodOptions(mashDataArgs, "MashDataOptions", "mashDataArgs")
    .assertMethodOptions(mashArgs, "MashOptions", "mashArgs")
}

# mashPriorCovariances() owns the cov_* chain, the supplied-prior bypass and
# the mash() weight fit; it takes the prior settings as plain arguments, so
# the bundle is unrolled here -- the one place that knows MashPriorParam()'s
# defaults for an unset field.
# @noRd
.mashPipelinePrior <- function(
    sumStatsList,
    alpha,
    vhat,
    prior,
    inputScale,
    mashDataArgs,
    mashArgs
) {
    # Normalized to a Param first: `prior` may be list() for "no settings",
    # and reading $engine off that would give NULL -- which now MEANS
    # "none" -- rather than the default engine.
    prior <- if (is(prior, "MashPriorParam")) prior else MashPriorParam()
    mashPriorCovariances(
        sumStatsList,
        alpha,
        vhat = vhat,
        components = prior$components,
        engine = prior$engine,
        priorCovariances = prior$priorCovariances,
        nPcs = prior$nPcs,
        inputScale = inputScale,
        mashDataArgs = mashDataArgs,
        mashArgs = mashArgs,
        setSeed = NULL
    )
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

# Identity ignores the null partition entirely, so a caller who assembled one
# and still landed on identity has gone to real trouble for nothing -- most
# likely they expected the null set to be used. Say so rather than quietly
# discarding it. Identity assumes conditions share no residual correlation,
# which overlapping-sample designs (the usual multi-context QTL case, where
# the same donors appear in every context) violate.
# @noRd
.mashNoteUnusedNull <- function(method, sumStatsList) {
    hasNull <- is_in("null", names(sumStatsList)) && !is.null(sumStatsList$null)
    if (method != "identity" || !hasNull) {
        return(invisible(FALSE))
    }
    inform(glue(
        "mashPipeline: `sumStatsList` carries a 'null' partition, but ",
        "residualCorrelationMethod = 'identity' does not use it -- the ",
        "residual correlation is taken to be the identity. Pass ",
        "'simple' (or 'simpleSpecific' / 'corshrink') to estimate it from ",
        "those null variants instead."
    ))
    invisible(TRUE)
}

# A prior covariance list reaches the consumers in one of two shapes, because
# the two producers disagree: mashCovarianceComponents() returns a bare Ulist
# (a named list of matrices) while mashPriorCovariances() wraps it as
# list(U, w, loglik). Accept either, so the output of either producer can be
# handed straight to any consumer. A real Ulist entry named "U" would be a
# matrix rather than a list, so the two shapes cannot be confused.
# @noRd
.mashAsUlist <- function(priorCovariances) {
    isWrapped <- is.list(priorCovariances) &&
        is_in("U", names(priorCovariances)) &&
        is.list(priorCovariances$U) &&
        !is.matrix(priorCovariances$U)
    if (isWrapped) {
        return(priorCovariances$U)
    }
    priorCovariances
}

# Vhat: the supplied residualCorrelation, else whatever
# `residualCorrelationMethod` names. The method used to be inferred from the
# partitions present ('simple' when a null set existed, else 'identity'),
# which meant a caller without a null set silently assumed zero residual
# correlation between conditions -- a statistical choice, made invisibly.
# It is now always the caller's, named up front.
# setSeed = NULL leaves the RNG stream (seeded by the caller) untouched, so
# delegated calls consume it in the original order.
# @noRd
.mashResolveVhat <- function(
    sumStatsList,
    alpha,
    inputScale,
    residualCorrelation,
    method,
    priorCovariances,
    mashDataArgs = MashDataOptions()
) {
    if (!is.null(residualCorrelation)) {
        return(residualCorrelation)
    }
    .mashNoteUnusedNull(method, sumStatsList)
    mashResidualCorrelation(
        sumStatsList,
        alpha,
        method = method,
        priorCovariances = priorCovariances,
        inputScale = inputScale,
        mashDataArgs = mashDataArgs,
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
#'   \code{"simpleSpecific"} = \code{Matrix::nearPD(cov(nullZ), corr = TRUE)};
#'   \code{"corshrink"} = \code{CorShrink::CorShrinkData()} adaptive-shrinkage
#'   correlation (needs the \pkg{CorShrink} package); \code{"mle"} =
#'   \code{mashr::mash_estimate_corr_em()} on a random subset (needs
#'   \code{"random"} + \code{priorCovariances}).
#' @param priorCovariances Prior \code{U} list (required by \code{method =
#'   "mle"}). Accepts a bare named list of covariance matrices or a
#'   \code{\link{mashPriorCovariances}} result, which is unwrapped to its
#'   \code{U}.
#' @param nSubset,maxIter \code{method = "mle"} controls (random-subset size and
#'   EM iterations).
#' @param inputScale SumStats -> matrix conversion scale (\code{"auto"} /
#'   \code{"beta"} / \code{"z"}).
#' @param mashDataArgs Extra arguments for
#'   \code{mashr::mash_set_data()}, built with
#'   \code{\link{MashDataOptions}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
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
    method = c("simple", "identity", "mle", "corshrink", "simpleSpecific"),
    priorCovariances = NULL,
    nSubset = 6000L,
    maxIter = 6L,
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = MashDataOptions(),
    setSeed = 999
) {
    # `method` takes a NAME or that estimator's own constructor, the same
    # character-or-constructor rule mashPriorCovariances(engine =) uses.
    chosen <- .resolveEngineChoice(
        if (is.character(method)) method[[1L]] else method,
        c("simple", "identity", "mle", "corshrink", "simpleSpecific"),
        "method"
    )
    method <- chosen$engine
    corArgs <- chosen$args
    .assertMethodOptions(mashDataArgs, "MashDataOptions", "mashDataArgs")
    inputScale <- arg_match(inputScale)
    if (!requireNamespace("mashr", quietly = TRUE)) {
        abort("Package 'mashr' is required for this function.")
    }
    if (methods::is(sumStatsList, "SimpleList")) {
        sumStatsList <- as.list(sumStatsList)
    }
    if (!is.null(setSeed)) {
        withr::local_seed(setSeed)
    }
    .mashResidCorDispatch(
        method,
        sumStatsList,
        alpha,
        inputScale,
        priorCovariances,
        nSubset,
        maxIter,
        mashDataArgs,
        corArgs
    )
}

# Run the chosen estimator. Split from mashResidualCorrelation() so that
# function is only the resolve-and-validate preamble: the five estimators
# share nothing but their return shape, so the choice reads as one table
# rather than as the tail of a longer function.
# @noRd
.mashResidCorDispatch <- function(
    method,
    sumStatsList,
    alpha,
    inputScale,
    priorCovariances,
    nSubset,
    maxIter,
    mashDataArgs,
    corArgs
) {
    if (method == "identity") {
        strongMats <- .mashSumStatsToMatrices(
            sumStatsList$strong,
            "strong",
            inputScale = inputScale
        )
        return(diag(rep(1, ncol(strongMats$b))))
    }
    if (method == "simple") {
        return(.mashResidCorSimple(
            sumStatsList,
            alpha,
            inputScale,
            mashDataArgs,
            corArgs
        ))
    }
    if (method == "mle") {
        return(.mashResidCorMle(
            sumStatsList,
            alpha,
            inputScale,
            priorCovariances,
            nSubset,
            maxIter,
            mashDataArgs,
            corArgs
        ))
    }
    .mashResidCorNullBased(sumStatsList, inputScale, method, corArgs)
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
#'   \code{"flashNonneg"}. Built in that fixed order.
#' @param nPcs PCs seeded into \code{cov_pca}. Default \code{ncol(Bhat) - 1}.
#' @param inputScale SumStats -> matrix conversion scale.
#' @param mashDataArgs Extra arguments for
#'   \code{mashr::mash_set_data()}, built with
#'   \code{\link{MashDataOptions}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
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
    components = c("canonical", "pca", "flash", "flashNonneg"),
    nPcs = NULL,
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = MashDataOptions(),
    setSeed = 999
) {
    .assertMethodOptions(mashDataArgs, "MashDataOptions", "mashDataArgs")
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
        inputScale,
        mashDataArgs
    )
    # This function's contract is one flat named list of covariance matrices,
    # so the role split that .mashBuildComponents keeps for the engine is
    # collapsed again here, canonical first as before.
    built <- .mashBuildComponents(
        mashData,
        components = components,
        nPcs = nPcs
    )
    c(built$canonical, built$dataDriven)
}

# `components` and `engine` each take a name or that choice's own
# constructor, so both resolve to a (name, options) pair the same way.
# @noRd
# Both arguments arrive already normalized -- a MashComponentParam and an
# Options record or NULL -- so this only splits them into the shapes the
# cov_* chain downstream expects. NULL engine is "none": no refinement.
# @noRd
.mashResolvePriorChoice <- function(components, engine) {
    entries <- slot(components, "components") %||% list()
    list(
        components = names(entries),
        componentArgs = entries,
        engine = if (is.null(engine)) "none" else metadata(engine)$engine,
        engineArgs = engine
    )
}

#' @title Estimate mash Prior Covariances and Mixture Weights
#' @description Build the mash prior covariance list \code{U} and mixture
#'   weights \code{w} from a \code{{strong[, random, null]}} SumStats
#'   collection. Split out of \code{\link{mashPipeline}} so the
#'   covariance-estimation chain lives in one place. The default builds every
#'   non-udr covariance component (\code{canonical + pca + flash +
#'   flashNonneg}) and refines them with mashr extreme deconvolution
#'   (\code{cov_ed}).
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) with at
#'   least \code{"strong"} (the discovery set the covariances are learned on).
#' @param alpha mash \code{alpha}.
#' @param vhat Residual correlation matrix (\code{V}); typically from
#'   \code{\link{mashResidualCorrelation}}. \code{NULL} -> identity.
#' @param components Data-driven covariance components, any of
#'   \code{"canonical"}, \code{"pca"}, \code{"flash"} (default
#'   \code{cov_flash}), \code{"flashNonneg"} (\code{cov_flash(factors =
#'   "nonneg")}). Built in that fixed order. May instead be a
#'   \code{\link{MashComponentParam}} record, which names the components and
#'   carries each one's options.
#' @param engine Covariance-refinement engine, either a name or the matching
#'   constructor carrying that engine's settings. \code{"covEd"} (default;
#'   \code{\link{CovEdOptions}} -- mashr's exported \code{cov_ed()} extreme
#'   deconvolution, whose default \code{algorithm = "bovy"} IS the Bovy et al.
#'   2011 method, weights from a final \code{mash()}); \code{"covUdr"}
#'   (\code{\link{CovUdrOptions}} -- \pkg{udr} ED / TED updates, returning
#'   weights directly -- OPT-IN, known numerical issues, so not the default;
#'   \code{UdFitOptions(unconstrained.update = "ted")} additionally needs
#'   i.i.d. (z-scale) data); \code{"none"} to skip refinement.
#' @param nPcs PCs seeded into \code{cov_pca}. Default \code{ncol(Bhat) - 1}.
#' @param priorCovariances Optional caller-supplied prior \code{U}: a non-empty
#'   named list of \code{nCond x nCond} matrices, or a
#'   \code{\link{mashPriorCovariances}} result, which is unwrapped to its
#'   \code{U}. When supplied, the covariance chain is bypassed and only the
#'   mixture-weight fit runs.
#' @param priorComponents Optional caller-supplied raw covariance components (a
#'   non-empty named list, e.g. the concatenated
#'   \code{\link{mashCovarianceComponents}} outputs, or a
#'   \code{\link{mashPriorCovariances}} result, which is unwrapped to its
#'   \code{U}). When supplied, these are
#'   refined by \code{engine} instead of being rebuilt internally -- the
#'   mixture-prior pipeline where separate steps built the components.
#'   \code{components} / \code{nPcs} are then ignored. Distinct from
#'   \code{priorCovariances}, which bypasses the engine entirely.
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
    components = c("canonical", "pca", "flash", "flashNonneg"),
    engine = c("covEd", "covUdr", "none"),
    nPcs = NULL,
    priorCovariances = NULL,
    priorComponents = NULL,
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = MashDataOptions(),
    mashArgs = MashOptions(),
    setSeed = 999
) {
    .assertMethodOptions(mashDataArgs, "MashDataOptions", "mashDataArgs")
    .assertMethodOptions(mashArgs, "MashOptions", "mashArgs")
    inputScale <- arg_match(inputScale)
    # Exported with the same convenience arguments as MashPriorParam(), so
    # it normalizes through the same helpers rather than resolving a second
    # time in its own way.
    choice <- .mashResolvePriorChoice(
        .mashNormalizeComponents(components),
        .mashNormalizeEngine(engine)
    )
    components <- choice$components
    componentArgs <- choice$componentArgs
    engine <- choice$engine
    engineArgs <- choice$engineArgs
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
        inputScale,
        mashDataArgs
    )
    .mashPriorFit(
        mashData,
        priorCovariances = priorCovariances,
        priorComponents = priorComponents,
        components = components,
        nPcs = nPcs,
        engine = engine,
        engineArgs = engineArgs,
        componentArgs = componentArgs,
        mashArgs = mashArgs
    )
}

# --- prior-covariance constructors ------------------------------------------
#
# mashr draws three roles that pecotmr had collapsed into one list. Canonical
# components are shape-driven -- fixed structural hypotheses built from the
# condition count. PCA and FLASH are data-driven GENERATORS: data in,
# covariances out, independent of each other. ED and udr are REFINERS, which
# consume a generator's output. mashr's own eQTL vignette refines only the
# data-driven components and passes canonical to mash() untouched.

# udr checks these names but not their VALUES, and a bad one fails deep and
# obscurely: udr::compute_penalty is two `if (update.type == ...)` branches
# with no else, so an unrecognised value assigns nothing and the function
# dies on `object 'log_penalty' not found` -- a variable the caller has
# never heard of. It is reached whenever `lambda` is non-zero, which
# .mashUdControl always sets.
#
# The vocabulary is transcribed, which this design normally refuses to do.
# The trade is three stable values against an error naming an internal
# variable: "ed" and "ted" are compute_penalty's own two branches, and
# "none" plus NA are resolved by udr::assign_prior_covariance_updates.
#
# ONLY this field. The sibling updates have their OWN vocabularies --
# scaled.update is "fa"/"none", rank1.update "ted"/"fa"/"none" -- and I have
# not traced their failure modes, so guessing at them would risk rejecting
# something udr accepts.
# @noRd
.udAssertUnconstrainedUpdate <- function(extra) {
    v <- extra[["unconstrained.update"]]
    if (is.null(v) || length(v) != 1L) {
        return(invisible(NULL))
    }
    if (is.logical(v) && is.na(v)) {
        return(invisible(NULL))
    }
    if (is.character(v) && is_in(v, c("ed", "ted", "none"))) {
        return(invisible(NULL))
    }
    abort(glue(
        "UdFitOptions: `unconstrained.update` must be 'ed', 'ted', 'none' ",
        "or NA (NA lets udr choose from the data); got ",
        "'{as.character(v)[[1L]]}'. udr accepts the name but not the ",
        "value, and fails later inside compute_penalty()."
    ))
}

# The generator constructors, keyed by the component name `components`
# accepts. flashNonneg shares cov_flash's options; pecotmr sets `factors`.
# @noRd
.mashComponentCtors <- function() {
    list(
        canonical = CovCanonicalOptions,
        pca = CovPcaOptions,
        flash = CovFlashOptions,
        flashNonneg = CovFlashOptions
    )
}

# The refiner constructors, keyed by the name `engine` accepts.
# @noRd
.mashEngineCtors <- function() {
    list(covEd = CovEdOptions, covUdr = CovUdrOptions)
}

# --- mashr fit / data / posterior settings ----------------------------------
#
# These three were the audit's top finding: the core mash fit had zero user
# control. Each names the ONE mashr function it configures, so the accepted
# set is that function's live formals.
#
# `seed` is refused by all of them: pecotmr owns run-to-run reproducibility
# through `setSeed`, and two seeds that disagree is worse than one that
# cannot be set.

# --- residual-correlation estimator settings --------------------------------
#
# `mashResidualCorrelation(method =)` already picks the estimator, so these
# follow the character-or-constructor rule the file already uses for
# `engine =`: the argument takes either the estimator's NAME or that
# estimator's constructor, which carries the identity in metadata().

#' @title Per-Component Arguments For mashPriorCovariances
#' @description Options for each prior-covariance component, keyed by
#'   component name. Each entry may be a plain list or the matching
#'   constructor -- a plain list is spliced into that constructor, so it gets
#'   the same defaults and the same checking either way.
#'
#'   Naming a component here also selects it, so \code{components} need not
#'   be given separately.
#' @param canonical Options for \code{\link{CovCanonicalOptions}}.
#' @param pca Options for \code{\link{CovPcaOptions}}.
#' @param flash,flashNonneg Options for \code{\link{CovFlashOptions}}.
#' @return A \code{MashComponentParam} object, a \code{\link{MethodParam}}.
#'   \code{names()} lists the components that were named, which is also the
#'   set that was selected.
#' @examples
#' MashComponentParam(
#'   pca = list(subset = 1:50),
#'   canonical = CovCanonicalOptions()
#' )
#' @export
MashComponentParam <- function(
    canonical = NULL,
    pca = NULL,
    flash = NULL,
    flashNonneg = NULL
) {
    # discard(is.null), not compact(): compact() also drops zero-length
    # elements, and an unconfigured constructor is a legitimately empty
    # record -- naming a component with default options must still select it.
    given <- discard(
        list(
            canonical = canonical,
            pca = pca,
            flash = flash,
            flashNonneg = flashNonneg
        ),
        is.null
    )
    # One named-list slot rather than one nullable slot per component, so a
    # record selecting pca shows that and not three NULLs. The formals stay
    # named, so R still rejects a misspelled component.
    entries <- .mashComponentEntries(given)
    new(
        "MashComponentParam",
        components = if (length(entries) == 0L) NULL else entries
    )
}

# Each entry validated and built the same way: a plain list is spliced into
# that component's constructor, an already-built record passes through, and
# one built by the WRONG constructor is refused.
#
# There is no "unknown component" check because it cannot fail any more: the
# four components are now formals, so R rejects a misspelled one as an
# unused argument before this is reached.
# @noRd
.mashComponentEntries <- function(given) {
    if (length(given) == 0L) {
        return(list())
    }
    ctors <- .mashComponentCtors()
    .nestedAssertEngines(given, ctors, NULL, "MashComponentParam")
    imap(given, .nestedElement, ctors = ctors, label = "MashComponentParam")
}

# Reject unknown prior-covariance component names.
# @noRd
# `components` may be a character vector or a MashComponentParam() record.
# Naming a component in the record selects it, so the two forms carry the
# same information and the record additionally carries per-component options.
# @noRd
.mashResolveComponentChoice <- function(components) {
    if (is(components, "MashComponentParam")) {
        entries <- slot(components, "components") %||% list()
        return(list(names = names(entries), args = entries))
    }
    if (!is.character(components)) {
        abort(glue(
            "mashPriorCovariances: `components` must be a character vector ",
            "or a MashComponentParam() record."
        ))
    }
    .mashValidateComponents(components)
    list(names = components, args = NULL)
}

.mashValidateComponents <- function(
    components,
    caller = "mashPriorCovariances"
) {
    valid <- c("canonical", "pca", "flash", "flashNonneg")
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
        abort("Package 'mashr' is required for this function.")
    }
    if (!requireNamespace("flashier", quietly = TRUE)) {
        abort("Package 'flashier' is required for this function.")
    }
}

# mashr mash_set_data over the STRONG effects (V defaults to identity).
# @noRd
.mashMakeMashData <- function(
    partition,
    label,
    vhat,
    alpha,
    inputScale,
    mashDataArgs
) {
    mats <- .mashSumStatsToMatrices(partition, label, inputScale = inputScale)
    if (is.null(vhat)) {
        vhat <- diag(rep(1, ncol(mats$b)))
    }
    .mashSetData(mats$b, mats$s, alpha, vhat, mashDataArgs)
}

# Caller-supplied prior covariance matrices (bypasses the cov_* chain; mashr
# sees only these). Validates they are a non-empty named list of nCond square
# matrices. Returns list(U, w = NULL, loglik = NULL).
# @noRd
.mashUserPriorCovariances <- function(priorCovariances, mashData) {
    priorCovariances <- .mashAsUlist(priorCovariances)
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
    engine,
    engineArgs = NULL,
    componentArgs = NULL
) {
    comps <- .mashResolveComponents(
        priorComponents,
        mashData,
        components,
        nPcs,
        componentArgs = componentArgs
    )
    # Canonical components are fixed structural hypotheses: they go to the
    # mash fit as they are. Only the data-driven set is refined -- the split
    # mashr's own eQTL vignette makes.
    refined <- if (engine == "none" || length(comps$dataDriven) == 0L) {
        list(U = comps$dataDriven, w = NULL, loglik = NULL)
    } else if (engine == "covEd") {
        .mashEngineCovEd(mashData, comps$dataDriven, engineArgs)
    } else {
        .mashEngineUd(
            mashData,
            comps$dataDriven,
            comps$canonical,
            engineArgs
        )
    }
    list(
        U = c(refined$U, comps$canonical),
        w = refined$w,
        loglik = refined$loglik
    )
}

# Raw covariance components: caller-supplied `priorComponents` (validated) or
# freshly built via .mashBuildComponents.
# @noRd
.mashResolveComponents <- function(
    priorComponents,
    mashData,
    components,
    nPcs,
    componentArgs = NULL
) {
    priorComponents <- .mashAsUlist(priorComponents)
    if (is.null(priorComponents)) {
        return(.mashBuildComponents(
            mashData,
            components = components,
            nPcs = nPcs,
            componentArgs = componentArgs
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
    # Caller-supplied components are covariances estimated from data, so they
    # are the engine's input rather than fixed hypotheses.
    list(canonical = NULL, dataDriven = priorComponents)
}

# udr control list for the ud / ud_ted engines: pecotmr's own settings,
# overlaid with whatever UdFitOptions() carried. No name translation --
# `fitOpts` is already spelled the way udr spells it, which is the point of
# splitting ud_init's arguments out of this bundle.
#
# `unconstrained.update` is left to udr where udr chooses well, and pinned
# where it does not. udr defaults it to NA and resolves that in
# assign_prior_covariance_updates() as
# `ifelse(is.matrix(fit$V), "ted", "none")`:
#
#   z scale (alpha = 1)    fit$V stays a single shared matrix, so udr picks
#                          TED -- its own and the better estimator. Omit the
#                          field and let it.
#   beta scale (alpha = 0) ud_init() expands the shared V into a per-variant
#                          one, so udr would pick "none" and leave the
#                          n_unconstrained components .mashEngineUd
#                          generated at their initialisation -- handing mash
#                          an unrefined prior with no error. Pin "ed", which
#                          does not need i.i.d. data.
#
# A caller overrides either way through UdFitOptions(), NA included.
# @noRd
.mashUdControl <- function(fitOpts, nCond, iid) {
    defaults <- list(
        scaled.update = "fa",
        resid.update = "none",
        lambda = nCond,
        penalty.type = "iw",
        maxiter = 1000L,
        tol = 1e-2,
        tol.lik = 1e-2
    )
    if (!iid) {
        defaults$unconstrained.update <- "ed"
    }
    list_modify(defaults, !!!compact(as.list(fitOpts)))
}

# One half of a covUdr selection, defaulting to that half's empty bundle:
# `engine = "covUdr"` as a bare string carries no options at all.
# @noRd
.mashUdHalf <- function(engineArgs, key, ctor) {
    if (is.null(engineArgs)) {
        return(ctor())
    }
    (engineArgs[[key]] %||% ctor())
}

# =============================================================================
# Mash model fit + posterior (mash_fit / mash_posterior notebooks)
# =============================================================================

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
#' @importFrom checkmate assertCharacter
#' @export
makePairwiseContrastCol <- function(pair, template) {
    assertCharacter(pair, len = 2L, any.missing = FALSE)
    replace(template, pair, c(1, -1))
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
#' @importFrom stringr str_remove_all fixed
#' @importFrom utils combn
#' @examples
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' fitMashContrast(1L, om, pm, pv)
#' @importFrom checkmate assertCount assertNumeric
#' @export
fitMashContrast <- function(
    index,
    origMean,
    posteriorMean,
    posteriorVcov,
    grouping = NULL
) {
    assertCount(index, positive = TRUE)
    assertNumeric(grouping, null.ok = TRUE)
    rawNames <- colnames(posteriorMean)
    populationNames <- if (is.null(rawNames)) {
        NULL
    } else {
        str_remove_all(rawNames, "BETA_")
    }
    origMeanVector <- set_names(origMean[index, ], populationNames)
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
        contrast <- replace(pairwiseVector, tested[1:2], c(1, -1))
        return(matrix(
            contrast,
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
    # Three cases, stated directly: conditions sharing a group split the
    # deviation weight between them, the diagonal carries the full weight,
    # and everything else contributes -1.
    sameAll <- outer(grouping, grouping, "==")
    # Group sizes are counted off `sameAll`, not looked up in
    # table(grouping) by as.character(grouping): that round-tripped every
    # grouping code through its string form to find its own count.
    groupSize <- as.integer(colSums(sameAll))
    sameGroup <- sameAll & matrix(grouping > 0, nPop, nPop)
    matrix(
        ifelse(
            sameGroup,
            matrix((nPop - 1) / groupSize, nPop, nPop),
            ifelse(diag(TRUE, nPop), nPop - 1, -1)
        ),
        nPop,
        nPop,
        dimnames = list(tested, str_c(tested, "_deviation"))
    )
}

# Pairwise (all-pairs) contrasts, with grouped conditions' contributions split
# evenly across the group.
# @noRd
# One grouped condition's share of a pairwise column: the group's matched
# contribution split evenly across its members.
# @noRd
.mashSplitGroupShare <- function(column, dg, grouping, groups, pwCol) {
    rowsInGroup <- names(grouping[grouping == dg])
    matchedRow <- rowsInGroup[is_in(rowsInGroup, groups)]
    if (length(matchedRow) == 0) {
        return(column)
    }
    replace(column, rowsInGroup, pwCol[matchedRow] / length(rowsInGroup))
}

# One pairwise column, with each grouped condition's contribution split
# across its group. A column whose two sides share a grouping (or that
# involves no grouped condition) is left as it is.
# @noRd
.mashAdjustPairwiseColumn <- function(col, pw, grouping) {
    column <- pw[, col]
    groups <- str_split(col, "_vs_")[[1]]
    groupValues <- grouping[is_in(names(grouping), groups)]
    relevant <- names(groupValues[groupValues > 0])
    if (n_distinct(groupValues) <= 1 || length(relevant) == 0) {
        return(column)
    }
    reduce(
        unique(groupValues[groupValues > 0]),
        .mashSplitGroupShare,
        grouping = grouping,
        groups = groups,
        pwCol = column,
        .init = column
    )
}

.mashPairwiseContrast <- function(tested, grouping, pairwiseVector) {
    twoCombn <- combn(tested, 2)
    pwNames <- apply(twoCombn, 2, str_flatten, collapse = "_vs_")
    pw <- `colnames<-`(
        apply(twoCombn, 2, makePairwiseContrastCol, pairwiseVector),
        pwNames
    )
    pwAdj <- matrix(
        unname(list_c(map(
            colnames(pw),
            .mashAdjustPairwiseColumn,
            pw = pw,
            grouping = grouping
        ))),
        nrow = nrow(pw),
        ncol = ncol(pw),
        dimnames = dimnames(pw)
    )
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
    fid <- rownames(posteriorMean)[index] %||% as.character(index)
    # unname: contrast* carry contrastDesign colnames; tibble (unlike
    # data.frame) preserves a named scalar's name on the column. Columns are
    # interleaved mean/se/p per contrast, which is the order the loop built.
    cols <- list_c(map(
        seq_along(cnames),
        .mashContrastColumns,
        cnames = cnames,
        contrastDiff = contrastDiff,
        contrastSe = contrastSe,
        contrastP = contrastP
    ))
    tibble(feature_id = fid, !!!cols)
}

# One contrast's three columns, named for it.
# @noRd
.mashContrastColumns <- function(
    i,
    cnames,
    contrastDiff,
    contrastSe,
    contrastP
) {
    set_names(
        list(
            unname(contrastDiff[i]),
            unname(contrastSe[i]),
            unname(contrastP[i])
        ),
        str_c(c("mean_contrast_", "se_contrast_", "p_contrast_"), cnames[i])
    )
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
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' mashPosteriorContrast(pm, pv, om)
#' @importFrom checkmate assertNumeric
#' @export
mashPosteriorContrast <- function(
    posteriorMean,
    posteriorVcov,
    origMean,
    grouping = NULL
) {
    assertNumeric(grouping, null.ok = TRUE)
    aligned <- origMean[, colnames(posteriorMean), drop = FALSE]
    origMean <- replace(aligned, is.nan(aligned), 0)

    parts <- compact(map(
        seq_len(nrow(posteriorMean)),
        fitMashContrast,
        origMean = origMean,
        posteriorMean = posteriorMean,
        posteriorVcov = posteriorVcov,
        grouping = grouping
    ))
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
#' model <- mashModelFit(ssl, alpha = 0, priorCovariances = prior,
#'   vhat = vhat)
#' updateMashModelCov(model, allSamples = conds, samples = conds[1:3])
#' @importFrom checkmate assertCharacter
#' @export
updateMashModelCov <- function(mashModel, allSamples, samples) {
    assertCharacter(allSamples, any.missing = FALSE)
    assertCharacter(samples, any.missing = FALSE)
    unwanted <- setdiff(allSamples, samples)
    cov <- mashModel$fitted_g$Ulist
    retained <- discard(names(cov), .mashCovIsDropped, unwanted = unwanted)
    resized <- set_names(
        map(retained, .mashResizeCov, cov = cov, samples = samples),
        retained
    )
    keptPi <- discard(
        names(mashModel$fitted_g$pi),
        .mashPiMentionsDropped,
        unwanted = unwanted
    )
    list_assign(
        mashModel,
        fitted_g = list_assign(
            mashModel$fitted_g,
            Ulist = resized,
            pi = mashModel$fitted_g$pi[keptPi]
        )
    )
}

# A covariance component belongs to a dropped condition, under either its
# bare name or its ED_ prefixed one.
# @noRd
.mashCovIsDropped <- function(d, unwanted) {
    is_in(d, unwanted) || is_in(d, str_c("ED_", unwanted))
}

# A mixture-weight name mentioning any dropped condition.
# @noRd
.mashPiMentionsDropped <- function(nm, unwanted) {
    any(map_lgl(unwanted, .mashNameMentions, nm = nm))
}

# @noRd
.mashNameMentions <- function(s, nm) {
    str_detect(nm, fixed(s))
}

# One covariance component resized to the retained conditions. A
# condition-specific component is a single 1 on its own diagonal entry;
# `identity` is a single 1 in the first cell; anything else is subset by
# name when it has one, else positionally.
# @noRd
.mashResizeCov <- function(d, cov, samples) {
    n <- length(samples)
    if (is_in(d, samples)) {
        at <- which(samples == d)
        return(replace(matrix(0, n, n), (at - 1L) * n + at, 1))
    }
    if (d == "identity") {
        return(replace(matrix(0, n, n), 1L, 1))
    }
    if (is.null(colnames(cov[[d]]))) {
        return(as.matrix(cov[[d]][seq_along(samples), seq_along(samples)]))
    }
    as.matrix(cov[[d]][samples, samples])
}

# One matrix sliced to (snps, samples), relabelled with the retained
# condition names. `snps` / `samples` are subscripts, so NULL is meaningful.
# @noRd
.mashSliceMatrix <- function(m, rows, cols) {
    `colnames<-`(as.matrix(m[rows, cols]), cols)
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
#'   dimnames = list(sprintf("chr1:%d:A:G", 100L * (1:p)), cond))
#' sbhat <- matrix(abs(rnorm(p * 3)) + 0.1, p, 3,
#'   dimnames = list(sprintf("chr1:%d:A:G", 100L * (1:p)), cond))
#' dat <- list(bhat = bhat, sbhat = sbhat, Z = bhat / sbhat,
#'   snp = sprintf("chr1:%d:A:G", 100L * (1:p)))
#' vhat <- diag(3)
#' dimnames(vhat) <- list(cond, cond)
#' sliceMashData(dat, vhat = vhat, snps = 1:4, samples = NULL)
#' @importFrom checkmate assertList assertCharacter
#' @export
sliceMashData <- function(data, vhat, snps, samples) {
    assertList(data)
    # `snps` and `samples` are SUBSCRIPTS -- `data$bhat[snps, samples]` -- so
    # character names, integer indices and NULL are all valid. No type
    # assertion is correct here (the @example passes snps = 1:4 and
    # samples = NULL).
    sliced <- list_assign(
        data,
        bhat = .mashSliceMatrix(data$bhat, snps, samples),
        sbhat = .mashSliceMatrix(data$sbhat, snps, samples),
        Z = .mashSliceMatrix(data$Z, snps, samples),
        snp = data$snp[is_in(data$snp, snps)]
    )
    list(data = sliced, vhat = .mashSliceMatrix(vhat, samples, samples))
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
#' @importFrom checkmate assertList
#' @export
sanitizeMashData <- function(data) {
    assertList(data)
    list_assign(
        data,
        bhat = replace(data$bhat, is.nan(data$bhat), 0),
        sbhat = replace(
            data$sbhat,
            is.nan(data$sbhat) | is.infinite(data$sbhat),
            1e3
        )
    )
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
#' @param metaArgs Extra arguments for \code{metafor::rma()}, built with
#'   \code{\link{RmaOptions}} -- \code{test = "knha"} in particular, the
#'   small-study correction.
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
    metaMethod = "DL",
    metaArgs = RmaOptions()
) {
    .assertMethodOptions(metaArgs, "RmaOptions", "metaArgs")
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
        metaMethod = metaMethod,
        metaArgs = metaArgs
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
    metaMethod,
    metaArgs = list()
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
        metaMethod = metaMethod,
        metaArgs = metaArgs
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
    metaMethod,
    metaArgs = list()
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
    ma <- .rmaMeta(es, se, method = metaMethod, metaArgs = metaArgs)
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
#' @param metaArgs Extra arguments for \code{metafor::rma()}, built with
#'   \code{\link{RmaOptions}} -- \code{test = "knha"} in particular, the
#'   small-study correction.
#' @return A \code{data.frame} with \code{condition} and \code{zScore}.
#' @seealso \code{\link{nSignificantScore}}, \code{\link{scoreFromCs}}
#' @examples
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' cr <- fitMashContrast(1L, om, pm, pv)
#' calculateFeatureScores(cr, metaMethod = "mean")
#' @importFrom checkmate assertString
#' @export
calculateFeatureScores <- function(
    contrastResult,
    metaMethod = "REML",
    metaArgs = RmaOptions()
) {
    assertString(metaMethod)
    .assertMethodOptions(metaArgs, "RmaOptions", "metaArgs")
    cr <- as_tibble(contrastResult)
    effCols <- names(cr)[str_detect(names(cr), "mean_contrast_.*deviation")]
    if (length(effCols) == 0L) {
        return(tibble(condition = character(0), zScore = numeric(0)))
    }
    scores <- map_dbl(
        effCols,
        .metaContrastZScore,
        cr = cr,
        metaMethod = metaMethod,
        metaArgs = metaArgs
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
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
#' dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
#' cr <- fitMashContrast(1L, om, pm, pv)
#' nSignificantScore(cr, pCutoff = 0.05)
#' @importFrom checkmate assertNumber
#' @export
nSignificantScore <- function(contrastResult, pCutoff = 1e-5) {
    assertNumber(pCutoff, lower = 0, upper = 1)
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
#' om <- matrix(c(0.1, 0.2, 0.3), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
#' pm <- matrix(c(0.5, 0.3, -0.2), 1, 3,
#'   dimnames = list("chr1:100:A:G", c("a", "b", "c")))
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
    cr <- as_tibble(contrastResults) |>
        filter(is_in(.data$feature_id, leadRows$variants))
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
    metaMethod,
    metaArgs = list()
) {
    .metaOneContrast(
        condition,
        condContrasts[i],
        abs(as.numeric(condEffects[, i])),
        as.numeric(condSes[, i]),
        seCutoff,
        metaMethod,
        metaArgs
    )
}

# Meta-analysis z-score (mean/se) for one mean-contrast column, NA when empty.
# @noRd
.metaContrastZScore <- function(ec, cr, metaMethod, metaArgs = list()) {
    seCol <- str_replace(ec, "^mean_contrast", "se_contrast")
    if (!is_in(seCol, names(cr))) {
        return(NA_real_)
    }
    esAll <- abs(as.numeric(cr[[ec]]))
    seAll <- as.numeric(cr[[seCol]])
    keep <- is.finite(esAll) & is.finite(seAll) & seAll > 0
    es <- esAll[keep]
    se <- seAll[keep]
    if (length(es) < 1L) {
        return(NA_real_)
    }
    ma <- .rmaMeta(es, se, method = metaMethod, metaArgs = metaArgs)
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


# =============================================================================
# Input preparation
# -----------------------------------------------------------------------------
# Assembling and filtering the matrices mash is fitted on. This is a
# preparatory step of the pipeline, not an engine interface -- nothing here
# calls mashr; the engine calls live in mashWrapper.R.
# =============================================================================

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
    m <- `storage.mode<-`(as.matrix(x), "double")
    replace(m, is.nan(m) | is.infinite(m) | is.na(m), replaceWith)
}

# Coerce z-scores to a matrix (NaN/Inf/NA -> 0) and, when a missing-rate
# threshold is given, drop rows falling below it.
# @noRd
.mashProcessZ <- function(zData, filterByMissingRate) {
    cleaned <- .mashReplaceValues(zData, 0)
    if (is.null(filterByMissingRate)) {
        return(cleaned)
    }
    proportionNonzero <- apply(cleaned, 1, .mashRowNonzeroRate)
    cleaned[proportionNonzero >= filterByMissingRate, , drop = FALSE]
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
#' @importFrom checkmate assertFlag assertList assertNumber
filterInvalidSummaryStat <- function(
    datList,
    bhat = NULL,
    sbhat = NULL,
    z = NULL,
    btoz = FALSE,
    sigPCutoff = 1E-6,
    filterByMissingRate = 0.2
) {
    assertList(datList)
    assertFlag(btoz)
    # NULL is how callers disable each filter.
    assertNumber(sigPCutoff, lower = 0, upper = 1, null.ok = TRUE)
    assertNumber(
        filterByMissingRate,
        lower = 0,
        upper = 1,
        null.ok = TRUE
    )
    reset <- if (
        !is.null(bhat) &&
            !is.null(sbhat) &&
            all(is_in(c(bhat, sbhat), names(datList)))
    ) {
        .mashFilterBhatSbhat(datList, bhat, sbhat, filterByMissingRate)
    } else {
        datList
    }
    withZ <- if (btoz) {
        .mashFilterBtoz(reset, bhat, sbhat, sigPCutoff)
    } else {
        reset
    }
    if (is.null(z)) {
        return(withZ)
    }
    .mashFilterZ(withZ, filterByMissingRate, sigPCutoff)
}

# Reset invalid bhat/sbhat cells (bhat -> 0, sbhat -> 1000) and, when a
# null/random partition is present, drop variants below `filterByMissingRate`
# non-missing.
# @noRd
.mashFilterBhatSbhat <- function(datList, bhat, sbhat, filterByMissingRate) {
    if (is.null(datList[[bhat]]) || is.null(datList[[sbhat]])) {
        return(datList)
    }
    reset <- list_assign(
        datList,
        !!!set_names(
            list(
                .mashReplaceValues(datList[[bhat]], 0),
                .mashReplaceValues(datList[[sbhat]], 1000)
            ),
            c(bhat, sbhat)
        )
    )
    hasNullOrRandom <- is_in("null.b", names(reset)) ||
        is_in("random.b", names(reset))
    if (!hasNullOrRandom || is.null(filterByMissingRate)) {
        return(reset)
    }
    proportionNonzero <- apply(reset[[bhat]], 1, .mashRowNonzeroRate)
    keep <- proportionNonzero >= filterByMissingRate
    list_assign(
        reset,
        !!!set_names(
            # `drop = FALSE`: one surviving variant would collapse these
            # to vectors, and they go back into the mash data list where
            # mashr expects an N x R matrix.
            list(
                reset[[bhat]][keep, , drop = FALSE],
                reset[[sbhat]][keep, , drop = FALSE]
            ),
            c(bhat, sbhat)
        )
    )
}

# Derive z = bhat / sbhat (into a `<condition>.z` or `z` slot) and apply the
# significance cutoff to strong signals.
# @noRd
.mashFilterBtoz <- function(datList, bhat, sbhat, sigPCutoff) {
    perCondition <- any(str_detect(bhat, "\\.b$")) ||
        any(str_detect(sbhat, "\\.s$"))
    zName <- if (perCondition) {
        str_c(str_remove(bhat, "\\.b$"), ".z")
    } else {
        "z"
    }
    # list(NULL) not NULL: the z slot must EXIST and be empty, where assigning
    # NULL would delete it.
    zValue <- if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
        list(as.matrix(datList[[bhat]] / datList[[sbhat]]))
    } else {
        list(NULL)
    }
    withZ <- list_assign(datList, !!!set_names(zValue, zName))
    if (!is_in("strong.z", names(withZ)) || is.null(sigPCutoff)) {
        return(withZ)
    }
    keepIndex <- filterBySignificance(withZ$strong.z, sigPCutoff)
    list_assign(
        withZ,
        # `drop = FALSE`: exactly one significant variant is a common
        # outcome, and a vector here would reach mashr as the strong set.
        strong.z = withZ$strong.z[keepIndex, , drop = FALSE],
        strong.b = withZ$strong.b[keepIndex, , drop = FALSE],
        strong.s = withZ$strong.s[keepIndex, , drop = FALSE]
    )
}

# Process each partition's z-matrix (missing-rate filter) and apply the
# significance cutoff to strong z-scores.
# @noRd
.mashFilterZ <- function(datList, filterByMissingRate, sigPCutoff) {
    # Only partitions that are present and carry a z get rewritten:
    # `list_assign()` would otherwise CREATE an absent component as NULL.
    components <- keep(
        intersect(c("strong", "random", "null"), names(datList)),
        .mashPartitionHasZ,
        datList = datList
    )
    processed <- if (length(components) == 0L) {
        datList
    } else {
        list_assign(
            datList,
            !!!set_names(
                map(
                    components,
                    .mashProcessPartition,
                    datList = datList,
                    filterByMissingRate = filterByMissingRate
                ),
                components
            )
        )
    }
    if (
        is.null(processed$strong) ||
            is.null(processed$strong$z) ||
            is.null(sigPCutoff)
    ) {
        return(processed)
    }
    keepIndex <- filterBySignificance(processed$strong$z, sigPCutoff)
    list_assign(
        processed,
        strong = list_assign(
            processed$strong,
            z = processed$strong$z[keepIndex, , drop = FALSE]
        )
    )
}

# @noRd
.mashPartitionHasZ <- function(comp, datList) {
    !is.null(datList[[comp]]) && !is.null(datList[[comp]]$z)
}

# One partition with its z-matrix missing-rate filtered.
# @noRd
.mashProcessPartition <- function(comp, datList, filterByMissingRate) {
    part <- datList[[comp]]
    list_assign(part, z = .mashProcessZ(part$z, filterByMissingRate))
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
#' @importFrom checkmate assertCharacter assertNumber
filterMixtureComponents <- function(
    conditionsToKeep,
    U,
    w = NULL,
    wCutoff = 1e-04
) {
    assertCharacter(conditionsToKeep, any.missing = FALSE)
    assertNumber(wCutoff, lower = 0, finite = TRUE)
    conditionsToFilter <- setdiff(colnames(U[[1]]), conditionsToKeep)
    sumW <- sum(w)
    subsetU <- .mashSubsetU(U, conditionsToKeep)
    # Drop all-zero matrices, then those below the weight cutoff.
    nonzero <- names(keep(subsetU, .mashMatrixNonzero))
    keepNames <- if (is.null(w)) {
        nonzero
    } else {
        intersect(nonzero, names(w[w >= wCutoff]))
    }
    # Also drop the U components driven by non-relevant contexts: the EM can
    # leave tiny non-zero diagonals, so all-zero removal alone won't drop them,
    # yet real diagonal signal must be kept.
    keptU <- subsetU[setdiff(keepNames, conditionsToFilter)]
    keptW <- w[keepNames]
    survivors <- keptW[!is_in(names(keptW), conditionsToFilter)]
    # Rescale the surviving weights back to the original total.
    rescaled <- (survivors / sum(survivors)) * sumW
    msg <- glue(
        "{length(keptU)} components of matrices remained after filtering."
    )
    inform(msg)
    list(U = keptU, w = rescaled)
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
        keys <- intersect(names(dat), c("z", "bhat", "sbhat"))
        dat <- list_assign(
            dat,
            !!!set_names(
                map(dat[keys], .mashDropConditions, drop = excludeCondition),
                keys
            )
        )
    }
    .mashExtractOneData(dat, nRandom, nNull)
}

# One matrix without the excluded condition columns.
# @noRd
.mashDropConditions <- function(m, drop) {
    m[, setdiff(colnames(m), drop), drop = FALSE]
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
#' @importFrom checkmate assertList
#' @export
mergeMashData <- function(resData, oneData) {
    assertList(resData, null.ok = TRUE)
    assertList(oneData, null.ok = TRUE)
    if (length(resData) == 0 || is.null(resData)) {
        return(oneData)
    }
    if (length(oneData) == 0 || is.null(oneData)) {
        return(resData)
    }

    set_names(
        map(
            names(oneData),
            .mashCombineDatum,
            oneData = oneData,
            resData = resData
        ),
        names(oneData)
    )
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
    rawMe <- getMarginalEffects(obj)
    rawCs <- getCs(obj, coverage = coverage)
    if (!all(is_in(c("variant_id", "context", "beta", "se"), names(rawMe)))) {
        msg <- glue(
            "mashInput: getMarginalEffects() must return variant_id/context/",
            "beta/se columns; a FineMappingResult with >= 2 contexts is ",
            "required."
        )
        abort(msg)
    }
    pinned <- .mashFmrMethodPin(rawMe, rawCs)
    me <- pinned$me
    cs <- pinned$cs
    contexts <- unique(me$context)
    variants <- unique(me$variant_id)
    empty <- matrix(
        NA_real_,
        length(variants),
        length(contexts),
        dimnames = list(variants, contexts)
    )
    # One (variant, context) cell per long-format row; the rest stay NA.
    cell <- cbind(match(me$variant_id, variants), match(me$context, contexts))
    list(
        b = replace(empty, cell, me$beta),
        s = replace(empty, cell, me$se),
        strongRows = .mashFmrStrongRows(cs, variants)
    )
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

# The max-PIP variant among one credible set's rows.
# @noRd
.mashLeadVariant <- function(rows, cs) {
    cs$variant_id[[rows[[which.max(cs$pip[rows])]]]]
}

# Row indices (into `variants`) of each credible set's lead (max PIP) variant.
# @noRd
.mashFmrStrongRows <- function(cs, variants) {
    csCol <- names(cs)[str_detect(names(cs), "^cs_")]
    strongVar <- if (
        nrow(cs) > 0L && length(csCol) > 0L && is_in("pip", names(cs))
    ) {
        grp <- interaction(cs$context, cs[[csCol[[1L]]]], drop = TRUE)
        map_chr(split(seq_len(nrow(cs)), grp), .mashLeadVariant, cs = cs)
    } else {
        character(0)
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
    hasStrong <- length(mats$strongRows) > 0L
    hasRandom <- !is.null(rn$random) && length(rn$random) > 0L
    hasNull <- !is.null(rn$null) && length(rn$null) > 0L
    compact(list(
        strong.b = if (hasStrong) {
            mats$b[mats$strongRows, keepCols, drop = FALSE]
        },
        strong.s = if (hasStrong) {
            mats$s[mats$strongRows, keepCols, drop = FALSE]
        },
        random.b = if (hasRandom) rn$random$bhat,
        random.s = if (hasRandom) rn$random$sbhat,
        null.b = if (hasNull) rn$null$bhat,
        null.s = if (hasNull) rn$null$sbhat
    ))
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
        return(set_names(objects, str_c("region", seq_along(objects))))
    }
    objects
}

# Extract + merge each object's strong/random/null partitions, disambiguating
# rownames by region before accumulating.
# @noRd
.mashCombinePartitions <- function(objects, cfg) {
    reduce(
        map(
            names(objects),
            .mashRegionPartitions,
            objects = objects,
            cfg = cfg
        ),
        mergeMashData,
        .init = list()
    )
}

# One region's partitions, with its rownames region-prefixed so the merge can
# tell same-named variants from different regions apart.
# @noRd
.mashRegionPartitions <- function(nm, objects, cfg) {
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
    map(part, .mashPrefixRownames, nm = nm)
}

# Coerce to data.frame, clean each partition + derive z (random/null before
# strong so the strong significance filter runs once), restore the strong
# 1-row matrix shape, add the strong XtX, and optionally drop b/s slots.
# @noRd
.mashFinalizeCombined <- function(combined, sigPCutoff, zOnly) {
    # Each condition's z derivation sees the frame the previous one produced,
    # so the sweep is a fold rather than a variable rewritten three times.
    withZ <- reduce(
        c("random", "null", "strong"),
        .mashDeriveZFor,
        sigPCutoff = sigPCutoff,
        .init = map(combined, .mashAsDataFrameOrNull)
    )
    shaped <- .mashAddXtX(.mashRestoreStrongShape(withZ))
    if (!zOnly) {
        return(shaped)
    }
    shaped[!str_detect(names(shaped), "\\.(b|s)$")]
}

# Derive z for one condition, when it carries both b and s.
# @noRd
.mashDeriveZFor <- function(combined, cond, sigPCutoff) {
    bKey <- str_c(cond, ".b")
    sKey <- str_c(cond, ".s")
    if (is.null(combined[[bKey]]) || is.null(combined[[sKey]])) {
        return(combined)
    }
    filterInvalidSummaryStat(
        combined,
        bhat = bKey,
        sbhat = sKey,
        btoz = TRUE,
        sigPCutoff = sigPCutoff
    )
}

# filterInvalidSummaryStat subsets strong without drop = FALSE, so a single
# surviving strong variant degrades to a vector; restore the 1-row matrix.
# @noRd
.mashRestoreStrongShape <- function(combined) {
    # Only keys that exist AND lost their dim are rewritten; `list_assign()`
    # would otherwise create an absent key as NULL.
    needs <- keep(
        intersect(c("strong.b", "strong.s", "strong.z"), names(combined)),
        .mashLostDim,
        combined = combined
    )
    if (length(needs) == 0L) {
        return(combined)
    }
    list_assign(
        combined,
        !!!set_names(map(combined[needs], .mashAsOneRow), needs)
    )
}

# @noRd
.mashLostDim <- function(k, combined) {
    !is.null(combined[[k]]) && is.null(dim(combined[[k]]))
}

# @noRd
.mashAsOneRow <- function(v) {
    matrix(v, nrow = 1L, dimnames = list(NULL, names(v)))
}

# Strong XtX cross-product (conditions x conditions), when strong.z is present.
# @noRd
.mashAddXtX <- function(combined) {
    if (
        is.null(combined$strong.z) || nrow(as.matrix(combined$strong.z)) == 0L
    ) {
        return(combined)
    }
    sz <- as.matrix(combined$strong.z)
    list_assign(combined, XtX = crossprod(sz) / nrow(sz))
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
    vids <- rownames(z) %||% str_c("var", seq_len(nrow(z)))
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
        mcolArgs = list(a1 = a1, a2 = a2, z = z, n = n)
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
    vids <- rownames(bhat) %||% str_c("var", seq_len(nrow(bhat)))
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
        mcolArgs = list(a1 = a1, a2 = a2, bhat = bhat, shat = shat, n = n)
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

# Decode chrom/pos from the variant ids, synthesising where they do not
# parse. An unparseable id still needs a placeable coordinate: the GRanges
# is keyed by the variant id, so the range only has to be unique and
# ordered, not correct.
# @noRd
#' @importFrom rlang try_fetch
.qszmCoords <- function(vids) {
    parsed <- try_fetch(
        suppressWarnings(parseVariantId(vids)),
        error = function(cnd) NULL
    )
    rawChrom <- if (!is.null(parsed)) {
        as.character(parsed$chrom)
    } else {
        rep(NA_character_, length(vids))
    }
    rawPos <- if (!is.null(parsed)) {
        suppressWarnings(as.integer(parsed$pos))
    } else {
        rep(NA_integer_, length(vids))
    }
    list(
        chrom = replace(
            rawChrom,
            is.na(rawChrom) | str_length(rawChrom) == 0L,
            "chr1"
        ),
        pos = replace(rawPos, is.na(rawPos), seq_along(rawPos)[is.na(rawPos)])
    )
}

# Internal: shared assembly for the z / beta matrix constructors. Recycles the
# context / trait labels, decodes chrom/pos from the variant ids (synthesising
# where they don't parse), builds one GRanges entry per condition with mcols
# from `mcolFn(j)`, and wraps the entries as a QtlSumStats.
# @noRd
#' @importFrom rlang try_fetch
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
    mcolArgs = list()
) {
    context <- .qszmRecycle(context, nCond, "context")
    trait <- .qszmRecycle(trait, nCond, "trait")
    coords <- .qszmCoords(vids)
    chrom <- coords$chrom
    pos <- coords$pos
    entries <- map(
        seq_len(nCond),
        .qszmEntry,
        chrom = chrom,
        pos = pos,
        vids = vids,
        mcolFn = mcolFn,
        mcolArgs = mcolArgs
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
    bhat <- exec(rbind, !!!blocks$bhat)
    shat <- exec(rbind, !!!blocks$shat)
    # bhat NA -> 0, shat NA / <= 0 -> 1000 (the mash_set_data
    # zero_Bhat_Shat_reset convention; missing-cell variants do not drive the
    # fit).
    list(
        b = replace(bhat, is.na(bhat), 0),
        s = replace(shat, is.na(shat) | shat <= 0, 1000)
    )
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
    rows <- map(
        rowsInBlock,
        .mashContextRow,
        x = x,
        setup = setup,
        requireCols = requireCols,
        resolvedScale = resolvedScale
    )
    contexts <- map_chr(rows, "context")
    # A later row overwrites an earlier one sharing a context, as the keyed
    # assignment did; the variant order is first-seen across all rows.
    lastPerContext <- !duplicated(contexts, fromLast = TRUE)
    list(
        variantOrder = unique(.mashConcatChr(map(rows, "snps"))),
        perContextB = set_names(
            map(rows[lastPerContext], "b"),
            contexts[lastPerContext]
        ),
        perContextSe = set_names(
            map(rows[lastPerContext], "se"),
            contexts[lastPerContext]
        )
    )
}

# @noRd
.mashConcatChr <- function(pieces) {
    if (length(pieces) == 0L) {
        return(character(0))
    }
    as.character(list_c(pieces))
}

# One row's effect / standard-error vectors for its context. On the z scale
# the standard errors are unit by construction.
# @noRd
.mashContextRow <- function(rIdx, x, setup, requireCols, resolvedScale) {
    df <- .mashRowDf(x, rIdx, setup, requireCols)
    snps <- df$variant_id
    onBeta <- resolvedScale == "beta"
    list(
        context = setup$contextCol[[rIdx]],
        snps = snps,
        b = set_names(if (onBeta) df$beta else df$z, snps),
        se = set_names(
            if (onBeta) df$se else rep(1, length(snps)),
            snps
        )
    )
}

# One row's sumstat data.frame (QtlSumStats keyed by study/context/trait;
# GwasSumStats by study).
# @noRd
.mashRowDf <- function(x, rIdx, setup, requireCols) {
    if (setup$isQtl) {
        getSumStatsDf(
            x,
            study = setup$studyCol[[rIdx]],
            context = setup$contextCol[[rIdx]],
            trait = setup$traitCol[[rIdx]],
            require = requireCols
        )
    } else {
        getSumStatsDf(x, study = setup$studyCol[[rIdx]], require = requireCols)
    }
}

# Assemble one block's (variant x context) Bhat / Shat matrices, disambiguating
# rownames by block key to avoid silent cross-block dedup.
# @noRd
.mashBlockMatrix <- function(bkey, x, setup, resolvedScale) {
    rowsInBlock <- which(setup$blockKeys == bkey)
    pc <- .mashBlockPerContext(x, rowsInBlock, setup, resolvedScale)
    # Every context owns one column, so each is built whole and the columns
    # are laid side by side -- no scatter into a preallocated matrix, and the
    # block-qualified rownames go on at construction.
    dims <- list(
        str_c(bkey, pc$variantOrder, sep = "::"),
        setup$columnLabels
    )
    list(
        b = .mashContextMatrix(pc$perContextB, pc$variantOrder, dims),
        s = .mashContextMatrix(pc$perContextSe, pc$variantOrder, dims)
    )
}

# One context's column, aligned to `variantOrder`. Indexing a named vector by
# a variant it lacks yields NA, which is the unfilled cell.
# @noRd
.mashContextColumn <- function(ctx, perContext, variantOrder) {
    v <- perContext[[ctx]]
    if (is.null(v)) {
        return(rep(NA_real_, length(variantOrder)))
    }
    unname(v[variantOrder])
}

# @noRd
.mashContextMatrix <- function(perContext, variantOrder, dims) {
    cols <- map(
        dims[[2L]],
        .mashContextColumn,
        perContext = perContext,
        variantOrder = variantOrder
    )
    matrix(
        unname(list_c(cols)),
        nrow = length(variantOrder),
        ncol = length(dims[[2L]]),
        dimnames = dims
    )
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
    joined <- bind_rows(as.data.frame(rd), as.data.frame(od))
    # NaN, not NA: mash reads a missing cell as NaN.
    combined <- replace(joined, is.na(joined), NaN)
    `rownames<-`(combined, c(rnRes, rnOne))
}

# Region-prefix one partition matrix's rownames (no-op for empty/NULL).
# @noRd
.mashPrefixRownames <- function(m, nm) {
    if (is.null(m) || nrow(m) == 0L) {
        return(m)
    }
    `rownames<-`(m, str_c(rownames(m), nm, sep = "_"))
}

# Coerce one partition to a data.frame (NULL passes through). Kept as a base
# data.frame (not a tibble) so the variant-id rownames survive to the output
# matrices -- see .mashCombineDatum.
# @noRd
.mashAsDataFrameOrNull <- function(m) {
    if (is.null(m)) {
        return(NULL)
    }
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
    S4Vectors::`mcols<-`(gr, value = exec(mcolFn, !!!mcolCallArgs))
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
