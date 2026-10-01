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
#'   \code{"flashNonneg"}, or a \code{\link{mashComponentConfig}} record to
#'   configure them. Ignored when \code{priorCovariances} is supplied.
#' @param engine How the data-driven components are refined:
#'   \code{"covEd"} (default), \code{"covUdr"} or \code{"none"}, or the
#'   matching constructor --- \code{\link{covEdConfig}} /
#'   \code{\link{covUdrConfig}} --- to configure it at the same time.
#'   Ignored when \code{priorCovariances} is supplied.
#' @param nPcs Optional integer; principal components seeded into
#'   \code{mashr::cov_pca()}. Defaults to \code{ncol} of the data. Read
#'   only when \code{components} includes \code{"pca"}.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashPriorConfig(components = c("canonical", "pca"), nPcs = 3)
#' @export
mashPriorConfig <- function(
    priorCovariances = NULL,
    components = c("canonical", "pca", "flash", "flashNonneg"),
    engine = c("covEd", "covUdr", "none"),
    nPcs = NULL
) {
    # components / engine are NOT arg_match()ed here: both take either a name
    # or that engine's own constructor, and mashPriorCovariances() resolves
    # the choice. Collapsing them to a string would drop the options.
    .newMethodConfig(
        NULL,
        defaults = list(
            priorCovariances = priorCovariances,
            components = components,
            engine = engine,
            nPcs = nPcs
        ),
        extra = list(),
        label = "mashPriorConfig",
        engine = "mashPrior"
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
#'   \code{\link{mashPriorConfig}}: \code{priorCovariances} supplies them
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
#'   \code{\link{mashDataConfig}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
#' @param mashArgs Extra arguments for \code{mashr::mash()}, built with
#'   \code{\link{mashConfig}} -- for example \code{nullweight},
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
#' mashPipeline(sumStatsList, alpha = 0, prior = mashPriorConfig(nPcs = 2L))
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
    prior = mashPriorConfig(),
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = mashDataConfig(),
    mashArgs = mashConfig(),
    setSeed = 999
) {
    inputScale <- arg_match(inputScale)
    residualCorrelationMethod <- arg_match(residualCorrelationMethod)
    .assertMethodConfig(prior, "mashPriorConfig", "prior")
    .assertMethodConfig(mashDataArgs, "mashDataConfig", "mashDataArgs")
    .assertMethodConfig(mashArgs, "mashConfig", "mashArgs")
    # The bundle's terminal: .mashResolveVhat and mashPriorCovariances below
    # take the four as plain arguments.
    priorCovariances <- prior$priorCovariances
    components <- prior$components %||%
        c("canonical", "pca", "flash", "flashNonneg")
    engine <- prior$engine %||% c("covEd", "covUdr", "none")
    nPcs <- prior$nPcs
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
    # mashPriorCovariances() owns the cov_* chain, the supplied-prior bypass,
    # and the mash() weight fit; mashPipeline just forwards its arguments.
    prior <- mashPriorCovariances(
        sumStatsList,
        alpha,
        vhat = vhat,
        components = components,
        engine = engine,
        priorCovariances = priorCovariances,
        nPcs = nPcs,
        inputScale = inputScale,
        mashDataArgs = mashDataArgs,
        mashArgs = mashArgs,
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
    mashDataArgs = mashDataConfig()
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
#'   \code{\link{mashDataConfig}} -- for example
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
    mashDataArgs = mashDataConfig(),
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
    .assertMethodConfig(mashDataArgs, "mashDataConfig", "mashDataArgs")
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
#'   \code{\link{mashDataConfig}} -- for example
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
    mashDataArgs = mashDataConfig(),
    setSeed = 999
) {
    .assertMethodConfig(mashDataArgs, "mashDataConfig", "mashDataArgs")
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
#'   \code{\link{mashComponentConfig}} record, which names the components and
#'   carries each one's options.
#' @param engine Covariance-refinement engine, either a name or the matching
#'   constructor carrying that engine's settings. \code{"covEd"} (default;
#'   \code{\link{covEdConfig}} -- mashr's exported \code{cov_ed()} extreme
#'   deconvolution, whose default \code{algorithm = "bovy"} IS the Bovy et al.
#'   2011 method, weights from a final \code{mash()}); \code{"covUdr"}
#'   (\code{\link{covUdrConfig}} -- \pkg{udr} ED / TED updates, returning
#'   weights directly -- OPT-IN, known numerical issues, so not the default;
#'   \code{unconstrainedUpdate = "ted"} additionally needs i.i.d. (z-scale)
#'   data); \code{"none"} to skip refinement.
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
#'   \code{\link{mashDataConfig}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
#' @param mashArgs Extra arguments for \code{mashr::mash()}, built with
#'   \code{\link{mashConfig}} -- for example \code{nullweight},
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
    mashDataArgs = mashDataConfig(),
    mashArgs = mashConfig(),
    setSeed = 999
) {
    .assertMethodConfig(mashDataArgs, "mashDataConfig", "mashDataArgs")
    .assertMethodConfig(mashArgs, "mashConfig", "mashArgs")
    inputScale <- arg_match(inputScale)
    resolvedComponents <- .mashResolveComponentChoice(components)
    components <- resolvedComponents$names
    componentArgs <- resolvedComponents$args
    chosenEngine <- .resolveEngineChoice(
        if (is.character(engine)) engine[[1L]] else engine,
        c(names(.mashEngineCtors()), "none"),
        "engine"
    )
    engine <- chosenEngine$engine
    engineArgs <- chosenEngine$args
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

# --- prior-covariance constructors ------------------------------------------
#
# mashr draws three roles that pecotmr had collapsed into one list. Canonical
# components are shape-driven -- fixed structural hypotheses built from the
# condition count. PCA and FLASH are data-driven GENERATORS: data in,
# covariances out, independent of each other. ED and udr are REFINERS, which
# consume a generator's output. mashr's own eQTL vignette refines only the
# data-driven components and passes canonical to mash() untouched.

#' @title Arguments For mashr's Canonical Covariance Components
#' @description Options for \code{mashr::cov_canonical}, which builds the
#'   fixed structural hypotheses (identity, singletons, equal effects, simple
#'   heterogeneity) from the condition count. These are hypotheses about
#'   shape, so they go to the mash fit unrefined.
#' @param cov_methods Which canonical components to build. \code{NULL}
#'   (default) leaves mashr's own selection in place.
#' @param ... Any other \code{mashr::cov_canonical} argument.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' covCanonicalConfig(cov_methods = c("identity", "equal_effects"))
#' @export
covCanonicalConfig <- function(cov_methods = NULL, ...) {
    .newMethodConfig(
        "mashr::cov_canonical",
        defaults = list(cov_methods = cov_methods),
        extra = list(...),
        label = "covCanonicalConfig",
        engine = "canonical"
    )
}

#' @title Arguments For mashr's PCA Covariance Components
#' @description Options for \code{mashr::cov_pca}, a data-driven generator.
#'   The number of components is not settable here: it comes from
#'   \code{mashPriorCovariances(nPcs = )}, which the caller also reports on.
#' @param subset Rows of the data to use. \code{NULL} (default) uses all.
#' @param ... Any other \code{mashr::cov_pca} argument.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' covPcaConfig(subset = 1:50)
#' @export
covPcaConfig <- function(subset = NULL, ...) {
    .newMethodConfig(
        "mashr::cov_pca",
        defaults = list(subset = subset),
        extra = list(...),
        label = "covPcaConfig",
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
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' covFlashConfig(remove_singleton = TRUE)
#' @export
covFlashConfig <- function(
    subset = NULL,
    remove_singleton = FALSE,
    tag = NULL,
    output_model = NULL,
    greedy_args = list(),
    backfit_args = list(),
    ...
) {
    .newMethodConfig(
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
        label = "covFlashConfig",
        engine = "flash"
    )
}

#' @title Arguments For Extreme Deconvolution Refinement
#' @description Options for \code{mashr::cov_ed}, which refines the
#'   data-driven components by Extreme Deconvolution. The components to refine
#'   are not settable here -- they are whatever \code{components} produced.
#'
#'   \code{mashr::cov_ed} takes \code{...}, so argument names \strong{cannot
#'   be checked} against it; see \code{\link{MethodConfig}}.
#' @param subset Rows of the data to use. \code{NULL} (default) uses all.
#' @param algorithm Deconvolution algorithm, \code{"bovy"} (default) or
#'   \code{"teem"}.
#' @param ... Any other argument \code{mashr::cov_ed} forwards.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' covEdConfig(algorithm = "teem")
#' @export
covEdConfig <- function(subset = NULL, algorithm = "bovy", ...) {
    .newMethodConfig(
        "mashr::cov_ed",
        defaults = list(subset = subset, algorithm = algorithm),
        extra = list(...),
        label = "covEdConfig",
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

#' @title Arguments For Unconstrained Deconvolution Refinement
#' @description Options for refining the data-driven components with
#'   \code{udr}. Names are checked against
#'   \code{udr::ud_fit_control_default()}, which is where udr's tunables
#'   live; \code{udr::ud_fit} itself takes an opaque \code{control} list.
#' @param unconstrainedUpdate How the unconstrained covariances are updated:
#'   \code{"ed"} (default) or \code{"ted"}. This was previously spelled as
#'   two separate engines, \code{"ud"} and \code{"ud_ted"}.
#' @param nUnconstrained Number of unconstrained covariances udr initialises
#'   when no data-driven components are supplied. Default \code{50}.
#' @param ... Any \code{udr} control field (\code{maxiter}, \code{tol},
#'   \code{tol.lik}, \code{penalty.type}, ...).
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' covUdrConfig(unconstrainedUpdate = "ted", maxiter = 500)
#' @export
covUdrConfig <- function(
    unconstrainedUpdate = c("ed", "ted"),
    nUnconstrained = 50L,
    ...
) {
    unconstrainedUpdate <- arg_match(unconstrainedUpdate)
    .newMethodConfig(
        NULL,
        defaults = list(
            unconstrainedUpdate = unconstrainedUpdate,
            nUnconstrained = nUnconstrained
        ),
        extra = list(...),
        label = "covUdrConfig",
        engine = "covUdr",
        accepted = c(
            "unconstrainedUpdate",
            "nUnconstrained",
            .udrControlNames()
        )
    )
}

# The generator constructors, keyed by the component name `components`
# accepts. flashNonneg shares cov_flash's options; pecotmr sets `factors`.
# @noRd
.mashComponentCtors <- function() {
    list(
        canonical = covCanonicalConfig,
        pca = covPcaConfig,
        flash = covFlashConfig,
        flashNonneg = covFlashConfig
    )
}

# The refiner constructors, keyed by the name `engine` accepts.
# @noRd
.mashEngineCtors <- function() {
    list(covEd = covEdConfig, covUdr = covUdrConfig)
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
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashConfig(nullweight = 10, optmethod = "mixSQP")
#' @export
mashConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "the sumstats pecotmr assembles",
            Ulist = "the prior covariances",
            outputlevel = "fixed by the calling entry point",
            seed = "the caller's `setSeed`"
        ),
        "mashConfig"
    )
    .newMethodConfig(
        "mashr::mash",
        defaults = list(),
        extra = extra,
        label = "mashConfig",
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
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashDataConfig(zero_Bhat_Shat_reset = 500)
#' @export
mashDataConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            Bhat = "the effect-size matrix pecotmr builds",
            Shat = "the standard-error matrix pecotmr builds",
            alpha = "the caller's `alpha`",
            V = "the caller's `vhat`"
        ),
        "mashDataConfig"
    )
    .newMethodConfig(
        "mashr::mash_set_data",
        defaults = list(zero_Bhat_Shat_reset = 1000),
        extra = extra,
        label = "mashDataConfig",
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
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashPosteriorConfig(pi_thresh = 1e-8)
#' @export
mashPosteriorConfig <- function(...) {
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
        "mashPosteriorConfig"
    )
    .newMethodConfig(
        "mashr::mash_compute_posterior_matrices",
        defaults = list(),
        extra = extra,
        label = "mashPosteriorConfig",
        engine = "mashPosterior"
    )
}

# --- residual-correlation estimator settings --------------------------------
#
# `mashResidualCorrelation(method =)` already picks the estimator, so these
# follow the character-or-constructor rule the file already uses for
# `engine =`: the argument takes either the estimator's NAME or that
# estimator's constructor, which carries the identity in metadata().

#' @title Settings For The Simple Null-Correlation Estimator
#' @description Options forwarded to
#'   \code{mashr::estimate_null_correlation_simple}. \code{data} is refused:
#'   it is the null-partition data object pecotmr assembles.
#' @param ... Arguments for the estimator, under its own names.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashCorSimpleConfig(z_thresh = 3)
#' @export
mashCorSimpleConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(data = "the null-partition data pecotmr assembles"),
        "mashCorSimpleConfig"
    )
    .newMethodConfig(
        "mashr::estimate_null_correlation_simple",
        defaults = list(),
        extra = extra,
        label = "mashCorSimpleConfig",
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
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashCorEmConfig(tol = 1e-5)
#' @export
mashCorEmConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            data = "the random-subset data pecotmr assembles",
            Ulist = "the prior covariances",
            max_iter = "the caller's `maxIter`"
        ),
        "mashCorEmConfig"
    )
    .newMethodConfig(
        "mashr::mash_estimate_corr_em",
        defaults = list(),
        extra = extra,
        label = "mashCorEmConfig",
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
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' corShrinkConfig(nboot = 100)
#' @export
corShrinkConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(data = "the null z-matrix pecotmr assembles"),
        "corShrinkConfig"
    )
    .newMethodConfig(
        "CorShrink::CorShrinkData",
        defaults = list(
            ash.control = list(mixcompdist = "halfuniform"),
            image = "null"
        ),
        extra = extra,
        label = "corShrinkConfig",
        engine = "corshrink"
    )
}

#' @title Per-Component Arguments For mashPriorCovariances
#' @description Options for each prior-covariance component, keyed by
#'   component name. Each entry may be a plain list or the matching
#'   constructor -- a plain list is spliced into that constructor, so it gets
#'   the same defaults and the same checking either way.
#'
#'   Naming a component here also selects it, so \code{components} need not
#'   be given separately.
#' @param canonical Options for \code{\link{covCanonicalConfig}}.
#' @param pca Options for \code{\link{covPcaConfig}}.
#' @param flash,flashNonneg Options for \code{\link{covFlashConfig}}.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mashComponentConfig(
#'   pca = list(subset = 1:50),
#'   canonical = covCanonicalConfig()
#' )
#' @export
mashComponentConfig <- function(
    canonical = NULL,
    pca = NULL,
    flash = NULL,
    flashNonneg = NULL
) {
    # discard(is.null), not compact(): compact() also drops zero-length
    # elements, and an unconfigured constructor is a legitimately empty
    # record -- naming a component with default options must still select it.
    .newNestedConfig(
        discard(
            list(
                canonical = canonical,
                pca = pca,
                flash = flash,
                flashNonneg = flashNonneg
            ),
            is.null
        ),
        .mashComponentCtors(),
        "mashComponentConfig"
    )
}

# Reject unknown prior-covariance component names.
# @noRd
# `components` may be a character vector or a mashComponentConfig() record.
# Naming a component in the record selects it, so the two forms carry the
# same information and the record additionally carries per-component options.
# @noRd
.mashResolveComponentChoice <- function(components) {
    if (.isMethodConfig(components)) {
        return(list(names = names(components), args = components))
    }
    if (!is.character(components)) {
        abort(glue(
            "mashPriorCovariances: `components` must be a character vector ",
            "or a mashComponentConfig() record."
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

# udr control list for the ud / ud_ted engines.
# @noRd
# pecotmr's udr control settings, overlaid with whatever covUdrConfig() carried.
# `unconstrainedUpdate` and `nUnconstrained` are pecotmr's spellings for
# udr's `unconstrained.update` and ud_init's `n_unconstrained`, so they are
# translated rather than forwarded.
# @noRd
.mashUdControl <- function(opts, nCond) {
    list_modify(
        list(
            unconstrained.update = opts$unconstrainedUpdate %||% "ed",
            scaled.update = "fa",
            resid.update = "none",
            lambda = nCond,
            penalty.type = "iw",
            maxiter = 1000L,
            tol = 1e-2,
            tol.lik = 1e-2
        ),
        !!!compact(opts[setdiff(
            names(opts),
            c("unconstrainedUpdate", "nUnconstrained")
        )])
    )
}

# ud_fit with a directed error for the ud_ted / non-i.i.d. (per-variant SE)
# incompatibility.
# @noRd
#' @importFrom rlang try_fetch
.mashUdFit <- function(fit0, opts, nCond) {
    control <- .mashUdControl(opts, nCond)
    try_fetch(
        udr::ud_fit(fit0, control = control, verbose = FALSE),
        error = function(cnd) {
            if (
                identical(opts$unconstrainedUpdate, "ted") &&
                    str_detect(conditionMessage(cnd), "i.i.d")
            ) {
                msg <- glue(
                    "mashPriorCovariances: the udr TED update needs ",
                    "i.i.d. data (a single shared V), which the beta scale ",
                    "does not provide (per-variant SE). Use ",
                    "covUdrConfig(unconstrainedUpdate = 'ed'), or a z-scale ",
                    "input."
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
    opts <- as.list(engineArgs %||% list())
    # udr distinguishes scaled components from unconstrained ones, which is
    # the same line mashr draws: the fixed hypotheses initialise `U_scaled`,
    # the data-driven components `U_unconstrained`. Previously the
    # data-driven set was discarded here and udr always started from
    # canonical alone.
    scaled <- canonical %||% mashr::cov_canonical(mashData)
    initArgs <- compact(list(
        U_scaled = scaled,
        U_unconstrained = if (length(dataDriven) > 0L) dataDriven,
        n_unconstrained = if (length(dataDriven) == 0L) {
            opts$nUnconstrained %||% 50L
        }
    ))
    fit0 <- exec(udr::ud_init, mashData, !!!initArgs)
    fit <- .mashUdFit(fit0, opts, ncol(mashData$Bhat))
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
#'   \code{\link{mashDataConfig}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
#' @param mashArgs Extra arguments for \code{mashr::mash()}, built with
#'   \code{\link{mashConfig}} -- for example \code{nullweight},
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
#' @export
mashModelFit <- function(
    sumStatsList,
    alpha,
    priorCovariances,
    vhat = NULL,
    fitOn = c("random", "strong"),
    outputLevel = 4L,
    inputScale = c("auto", "beta", "z"),
    mashDataArgs = mashDataConfig(),
    mashArgs = mashConfig(),
    setSeed = 999
) {
    .assertMethodConfig(mashDataArgs, "mashDataConfig", "mashDataArgs")
    .assertMethodConfig(mashArgs, "mashConfig", "mashArgs")
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
#' @param mashDataArgs Extra arguments for
#'   \code{mashr::mash_set_data()}, built with
#'   \code{\link{mashDataConfig}} -- for example
#'   \code{zero_Bhat_Shat_reset} or \code{zero_Shat_reset}. \code{Bhat},
#'   \code{Shat}, \code{alpha} and \code{V} are supplied by pecotmr and
#'   are refused by the constructor.
#' @param posteriorArgs Extra arguments for
#'   \code{mashr::mash_compute_posterior_matrices()}, built with
#'   \code{\link{mashPosteriorConfig}} -- for example
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
    mashDataArgs = mashDataConfig(),
    posteriorArgs = mashPosteriorConfig()
) {
    .assertMethodConfig(mashDataArgs, "mashDataConfig", "mashDataArgs")
    .assertMethodConfig(posteriorArgs, "mashPosteriorConfig", "posteriorArgs")
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
    groupSize <- as.integer(table(grouping)[as.character(grouping)])
    sameGroup <- outer(grouping, grouping, "==") &
        matrix(grouping > 0, nPop, nPop)
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
#'   \code{\link{rmaConfig}} -- \code{test = "knha"} in particular, the
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
    metaArgs = rmaConfig()
) {
    .assertMethodConfig(metaArgs, "rmaConfig", "metaArgs")
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
#'   \code{\link{rmaConfig}} -- \code{test = "knha"} in particular, the
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
    metaArgs = rmaConfig()
) {
    assertString(metaMethod)
    .assertMethodConfig(metaArgs, "rmaConfig", "metaArgs")
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
