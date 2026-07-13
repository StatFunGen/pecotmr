#' @title Run mashr Across Multi-Context QTL or GWAS Summary Statistics
#' @description End-to-end driver: from `{strong, random, null}` sumstats
#'   collections, builds the variant × context Bhat / Shat matrices,
#'   estimates the residual correlation (\code{Vhat}), and fits the mash
#'   model with canonical + PCA + flash + ED covariance components,
#'   returning the fitted covariance list and the estimated mixture
#'   weights.
#' @param sumStatsList Named list (or \code{SimpleList}) of
#'   \code{\link{QtlSumStats}} or \code{\link{GwasSumStats}} objects.
#'   Required names: \code{"strong"} (discovery variants), \code{"random"}
#'   (random background). Optional: \code{"null"} (for residual
#'   correlation estimation).
#' @param alpha Numeric (length 1). Variance-stabilising-transform
#'   exponent forwarded to \code{mashr::mash_set_data()}. Use
#'   \code{alpha = 0} on the BETA scale, \code{alpha = 1} on the Z scale.
#' @param residualCorrelation Optional pre-computed residual correlation
#'   matrix (\code{Vhat}). When supplied, replaces the inline
#'   \code{mashr::estimate_null_correlation_simple()} call entirely and
#'   the \code{"random"} slot of \code{sumStatsList} becomes optional
#'   (the function does not need it for anything else). Useful when
#'   \code{Vhat} was estimated previously on a larger reference and
#'   shipped as a static artefact (the legacy MWE pattern).
#' @param priorCovariances Optional named list of square covariance
#'   matrices (the \code{Ulist} \code{mashr::mash()} consumes). When
#'   supplied, replaces the canonical + PCA + flash + ED chain
#'   (\code{cov_canonical} / \code{cov_pca} / \code{cov_flash} /
#'   \code{cov_ed}) entirely; mash sees only the supplied matrices.
#'   Every entry must be a \code{ncol(Bhat) x ncol(Bhat)} matrix.
#'   Useful when \code{U} was learnt on a larger reference and shipped
#'   as a static artefact (the legacy MWE pattern).
#' @param nPcs Optional integer; number of principal components seeded
#'   into \code{mashr::cov_pca()}. Defaults to \code{ncol(Bhat) - 1}.
#'   Ignored when \code{priorCovariances} is supplied.
#' @param inputScale One of \code{"auto"} (default), \code{"beta"},
#'   \code{"z"}. Controls which (Bhat, Shat) pair is extracted from each
#'   sumstats entry:
#'   \describe{
#'     \item{\code{"beta"}}{Bhat = BETA, Shat = SE — the standard
#'       effect-size scale mashr was designed around. Requires every
#'       entry to carry BETA + SE mcols.}
#'     \item{\code{"z"}}{Bhat = Z, Shat = 1 — z-score scale. Requires Z.}
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
#'   canonical + PCA + flash + ED) and \code{w} (the estimated mixture
#'   weights).
#' @export
mashPipeline <- function(sumStatsList, alpha,
                          residualCorrelation = NULL,
                          priorCovariances = NULL,
                          nPcs = NULL,
                          inputScale = c("auto", "beta", "z"),
                          setSeed = 999) {
  inputScale <- match.arg(inputScale)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
    # nocov end
  }
  if (!requireNamespace("flashier", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install flashier: ",
         "https://github.com/willwerscheid/flashier")
    # nocov end
  }

  # Accept either a base list or a S4Vectors::SimpleList.
  if (methods::is(sumStatsList, "SimpleList")) {
    sumStatsList <- as.list(sumStatsList)
  }
  if (!is.list(sumStatsList) || is.null(names(sumStatsList))) {
    stop("mashPipeline: `sumStatsList` must be a named list (or SimpleList) ",
         "of QtlSumStats / GwasSumStats objects, named with at least ",
         "'strong' and 'random' (optionally 'null').")
  }
  # `random` is required only when we need to derive vhat from data; when
  # the caller supplies a pre-computed `residualCorrelation`, the random
  # subset is no longer needed.
  required <- if (is.null(residualCorrelation)) c("strong", "random") else "strong"
  missingNames <- setdiff(required, names(sumStatsList))
  if (length(missingNames) > 0L) {
    stop("mashPipeline: `sumStatsList` is missing required entr",
         if (length(missingNames) == 1L) "y: " else "ies: ",
         paste(shQuote(missingNames), collapse = ", "), ".")
  }
  extraNames <- setdiff(names(sumStatsList), c("strong", "random", "null"))
  if (length(extraNames) > 0L) {
    stop("mashPipeline: `sumStatsList` has unrecognised entries: ",
         paste(shQuote(extraNames), collapse = ", "),
         ". Only 'strong', 'random', and 'null' are accepted.")
  }

  set.seed(setSeed)

  # Vhat: use the supplied residual correlation if given; otherwise derive it
  # from the null set (estimate_null_correlation_simple) when a `null` entry is
  # present, else fall back to identity. The estimation now lives in exactly one
  # place -- mashResidualCorrelation(). `setSeed = NULL` there leaves the RNG
  # stream seeded above untouched, so the delegated calls consume it in the same
  # order the inline code used to (identity / null-simple consume no RNG before
  # cov_flash / cov_ed run downstream).
  vhat <- if (!is.null(residualCorrelation)) {
    residualCorrelation
  } else if ("null" %in% names(sumStatsList) && !is.null(sumStatsList$null)) {
    mashResidualCorrelation(sumStatsList, alpha, method = "simple",
                            inputScale = inputScale, setSeed = NULL)
  } else {
    mashResidualCorrelation(sumStatsList, alpha, method = "identity",
                            inputScale = inputScale, setSeed = NULL)
  }

  # Prior covariances + mixture weights. mashPriorCovariances() owns the
  # cov_canonical / cov_pca / cov_flash / cov_ed chain, the supplied-prior
  # bypass, and the mash() weight fit; mashPipeline just forwards its arguments.
  prior <- mashPriorCovariances(sumStatsList, alpha, vhat = vhat,
                                priorCovariances = priorCovariances,
                                nPcs = nPcs, inputScale = inputScale,
                                setSeed = NULL)

  list(U = prior$U, w = prior$w)
}

#' @title Estimate the mash Residual Correlation Matrix (Vhat)
#' @description Estimate the residual (null) correlation matrix \eqn{\hat V} that
#'   \code{mashr::mash_set_data()} consumes as \code{V}. Split out of
#'   \code{\link{mashPipeline}} so the estimation logic lives in one place and is
#'   reusable by the mixture-prior workflows.
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
#' @export
mashResidualCorrelation <- function(sumStatsList, alpha,
                                    method = c("simple", "identity", "mle",
                                               "corshrink", "simple_specific"),
                                    priorCovariances = NULL,
                                    nSubset = 6000L, maxIter = 6L,
                                    inputScale = c("auto", "beta", "z"),
                                    setSeed = 999) {
  method <- match.arg(method)
  inputScale <- match.arg(inputScale)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
    # nocov end
  }
  if (methods::is(sumStatsList, "SimpleList")) {
    sumStatsList <- as.list(sumStatsList)
  }
  if (!is.null(setSeed)) set.seed(setSeed)

  if (method == "identity") {
    strongMats <- .mashSumStatsToMatrices(sumStatsList$strong, "strong",
                                          inputScale = inputScale)
    return(diag(rep(1, ncol(strongMats$b))))
  }

  if (method == "simple") {
    if (is.null(sumStatsList$null)) {
      stop("mashResidualCorrelation: method 'simple' requires a 'null' entry ",
           "in `sumStatsList`.")
    }
    nullMats <- .mashSumStatsToMatrices(sumStatsList$null, "null",
                                        inputScale = inputScale)
    return(mashr::estimate_null_correlation_simple(
      mashr::mash_set_data(nullMats$b, Shat = nullMats$s,
                           alpha, zero_Bhat_Shat_reset = 1000)))
  }

  if (method == "mle") {
    if (is.null(sumStatsList$random)) {
      stop("mashResidualCorrelation: method 'mle' requires a 'random' entry ",
           "in `sumStatsList`.")
    }
    if (is.null(priorCovariances)) {
      stop("mashResidualCorrelation: method 'mle' requires `priorCovariances` ",
           "(the prior U the EM refines V against).")
    }
    randomMats <- .mashSumStatsToMatrices(sumStatsList$random, "random",
                                          inputScale = inputScale)
    n <- nrow(randomMats$b)
    idx <- sample(seq_len(n), min(nSubset, n))
    dsub <- mashr::mash_set_data(randomMats$b[idx, , drop = FALSE],
                                 Shat = randomMats$s[idx, , drop = FALSE],
                                 alpha, zero_Bhat_Shat_reset = 1000)
    fit <- mashr::mash_estimate_corr_em(dsub, priorCovariances,
                                        max_iter = maxIter, details = TRUE)
    return(fit$V)
  }

  # method %in% c("corshrink", "simple_specific"): estimate V on the null set.
  # In this pipeline the `null` partition is already the null variants
  # (mashRandNullSample kept max|z| < 2), so no re-thresholding is needed --
  # the GTEx-specific per-gene `GetSS()` fetch of the legacy notebook is
  # replaced by our own null z-matrix.
  if (is.null(sumStatsList$null)) {
    stop("mashResidualCorrelation: method '", method, "' requires a 'null' ",
         "entry in `sumStatsList` (the null variants V is estimated on).")
  }
  nullMats <- .mashSumStatsToMatrices(sumStatsList$null, "null",
                                      inputScale = inputScale)
  # Null z-scores (z-scale: b already = Z, s = 1; beta-scale: b/s = beta/se).
  nullZ <- nullMats$b / nullMats$s

  if (method == "simple_specific") {
    return(as.matrix(Matrix::nearPD(stats::cov(nullZ), conv.tol = 1e-06,
                                    doSym = TRUE, corr = TRUE)$mat))
  }

  # method == "corshrink": adaptive-shrinkage null correlation (ashr on the
  # Fisher-z of the pairwise sample correlations); image = "null" suppresses
  # CorShrink's plotting.
  if (!requireNamespace("CorShrink", quietly = TRUE)) {
    # nocov start  (optional-package guard; CorShrink is Suggests-only)
    stop("mashResidualCorrelation: method 'corshrink' needs the CorShrink ",
         "package. Install it, or use 'simple' / 'simple_specific'.")
    # nocov end
  }
  as.matrix(CorShrink::CorShrinkData(
    nullZ, ash.control = list(mixcompdist = "halfuniform"),
    image = "null")$cor)
}

# Internal: build the requested data-driven covariance components off a prepared
# mashr data object, in a fixed order (canonical, pca, flash, flash_nonneg) so a
# given `components` set reproduces the same ordering everywhere. RNG is consumed
# only by cov_flash / cov_flash(nonneg). Shared by mashPriorCovariances and the
# exported mashCovarianceComponents.
# @noRd
.mashBuildComponents <- function(mashData, components, nPcs = NULL) {
  comps <- list()
  if ("canonical" %in% components) {
    comps <- c(comps, mashr::cov_canonical(mashData))
  }
  if ("pca" %in% components) {
    if (is.null(nPcs)) nPcs <- ncol(mashData$Bhat) - 1
    comps <- c(comps, mashr::cov_pca(mashData, npc = nPcs))
  }
  if ("flash" %in% components) {
    comps <- c(comps, mashr::cov_flash(mashData))
  }
  if ("flash_nonneg" %in% components) {
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
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) with at least
#'   \code{"strong"} (the discovery set the covariances are learned on).
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
#' @export
mashCovarianceComponents <- function(sumStatsList, alpha, vhat = NULL,
                                     components = c("canonical", "pca", "flash",
                                                   "flash_nonneg"),
                                     nPcs = NULL,
                                     inputScale = c("auto", "beta", "z"),
                                     setSeed = 999) {
  inputScale <- match.arg(inputScale)
  components <- as.character(components)
  validComponents <- c("canonical", "pca", "flash", "flash_nonneg")
  if (!all(components %in% validComponents)) {
    stop("mashCovarianceComponents: unknown component(s): ",
         paste(setdiff(components, validComponents), collapse = ", "),
         ". Valid: ", paste(validComponents, collapse = ", "), ".")
  }
  if (!requireNamespace("mashr", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
    # nocov end
  }
  if (!requireNamespace("flashier", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install flashier: ",
         "https://github.com/willwerscheid/flashier")
    # nocov end
  }
  if (methods::is(sumStatsList, "SimpleList")) {
    sumStatsList <- as.list(sumStatsList)
  }
  if (!is.null(setSeed)) set.seed(setSeed)
  strongMats <- .mashSumStatsToMatrices(sumStatsList$strong, "strong",
                                        inputScale = inputScale)
  if (is.null(vhat)) vhat <- diag(rep(1, ncol(strongMats$b)))
  mashData <- mashr::mash_set_data(strongMats$b, Shat = strongMats$s,
                                   V = vhat, alpha, zero_Bhat_Shat_reset = 1000)
  .mashBuildComponents(mashData, components = components, nPcs = nPcs)
}

#' @title Estimate mash Prior Covariances and Mixture Weights
#' @description Build the mash prior covariance list \code{U} and mixture weights
#'   \code{w} from a \code{{strong[, random, null]}} SumStats collection. Split
#'   out of \code{\link{mashPipeline}} so the covariance-estimation chain lives in
#'   one place. The default builds every non-udr covariance component
#'   (\code{canonical + pca + flash + flash_nonneg}) and refines them with mashr
#'   extreme deconvolution (\code{cov_ed}).
#' @param sumStatsList Named list (or \code{S4Vectors::SimpleList}) with at least
#'   \code{"strong"} (the discovery set the covariances are learned on).
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
#'   non-empty named list, e.g. the concatenated \code{\link{mashCovarianceComponents}}
#'   outputs). When supplied, these are refined by \code{engine} instead of being
#'   rebuilt internally -- the mixture-prior pipeline where separate steps built
#'   the components. \code{components} / \code{nPcs} are then ignored. Distinct
#'   from \code{priorCovariances}, which bypasses the engine entirely.
#' @param udControl Named list overriding the \pkg{udr} controls for the
#'   \code{"ud"} / \code{"ud_ted"} engines: \code{n_unconstrained} (data-driven
#'   matrices to fit; default 50 -- the dominant cost, so reduce it for
#'   few-condition data), \code{maxiter} (default 1000), \code{tol},
#'   \code{tol.lik}. Ignored by \code{cov_ed}.
#' @param inputScale SumStats -> matrix conversion scale.
#' @param setSeed Integer seed, or \code{NULL} to leave the ambient RNG untouched.
#' @return \code{list(U, w, loglik)}: the covariance list, the
#'   \code{mashr::get_estimated_pi()} mixture weights, and the fit
#'   log-likelihood (\code{NULL} for the \code{cov_ed} engine).
#' @seealso \code{\link{mashPipeline}}, \code{\link{mashResidualCorrelation}}
#' @export
mashPriorCovariances <- function(sumStatsList, alpha, vhat = NULL,
                                 components = c("canonical", "pca", "flash",
                                               "flash_nonneg"),
                                 engine = c("cov_ed", "ud", "ud_ted"),
                                 nPcs = NULL,
                                 priorCovariances = NULL,
                                 priorComponents = NULL,
                                 udControl = list(),
                                 inputScale = c("auto", "beta", "z"),
                                 setSeed = 999) {
  engine <- match.arg(engine)
  inputScale <- match.arg(inputScale)
  # udr controls (engine = "ud"/"ud_ted" only). Defaults are the legacy
  # mixture-prior values, but `n_unconstrained` dominates the udr cost -- it
  # generates that many full data-driven covariances, so 50 (sized for
  # many-condition data) massively overparameterises small problems. Reduce it
  # for few-condition data.
  udControl <- utils::modifyList(
    list(n_unconstrained = 50L, maxiter = 1000L, tol = 1e-2, tol.lik = 1e-2),
    udControl)
  components <- as.character(components)
  validComponents <- c("canonical", "pca", "flash", "flash_nonneg")
  if (!all(components %in% validComponents)) {
    stop("mashPriorCovariances: unknown component(s): ",
         paste(setdiff(components, validComponents), collapse = ", "),
         ". Valid: ", paste(validComponents, collapse = ", "), ".")
  }
  if (!requireNamespace("mashr", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
    # nocov end
  }
  if (!requireNamespace("flashier", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install flashier: ",
         "https://github.com/willwerscheid/flashier")
    # nocov end
  }
  if (methods::is(sumStatsList, "SimpleList")) {
    sumStatsList <- as.list(sumStatsList)
  }
  if (!is.null(setSeed)) set.seed(setSeed)

  strongMats <- .mashSumStatsToMatrices(sumStatsList$strong, "strong",
                                        inputScale = inputScale)
  if (is.null(vhat)) vhat <- diag(rep(1, ncol(strongMats$b)))
  mashData <- mashr::mash_set_data(strongMats$b, Shat = strongMats$s,
                                   V = vhat, alpha, zero_Bhat_Shat_reset = 1000)

  if (!is.null(priorCovariances)) {
    # Caller-supplied prior covariance matrices (e.g. pre-baked U from a
    # reference dataset). Bypasses the cov_* chain; mashr sees only these.
    if (!is.list(priorCovariances) || length(priorCovariances) == 0L ||
        is.null(names(priorCovariances)) || any(names(priorCovariances) == "")) {
      stop("mashPriorCovariances: `priorCovariances` must be a non-empty named ",
           "list of square covariance matrices.")
    }
    nCond <- ncol(mashData$Bhat)
    bad <- !vapply(priorCovariances, function(M) {
      is.matrix(M) && nrow(M) == nCond && ncol(M) == nCond
    }, logical(1L))
    if (any(bad)) {
      stop(sprintf(paste0("mashPriorCovariances: every `priorCovariances` entry ",
                          "must be a %d x %d matrix; offenders: %s."),
                   nCond, nCond,
                   paste(names(priorCovariances)[bad], collapse = ", ")))
    }
    U.all <- priorCovariances
    w <- NULL
    loglik <- NULL
  } else {
    # Data-driven covariance chain. Use caller-supplied raw components when given
    # -- the mixture-prior pipeline, where the per-method component steps built
    # them (via mashCovarianceComponents) and this step only refines + weights.
    # Otherwise build them here with the same builder mashCovarianceComponents
    # surfaces. Either way the chosen engine refines + weights.
    if (!is.null(priorComponents)) {
      if (!is.list(priorComponents) || length(priorComponents) == 0L ||
          is.null(names(priorComponents)) || any(names(priorComponents) == "")) {
        stop("mashPriorCovariances: `priorComponents` must be a non-empty named ",
             "list of covariance matrices (e.g. mashCovarianceComponents() output).")
      }
      comps <- priorComponents
    } else {
      comps <- .mashBuildComponents(mashData, components = components, nPcs = nPcs)
    }

    if (engine == "cov_ed") {
      # cov_ed() is mashr's exported extreme-deconvolution wrapper; its default
      # algorithm = "bovy" IS the internal bovy_wrapper (Bovy et al. 2011), so
      # this is the exported route to bovy ED -- no internal reach-in needed.
      # It refines the components; mash() (below) then estimates the mixture
      # weights over the components + their ED versions (bovy_wrapper's own $pi
      # is not needed -- get_estimated_pi over the full set is the proper weight).
      U.ed <- mashr::cov_ed(mashData, Ulist_init = comps)
      U.all <- c(comps, U.ed)
      w <- NULL
      loglik <- NULL
    } else {
      # engine %in% c("ud", "ud_ted") -- udr (OPT-IN; known numerical issues, so
      # not the default). udr seeds canonical as its scaled prior and generates
      # n_unconstrained data-driven matrices, returning U + w + loglik directly;
      # `components` is not used on this path.
      if (!requireNamespace("udr", quietly = TRUE)) {
        # nocov start  (optional-package guard; udr is Suggests-only)
        stop("mashPriorCovariances: engine '", engine, "' needs the udr ",
             "package. Install it, or use the default 'cov_ed'.")
        # nocov end
      }
      U.can <- mashr::cov_canonical(mashData)
      fit0 <- udr::ud_init(mashData,
                           n_unconstrained = udControl$n_unconstrained,
                           U_scaled = U.can)
      fit <- tryCatch(
        udr::ud_fit(fit0, control = list(
          unconstrained.update = if (engine == "ud_ted") "ted" else "ed",
          scaled.update = "fa", resid.update = "none",
          lambda = ncol(mashData$Bhat), penalty.type = "iw",
          maxiter = udControl$maxiter, tol = udControl$tol,
          tol.lik = udControl$tol.lik), verbose = FALSE),
        error = function(e) {
          # udr's TED update only supports i.i.d. data (a single shared V). On
          # the beta scale each variant carries its own SE, so V differs per
          # point and TED errors -- surface a directed message.
          if (engine == "ud_ted" && grepl("i.i.d", conditionMessage(e))) {
            stop("mashPriorCovariances: engine 'ud_ted' (udr TED update) needs ",
                 "i.i.d. data (a single shared V), which the beta scale does not ",
                 "provide (per-variant SE). Use engine 'ud' (ED update), or a ",
                 "z-scale input.", call. = FALSE)
          }
          stop(e)
        })
      U.all <- lapply(fit$U, function(e) e$mat)
      w <- fit$w
      loglik <- fit$loglik
    }
  }

  # cov_ed and the supplied-prior bypass learn weights with a final mash() fit;
  # the udr engines already produced them.
  if (is.null(w)) {
    m <- mashr::mash(mashData, Ulist = U.all, outputlevel = 1)
    w <- mashr::get_estimated_pi(m)
  }
  list(U = U.all, w = w, loglik = loglik)
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
#' @param setSeed Integer seed, or \code{NULL} to leave the ambient RNG untouched.
#' @return The fitted \pkg{mashr} model (the \code{mashr::mash()} object).
#' @seealso \code{\link{mashPosterior}}, \code{\link{mashPriorCovariances}}
#' @export
mashModelFit <- function(sumStatsList, alpha, priorCovariances, vhat = NULL,
                         fitOn = c("random", "strong"), outputLevel = 4L,
                         inputScale = c("auto", "beta", "z"), setSeed = 999) {
  fitOn <- match.arg(fitOn)
  inputScale <- match.arg(inputScale)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
    # nocov end
  }
  if (methods::is(sumStatsList, "SimpleList")) {
    sumStatsList <- as.list(sumStatsList)
  }
  if (is.null(priorCovariances) || !is.list(priorCovariances) ||
      length(priorCovariances) == 0L || is.null(names(priorCovariances))) {
    stop("mashModelFit: `priorCovariances` must be a non-empty named list of ",
         "covariance matrices (e.g. mashPriorCovariances()$U).")
  }
  if (is.null(sumStatsList[[fitOn]])) {
    stop("mashModelFit: `sumStatsList` has no '", fitOn, "' entry to fit on.")
  }
  if (!is.null(setSeed)) set.seed(setSeed)
  mats <- .mashSumStatsToMatrices(sumStatsList[[fitOn]], fitOn,
                                  inputScale = inputScale)
  if (is.null(vhat)) vhat <- diag(rep(1, ncol(mats$b)))
  mashData <- mashr::mash_set_data(mats$b, Shat = mats$s, V = vhat, alpha,
                                   zero_Bhat_Shat_reset = 1000)
  mashr::mash(mashData, Ulist = priorCovariances, outputlevel = outputLevel)
}

#' @title Compute mash Posterior Matrices for a Target Set
#' @description Apply a fitted \pkg{mashr} model (from \code{\link{mashModelFit}})
#'   to a target SumStats set, returning the posterior matrices
#'   (\code{PosteriorMean}, \code{PosteriorSD}, \code{lfsr}, ...). This is the
#'   posterior step of the mash workflow -- \code{mash_fit}'s second step and
#'   \code{mash_posterior}'s per-analysis-unit computation.
#' @param model A fitted \pkg{mashr} model (the \code{\link{mashModelFit}}
#'   output).
#' @param sumStats A single \code{\link{QtlSumStats}} / \code{\link{GwasSumStats}}
#'   -- the target (e.g. strong) set to compute posteriors on.
#' @param alpha mash \code{alpha} (match the fit).
#' @param vhat Residual correlation matrix (\code{V}); \code{NULL} -> identity.
#' @param excludeCondition Character vector of condition (column) names to drop
#'   from the target AND the model before computing posteriors (the model's
#'   covariances are resized via \code{\link{updateMashModelCov}}). Default none.
#' @param outputPosteriorCov Return the full posterior covariance array (needed
#'   by \code{\link{fitMashContrast}}). Default \code{TRUE}.
#' @param inputScale SumStats -> matrix conversion scale.
#' @return The \code{mashr::mash_compute_posterior_matrices()} result: a list of
#'   \code{PosteriorMean} / \code{PosteriorSD} / \code{lfsr} / \code{NegativeProb}
#'   (+ \code{PosteriorCov} when \code{outputPosteriorCov = TRUE}).
#' @seealso \code{\link{mashModelFit}}, \code{\link{fitMashContrast}}
#' @export
mashPosterior <- function(model, sumStats, alpha, vhat = NULL,
                          excludeCondition = character(0),
                          outputPosteriorCov = TRUE,
                          inputScale = c("auto", "beta", "z")) {
  inputScale <- match.arg(inputScale)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    # nocov start
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
    # nocov end
  }
  mats <- .mashSumStatsToMatrices(sumStats, "target", inputScale = inputScale)
  b <- mats$b
  s <- mats$s
  allConditions <- colnames(b)
  excludeCondition <- as.character(excludeCondition)
  if (length(excludeCondition) > 0L) {
    bad <- setdiff(excludeCondition, allConditions)
    if (length(bad) > 0L) {
      stop("mashPosterior: excludeCondition not found in the target: ",
           paste(bad, collapse = ", "), ".")
    }
    keep <- setdiff(allConditions, excludeCondition)
    if (length(keep) == 0L) {
      stop("mashPosterior: excludeCondition drops every condition.")
    }
    # Index positionally: b/s carry condition colnames, but a supplied `vhat`
    # (e.g. diag(K)) may have no dimnames, so a name subscript would error.
    keepIdx <- match(keep, allConditions)
    b <- b[, keepIdx, drop = FALSE]
    s <- s[, keepIdx, drop = FALSE]
    if (!is.null(vhat)) vhat <- vhat[keepIdx, keepIdx, drop = FALSE]
    model <- updateMashModelCov(model, allConditions, keep)
  }
  if (is.null(vhat)) vhat <- diag(rep(1, ncol(b)))
  mashData <- mashr::mash_set_data(b, Shat = s, V = vhat, alpha,
                                   zero_Bhat_Shat_reset = 1000)
  mashr::mash_compute_posterior_matrices(
    model, mashData, output_posterior_cov = outputPosteriorCov)
}

# =============================================================================
# Mash pairwise contrast functions
# =============================================================================

#' Create a pairwise contrast column
#'
#' Sets +1 for the first condition and -1 for the second in a zero vector.
#' Used as a building block for contrast design matrices.
#'
#' @param pair A length-2 character vector naming the two conditions to contrast.
#' @param template A named numeric vector of zeros with names matching all
#'   conditions.
#' @return The template vector with +1 at \code{pair[1]} and -1 at
#'   \code{pair[2]}.
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
#' @param origMean Matrix of original effect sizes (variants x conditions).
#'   Used to determine which conditions are "tested" (non-zero).
#' @param posteriorMean Matrix of mash posterior means (variants x conditions).
#' @param posteriorVcov 3D array of posterior covariance matrices
#'   (conditions x conditions x variants).
#' @param grouping Named integer vector mapping condition names to group IDs.
#'   Conditions with the same positive group ID are treated as replicates
#'   (e.g., multiple datasets for the same cell type). Use 0 for ungrouped.
#'   If NULL (default), all conditions are treated independently.
#' @return A single-row data.frame with columns
#'   \code{mean_contrast_*}, \code{se_contrast_*}, \code{p_contrast_*} for
#'   both deviation and pairwise contrasts. Returns NULL if fewer than 2
#'   tested conditions.
#' @importFrom stringr str_remove_all
#' @importFrom utils combn
#' @export
fitMashContrast <- function(index, origMean, posteriorMean, posteriorVcov,
                               grouping = NULL) {
  populationNames <- colnames(posteriorMean)
  if (!is.null(populationNames))
    populationNames <- str_remove_all(populationNames, "BETA_")

  origMeanVector <- origMean[index, ]
  names(origMeanVector) <- populationNames
  tested <- names(origMeanVector[origMeanVector != 0])

  if (length(tested) < 2) return(NULL)

  nPop <- length(tested)
  pairwiseVector <- setNames(rep(0, nPop), tested)

  # Default grouping: all independent

  if (is.null(grouping)) {
    grouping <- setNames(rep(0L, nPop), tested)
  } else {
    grouping <- grouping[tested]
  }

  if (nPop > 2) {
    # 1. Deviation contrasts
    dev <- matrix(-1, nPop, nPop, dimnames = list(tested, tested))
    diag(dev) <- nPop - 1

    # Adjust for grouped conditions
    uniqueGroups <- unique(grouping)
    for (grp in uniqueGroups[uniqueGroups > 0]) {
      grpMask <- grouping == grp
      grpSize <- sum(grpMask)
      diag(dev)[grpMask] <- (nPop - 1) / grpSize
      dev[grpMask, grpMask] <- (nPop - 1) / grpSize
    }
    colnames(dev) <- paste0(tested, "_deviation")

    # 2. Pairwise contrasts
    twoCombn <- combn(tested, 2)
    pwNames <- apply(twoCombn, 2, paste, collapse = "_vs_")
    pw <- apply(twoCombn, 2, makePairwiseContrastCol, pairwiseVector)
    colnames(pw) <- pwNames

    # Adjust pairwise contrasts for grouped conditions
    pwAdj <- pw
    for (col in colnames(pw)) {
      groups <- strsplit(col, "_vs_")[[1]]
      groupValues <- grouping[names(grouping) %in% groups]
      relevant <- names(groupValues[groupValues > 0])
      if (length(unique(groupValues)) > 1 && length(relevant) > 0) {
        for (dg in unique(groupValues[groupValues > 0])) {
          rowsInGroup <- names(grouping[grouping == dg])
          matchedRow <- rowsInGroup[rowsInGroup %in% groups]
          if (length(matchedRow) > 0)
            pwAdj[rowsInGroup, col] <- pw[matchedRow, col] / length(rowsInGroup)
        }
      }
    }

    contrastDesign <- cbind(dev / (nPop - 1), pwAdj)
  } else {
    pairwiseVector[tested[1]] <- 1
    pairwiseVector[tested[2]] <- -1
    contrastDesign <- matrix(pairwiseVector, ncol = 1,
                              dimnames = list(tested, paste0(tested[1], "_vs_", tested[2])))
  }

  # Subset posterior to tested conditions
  pm <- posteriorMean[index, tested]
  pv <- posteriorVcov[tested, tested, index]

  # Compute contrasts
  contrastDiff <- drop(t(contrastDesign) %*% pm)
  contrastVcov <- t(contrastDesign) %*% pv %*% contrastDesign
  contrastSe <- sqrt(diag(contrastVcov))
  contrastP <- .zToPvalue(contrastDiff / contrastSe)

  # Build output data.frame
  cnames <- colnames(contrastDesign)
  df <- data.frame(
    row.names = rownames(posteriorMean)[index],
    stringsAsFactors = FALSE)
  for (i in seq_along(cnames)) {
    df[[paste0("mean_contrast_", cnames[i])]] <- contrastDiff[i]
    df[[paste0("se_contrast_", cnames[i])]] <- contrastSe[i]
    df[[paste0("p_contrast_", cnames[i])]] <- contrastP[i]
  }
  df
}

#' Posterior contrast table over an entire mash posterior
#'
#' Orchestrates \code{\link{fitMashContrast}} across every feature (variant)
#' of a mash posterior: aligns the original effect matrix to the posterior
#' columns, runs the per-feature deviation + pairwise contrast, row-binds the
#' results (aligning the union of contrast columns), and orders them
#' (\code{mean}/\code{se}/\code{p}, deviation before pairwise). Features with
#' fewer than two tested conditions are dropped.
#'
#' @param posteriorMean Numeric matrix (features x conditions) of posterior
#'   means (\code{PosteriorMean} from \code{\link{mashPosterior}}).
#' @param posteriorVcov Numeric array (conditions x conditions x features) of
#'   posterior covariances (\code{PosteriorCov}).
#' @param origMean Numeric matrix (features x conditions) of the original
#'   effect estimates (e.g. \code{bhat}); used to decide which conditions were
#'   tested per feature. Aligned to \code{posteriorMean}'s columns by name;
#'   \code{NaN} is treated as 0 (untested).
#' @param grouping Optional named integer vector assigning conditions to groups
#'   (0 = independent); forwarded to \code{\link{fitMashContrast}} so replicate
#'   populations of one cell type share weight. Names are condition labels.
#' @return A \code{data.frame} (features x contrasts) with
#'   \code{mean_contrast_*}, \code{se_contrast_*}, \code{p_contrast_*} columns
#'   and feature ids as row names. Empty when nothing is testable.
#' @seealso \code{\link{fitMashContrast}}, \code{\link{metaAnalysisPerCondition}}
#' @importFrom dplyr bind_rows select matches
#' @export
mashPosteriorContrast <- function(posteriorMean, posteriorVcov, origMean,
                                  grouping = NULL) {
  origMean <- origMean[, colnames(posteriorMean), drop = FALSE]
  origMean[is.nan(origMean)] <- 0

  parts <- lapply(seq_len(nrow(posteriorMean)), function(i)
    fitMashContrast(i, origMean, posteriorMean, posteriorVcov,
                    grouping = grouping))
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) return(data.frame())

  # fitMashContrast sets each 1-row frame's rowname to the feature id, but
  # bind_rows drops row names -- capture them first, restore after.
  featureIds <- unlist(lapply(parts, rownames), use.names = FALSE)
  res <- as.data.frame(dplyr::bind_rows(parts))
  res <- dplyr::select(res,
    dplyr::matches("mean_contrast.*deviation"), dplyr::matches("mean_contrast.*_vs_"),
    dplyr::matches("se_contrast.*deviation"),   dplyr::matches("se_contrast.*_vs_"),
    dplyr::matches("p_contrast.*deviation"),    dplyr::matches("p_contrast.*_vs_"))
  rownames(res) <- featureIds
  res
}

# =============================================================================
# Mash model subsetting functions
# =============================================================================

#' Subset a fitted mash model to a subset of conditions
#'
#' Updates the prior covariance matrices (\code{Ulist}) and mixture weights
#' (\code{pi}) in a fitted \code{mashr} model to match a reduced set of
#' conditions. Handles condition-specific, identity, and data-driven
#' covariance components.
#'
#' @param mashModel A fitted mash model object (from \code{mashr::mash}).
#' @param allSamples Character vector of all original condition names.
#' @param samples Character vector of the conditions to retain.
#' @return The updated mash model with resized covariance matrices and
#'   pruned mixture weights.
#' @export
updateMashModelCov <- function(mashModel, allSamples, samples) {
  cov <- mashModel$fitted_g$Ulist

  # Remove matrices for dropped conditions
  unwanted <- setdiff(allSamples, samples)
  for (d in names(cov)) {
    if (d %in% unwanted || d %in% paste0("ED_", unwanted))
      cov[[d]] <- NULL
  }

  # Resize remaining matrices to match retained conditions
  for (d in names(cov)) {
    if (d %in% samples) {
      # Condition-specific: single 1 on diagonal
      m <- matrix(0, length(samples), length(samples))
      m[which(samples == d), which(samples == d)] <- 1
      cov[[d]] <- m
    } else if (d == "identity") {
      m <- matrix(0, length(samples), length(samples))
      m[1, 1] <- 1
      cov[[d]] <- m
    } else if (is.null(colnames(cov[[d]]))) {
      cov[[d]] <- cov[[d]][seq_len(length(samples)), seq_len(length(samples))]
    } else {
      cov[[d]] <- cov[[d]][samples, samples]
    }
    cov[[d]] <- as.matrix(cov[[d]])
  }

  mashModel$fitted_g$Ulist <- cov

  # Prune mixture weights for removed conditions
  for (s in unwanted) {
    dropIdx <- grep(s, names(mashModel$fitted_g$pi), fixed = TRUE)
    if (length(dropIdx) > 0)
      mashModel$fitted_g$pi <- mashModel$fitted_g$pi[-dropIdx]
  }

  mashModel
}

#' Subset mash data matrices to specific SNPs and conditions
#'
#' Slices the \code{bhat}, \code{sbhat}, and \code{Z} matrices by row (SNPs)
#' and column (samples/conditions), and correspondingly subsets the \code{vhat}
#' covariance matrix.
#'
#' @param data A mash data list with elements \code{bhat}, \code{sbhat},
#'   \code{Z} (matrices), and \code{snp} (character vector).
#' @param vhat A square covariance matrix (conditions x conditions).
#' @param snps Character vector of SNP IDs to retain (row names).
#' @param samples Character vector of condition names to retain (column names).
#' @return A list with \code{data} (sliced data list) and \code{vhat}
#'   (sliced covariance matrix).
#' @export
sliceMashData <- function(data, vhat, snps, samples) {
  data$bhat <- as.matrix(data$bhat[snps, samples])
  data$sbhat <- as.matrix(data$sbhat[snps, samples])
  data$Z <- as.matrix(data$Z[snps, samples])
  vhat <- as.matrix(vhat[samples, samples])
  data$snp <- data$snp[data$snp %in% snps]
  colnames(data$bhat) <- colnames(data$sbhat) <- colnames(data$Z) <- colnames(vhat) <- samples
  list(data = data, vhat = vhat)
}

#' Sanitize NaN/Inf values in mash data
#'
#' Replaces NaN in \code{bhat} with 0 and NaN/Inf in \code{sbhat} with 1e3
#' (indicating high uncertainty).
#'
#' @param data A mash data list with \code{bhat} and \code{sbhat} matrices.
#' @return The data list with sanitized values.
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
#' @param effectSizes Numeric matrix (features x contrasts) of contrast
#'   effect sizes. Column names must follow the pattern
#'   \code{mean_contrast_<conditionA>_vs_<conditionB>}.
#' @param seValues Numeric matrix (features x contrasts) of contrast
#'   standard errors. Must have the same dimensions and column names as
#'   \code{effectSizes}.
#' @param seCutoff Numeric; minimum SE below which a contrast is excluded
#'   from the meta-analysis for a given feature (default 0).
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
#' @export
metaAnalysisPerCondition <- function(effectSizes, seValues,
                                     seCutoff = 0, metaMethod = "DL") {
  stopifnot(identical(dim(effectSizes), dim(seValues)))
  stopifnot(identical(colnames(effectSizes), colnames(seValues)))

  contrasts <- sub("^mean_contrast_", "", colnames(effectSizes))
  conditions <- unique(c(sub("_vs_.*", "", contrasts),
                         sub(".*_vs_", "", contrasts)))

  results <- list()
  for (condition in conditions) {
    # Contrast columns involving this condition
    idx <- grep(condition, colnames(effectSizes))
    if (length(idx) == 0) next

    condEffects <- effectSizes[, idx, drop = FALSE]
    condSes <- seValues[, idx, drop = FALSE]
    condContrasts <- contrasts[idx]

    for (i in seq_along(condContrasts)) {
      es <- abs(as.numeric(condEffects[, i]))
      se <- as.numeric(condSes[, i])

      # Filter by SE cutoff
      keep <- se > seCutoff & is.finite(es) & is.finite(se)
      es <- es[keep]
      se <- se[keep]

      if (length(es) < 2) {
        results[[length(results) + 1]] <- tibble(
          condition = condition,
          contrast = condContrasts[i],
          meta_pvalue = if (length(es) == 1) {
            .zToPvalue(es / se)
          } else NA_real_,
          meta_effect = if (length(es) == 1) es else NA_real_,
          meta_se = if (length(es) == 1) se else NA_real_,
          tau2 = NA_real_,
          I2 = NA_real_
        )
        next
      }

      ma <- .rmaMeta(es, se, method = metaMethod)
      z <- ma$mean / ma$se
      results[[length(results) + 1]] <- tibble(
        condition = condition,
        contrast = condContrasts[i],
        meta_pvalue = .zToPvalue(z),
        meta_effect = ma$mean,
        meta_se = ma$se,
        tau2 = ma$tau2,
        I2 = ma$I2
      )
    }
  }

  bind_rows(results)
}

#' Feature score from deviation contrasts (random-effects meta per condition)
#'
#' For each condition's deviation contrast, meta-analyzes the per-variant
#' absolute effect sizes (random-effects, via \code{metafor::rma}) and returns
#' the pooled Z-score (\eqn{\hat\mu / \mathrm{se}}). One score per condition --
#' the "meta" feature score of \code{mash_posterior.ipynb}.
#'
#' @param contrastResult A contrast table from \code{\link{mashPosteriorContrast}}
#'   (variants x contrasts) carrying \code{mean_contrast_*_deviation} and
#'   \code{se_contrast_*_deviation} columns.
#' @param metaMethod Between-study variance estimator forwarded to
#'   \code{metafor::rma} (default \code{"REML"}).
#' @return A \code{data.frame} with \code{condition} and \code{zScore}.
#' @seealso \code{\link{nSignificantScore}}, \code{\link{scoreFromCs}}
#' @export
calculateFeatureScores <- function(contrastResult, metaMethod = "REML") {
  cr <- as.data.frame(contrastResult, stringsAsFactors = FALSE)
  effCols <- grep("mean_contrast_.*deviation", names(cr), value = TRUE)
  if (length(effCols) == 0L)
    return(data.frame(condition = character(0), zScore = numeric(0),
                      stringsAsFactors = FALSE))
  scores <- vapply(effCols, function(ec) {
    seCol <- sub("^mean_contrast", "se_contrast", ec)
    if (!seCol %in% names(cr)) return(NA_real_)
    es <- abs(as.numeric(cr[[ec]]))
    se <- as.numeric(cr[[seCol]])
    keep <- is.finite(es) & is.finite(se) & se > 0
    es <- es[keep]; se <- se[keep]
    if (length(es) < 1L) return(NA_real_)
    ma <- .rmaMeta(es, se, method = metaMethod)
    ma$mean / ma$se
  }, numeric(1))
  data.frame(
    condition = sub("_deviation$", "", sub("^mean_contrast_", "", effCols)),
    zScore = as.numeric(scores), stringsAsFactors = FALSE)
}

#' Feature score from the fraction of significant deviation contrasts
#'
#' For each condition's deviation contrast, the proportion of variants whose
#' contrast p-value falls below \code{pCutoff} -- the "n-significant" feature
#' score of \code{mash_posterior.ipynb}. No meta-analysis is involved.
#'
#' @param contrastResult A contrast table from \code{\link{mashPosteriorContrast}}
#'   carrying \code{p_contrast_*_deviation} columns.
#' @param pCutoff Significance threshold (default 1e-5).
#' @return A \code{data.frame} with \code{condition} and \code{ratio}
#'   (\eqn{n_{sig} / n}).
#' @export
nSignificantScore <- function(contrastResult, pCutoff = 1e-5) {
  cr <- as.data.frame(contrastResult, stringsAsFactors = FALSE)
  pCols <- grep("p_contrast_.*deviation", names(cr), value = TRUE)
  if (length(pCols) == 0L)
    return(data.frame(condition = character(0), ratio = numeric(0),
                      stringsAsFactors = FALSE))
  ratios <- vapply(pCols, function(pc) {
    p <- as.numeric(cr[[pc]])
    nTot <- sum(!is.na(p))
    if (nTot == 0L) NA_real_ else sum(p < pCutoff, na.rm = TRUE) / nTot
  }, numeric(1))
  data.frame(
    condition = sub("_deviation$", "", sub("^p_contrast_", "", pCols)),
    ratio = as.numeric(ratios), stringsAsFactors = FALSE)
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
#'   \code{pip} and \code{variants} columns (credible-set index 0 = not in a CS).
#' @param contrastResults A contrast table from
#'   \code{\link{mashPosteriorContrast}} with variant row names.
#' @param condition Condition label whose deviation p-value column is used; if
#'   that column is absent and exactly one pairwise contrast exists, the
#'   pairwise p-value is used instead.
#' @return A single numeric score, or \code{NA} when nothing overlaps.
#' @export
scoreFromCs <- function(fineMapping, contrastResults, condition) {
  css <- setdiff(unique(fineMapping$cs_order), 0)
  if (length(css) == 0L) return(NA_real_)
  leadRows <- do.call(rbind, lapply(css, function(cs) {
    tmp <- fineMapping[fineMapping$cs_order == cs, , drop = FALSE]
    tmp[tmp$pip == max(tmp$pip), , drop = FALSE]
  }))
  cr <- as.data.frame(contrastResults, stringsAsFactors = FALSE)
  snps <- intersect(rownames(cr), leadRows$variants)
  if (length(snps) == 0L) return(NA_real_)

  pDevCol <- grep(paste0("p_contrast_", condition, "_deviation"),
                  names(cr), value = TRUE)
  if (length(pDevCol) > 0L) {
    contrastP <- cr[snps, pDevCol[1L], drop = FALSE]
  } else {
    pvCols <- grep("p_contrast.*_vs_", names(cr), value = TRUE)
    if (length(pvCols) != 1L) return(NA_real_)
    contrastP <- cr[snps, pvCols, drop = FALSE]
  }
  maxSnp <- rownames(contrastP)[which(contrastP[, 1L] ==
                                       max(contrastP[, 1L], na.rm = TRUE))]
  meanVs <- as.numeric(as.matrix(cr[maxSnp,
              grep("mean_contrast.*_vs_", names(cr)), drop = FALSE]))
  seVs <- as.numeric(as.matrix(cr[maxSnp,
              grep("se_contrast.*_vs_", names(cr)), drop = FALSE]))
  max(abs(meanVs / seVs), na.rm = TRUE)
}
