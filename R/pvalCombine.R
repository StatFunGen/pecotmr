#' @title P-value Combination Methods
#' @description Functions for combining p-values across multiple tests or
#'   methods: Cauchy combination (ACAT), harmonic mean p-value (HMP),
#'   poolr, GBJ, and aSPU methods. Also includes the null correlation
#'   matrix for TWAS z-scores used by correlation-adjusted methods.
#' @name pecotmr-pval-combine
#' @keywords internal
#' @importFrom magrittr %>%
NULL

#' @export
waldTestPval <- function(beta, se, n) {
  # Calculate the t statistic
  tValue <- beta / se
  # Degrees of freedom
  df <- n - 2
  # Calculate two-tailed p-value
  pValue <- 2 * pt(-abs(tValue), df = df, lower.tail = TRUE)

  return(pValue)
}

pvalAcat <- function(pvals, naRm = TRUE) {
  # ACAT (Aggregated Cauchy Association Test) — Liu & Xie (2020).
  # T = mean(tan(pi * (0.5 - p_i)))
  #
  # Robustness handling merged from the former pvalCauchy implementation:
  #   - naRm: optionally drop NAs
  #   - clip p > 0.99 to 0.99 to bound the contribution of near-1 p-values
  #   - small-p asymptotic: tan(pi*(0.5 - p)) ~ 1/(pi*p) for p < 1e-15 to
  #     avoid Inf from floating-point precision loss in pi*0.5
  #   - large-stat asymptotic: when the mean Cauchy variate is > 1e15 the
  #     CDF tail collapses to (1/T) / pi (Cauchy survival expansion)
  if (naRm) pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0L) return(NA_real_)
  if (length(pvals) == 1L) return(pvals[[1]])
  pvals <- pmin(pvals, 0.99)
  cauchyVals <- ifelse(pvals < 1e-15,
                       1 / (pvals * pi),
                       tan(pi * (0.5 - pvals)))
  stat <- mean(cauchyVals)
  if (!is.finite(stat)) return(NA_real_)
  if (stat > 1e15) return((1 / stat) / pi)
  pcauchy(stat, lower.tail = FALSE)
}

pvalHmp <- function(pvals) {
  # Make sure harmonicmeanp is installed
  if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
    stop("To use this function, please install harmonicmeanp: https://cran.r-project.org/web/packages/harmonicmeanp/index.html")
  }
  # https://search.r-project.org/CRAN/refmans/harmonicmeanp/html/pLandau.html
  L <- length(pvals)
  HMP <- L / sum(pvals^-1)

  LOC_L1 <- 0.874367040387922
  SCALE <- 1.5707963267949

  return(harmonicmeanp::pLandau(1 / HMP, mu = log(L) + LOC_L1, sigma = SCALE, lower.tail = FALSE))
}

pvalPoolr <- function(pvals, method, R) {
  if (!requireNamespace("poolr", quietly = TRUE)) {
    stop("To use this method, please install poolr: install.packages('poolr')")
  }
  fn <- switch(method,
    fisher = poolr::fisher,
    stouffer = poolr::stouffer,
    invchisq = poolr::invchisq,
    stop(sprintf("Unknown poolr method: '%s'", method))
  )
  fn(pvals, adjust = "generalized", R = R)$p
}

pvalGbj <- function(zScores, R, method) {
  if (!requireNamespace("GBJ", quietly = TRUE)) {
    stop("To use this method, please install GBJ: install.packages('GBJ')")
  }
  result <- switch(method,
    gbj = GBJ::GBJ(test_stats = zScores, cor_mat = R),
    bj = GBJ::BJ(test_stats = zScores, cor_mat = R),
    hc = GBJ::HC(test_stats = zScores, cor_mat = R),
    ghc = GBJ::GHC(test_stats = zScores, cor_mat = R),
    minp = GBJ::minP(test_stats = zScores, cor_mat = R),
    gbj_omni = GBJ::OMNI_ss(test_stats = zScores, cor_mat = R),
    stop(sprintf("Unknown GBJ method: '%s'", method))
  )
  pvalName <- switch(method,
    gbj = "GBJ_pvalue",
    bj = "BJ_pvalue",
    hc = "HC_pvalue",
    ghc = "GHC_pvalue",
    minp = "minP_pvalue",
    gbj_omni = "OMNI_pvalue"
  )
  result[[pvalName]]
}

pvalAspu <- function(zScores = NULL, pvals = NULL, R, method) {
  if (!requireNamespace("aSPU", quietly = TRUE)) {
    stop("To use this method, please install aSPU: install.packages('aSPU')")
  }
  switch(method,
    aspu = {
      result <- aSPU::aSPUs(Zs = zScores, corSNP = R)
      result$pvs["aSPUs"]
    },
    gates = {
      result <- aSPU::GATES2(ldmatrix = R, p = pvals)
      result[["Pg"]]
    },
    stop(sprintf("Unknown aSPU method: '%s'", method))
  )
}

# =============================================================================
# combinePValues — unified dispatcher
# =============================================================================

# Methods that require an R correlation matrix.
.combinePvalMethodsNeedingR <- c(
  "fisher", "stouffer", "invchisq",
  "gbj", "bj", "hc", "ghc", "minp", "gbj_omni",
  "aspu", "gates")

# Methods that require signed zScores (and therefore cannot be derived from
# pvals alone). All other methods can work from pvals.
.combinePvalMethodsNeedingZ <- c(
  "gbj", "bj", "hc", "ghc", "minp", "gbj_omni",
  "aspu")

.combinePvalKnownMethods <- c(
  "acat", "hmp", "bonferroni",
  "fisher", "stouffer", "invchisq",
  "gbj", "bj", "hc", "ghc", "minp", "gbj_omni",
  "aspu", "gates")

# Internal: align an R correlation matrix to a target order. If R has
# rownames/colnames, reorder to match `targetNames`; require every target
# name to be present. If R is unnamed, only length check.
.combinePvalAlignR <- function(R, targetNames) {
  if (is.null(R)) return(NULL)
  if (!is.matrix(R)) stop("`R` must be a matrix.")
  if (nrow(R) != ncol(R)) stop("`R` must be square.")
  rNames <- rownames(R)
  cNames <- colnames(R)
  hasNames <- !is.null(rNames) && !is.null(cNames)
  if (hasNames) {
    if (!identical(rNames, cNames))
      stop("`R` rownames and colnames must be identical.")
    missing <- setdiff(targetNames, rNames)
    if (length(missing) > 0L)
      stop("`R` is missing entries for: ",
           paste(utils::head(missing, 5), collapse = ", "),
           if (length(missing) > 5L) sprintf(" (and %d more)",
                                              length(missing) - 5L)
           else "")
    R <- R[targetNames, targetNames, drop = FALSE]
  } else {
    if (nrow(R) != length(targetNames))
      stop("Unnamed `R` must have nrow = length(pvals); got nrow(R) = ",
           nrow(R), ", length(pvals) = ", length(targetNames), ".")
  }
  R
}

# Internal: compute one method's combined p-value. Assumes inputs have
# already been validity-filtered (positive, finite, < 1) and aligned with R.
.combinePvalSingle <- function(method, pvals, zScores, R) {
  switch(method,
    acat       = pvalAcat(pvals),
    hmp        = pvalHmp(pvals),
    bonferroni = min(length(pvals) * min(pvals), 1.0),
    fisher     = ,
    stouffer   = ,
    invchisq   = pvalPoolr(pvals, method = method, R = R),
    gbj        = ,
    bj         = ,
    hc         = ,
    ghc        = ,
    minp       = ,
    gbj_omni   = pvalGbj(zScores, R = R, method = method),
    aspu       = pvalAspu(zScores = zScores, R = R, method = "aspu"),
    gates      = pvalAspu(pvals = pvals, R = R, method = "gates"),
    stop(sprintf("Unknown combination method: '%s'", method))
  )
}

#' Combine P-values via Any of a Menu of Methods
#'
#' Unified dispatcher for combining a vector of p-values (and/or z-scores)
#' into a single combined p-value. Supports independent-test methods
#' (ACAT, HMP, Bonferroni) and correlation-adjusted methods (Fisher /
#' Stouffer / inverse-chi-square via \code{poolr}; GBJ / BJ / HC / GHC /
#' minP / GBJ-omnibus via \code{GBJ}; aSPU / GATES via \code{aSPU}).
#' Multiple methods may be requested in a single call; the function
#' returns a per-method result list keyed by method name.
#'
#' Either \code{pvals} or \code{zScores} may be supplied. If only
#' \code{zScores} is given, two-sided p-values are derived as
#' \code{p = 2 * (1 - pnorm(|z|))}. Methods that require signed
#' \code{zScores} (\code{gbj}, \code{bj}, \code{hc}, \code{ghc},
#' \code{minp}, \code{gbj_omni}, \code{aspu}) cannot be derived from
#' p-values alone and error if \code{zScores} is missing.
#'
#' Methods that require a correlation matrix \code{R}
#' (\code{fisher}, \code{stouffer}, \code{invchisq}, \code{gbj},
#' \code{bj}, \code{hc}, \code{ghc}, \code{minp}, \code{gbj_omni},
#' \code{aspu}, \code{gates}) error if \code{R} is missing. Methods that
#' do not use \code{R} silently ignore it. When \code{R} is named, it is
#' realigned to match the order of \code{pvals} / \code{zScores}, with a
#' hard error if any entry is missing from \code{R}'s names.
#'
#' An internal validity filter drops entries where
#' \code{!is.finite(pvals) | pvals <= 0 | pvals >= 1}; a warning is
#' emitted when any are dropped.
#'
#' @param pvals Optional numeric vector of p-values. Required for
#'   methods that work on p-values; derivable from \code{zScores}.
#' @param zScores Optional numeric vector of signed z-scores.
#' @param methods Character vector of combination method names; see
#'   above for the menu (lowercase).
#' @param R Optional correlation matrix aligned to \code{pvals} /
#'   \code{zScores}. Required for the correlation-adjusted methods.
#' @param naRm Logical; if \code{TRUE} (default), drop NA p-values
#'   before combination.
#' @return A list with two elements:
#'   \describe{
#'     \item{input}{Summary of the call: \code{nPvalsIn},
#'       \code{nZScoresIn}, \code{nValid}, and the aligned
#'       \code{Raligned} matrix (or \code{NULL}).}
#'     \item{results}{Named list keyed by method, each element a list
#'       with \code{method} and \code{pval}.}
#'   }
#' @export
combinePValues <- function(pvals = NULL, zScores = NULL,
                           methods, R = NULL, naRm = TRUE) {
  if (missing(methods) || length(methods) == 0L)
    stop("`methods` is required (one or more of: ",
         paste(.combinePvalKnownMethods, collapse = ", "), ").")
  methods <- as.character(methods)
  unknown <- setdiff(methods, .combinePvalKnownMethods)
  if (length(unknown) > 0L)
    stop("Unknown method(s): ", paste(unknown, collapse = ", "),
         ". Known: ", paste(.combinePvalKnownMethods, collapse = ", "))

  nPvalsIn   <- if (is.null(pvals))   0L else length(pvals)
  nZScoresIn <- if (is.null(zScores)) 0L else length(zScores)

  # Method-level prerequisite checks.
  needZ <- intersect(methods, .combinePvalMethodsNeedingZ)
  if (length(needZ) > 0L && is.null(zScores))
    stop("Method(s) ", paste(needZ, collapse = ", "),
         " require `zScores`; supplied input only has pvals. Signed ",
         "z-scores cannot be recovered from p-values alone.")
  needR <- intersect(methods, .combinePvalMethodsNeedingR)
  if (length(needR) > 0L && is.null(R))
    stop("Method(s) ", paste(needR, collapse = ", "),
         " require an `R` correlation matrix; got NULL.")

  # Derive missing input where possible: pvals from zScores via two-sided.
  if (is.null(pvals) && !is.null(zScores))
    pvals <- 2 * stats::pnorm(-abs(as.numeric(zScores)))
  if (is.null(pvals))
    stop("Either `pvals` or `zScores` must be supplied.")

  if (!is.null(zScores) && length(zScores) != length(pvals))
    stop("`pvals` and `zScores` must have the same length when both supplied.")

  # Internal validity filter + optional NA drop.
  naMask <- is.na(pvals) | if (is.null(zScores)) FALSE else is.na(zScores)
  invalidMask <- !naMask & (!is.finite(pvals) | pvals <= 0 | pvals >= 1)
  dropMask <- (naRm & naMask) | invalidMask
  if (any(dropMask)) {
    warning(sprintf(
      "combinePValues: dropped %d entry/entries (%d NA, %d invalid).",
      sum(dropMask), sum(naMask & dropMask), sum(invalidMask)))
  }
  keep <- !dropMask
  pvalsK   <- pvals[keep]
  zScoresK <- if (is.null(zScores)) NULL else zScores[keep]

  # Build a stable target-name vector. If pvals is named, use those;
  # otherwise use the corresponding R names when available;
  # otherwise positional integer labels.
  targetNames <- if (!is.null(names(pvals))) names(pvals)[keep]
                 else if (!is.null(zScores) && !is.null(names(zScores)))
                   names(zScores)[keep]
                 else if (!is.null(R) && !is.null(rownames(R)) &&
                          nrow(R) == length(pvals))
                   rownames(R)[keep]
                 else as.character(seq_along(pvalsK))

  Raligned <- .combinePvalAlignR(R, targetNames)

  nValid <- length(pvalsK)
  if (nValid < 1L) {
    perMethod <- lapply(methods, function(m) {
      list(method = m, pval = NA_real_)
    })
    names(perMethod) <- methods
    return(list(
      input   = list(nPvalsIn = nPvalsIn, nZScoresIn = nZScoresIn,
                     nValid = nValid, Raligned = Raligned),
      results = perMethod))
  }

  perMethod <- lapply(methods, function(m) {
    p <- tryCatch(
      .combinePvalSingle(m, pvals = pvalsK, zScores = zScoresK,
                         R = Raligned),
      error = function(e) {
        warning(sprintf("combinePValues: method '%s' failed: %s",
                        m, conditionMessage(e)))
        NA_real_
      })
    list(method = m, pval = as.numeric(p))
  })
  names(perMethod) <- methods

  list(
    input   = list(nPvalsIn = nPvalsIn, nZScoresIn = nZScoresIn,
                   nValid = nValid, Raligned = Raligned),
    results = perMethod)
}
