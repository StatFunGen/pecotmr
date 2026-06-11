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

pvalAcat <- function(pvals) {
  if (length(pvals) == 1) {
    return(pvals[1])
  }
  # ACAT statistic: T = mean(tan(pi*(0.5 - p_i)))
  # Liu & Xie (2020) "Cauchy combination test"
  #
  # For very small p, tan(pi*(0.5-p)) overflows due to floating-point
  # precision loss in pi*0.5. Use the asymptotic approximation
  # tan(pi*(0.5-p)) ~ 1/(pi*p) for p < 1e-15 to avoid Inf/NaN.
  cauchyVals <- ifelse(pvals < 1e-15,
                       1 / (pvals * pi),
                       tan(pi * (0.5 - pvals)))
  stat <- mean(cauchyVals)
  return(pcauchy(stat, lower.tail = FALSE))
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

pvalGlobal <- function(pvals, combMethod = "HMP", naive = FALSE) {
  # assuming sstats has tissues as columns and rows as pvals
  minPval <- min(pvals)
  nTotalTests <- pvals %>%
    unique() %>%
    length() # There should be one unique pval per tissue
  globalPval <- if (combMethod == "HMP") pvalHmp(pvals) else pvalAcat(pvals) # pval vector
  naivePval <- min(nTotalTests * minPval, 1.0)
  return(if (naive) naivePval else globalPval) # globalPval and naivePval
}

pvalCauchy <- function(p, na.rm = TRUE) {
  if (na.rm) {
    if (sum(is.na(p))) {
      p <- p[!is.na(p)]
    }
  }
  p[p > 0.99] <- 0.99
  is.small <- (p < 1e-16) & !is.na(p)
  is.regular <- (p >= 1e-16) & !is.na(p)
  temp <- rep(NA, length(p))
  temp[is.small] <- 1 / p[is.small] / pi
  temp[is.regular] <- as.numeric(tan((0.5 - p[is.regular]) * pi))

  cct.stat <- mean(temp, na.rm = TRUE)
  if (is.na(cct.stat)) {
    return(NA)
  }
  if (cct.stat > 1e+15) {
    return((1 / cct.stat) / pi)
  } else {
    return(1 - pcauchy(cct.stat))
  }
}

# Compute null correlation matrix of TWAS z-scores across methods.
# z_k = w_k' z / sqrt(w_k' R w_k), so under H0 (z ~ N(0, R)):
# cor(z_k, z_j) = w_k' R w_j / sqrt(w_k' R w_k * w_j' R w_j)
twasMethodCor <- function(weightsList, LD) {
  K <- length(weightsList)
  corMat <- diag(K)
  Rw <- lapply(weightsList, function(w) as.numeric(LD %*% w))
  varW <- vapply(seq_len(K), function(i) sum(weightsList[[i]] * Rw[[i]]), numeric(1))
  for (i in seq_len(K - 1)) {
    for (j in (i + 1):K) {
      if (varW[i] > 0 && varW[j] > 0) {
        rho <- sum(weightsList[[i]] * Rw[[j]]) / sqrt(varW[i] * varW[j])
        corMat[i, j] <- rho
        corMat[j, i] <- rho
      }
    }
  }
  corMat
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
