#' @title P-value Combination Methods
#' @description Functions for combining p-values across multiple tests or
#'   methods: Cauchy combination (ACAT), harmonic mean p-value (HMP),
#'   poolr, GBJ, and aSPU methods. Also includes the null correlation
#'   matrix for TWAS z-scores used by correlation-adjusted methods.
#' @importFrom magrittr %>%
NULL

#' @export
wald_test_pval <- function(beta, se, n) {
  # Calculate the t statistic
  t_value <- beta / se
  # Degrees of freedom
  df <- n - 2
  # Calculate two-tailed p-value
  p_value <- 2 * pt(-abs(t_value), df = df, lower.tail = TRUE)

  return(p_value)
}

pval_acat <- function(pvals) {
  if (length(pvals) == 1) {
    return(pvals[1])
  }
  # ACAT statistic: T = mean(tan(pi*(0.5 - p_i)))
  # Liu & Xie (2020) "Cauchy combination test"
  #
  # For very small p, tan(pi*(0.5-p)) overflows due to floating-point
  # precision loss in pi*0.5. Use the asymptotic approximation
  # tan(pi*(0.5-p)) ~ 1/(pi*p) for p < 1e-15 to avoid Inf/NaN.
  cauchy_vals <- ifelse(pvals < 1e-15,
                        1 / (pvals * pi),
                        tan(pi * (0.5 - pvals)))
  stat <- mean(cauchy_vals)
  return(pcauchy(stat, lower.tail = FALSE))
}

pval_hmp <- function(pvals) {
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

pval_global <- function(pvals, comb_method = "HMP", naive = FALSE) {
  # assuming sstats has tissues as columns and rows as pvals
  min_pval <- min(pvals)
  n_total_tests <- pvals %>%
    unique() %>%
    length() # There should be one unique pval per tissue
  global_pval <- if (comb_method == "HMP") pval_hmp(pvals) else pval_acat(pvals) # pval vector
  naive_pval <- min(n_total_tests * min_pval, 1.0)
  return(if (naive) naive_pval else global_pval) # global_pval and naive_pval
}

pval_cauchy <- function(p, na.rm = TRUE) {
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
twas_method_cor <- function(weights_list, LD) {
  K <- length(weights_list)
  cor_mat <- diag(K)
  Rw <- lapply(weights_list, function(w) as.numeric(LD %*% w))
  var_w <- vapply(seq_len(K), function(i) sum(weights_list[[i]] * Rw[[i]]), numeric(1))
  for (i in seq_len(K - 1)) {
    for (j in (i + 1):K) {
      if (var_w[i] > 0 && var_w[j] > 0) {
        rho <- sum(weights_list[[i]] * Rw[[j]]) / sqrt(var_w[i] * var_w[j])
        cor_mat[i, j] <- rho
        cor_mat[j, i] <- rho
      }
    }
  }
  cor_mat
}

pval_poolr <- function(pvals, method, R) {
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

pval_gbj <- function(z_scores, R, method) {
  if (!requireNamespace("GBJ", quietly = TRUE)) {
    stop("To use this method, please install GBJ: install.packages('GBJ')")
  }
  result <- switch(method,
    gbj = GBJ::GBJ(test_stats = z_scores, cor_mat = R),
    bj = GBJ::BJ(test_stats = z_scores, cor_mat = R),
    hc = GBJ::HC(test_stats = z_scores, cor_mat = R),
    ghc = GBJ::GHC(test_stats = z_scores, cor_mat = R),
    minp = GBJ::minP(test_stats = z_scores, cor_mat = R),
    gbj_omni = GBJ::OMNI_ss(test_stats = z_scores, cor_mat = R),
    stop(sprintf("Unknown GBJ method: '%s'", method))
  )
  pval_name <- switch(method,
    gbj = "GBJ_pvalue",
    bj = "BJ_pvalue",
    hc = "HC_pvalue",
    ghc = "GHC_pvalue",
    minp = "minP_pvalue",
    gbj_omni = "OMNI_pvalue"
  )
  result[[pval_name]]
}

pval_aspu <- function(z_scores = NULL, pvals = NULL, R, method) {
  if (!requireNamespace("aSPU", quietly = TRUE)) {
    stop("To use this method, please install aSPU: install.packages('aSPU')")
  }
  switch(method,
    aspu = {
      result <- aSPU::aSPUs(Zs = z_scores, corSNP = R)
      result$pvs["aSPUs"]
    },
    gates = {
      result <- aSPU::GATES2(ldmatrix = R, p = pvals)
      result[["Pg"]]
    },
    stop(sprintf("Unknown aSPU method: '%s'", method))
  )
}
