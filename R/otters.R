#' Train eQTL weights using multiple RSS methods (OTTERS Stage I)
#'
#' Implements the training stage of the OTTERS framework (Omnibus Transcriptome
#' Test using Expression Reference Summary data, Zhang et al. 2024). Trains
#' eQTL effect size weights for a gene region using multiple summary-statistics-based
#' methods in parallel, enabling downstream omnibus TWAS testing.
#'
#' Methods are dispatched dynamically via \code{do.call(paste0(method, "_weights"), ...)},
#' so any function following the \code{*_weights(stat, LD, ...)} convention can be used
#' (e.g., \code{lassosum_rss_weights}, \code{prs_cs_weights}, \code{sdpr_weights},
#' \code{mr_ash_rss_weights}).
#'
#' P+T (pruning and thresholding) is handled internally: for each threshold, SNPs with
#' eQTL p-value below the threshold are selected, and their marginal z-scores (scaled
#' to correlation units: \code{z / sqrt(n)}) are used as weights.
#'
#' @param sumstats A data.frame of eQTL summary statistics. Must contain column \code{z}
#'   (z-scores). If \code{z} is absent but \code{beta} and \code{se} are present,
#'   z-scores are computed as \code{beta / se}.
#' @param LD LD correlation matrix R for the gene region (single matrix, not a list).
#'   Should have row/column names matching variant identifiers if variant alignment
#'   is desired.
#' @param n eQTL study sample size (scalar).
#' @param methods Named list of RSS methods and their extra arguments. Each element
#'   name must correspond to a \code{*_weights} function in pecotmr (without the
#'   \code{_weights} suffix). Defaults match the original OTTERS pipeline
#'   (Zhang et al. 2024):
#'   \itemize{
#'     \item \code{lassosum_rss}: s grid = c(0.2, 0.5, 0.9, 1.0), lambda from
#'       0.0001 to 0.1 (20 values on log scale)
#'     \item \code{prs_cs}: phi = 1e-4 (fixed, not learned), 1000 iterations,
#'       500 burn-in, thin = 5
#'     \item \code{sdpr}: 1000 iterations, 200 burn-in, thin = 1 (no thinning)
#'   }
#'   To add learners (e.g., \code{mr_ash_rss}), simply append to this list.
#' @param p_thresholds Numeric vector of p-value thresholds for P+T. Set to
#'   \code{NULL} to skip P+T. Default: \code{c(0.001, 0.05)}.
#' @param check_ld_method LD quality check method passed to \code{\link{check_ld}}.
#'   Default \code{"eigenfix"} sets negative eigenvalues to zero (required for
#'   PRS-CS Cholesky, matching OTTERS' SVD-based PD forcing). Set to \code{NULL}
#'   to skip checking.
#'
#' @return A named list of weight vectors (one per method). Each vector has length
#'   equal to \code{nrow(sumstats)}. P+T results are named \code{PT_<threshold>}.
#'
#' @examples
#' set.seed(42)
#' n <- 500; p <- 20
#' z <- rnorm(p, sd = 2)
#' R <- diag(p)
#' sumstats <- data.frame(z = z)
#' weights <- otters_weights(sumstats, R, n,
#'   methods = list(lassosum_rss = list()),
#'   p_thresholds = c(0.05))
#'
#' @export
otters_weights <- function(sumstats, LD, n,
                           methods = list(
                             lassosum_rss = list(),
                             prs_cs = list(phi = 1e-4,
                                           n_iter = 1000, n_burnin = 500, thin = 5),
                             sdpr = list(iter = 1000, burn = 200, thin = 1, verbose = FALSE)
                           ),
                           p_thresholds = c(0.001, 0.05),
                           check_ld_method = "eigenfix") {
  # Check and optionally repair LD matrix quality
  # PRS-CS requires positive-definite LD for Cholesky; OTTERS forces PD via SVD.
  # Default "eigenfix" sets negative eigenvalues to 0 (susieR approach).
  # Set to NULL to skip (e.g., if LD is known to be clean).
  if (!is.null(check_ld_method)) {
    ld_check <- check_ld(LD, method = check_ld_method)
    if (ld_check$method_applied != "none") {
      message(sprintf("check_ld: repaired LD via '%s' (min eigenvalue was %.2e, %d negative).",
                      ld_check$method_applied, ld_check$min_eigenvalue, ld_check$n_negative))
    }
    LD <- ld_check$R
  }

  # Compute z-scores if not present
  if (is.null(sumstats$z)) {
    if (!is.null(sumstats$beta) && !is.null(sumstats$se)) {
      sumstats$z <- sumstats$beta / sumstats$se
    } else {
      stop("sumstats must have 'z' or ('beta' and 'se') columns.")
    }
  }

  p <- nrow(sumstats)
  z <- sumstats$z

  # Build stat object for _weights() convention
  b <- z / sqrt(n)
  stat <- list(b = b, cor = b, z = z, n = rep(n, p))

  results <- list()

  # --- P+T (Pruning and Thresholding) ---
  if (!is.null(p_thresholds)) {
    pvals <- pchisq(z^2, df = 1, lower.tail = FALSE)
    for (thr in p_thresholds) {
      selected <- pvals < thr
      # Weights = clamped marginal correlation (stat$b) for selected SNPs
      w <- ifelse(selected, stat$b, 0)
      results[[paste0("PT_", thr)]] <- w
    }
  }

  # --- RSS methods ---
  for (method_name in names(methods)) {
    fn_name <- paste0(method_name, "_weights")
    if (!exists(fn_name, mode = "function")) {
      warning(sprintf("Method '%s' not found (looking for function '%s'). Skipping.",
                      method_name, fn_name))
      next
    }
    tryCatch({
      w <- do.call(fn_name, c(list(stat = stat, LD = LD), methods[[method_name]]))
      results[[method_name]] <- as.numeric(w)
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", method_name, e$message))
      results[[method_name]] <<- rep(0, p)
    })
  }

  results
}


#' TWAS association testing with omnibus combination (OTTERS Stage II)
#'
#' Computes per-method TWAS z-scores using the FUSION formula and combines
#' p-values across methods via ACAT (Aggregated Cauchy Association Test) or
#' HMP (Harmonic Mean P-value).
#'
#' The FUSION TWAS statistic (Gusev et al. 2016) is:
#' \deqn{Z_{TWAS} = \frac{w^T z}{\sqrt{w^T R w}}}
#' where \eqn{w} are eQTL weights, \eqn{z} are GWAS z-scores, and \eqn{R}
#' is the LD correlation matrix.
#'
#' @param weights Named list of weight vectors (output from \code{\link{otters_weights}}
#'   or any named list of numeric vectors).
#' @param gwas_z Numeric vector of GWAS z-scores, same length and order as the
#'   weights vectors. Must be aligned to the same variants and allele orientation
#'   as the weights and LD matrix. Use \code{\link{allele_qc}} or
#'   \code{\link{rss_basic_qc}} for harmonization before calling this function.
#' @param LD LD correlation matrix R, aligned to the same variants as weights
#'   and gwas_z.
#' @param combine_method Method to combine p-values across methods.
#'   Correlation-free (valid under arbitrary dependence):
#'   \code{"acat"} (default), \code{"hmp"}.
#'   Correlation-adjusted via poolr (generalized multivariate theory):
#'   \code{"fisher"} (Brown's method), \code{"stouffer"} (Strube's method),
#'   \code{"invchisq"}.
#'   Set-based tests via GBJ (uses TWAS z-scores and inter-method correlation):
#'   \code{"gbj"}, \code{"bj"}, \code{"hc"}, \code{"ghc"}, \code{"minp"},
#'   \code{"gbj_omni"}.
#'   Adaptive and Simes-type tests via aSPU:
#'   \code{"aspu"} (adaptive sum of powered scores),
#'   \code{"gates"} (extended Simes / GATES).
#'   The poolr, GBJ, and aSPU methods automatically compute the inter-method
#'   TWAS z-score correlation from the weight vectors and LD matrix.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{method}{Method name (per-method rows plus a combined row).}
#'   \item{twas_z}{TWAS z-score (\code{NA} for combined row).}
#'   \item{twas_pval}{TWAS p-value.}
#'   \item{n_snps}{Number of non-zero weight SNPs used.}
#' }
#'
#' @examples
#' set.seed(42)
#' p <- 20
#' gwas_z <- rnorm(p)
#' R <- diag(p)
#' weights <- list(method1 = rnorm(p, sd = 0.01), method2 = rnorm(p, sd = 0.01))
#' otters_association(weights, gwas_z, R)
#'
#' @export
otters_association <- function(weights, gwas_z, LD,
                               combine_method = c("acat", "hmp",
                                                   "fisher", "stouffer", "invchisq",
                                                   "gbj", "bj", "hc", "ghc",
                                                   "minp", "gbj_omni",
                                                   "aspu", "gates")) {
  combine_method <- match.arg(combine_method)

  # Validate dimensions
  p <- length(gwas_z)
  if (nrow(LD) != p || ncol(LD) != p) {
    stop(sprintf("LD dimensions (%d x %d) do not match gwas_z length (%d).",
                 nrow(LD), ncol(LD), p))
  }
  for (nm in names(weights)) {
    if (length(weights[[nm]]) != p) {
      stop(sprintf("Weight vector '%s' has length %d but gwas_z has length %d.",
                   nm, length(weights[[nm]]), p))
    }
  }

  results <- data.frame(
    method = character(),
    twas_z = numeric(),
    twas_pval = numeric(),
    n_snps = integer(),
    stringsAsFactors = FALSE
  )

  valid_pvals <- c()
  valid_zscores <- c()
  valid_weights <- list()

  for (method_name in names(weights)) {
    w <- weights[[method_name]]

    # Skip all-zero weights
    if (all(w == 0)) {
      results <- rbind(results, data.frame(
        method = method_name, twas_z = NA_real_,
        twas_pval = NA_real_, n_snps = 0L,
        stringsAsFactors = FALSE
      ))
      next
    }

    # Use non-zero SNPs
    nz <- which(w != 0)
    n_snps <- length(nz)

    # Compute TWAS z-score via twas_z()
    res <- twas_z(weights = w[nz], z = gwas_z[nz], R = LD[nz, nz, drop = FALSE])

    z_val <- as.numeric(res$z)
    p_val <- as.numeric(res$pval)

    results <- rbind(results, data.frame(
      method = method_name, twas_z = z_val,
      twas_pval = p_val, n_snps = n_snps,
      stringsAsFactors = FALSE
    ))

    if (!is.na(p_val) && is.finite(p_val) && p_val > 0 && p_val < 1) {
      valid_pvals <- c(valid_pvals, p_val)
      valid_zscores <- c(valid_zscores, z_val)
      valid_weights[[length(valid_weights) + 1]] <- w
    }
  }

  # Combine p-values across methods
  if (length(valid_pvals) >= 2) {
    poolr_methods <- c("fisher", "stouffer", "invchisq")
    gbj_methods <- c("gbj", "bj", "hc", "ghc", "minp", "gbj_omni")
    aspu_methods <- c("aspu", "gates")
    needs_cor <- combine_method %in% c(poolr_methods, gbj_methods, aspu_methods)

    method_cor <- NULL
    if (needs_cor) {
      method_cor <- twas_method_cor(valid_weights, LD)
    }

    combined_pval <- if (combine_method == "acat") {
      pval_acat(valid_pvals)
    } else if (combine_method == "hmp") {
      pval_hmp(valid_pvals)
    } else if (combine_method %in% poolr_methods) {
      pval_poolr(valid_pvals, combine_method, R = method_cor)
    } else if (combine_method %in% gbj_methods) {
      pval_gbj(valid_zscores, method_cor, combine_method)
    } else if (combine_method %in% aspu_methods) {
      pval_aspu(z_scores = valid_zscores, pvals = valid_pvals,
                R = method_cor, method = combine_method)
    }

    results <- rbind(results, data.frame(
      method = paste0(toupper(combine_method), "_combined"),
      twas_z = NA_real_,
      twas_pval = as.numeric(combined_pval),
      n_snps = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  results
}
