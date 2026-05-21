#' Convert Log Bayes Factors to Single Effects PIP
#'
#' This function converts log Bayes factors (LBF) to alpha values, optionally
#' using prior weights. It handles numerical stability by adjusting with the
#' maximum LBF value.
#'
#' @param lbf Numeric vector of log Bayes factors.
#' @param prior_weights Optional numeric vector of prior weights for each element in lbf.
#' @return A named numeric vector of alpha values corresponding to the input LBF.
#' @examples
#' lbf <- c(-0.5, 1.2, 0.3)
#' alpha <- lbf_to_alpha_vector(lbf)
#' print(alpha)
#' @noRd
lbf_to_alpha_vector <- function(lbf, prior_weights = NULL) {
  if (is.null(prior_weights)) prior_weights <- rep(1 / length(lbf), length(lbf))
  maxlbf <- max(lbf)

  # If maxlbf is 0, return a vector of zeros
  if (maxlbf == 0) {
    return(setNames(rep(0, length(lbf)), names(lbf)))
  }

  # w is proportional to BF, subtract max for numerical stability
  w <- exp(lbf - maxlbf)

  # Posterior prob for each SNP
  w_weighted <- w * prior_weights
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w

  return(alpha)
}

#' Applies the 'lbf_to_alpha_vector' function row-wise to a matrix of log Bayes factors
#' to convert them to Single Effect PIP values.
#'
#' @param lbf Matrix of log Bayes factors.
#' @return A matrix of alpha values with the same dimensions as the input LBF matrix.
#' @examples
#' lbf_matrix <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
#' alpha_matrix <- lbf_to_alpha(lbf_matrix)
#' print(alpha_matrix)
#' @export
lbf_to_alpha <- function(lbf) {
  alpha_matrix <- t(apply(as.matrix(lbf), 1, lbf_to_alpha_vector))
  if (ncol(lbf) == 1) alpha_matrix <- matrix(alpha_matrix, ncol = 1, dimnames = list(NULL, colnames(lbf)))
  return(alpha_matrix)
}

format_pip_column <- function(method) {
  paste0("pip_", method)
}

resolve_pip_column <- function(top_loci, method = NULL) {
  if (is.null(top_loci) || nrow(top_loci) == 0) return(NULL)
  if (!is.null(method)) {
    pip_col <- format_pip_column(method)
    if (pip_col %in% names(top_loci)) return(pip_col)
  }
  if ("pip" %in% names(top_loci)) return("pip")
  pip_cols <- grep("^pip_", names(top_loci), value = TRUE)
  if (length(pip_cols) == 1) return(pip_cols)
  NULL
}

format_cs_column <- function(coverage, method) {
  pct <- as.numeric(coverage) * 100
  if (is.na(pct)) stop("coverage must be numeric.")
  label <- if (abs(pct - round(pct)) < 1e-8) {
    as.character(as.integer(round(pct)))
  } else {
    gsub("\\.", "_", format(pct, scientific = FALSE, trim = TRUE))
  }
  paste0("CS_", label, "_", method)
}

.translate_legacy_cs_column_name <- function(coverage) {
  if (is.null(coverage)) return(NULL)
  vapply(coverage, function(x) {
    x <- as.character(x)
    old_match <- regexec("^cs_coverage_([0-9.]+)$", x, ignore.case = TRUE)
    old_parts <- regmatches(x, old_match)[[1]]
    if (length(old_parts) == 2) return(format_cs_column(as.numeric(old_parts[[2]]), "susie"))
    x
  }, character(1), USE.NAMES = FALSE)
}

.translate_legacy_top_loci_cs_columns <- function(top_loci) {
  if (!is.data.frame(top_loci)) return(top_loci)
  names(top_loci) <- .translate_legacy_cs_column_name(names(top_loci))
  if ("pip_susie" %in% names(top_loci) && !"pip" %in% names(top_loci)) {
    names(top_loci)[names(top_loci) == "pip_susie"] <- "pip"
  }
  top_loci
}

.set_finemapping_fit_class <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  method_class <- switch(method,
    susie = "susie",
    susie_inf = "susie_inf",
    susie_rss = "susie_rss",
    single_effect = "susie_rss",
    bayesian_conditional_regression = "susie_rss",
    fsusie = "susiF",
    mvsusie = "mvsusie",
    NULL
  )
  if (!is.null(method_class)) class(fit) <- unique(c(method_class, class(fit)))
  fit
}

prepare_susie_from_inf_args <- function(args, susie_inf_fit, refine_default = NULL) {
  L <- args[["L"]]
  if (is.null(L)) L <- length(susie_inf_fit$V)
  if (is.null(args[["refine"]]) && !is.null(refine_default)) args[["refine"]] <- refine_default
  args[["unmappable_effects"]] <- "none"
  args[["model_init"]] <- susie_inf_fit
  if (!is.null(args[["L_greedy"]])) args[["L_greedy"]] <- min(length(susie_inf_fit$V), L)
  args
}

fit_susie_inf_then_susie <- function(X, y, args = list(),
                                     susie_inf_args = list(),
                                     susie_args = list(),
                                     fitted_models = NULL) {
  if (is.null(fitted_models)) fitted_models <- list()
  susie_inf_fit <- fitted_models[["susie_inf"]]
  susie_fit <- fitted_models[["susie"]]

  if (is.null(susie_inf_fit)) {
    fit_args <- modifyList(args, susie_inf_args)
    fit_args <- modifyList(fit_args, list(
      X = X, y = y, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susie_inf_fit <- do.call(susie, fit_args)
  }
  susie_inf_fit <- .set_finemapping_fit_class(susie_inf_fit, "susie_inf")

  if (is.null(susie_fit)) {
    fit_args <- prepare_susie_from_inf_args(modifyList(args, susie_args), susie_inf_fit, refine_default = TRUE)
    susie_fit <- do.call(susie, c(list(X = X, y = y), fit_args))
  }
  susie_fit <- .set_finemapping_fit_class(susie_fit, "susie")

  list(susie = susie_fit, susie_inf = susie_inf_fit)
}

#' Post-process Fine-mapping Fits
#'
#' Applies method-aware post-processing to one or more SuSiE-family fits and
#' builds both a method-specific result list and shared top-loci tables.
#'
#' @param fits Named list of fine-mapping fits. Names define method identity,
#'   for example \code{susie}, \code{susie_inf}, \code{susie_rss},
#'   \code{mvsusie}, or \code{fsusie}.
#' @param data_x Genotype matrix, LD/correlation matrix, or other method-specific
#'   input used for credible-set purity and correlations.
#' @param data_y Phenotype vector/matrix or summary statistics. Default NULL.
#' @param X_scalar Scaling factor for genotype effects. Default 1.
#' @param y_scalar Scaling factor for phenotype effects. Default 1.
#' @param maf Minor allele frequencies. Default NULL.
#' @param coverage Primary credible-set coverage.
#' @param secondary_coverage Additional credible-set coverages.
#' @param signal_cutoff PIP cutoff for including non-CS variants in top loci.
#' @param other_quantities Optional list carried into each method result.
#' @param prior_eff_tol Tolerance for retaining effects by prior variance.
#' @param min_abs_corr Minimum absolute correlation for credible-set purity.
#' @return A list with \code{finemapping_results}, \code{top_loci_long}, and
#'   \code{top_loci}. The long table is lossless, with one row per
#'   variant-method-coverage-CS membership. The wide table stores one row per
#'   variant, method-specific \code{pip_<method>} columns, method-specific
#'   \code{CS_<coverage>_<method>} columns, and \code{model_source}.
#' @export
postprocess_finemapping_fits <- function(fits, data_x, data_y = NULL,
                                         X_scalar = 1, y_scalar = 1,
                                         maf = NULL, coverage = NULL,
                                         secondary_coverage = c(0.7, 0.5),
                                         signal_cutoff = 0.1,
                                         other_quantities = NULL,
                                         prior_eff_tol = 1e-9,
                                         min_abs_corr = 0.8) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) == 0) stop("At least one fine-mapping fit must be supplied.")
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("fits must be a named list; names define method identity.")
  }

  posts <- lapply(names(fits), function(method) {
    fit <- .set_finemapping_fit_class(fits[[method]], method)
    postprocess_finemapping_fit(
      fit, method = method, data_x = data_x, data_y = data_y,
      X_scalar = X_scalar, y_scalar = y_scalar, maf = maf,
      coverage = coverage, secondary_coverage = secondary_coverage,
      signal_cutoff = signal_cutoff, other_quantities = other_quantities,
      prior_eff_tol = prior_eff_tol, min_abs_corr = min_abs_corr
    )
  })
  names(posts) <- names(fits)

  top_loci_long <- bind_rows(lapply(posts, function(x) x$top_loci_long))
  posts <- lapply(posts, function(x) {
    x$top_loci_long <- NULL
    x
  })
  top_loci <- build_top_loci_wide(top_loci_long, posts)

  list(
    finemapping_results = posts,
    top_loci_long = if (nrow(top_loci_long) > 0) top_loci_long else NULL,
    top_loci = top_loci
  )
}

postprocess_finemapping_fit <- function(fit, ...) {
  UseMethod("postprocess_finemapping_fit")
}

#' @exportS3Method
postprocess_finemapping_fit.susie <- function(fit, method = "susie", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "X", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susie_inf <- function(fit, method = "susie_inf", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "X", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susie_rss <- function(fit, method = "susie_rss", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "Xcorr", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.mvsusie <- function(fit, method = "mvsusie", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "X", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susiF <- function(fit, method = "fsusie", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "fsusie", ...)
}

.postprocess_finemapping_fit_common <- function(fit, method, data_x, data_y = NULL,
                                                X_scalar = 1, y_scalar = 1,
                                                maf = NULL, coverage = NULL,
                                                secondary_coverage = c(0.7, 0.5),
                                                signal_cutoff = 0.1,
                                                other_quantities = NULL,
                                                prior_eff_tol = 1e-9,
                                                min_abs_corr = 0.8,
                                                cs_input = c("X", "Xcorr", "fsusie")) {
  cs_input <- match.arg(cs_input)
  variant_names <- extract_variant_names(fit)
  sumstats <- extract_sumstats(fit, data_x, data_y, X_scalar, y_scalar, method)
  effect_idx <- select_effects(fit, prior_eff_tol)
  cs_tables <- compute_cs_tables(
    fit, data_x = data_x, coverage = coverage,
    secondary_coverage = secondary_coverage, method = method,
    cs_input = cs_input, min_abs_corr = min_abs_corr
  )
  top_loci_long <- build_top_loci_long(
    fit, cs_tables, variant_names = variant_names, sumstats = sumstats,
    maf = maf, method = method, signal_cutoff = signal_cutoff
  )

  trimmed <- trim_finemapping_fit(fit, effect_idx, method, cs_tables)

  # Build FineMappingResult S4 object
  fm_result <- FineMappingResult(
    variant_names = variant_names,
    trimmed_fit = trimmed,
    top_loci = if (is.null(top_loci_long) || nrow(top_loci_long) == 0) {
      data.frame(variant_id = character(0), method = character(0))
    } else {
      top_loci_long
    },
    method = method,
    sumstats = sumstats
  )

  # Also return as list for backwards compatibility with existing consumers
  res <- list(
    variant_names = variant_names,
    result_trimmed = trimmed,
    top_loci_long = top_loci_long,
    finemapping_result = fm_result
  )
  if (!is.null(sumstats)) res$sumstats <- sumstats
  sample_names <- .sample_names_from_data_y(data_y)
  if (!is.null(sample_names)) res$sample_names <- sample_names
  if (method == "mvsusie" && !is.null(fit$outcome_names)) res$context_names <- fit$outcome_names
  analysis_script <- load_script()
  if (analysis_script != "") res$analysis_script <- analysis_script
  if (!is.null(other_quantities)) res$other_quantities <- other_quantities
  res
}

extract_variant_names <- function(fit) {
  variant_names <- names(fit$pip)
  if (is.null(variant_names)) variant_names <- colnames(fit$alpha)
  if (is.null(variant_names)) variant_names <- paste0("variant_", seq_along(fit$pip))
  tryCatch(normalize_variant_id(variant_names), error = function(e) variant_names)
}

extract_sumstats <- function(fit, data_x, data_y, X_scalar = 1, y_scalar = 1, method = "susie") {
  if (is.null(data_y)) return(NULL)
  if (method == "susie_rss") return(data_y)
  if (is.list(data_y) && !is.data.frame(data_y) &&
      any(c("betahat", "sebetahat", "z") %in% names(data_y))) {
    return(data_y)
  }
  if (is.null(data_x)) return(NULL)
  if (is.matrix(data_y) || is.data.frame(data_y)) {
    if (ncol(as.matrix(data_y)) != 1) return(NULL)
  }
  sumstats <- univariate_regression(data_x, data_y)
  y_scalar <- if (is.null(y_scalar) || all(y_scalar == 1)) 1 else y_scalar
  X_scalar <- if (is.null(X_scalar) || all(X_scalar == 1)) 1 else X_scalar
  sumstats$betahat <- sumstats$betahat * y_scalar / X_scalar
  sumstats$sebetahat <- sumstats$sebetahat * y_scalar / X_scalar
  sumstats
}

.sample_names_from_data_y <- function(data_y) {
  if (is.null(data_y) || is.list(data_y)) return(NULL)
  rownames(as.matrix(data_y))
}

select_effects <- function(fit, prior_eff_tol = 1e-9) {
  alpha <- .as_effect_matrix(fit$alpha)
  n_effects <- nrow(alpha)
  if (n_effects == 0) return(integer(0))
  if (!is.null(fit$V)) {
    which(fit$V > prior_eff_tol)
  } else {
    seq_len(n_effects)
  }
}

.as_effect_matrix <- function(x) {
  if (is.null(x)) return(matrix(numeric(0), nrow = 0))
  if (is.list(x) && !is.data.frame(x)) return(do.call(rbind, x))
  as.matrix(x)
}

.as_lbf_matrix <- function(fit) {
  if (!is.null(fit$lbf_variable)) return(.as_effect_matrix(fit$lbf_variable))
  if (!is.null(fit$lBF)) return(.as_effect_matrix(fit$lBF))
  NULL
}

#' @importFrom susieR get_cs_correlation
#' @noRd
compute_cs_tables <- function(fit, data_x, coverage = NULL,
                              secondary_coverage = c(0.7, 0.5),
                              method = "susie", cs_input = c("X", "Xcorr", "fsusie"),
                              min_abs_corr = 0.8) {
  cs_input <- match.arg(cs_input)
  primary_coverage <- coverage
  if (is.null(primary_coverage)) primary_coverage <- fit$sets$requested_coverage
  if (is.null(primary_coverage)) primary_coverage <- 0.95
  coverages <- unique(c(primary_coverage, secondary_coverage))
  coverages <- coverages[!is.na(coverages)]

  tables <- lapply(coverages, function(cov) {
    compute_cs_table(fit, data_x, coverage = cov, cs_input = cs_input, min_abs_corr = min_abs_corr)
  })
  names(tables) <- vapply(coverages, format_cs_column, character(1), method = method)
  attr(tables, "coverage") <- coverages
  tables
}

compute_cs_table <- function(fit, data_x, coverage, cs_input = c("X", "Xcorr", "fsusie"),
                             min_abs_corr = 0.8) {
  cs_input <- match.arg(cs_input)
  if (cs_input == "fsusie") {
    sets <- tryCatch(
      fsusie_get_cs(fit, data_x, requested_coverage = coverage),
      error = function(e) list(cs = list(), requested_coverage = coverage)
    )
    if (is.null(sets$cs) || length(sets$cs) == 0 || all(vapply(sets$cs, is.null, logical(1)))) {
      sets$cs <- list()
      return(list(sets = sets, cs_corr = NULL, pip = fit$pip))
    }
    tmp <- fit
    tmp$sets <- sets
    cs_corr <- if (requireNamespace("fsusieR", quietly = TRUE)) {
      tryCatch(fsusieR::cal_cor_cs(tmp, data_x), error = function(e) NULL)
    } else {
      NULL
    }
    return(list(sets = sets, cs_corr = cs_corr, pip = fit$pip))
  }

  if (cs_input == "X") {
    sets <- susie_get_cs(fit, X = data_x, coverage = coverage, min_abs_corr = min_abs_corr)
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, X = data_x)
  } else {
    sets <- susie_get_cs(fit, Xcorr = data_x, coverage = coverage, min_abs_corr = min_abs_corr)
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, Xcorr = data_x)
  }
  out
}

build_top_loci_long <- function(fit, cs_tables, variant_names, sumstats = NULL,
                                maf = NULL, method, signal_cutoff = 0.1) {
  if (length(cs_tables) == 0) return(.empty_top_loci_long())
  coverage_values <- attr(cs_tables, "coverage")
  rows <- lapply(seq_along(cs_tables), function(i) {
    cs_table <- cs_tables[[i]]
    cov <- coverage_values[[i]]
    top_variants_idx <- get_top_variants_idx(cs_table, signal_cutoff)
    cs_info <- get_cs_info(cs_table$sets$cs, top_variants_idx)
    if (is.null(cs_info) || nrow(cs_info) == 0) return(NULL)
    idx <- cs_info$variant_idx
    optional_cols <- .top_loci_optional_columns(idx, sumstats, maf)
    base <- data.frame(
      variant_id = variant_names[idx],
      method = method,
      coverage = cov,
      cs = as.integer(cs_info$cs_idx),
      pip = as.numeric(fit$pip[idx]),
      stringsAsFactors = FALSE
    )
    if (ncol(optional_cols) > 0) cbind(base, optional_cols) else base
  })
  out <- bind_rows(rows)
  if (nrow(out) == 0) .empty_top_loci_long() else out
}

.empty_top_loci_long <- function() {
  data.frame(
    variant_id = character(), method = character(), coverage = numeric(),
    cs = integer(), pip = numeric(), stringsAsFactors = FALSE
  )
}

.top_loci_optional_columns <- function(idx, sumstats = NULL, maf = NULL) {
  optional_cols <- list(
    betahat = if (!is.null(sumstats$betahat)) sumstats$betahat[idx] else NULL,
    sebetahat = if (!is.null(sumstats$sebetahat)) sumstats$sebetahat[idx] else NULL,
    z = if (!is.null(sumstats$z)) sumstats$z[idx] else NULL,
    maf = if (!is.null(maf)) maf[idx] else NULL
  )
  as.data.frame(Filter(Negate(is.null), optional_cols))
}

build_top_loci_wide <- function(top_loci_long, posts) {
  if (is.null(top_loci_long) || nrow(top_loci_long) == 0) return(NULL)
  ids <- unique(top_loci_long$variant_id)
  out <- data.frame(variant_id = ids, stringsAsFactors = FALSE)
  for (column in c("betahat", "sebetahat", "z", "maf")) {
    if (column %in% names(top_loci_long)) {
      out[[column]] <- vapply(ids, function(id) {
        values <- top_loci_long[[column]][top_loci_long$variant_id == id]
        values <- values[!is.na(values)]
        if (length(values) == 0) NA_real_ else values[[1]]
      }, numeric(1))
    }
  }

  methods <- names(posts)
  for (method in methods) {
    post <- posts[[method]]
    pip_col <- format_pip_column(method)
    pip <- post$result_trimmed$pip
    names(pip) <- post$variant_names
    out[[pip_col]] <- as.numeric(pip[ids])

    method_rows <- top_loci_long[top_loci_long$method == method, , drop = FALSE]
    for (cov in sort(unique(top_loci_long$coverage), decreasing = TRUE)) {
      cs_col <- format_cs_column(cov, method)
      out[[cs_col]] <- vapply(ids, function(id) {
        cs <- method_rows$cs[method_rows$variant_id == id & method_rows$coverage == cov]
        if (length(cs) == 0) return(0L)
        min(cs)
      }, integer(1))
    }
  }

  out$model_source <- vapply(ids, function(id) {
    selected_methods <- unique(top_loci_long$method[top_loci_long$variant_id == id])
    paste(selected_methods[selected_methods %in% methods], collapse = ";")
  }, character(1))
  rownames(out) <- NULL
  out
}

trim_finemapping_fit <- function(fit, effect_idx, method, cs_tables) {
  alpha <- .as_effect_matrix(fit$alpha)
  lbf_variable <- .as_lbf_matrix(fit)
  primary <- cs_tables[[1]]
  secondary <- if (length(cs_tables) > 1) {
    lapply(cs_tables[-1], function(x) x[names(x) != "pip"])
  } else {
    NULL
  }

  trimmed <- list(
    pip = as.numeric(fit$pip),
    sets = primary$sets,
    cs_corr = primary$cs_corr,
    sets_secondary = secondary,
    alpha = alpha[effect_idx, , drop = FALSE],
    lbf_variable = if (!is.null(lbf_variable)) lbf_variable[effect_idx, , drop = FALSE] else NULL,
    V = if (!is.null(fit$V)) fit$V[effect_idx] else NULL,
    niter = fit$niter,
    max_L = nrow(alpha),
    n_effects = nrow(alpha)
  )

  if (!is.null(fit$X_column_scale_factors)) trimmed$X_column_scale_factors <- fit$X_column_scale_factors
  if (!is.null(fit$mu)) {
    trimmed$mu <- if (length(dim(fit$mu)) == 3) fit$mu[effect_idx, , , drop = FALSE] else fit$mu[effect_idx, , drop = FALSE]
  }
  if (!is.null(fit$mu2)) {
    # mu2 is L x p for univariate susie and L x p x R for multivariate (mvsusie).
    # Match the shape handling used for mu just above.
    trimmed$mu2 <- if (length(dim(fit$mu2)) == 3) fit$mu2[effect_idx, , , drop = FALSE] else fit$mu2[effect_idx, , drop = FALSE]
  }
  if (!is.null(fit$theta)) trimmed$theta <- fit$theta
  if (!is.null(fit$omega_weights)) trimmed$omega_weights <- fit$omega_weights

  if (method == "mvsusie") {
    if (!is.null(fit$mu2_diag)) trimmed$mu2_diag <- fit$mu2_diag[effect_idx, , , drop = FALSE]
    if (requireNamespace("mvsusieR", quietly = TRUE)) {
      trimmed$coef <- mvsusieR::coef.mvsusie(fit)[-1, , drop = FALSE]
    }
    if (!is.null(fit$conditional_lfsr)) trimmed$clfsr <- fit$conditional_lfsr[effect_idx, , , drop = FALSE]
  }

  class(trimmed) <- unique(c(method, "susie"))
  trimmed
}

#' Format Fine-mapping Post-processing for Protocol Output
#'
#' Converts method-aware fine-mapping post-processing output into the root-level
#' fields consumed by protocol RDS files.
#'
#' @param post Output from \code{\link{postprocess_finemapping_fits}}.
#' @param primary_method Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{variant_names},
#'   \code{susie_result_trimmed}, \code{top_loci_long}, and \code{top_loci}.
#' @export
format_finemapping_output <- function(post, primary_method) {
  method_post <- post$finemapping_results[[primary_method]]
  if (is.null(method_post)) {
    stop("primary_method was not found in finemapping_results: ", primary_method)
  }
  keep_names <- setdiff(names(method_post), c("result_trimmed", "top_loci_long"))
  c(
    method_post[keep_names],
    list(
      susie_result_trimmed = method_post$result_trimmed,
      top_loci_long = post$top_loci_long,
      top_loci = post$top_loci
    )
  )
}

#' Adjust SuSiE Weights
#'
#' Adjusts SuSiE TWAS weights by subsetting to intersected variants and
#' optionally running allele QC against LD reference variants.
#'
#' @param twas_weights_results A list containing TWAS weight data (nested structure).
#' @param keep_variants Vector of variant names to keep.
#' @param run_allele_qc Whether to run allele_qc to align alleles. Default TRUE.
#' @param variable_name_obj Path to variant names in the nested list.
#' @param susie_obj Path to susie result in the nested list.
#' @param twas_weights_table Path to weights table in the nested list.
#' @param LD_variants Vector of LD reference variant IDs for allele QC.
#' @param match_min_prop Minimum proportion of matched variants. Default 0.2.
#' @return A list with adjusted_susie_weights and remained_variants_ids.
#' @export
adjust_susie_weights <- function(twas_weights_results, keep_variants, run_allele_qc = TRUE,
                                 variable_name_obj = c("susie_results", context, "variant_names"),
                                 susie_obj = c("susie_results", context, "susie_result_trimmed"),
                                 twas_weights_table = c("weights", context), LD_variants, match_min_prop = 0.2) {
  # Intersect the rownames of weights with keep_variants
  twas_weights_variants <- get_nested_element(twas_weights_results, variable_name_obj)
  # Normalize to canonical format (with chr prefix)
  twas_weights_variants <- normalize_variant_id(twas_weights_variants)
  # allele flip twas weights matrix variants name
  if (run_allele_qc) {
    weights_matrix <- get_nested_element(twas_weights_results, twas_weights_table)
    if (!all(c("chrom", "pos", "A2", "A1") %in% colnames(weights_matrix))) {
      weights_matrix <- cbind(parse_variant_id(twas_weights_variants), weights_matrix)
    }
    weights_matrix_qced <- match_ref_panel(weights_matrix, LD_variants, colnames(weights_matrix)[!colnames(weights_matrix) %in% c(
      "chrom",
      "pos", "A2", "A1"
    )], match_min_prop = match_min_prop)
    # match_ref_panel outputs canonical variant_ids (with chr prefix)
    original_idx <- match(weights_matrix_qced$qc_summary$variants_id_original, twas_weights_variants)
    intersected_indices <- original_idx[weights_matrix_qced$qc_summary$keep == TRUE]
  } else {
    # Normalize keep_variants to canonical format for matching
    keep_variants_normalized <- normalize_variant_id(keep_variants)
    intersected_variants <- intersect(twas_weights_variants, keep_variants_normalized)
    intersected_indices <- match(intersected_variants, twas_weights_variants)
  }
  if (length(intersected_indices) == 0) {
    stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
  }
  # Subset lbf_matrix, mu, and x_column_scale_factors
  lbf_matrix <- get_nested_element(twas_weights_results, c(susie_obj, "lbf_variable"))
  mu <- get_nested_element(twas_weights_results, c(susie_obj, "mu"))
  x_column_scal_factors <- get_nested_element(twas_weights_results, c(susie_obj, "X_column_scale_factors"))

  lbf_matrix_subset <- lbf_matrix[, intersected_indices, drop = FALSE]
  mu_subset <- mu[, intersected_indices, drop = FALSE]
  x_column_scal_factors_subset <- x_column_scal_factors[intersected_indices]

  # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
  adjusted_xqtl_alpha <- lbf_to_alpha(lbf_matrix_subset)
  adjusted_xqtl_coef <- colSums(adjusted_xqtl_alpha * mu_subset) / x_column_scal_factors_subset
  # allele_qc now outputs canonical variant_ids (with chr prefix) -- no need to add chr
  remained_variants_ids <- if (run_allele_qc) {
    weights_matrix_qced$target_data_qced$variant_id
  } else {
    intersected_variants
  }
  return(list(adjusted_susie_weights = adjusted_xqtl_coef, remained_variants_ids = remained_variants_ids))
}

#' Run the SuSiE RSS pipeline
#'
#' Runs SuSiE RSS analysis with the specified method. Supports both z+R
#' (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstats Data frame with 'z' or ('beta' and 'se') columns.
#' @param LD_mat LD correlation matrix. Mutually exclusive with X_mat.
#' @param X_mat Genotype matrix (samples x variants). Mutually exclusive with LD_mat.
#' @param n Sample size.
#' @param L Maximum number of causal configurations (default: 30).
#' @param L_greedy Initial greedy number of causal configurations (default: 5).
#' @param analysis_method One of "susie_rss", "single_effect", "bayesian_conditional_regression".
#' @param coverage Coverage level (default: 0.95).
#' @param secondary_coverage Secondary coverage levels (default: c(0.7, 0.5)).
#' @param signal_cutoff PIP cutoff for selecting top loci (default: 0.1).
#' @param min_abs_corr Minimum absolute correlation for CS purity (default: 0.8).
#' @param R_finite Controls variance inflation to account for estimating
#'   the R matrix from a finite reference panel. NULL (default): no
#'   variance inflation. Passed directly to susie_rss.
#' @param R_mismatch LD mismatch correction method passed directly to susie_rss.
#'   Default NULL disables mismatch correction.
#' @param ... Additional parameters passed to susie_rss (e.g., var_y).
#' @return A list with post-processed SuSiE RSS results. Method-specific
#'   \code{top_loci} columns use the selected \code{analysis_method}, for
#'   example \code{pip_susie_rss} or \code{pip_single_effect}.
#' @importFrom susieR susie_rss
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange select
#' @export
susie_rss_pipeline <- function(sumstats, LD_mat = NULL, X_mat = NULL, n = NULL,
                               L = 30, L_greedy = 5,
                               analysis_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
                               coverage = 0.95,
                               secondary_coverage = c(0.7, 0.5),
                               signal_cutoff = 0.1,
                               min_abs_corr = 0.8,
                               R_finite = NULL, R_mismatch = NULL, ...) {
  analysis_method <- match.arg(analysis_method)
  if (is.null(LD_mat) && is.null(X_mat)) stop("Either LD_mat or X_mat must be provided.")
  if (!is.null(LD_mat) && !is.null(X_mat)) stop("Only one of LD_mat or X_mat should be provided, not both.")
  if (!is.null(L_greedy)) L_greedy <- min(L_greedy, L)

  if (!is.null(sumstats$z)) {
    z <- sumstats$z
  } else if (!is.null(sumstats$beta) && !is.null(sumstats$se)) {
    z <- sumstats$beta / sumstats$se
  } else {
    stop("sumstats must have 'z' or ('beta' and 'se') columns.")
  }
  if (is.null(names(z)) && !is.null(sumstats$variant_id) && length(sumstats$variant_id) == length(z)) {
    names(z) <- sumstats$variant_id
  }
  if (is.null(names(z)) && !is.null(rownames(sumstats)) && length(rownames(sumstats)) == length(z)) {
    names(z) <- rownames(sumstats)
  }

  common <- list(z = z, n = n, coverage = coverage,
                 R_finite = R_finite, R_mismatch = R_mismatch, ...)
  if (!is.null(X_mat)) common$X <- X_mat else common$R <- LD_mat

  if (analysis_method == "single_effect") {
    res <- do.call(susie_rss, c(common, list(L = 1, L_greedy = NULL, max_iter = 1)))
  } else if (analysis_method == "bayesian_conditional_regression") {
    res <- do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy, max_iter = 1)))
  } else {
    res <- do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy)))
  }

  # For post-processing, need a square matrix (R or computed from X).
  # For mixture panels (list of X), use the first panel to compute R.
  if (!is.null(LD_mat)) {
    data_x <- LD_mat
  } else if (is.list(X_mat) && !is.matrix(X_mat)) {
    data_x <- compute_LD(X_mat[[1]][, seq_along(z), drop = FALSE], method = "sample")
  } else {
    data_x <- compute_LD(X_mat[, seq_along(z), drop = FALSE], method = "sample")
  }

  rss_method <- analysis_method
  rss_fit <- .set_finemapping_fit_class(res, rss_method)
  post <- postprocess_finemapping_fits(
    fits = setNames(list(rss_fit), rss_method),
    data_x = data_x,
    data_y = list(z = z),
    coverage = coverage,
    secondary_coverage = secondary_coverage,
    signal_cutoff = signal_cutoff,
    min_abs_corr = min_abs_corr
  )
  format_finemapping_output(post, primary_method = rss_method)
}

#' @noRd
get_cs_index <- function(snps_idx, susie_cs) {
  # Return ALL CS indices that contain this variant (not just one)
  idx <- which(vapply(susie_cs, function(x) snps_idx %in% x, logical(1)))
  if (length(idx) == 0) return(NA_integer_)
  return(idx)
}
#' @noRd
get_top_variants_idx <- function(susie_output, signal_cutoff) {
  c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>%
    unique() %>%
    sort()
}
# Returns a data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair.
# Variants in multiple CSs get multiple rows.
#' @importFrom stringr str_replace
#' @noRd
get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
  cs_names <- names(susie_output_sets_cs)
  rows <- lapply(top_variants_idx, function(vi) {
    idx <- get_cs_index(vi, susie_output_sets_cs)
    if (length(idx) == 1 && is.na(idx)) {
      data.frame(variant_idx = vi, cs_idx = 0L, stringsAsFactors = FALSE)
    } else {
      cs_nums <- as.integer(str_replace(cs_names[idx], "L", ""))
      data.frame(variant_idx = rep(vi, length(cs_nums)), cs_idx = cs_nums, stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}
