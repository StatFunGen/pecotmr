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

# Build the argument list for a SuSiE / SuSiE-ash fit initialised from a
# prior SuSiE-inf fit. `unmappable_effects` controls which branch the
# downstream fit takes: "none" yields the standard SuSiE-inf-initialised
# SuSiE; "ash" yields SuSiE-ash with the SuSiE-inf warm start.
prepare_susie_from_inf_args <- function(args, susie_inf_fit, refine_default = NULL,
                                        unmappable_effects = c("none", "ash")) {
  unmappable_effects <- match.arg(unmappable_effects)
  L <- args[["L"]]
  if (is.null(L)) L <- length(susie_inf_fit$V)
  if (is.null(args[["refine"]]) && !is.null(refine_default)) args[["refine"]] <- refine_default
  args[["unmappable_effects"]] <- unmappable_effects
  args[["model_init"]] <- susie_inf_fit
  if (unmappable_effects == "ash") {
    args[["convergence_method"]] <- args[["convergence_method"]] %||% "pip"
  }
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

#' Two-stage SuSiE-RSS Fine-mapping
#'
#' RSS analog of \code{fit_susie_inf_then_susie}. Fits SuSiE-inf via
#' \code{susie_rss} first, then initialises standard SuSiE-RSS from
#' the SuSiE-inf result. The single pair of fits can be used both for
#' fine-mapping post-processing and TWAS weight extraction.
#'
#' @param z Numeric vector of z-scores.
#' @param R LD correlation matrix.
#' @param n Sample size (scalar).
#' @param args Default arguments forwarded to both fits.
#' @param susie_inf_args SuSiE-inf-specific overrides.
#' @param susie_args Standard SuSiE-RSS-specific overrides.
#' @param fitted_models Optional list with pre-fitted \code{$susie} and/or
#'   \code{$susie_inf} objects to skip re-fitting.
#' @return A list with \code{susie} and \code{susie_inf} fit objects.
#' @importFrom susieR susie_rss
#' @export
fit_susie_inf_then_susie_rss <- function(z, R, n, args = list(),
                                         susie_inf_args = list(),
                                         susie_args = list(),
                                         fitted_models = NULL) {
  if (is.null(fitted_models)) fitted_models <- list()
  susie_inf_fit <- fitted_models[["susie_inf"]]
  susie_fit <- fitted_models[["susie"]]

  if (is.null(susie_inf_fit)) {
    fit_args <- modifyList(args, susie_inf_args)
    fit_args <- modifyList(fit_args, list(
      z = z, R = R, n = n, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susie_inf_fit <- do.call(susie_rss, fit_args)
  }
  susie_inf_fit <- .set_finemapping_fit_class(susie_inf_fit, "susie_inf")

  if (is.null(susie_fit)) {
    fit_args <- prepare_susie_from_inf_args(modifyList(args, susie_args), susie_inf_fit, refine_default = TRUE)
    susie_fit <- do.call(susie_rss, c(list(z = z, R = R, n = n), fit_args))
  }
  susie_fit <- .set_finemapping_fit_class(susie_fit, "susie_rss")

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
#' @return A list with \code{finemapping_results} (per-method post-processed
#'   objects, each carrying a trimmed fit and method-specific intermediates)
#'   and a single unified \code{top_loci} table in the fixed 22-column shape
#'   (see \code{\link{build_top_loci}}). Per-method contributions are
#'   row-bound into \code{top_loci} by an outer method for-loop.
#' @export
postprocess_finemapping_fits <- function(fits, data_x, data_y = NULL,
                                         X_scalar = 1, y_scalar = 1,
                                         maf = NULL, coverage = NULL,
                                         secondary_coverage = c(0.7, 0.5),
                                         signal_cutoff = 0.1,
                                         other_quantities = NULL,
                                         region = NULL,
                                         prior_eff_tol = 1e-9,
                                         min_abs_corr = 0.8,
                                         cs_input = NULL) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) == 0) stop("At least one fine-mapping fit must be supplied.")
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("fits must be a named list; names define method identity.")
  }

  # One method for-loop: each method calls build_top_loci() once per fit; the
  # per-method 22-column contributions are row-bound below into the single
  # final `top_loci` table. There is no separately exposed long or wide table.
  posts <- lapply(names(fits), function(method) {
    fit <- .set_finemapping_fit_class(fits[[method]], method)
    postprocess_finemapping_fit(
      fit, method = method, data_x = data_x, data_y = data_y,
      X_scalar = X_scalar, y_scalar = y_scalar, maf = maf,
      coverage = coverage, secondary_coverage = secondary_coverage,
      signal_cutoff = signal_cutoff, other_quantities = other_quantities,
      region = region,
      prior_eff_tol = prior_eff_tol, min_abs_corr = min_abs_corr,
      cs_input = cs_input
    )
  })
  names(posts) <- names(fits)

  per_method <- lapply(posts, function(x) x$top_loci)
  per_method <- per_method[!vapply(per_method, is.null, logical(1))]
  top_loci <- if (length(per_method) == 0L) {
    .empty_top_loci()
  } else {
    do.call(rbind, per_method)
  }
  rownames(top_loci) <- NULL
  posts <- lapply(posts, function(x) {
    x$top_loci <- NULL
    x
  })

  list(
    finemapping_results = posts,
    top_loci = top_loci
  )
}

postprocess_finemapping_fit <- function(fit, ...) {
  UseMethod("postprocess_finemapping_fit")
}

#' @exportS3Method
postprocess_finemapping_fit.susie <- function(fit, method = "susie", cs_input = NULL, ...) {
  if (is.null(cs_input)) cs_input <- "X"
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = cs_input, ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susie_inf <- function(fit, method = "susie_inf", cs_input = NULL, ...) {
  if (is.null(cs_input)) cs_input <- "X"
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = cs_input, ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susie_rss <- function(fit, method = "susie_rss", cs_input = NULL, ...) {
  if (is.null(cs_input)) cs_input <- "Xcorr"
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = cs_input, ...)
}

#' @exportS3Method
postprocess_finemapping_fit.mvsusie <- function(fit, method = "mvsusie", cs_input = NULL, ...) {
  if (is.null(cs_input)) cs_input <- "X"
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = cs_input, ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susiF <- function(fit, method = "fsusie", cs_input = NULL, ...) {
  if (is.null(cs_input)) cs_input <- "fsusie"
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = cs_input, ...)
}

.postprocess_finemapping_fit_common <- function(fit, method, data_x, data_y = NULL,
                                                X_scalar = 1, y_scalar = 1,
                                                maf = NULL, coverage = NULL,
                                                secondary_coverage = c(0.7, 0.5),
                                                signal_cutoff = 0.1,
                                                other_quantities = NULL,
                                                region = NULL,
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
  top_loci <- build_top_loci(
    fit, cs_tables, variant_names = variant_names, sumstats = sumstats,
    maf = maf, method = method, signal_cutoff = signal_cutoff,
    data_x = data_x, data_y = data_y, other_quantities = other_quantities,
    region = region
  )

  trimmed <- trim_finemapping_fit(fit, effect_idx, method, cs_tables)

  # Build FineMappingResult S4 object. The S4 contract (validity check,
  # vcf_writer, getPIP, getCS) still expects `variant_id`, `pip`, and an
  # integer `cs` column on the slot. To avoid rippling renames into
  # AllClasses / AllMethods / vcf_writer for this change, we project the
  # new 22-column `top_loci` into the legacy slot shape here, in
  # susie_wrapper only. The wrapper-facing `top_loci` returned to callers
  # is unchanged.
  s4_top_loci <- .top_loci_for_s4_slot(top_loci)
  fm_result <- FineMappingResult(
    variant_names = variant_names,
    trimmed_fit = trimmed,
    top_loci = s4_top_loci,
    method = method,
    sumstats = sumstats
  )

  res <- list(
    top_loci = top_loci,
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

#' Build the unified top-loci table for one fit and one method.
#'
#' Returns the per-fit, per-method contribution to the unified \code{top_loci}
#' table in the fixed 22-column shape. \code{postprocess_finemapping_fits()}
#' calls this once per method per fit and row-binds the results into the
#' single \code{top_loci} returned by \code{format_finemapping_output()}.
#'
#' Output columns, in order: \code{#chr}, \code{start}, \code{end}, \code{a1},
#' \code{a2}, \code{variant}, \code{gene}, \code{event}, \code{n}, \code{maf},
#' \code{beta}, \code{se}, \code{pip}, \code{posterior_effect_mean},
#' \code{posterior_effect_se}, \code{cs_95}, \code{cs_70}, \code{cs_50},
#' \code{cs_95_purity}, \code{method}, \code{grange_start}, \code{grange_end}.
#'
#' \code{cs_95} / \code{cs_70} / \code{cs_50} are character strings of the
#' form \code{"<method>_<cs_index>"} where each method numbers credible sets
#' independently from 1. Variants retained by the PIP cutoff but not in any
#' credible set at a coverage carry \code{"<method>_0"}. \code{cs_95_purity}
#' is the 0.95-coverage purity for the row's \code{(method, cs_95)}; rows
#' whose \code{cs_95} is \code{"<method>_0"} carry \code{0}.
#'
#' Row uniqueness is \code{(variant, gene, cs_membership)} at the given
#' \code{method}; overlapping CS within the same method produces one row per
#' CS.
#'
#' @param fit Fitted SuSiE-family object (must expose \code{alpha},
#'   \code{mu}, \code{mu2}, \code{pip}).
#' @param cs_tables List of CS tables (one per coverage) from
#'   \code{compute_cs_tables()}.
#' @param variant_names Character vector of variant IDs
#'   (\code{chr:pos:A2:A1}).
#' @param sumstats Optional marginal-association summary (\code{betahat},
#'   \code{sebetahat}) filling \code{beta} / \code{se}.
#' @param maf Optional numeric vector of minor-allele frequencies.
#' @param method Method name (e.g. \code{"susie"}, \code{"susie_inf"}). Required.
#' @param signal_cutoff PIP cutoff for retaining PIP-only (non-CS) variants.
#' @param data_x Optional regional genotype matrix.
#' @param data_y Optional regional phenotype matrix; \code{nrow(data_y)} fills
#'   \code{n}, \code{colnames(data_y)[1]} fills \code{gene}.
#' @param other_quantities Optional list. Default is NULL.
#' @param region Optional \code{"chr:start-end"} string. Default is NULL.
#' @return A data frame in the fixed 22-column shape for this fit and method,
#'   or an empty data frame if nothing is retained.
#' @export
build_top_loci <- function(fit, cs_tables, variant_names, sumstats = NULL,
                           maf = NULL, method, signal_cutoff = 0.1,
                           data_x = NULL, data_y = NULL,
                           other_quantities = NULL,
                           region = NULL) {
  if (missing(method) || is.null(method) ||
      length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("build_top_loci: `method` is required (e.g. \"susie\", \"susie_inf\").")
  }
  if (length(cs_tables) == 0) return(.empty_top_loci())
  coverage_values <- attr(cs_tables, "coverage")
  if (is.null(coverage_values)) coverage_values <- rep(NA_real_, length(cs_tables))

  # Per-fit constants.
  data_y_mat <- if (!is.null(data_y)) as.matrix(data_y) else NULL
  fit_n    <- if (is.null(data_y_mat)) NA_integer_ else as.integer(nrow(data_y_mat))
  fit_gene <- if (!is.null(data_y_mat) && !is.null(colnames(data_y_mat))) {
    colnames(data_y_mat)[1]
  } else NA_character_
  fit_event <- if (!is.null(other_quantities$condition_id) &&
                   !is.na(fit_gene) && nzchar(fit_gene)) {
    paste(other_quantities$condition_id, fit_gene, sep = "_")
  } else NA_character_
  grange <- .parse_grange(region)

  # Per-variant posterior effect / SE, computed once across all variants.
  alpha <- as.matrix(fit$alpha)
  mu    <- if (!is.null(fit$mu))  as.matrix(fit$mu)  else NULL
  mu2   <- if (!is.null(fit$mu2)) as.matrix(fit$mu2) else NULL
  post_mean <- if (!is.null(mu) && all(dim(alpha) == dim(mu))) {
    colSums(alpha * mu)
  } else rep(NA_real_, length(variant_names))
  post_se <- if (!is.null(mu2) && all(dim(alpha) == dim(mu2))) {
    sqrt(pmax(colSums(alpha * mu2) - post_mean^2, 0))
  } else rep(NA_real_, length(variant_names))

  # Collect CS-membership records (variant_idx, cs_idx, coverage) across all
  # requested coverages. This is the only intermediate; the 22-column shape
  # is projected from it below.
  cs_records <- do.call(rbind, lapply(seq_along(cs_tables), function(i) {
    ct <- cs_tables[[i]]
    info <- get_cs_info(ct$sets$cs, get_top_variants_idx(ct, signal_cutoff))
    if (is.null(info) || nrow(info) == 0) return(NULL)
    data.frame(variant_idx = as.integer(info$variant_idx),
               cs_idx      = as.integer(info$cs_idx),
               coverage    = as.numeric(coverage_values[[i]]),
               stringsAsFactors = FALSE)
  }))
  if (is.null(cs_records) || nrow(cs_records) == 0) return(.empty_top_loci())

  # Key grid: one row per (variant_idx, cs_idx). Overlapping CS membership
  # within this method is preserved as separate keys.
  key_grid <- unique(cs_records[, c("variant_idx", "cs_idx"), drop = FALSE])
  rownames(key_grid) <- NULL
  n_keys  <- nrow(key_grid)
  key_str <- paste(key_grid$variant_idx, key_grid$cs_idx, sep = ":")

  # For each requested coverage, which keys appear in cs_records at that
  # coverage? Returns the key's cs_idx if present, else 0L.
  idx_at <- function(cov) {
    at <- cs_records[abs(cs_records$coverage - cov) < 1e-12, , drop = FALSE]
    hits <- paste(at$variant_idx, at$cs_idx, sep = ":")
    ifelse(key_str %in% hits, key_grid$cs_idx, 0L)
  }
  idx95 <- idx_at(0.95); idx70 <- idx_at(0.70); idx50 <- idx_at(0.50)

  # Per-coverage CS purity vectors (indexed by 1-based CS index). Only the
  # 0.95-coverage purity is currently exported (as cs_95_purity); per-CS
  # purities for the other coverages are kept here for downstream / future
  # use even though they are not part of the 22-column output.
  purity_per_cov <- lapply(cs_tables, .cs_purity_vec)
  cov95          <- which(abs(coverage_values - 0.95) < 1e-12)
  purity_95      <- if (length(cov95) > 0L) purity_per_cov[[cov95[1]]] else numeric()
  cs_95_purity   <- vapply(idx95, function(i) {
    if (i <= 0L || i > length(purity_95)) return(0)
    v <- purity_95[i]; if (is.na(v)) 0 else as.numeric(v)
  }, numeric(1))

  v_idx          <- key_grid$variant_idx
  variant_id_vec <- variant_names[v_idx]
  parsed <- tryCatch(
    suppressWarnings(parse_variant_id(variant_id_vec)),
    error = function(e) stop("build_top_loci: parse_variant_id failed: ",
                             conditionMessage(e))
  )
  if (is.null(parsed) || nrow(parsed) != length(variant_id_vec)) {
    stop("build_top_loci: parse_variant_id did not return one row per variant.")
  }
  invalid <- is.na(parsed$chrom) | is.na(parsed$pos) |
    is.na(parsed$A1) | !nzchar(parsed$A1) |
    is.na(parsed$A2) | !nzchar(parsed$A2)
  if (any(invalid)) {
    stop("build_top_loci: parse_variant_id produced invalid coordinates ",
         "for variant_id: ", variant_id_vec[which(invalid)[[1]]])
  }
  pick <- function(x) if (is.null(x)) rep(NA_real_, n_keys) else x[v_idx]

  out <- data.frame(
    "#chr"                = parsed$chrom,
    start                 = as.integer(parsed$pos) - 1L,
    end                   = as.integer(parsed$pos),
    a1                    = parsed$A1,
    a2                    = parsed$A2,
    variant               = variant_id_vec,
    gene                  = rep(fit_gene, n_keys),
    event                 = rep(fit_event, n_keys),
    n                     = rep(fit_n, n_keys),
    maf                   = pick(maf),
    beta                  = pick(sumstats$betahat),
    se                    = pick(sumstats$sebetahat),
    pip                   = as.numeric(fit$pip[v_idx]),
    posterior_effect_mean = post_mean[v_idx],
    posterior_effect_se   = post_se[v_idx],
    cs_95                 = paste0(method, "_", idx95),
    cs_70                 = paste0(method, "_", idx70),
    cs_50                 = paste0(method, "_", idx50),
    cs_95_purity          = cs_95_purity,
    method                = rep(method, n_keys),
    grange_start          = rep(grange[["start"]], n_keys),
    grange_end            = rep(grange[["end"]],   n_keys),
    stringsAsFactors      = FALSE,
    check.names           = FALSE
  )
  rownames(out) <- NULL
  out
}

# Per-CS purity from one cs_table: prefer susieR's sets$purity$min.abs.corr;
# fall back to cs_corr when purity is unavailable.
.cs_purity_vec <- function(ct) {
  sp <- ct$sets$purity
  if (!is.null(sp) && "min.abs.corr" %in% names(sp)) {
    return(as.numeric(sp$min.abs.corr))
  }
  if (!is.null(ct$cs_corr)) {
    return(vapply(ct$cs_corr, function(m) {
      if (is.null(m)) return(NA_real_)
      if (!is.matrix(m) || nrow(m) <= 1) return(1)
      min(abs(m[upper.tri(m)]))
    }, numeric(1)))
  }
  rep(NA_real_, length(ct$sets$cs))
}

.empty_top_loci <- function() {
  data.frame(
    "#chr"                = character(),
    start                 = integer(),
    end                   = integer(),
    a1                    = character(),
    a2                    = character(),
    variant               = character(),
    gene                  = character(),
    event                 = character(),
    n                     = integer(),
    maf                   = numeric(),
    beta                  = numeric(),
    se                    = numeric(),
    pip                   = numeric(),
    posterior_effect_mean = numeric(),
    posterior_effect_se   = numeric(),
    cs_95                 = character(),
    cs_70                 = character(),
    cs_50                 = character(),
    cs_95_purity          = numeric(),
    method                = character(),
    grange_start          = integer(),
    grange_end            = integer(),
    stringsAsFactors      = FALSE,
    check.names           = FALSE
  )
}

.parse_grange <- function(region_str) {
  if (is.null(region_str) || length(region_str) == 0L ||
      is.na(region_str) || !nzchar(as.character(region_str))) {
    return(c(start = NA_integer_, end = NA_integer_))
  }
  pr <- tryCatch(parse_region(as.character(region_str)),
                 error = function(e) NULL)
  if (is.null(pr) || !is.data.frame(pr)) {
    return(c(start = NA_integer_, end = NA_integer_))
  }
  c(start = as.integer(pr$start), end = as.integer(pr$end))
}

# Project the new 22-column `top_loci` into the legacy shape expected by the
# FineMappingResult S4 slot, vcf_writer, getPIP, and getCS. We add backward-
# compatible aliases without renaming any column in the wrapper-facing
# `top_loci`:
#
#   * `variant_id` — copy of `variant`
#   * `cs`         — integer credible-set index derived from `cs_95` strings of
#                    the form `<method>_<idx>` (PIP-only `<method>_0` -> 0L)
#
# This isolates the schema change to susie_wrapper.R so AllClasses.R,
# AllMethods.R, and vcf_writer.R do not have to change.
.top_loci_for_s4_slot <- function(top_loci) {
  if (is.null(top_loci) || nrow(top_loci) == 0) {
    return(data.frame(variant_id = character(0),
                      method     = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- top_loci
  if ("variant" %in% names(out) && !"variant_id" %in% names(out)) {
    out$variant_id <- out$variant
  }
  if ("cs_95" %in% names(out) && !"cs" %in% names(out)) {
    out$cs <- vapply(out$cs_95, function(s) {
      if (is.na(s) || !nzchar(s)) return(0L)
      tail_str <- sub("^.*_", "", s)
      suppressWarnings(as.integer(tail_str))
    }, integer(1))
    out$cs[is.na(out$cs)] <- 0L
  }
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
#' fields consumed by protocol RDS files. The primary method's
#' \code{FineMappingResult} S4 object is promoted to the \code{finemapping_result}
#' field; use its accessors (\code{getTrimmedFit}, \code{getVariantNames},
#' \code{getTopLoci}, etc.) instead of legacy list keys.
#'
#' @param post Output from \code{\link{postprocess_finemapping_fits}}.
#' @param primary_method Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{finemapping_result}
#'   and \code{top_loci}.
#' @export
format_finemapping_output <- function(post, primary_method) {
  method_post <- post$finemapping_results[[primary_method]]
  if (is.null(method_post)) {
    stop("primary_method was not found in finemapping_results: ", primary_method)
  }
  c(
    method_post,
    list(
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
    qc_summary_df <- getQCSummary(weights_matrix_qced)
    original_idx <- match(qc_summary_df$variants_id_original, twas_weights_variants)
    intersected_indices <- original_idx[qc_summary_df$keep == TRUE]
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
    getHarmonizedData(weights_matrix_qced)$variant_id
  } else {
    intersected_variants
  }
  return(list(adjusted_susie_weights = adjusted_xqtl_coef, remained_variants_ids = remained_variants_ids))
}

#' Run the SuSiE RSS pipeline
#'
#' Runs SuSiE RSS analysis with one or more SuSiE-family variants. Supports
#' both z+R (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstats Data frame with 'z' or ('beta' and 'se') columns.
#' @param LD_mat LD correlation matrix. Mutually exclusive with X_mat.
#' @param X_mat Genotype matrix (samples x variants). Mutually exclusive with LD_mat.
#' @param n Sample size.
#' @param L Maximum number of causal configurations (default: 30).
#' @param L_greedy Initial greedy number of causal configurations (default: 5).
#' @param analysis_method Iteration mode for the \code{"susie_rss"} fit:
#'   \code{"susie_rss"} (default, normal IBSS), \code{"single_effect"} (L=1,
#'   single iteration), or \code{"bayesian_conditional_regression"}
#'   (full L, single iteration). Only affects the \code{"susie_rss"}
#'   method; ignored for \code{"susie_inf_rss"} and \code{"susie_ash_rss"}.
#' @param methods Optional character vector selecting which RSS variants to
#'   fit. Any subset of \code{c("susie_rss", "susie_inf_rss",
#'   "susie_ash_rss")}. Default \code{NULL} falls back to a single-method fit
#'   driven by \code{analysis_method} (backward-compatible behavior). When
#'   \code{methods} is passed explicitly, each requested method is fitted;
#'   if \code{"susie_inf_rss"} is paired with \code{"susie_rss"} or
#'   \code{"susie_ash_rss"} (or both) and \code{add_susie_inf = TRUE}, the
#'   SuSiE-inf-RSS fit initialises the downstream method. This exposes five
#'   distinct fitting modes mirroring the individual-level pipeline.
#' @param add_susie_inf Logical. When \code{methods} contains
#'   \code{"susie_inf_rss"} alongside \code{"susie_rss"} and/or
#'   \code{"susie_ash_rss"}, controls whether SuSiE-inf-RSS is chained into
#'   the downstream method(s) as initialisation. Default \code{TRUE}.
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
#' @return A list with post-processed SuSiE RSS results. The unified
#'   \code{top_loci} table contains rows from every requested method,
#'   distinguished by the \code{method} column.
#' @importFrom susieR susie_rss
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange select
#' @export
susie_rss_pipeline <- function(sumstats, LD_mat = NULL, X_mat = NULL, n = NULL,
                               L = 30, L_greedy = 5,
                               analysis_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
                               methods = NULL,
                               add_susie_inf = TRUE,
                               coverage = 0.95,
                               secondary_coverage = c(0.7, 0.5),
                               signal_cutoff = 0.1,
                               min_abs_corr = 0.8,
                               R_finite = NULL, R_mismatch = NULL, ...) {
  analysis_method <- match.arg(analysis_method)
  if (is.null(LD_mat) && is.null(X_mat)) stop("Either LD_mat or X_mat must be provided.")
  if (!is.null(LD_mat) && !is.null(X_mat)) stop("Only one of LD_mat or X_mat should be provided, not both.")
  if (!is.null(L_greedy)) L_greedy <- min(L_greedy, L)

  # Resolve effective methods. NULL => legacy single-method via analysis_method.
  valid_rss_methods <- c("susie_rss", "susie_inf_rss", "susie_ash_rss")
  if (is.null(methods)) {
    # Backward-compatible: single fit using analysis_method, labeled accordingly.
    fit_methods <- analysis_method
  } else {
    if (!is.character(methods) || length(methods) == 0L) {
      stop("methods must be a non-empty character vector of method names.")
    }
    bad <- setdiff(methods, valid_rss_methods)
    if (length(bad) > 0) {
      stop("Unknown RSS method(s): ", paste(bad, collapse = ", "),
           ". Valid options: ", paste(valid_rss_methods, collapse = ", "))
    }
    fit_methods <- unique(methods)
  }
  chain_inf_to_susie_rss     <- isTRUE(add_susie_inf) &&
    all(c("susie_inf_rss", "susie_rss") %in% fit_methods)
  chain_inf_to_susie_ash_rss <- isTRUE(add_susie_inf) &&
    all(c("susie_inf_rss", "susie_ash_rss") %in% fit_methods)
  any_chained_init_rss <- chain_inf_to_susie_rss || chain_inf_to_susie_ash_rss

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

  fit_one_susie_rss <- function() {
    if (analysis_method == "single_effect") {
      do.call(susie_rss, c(common, list(L = 1, L_greedy = NULL, max_iter = 1)))
    } else if (analysis_method == "bayesian_conditional_regression") {
      do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy, max_iter = 1)))
    } else {
      do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy)))
    }
  }
  fit_one_susie_inf_rss <- function() {
    do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy,
                                       unmappable_effects = "inf",
                                       convergence_method = "pip",
                                       refine = FALSE, model_init = NULL)))
  }
  fit_one_susie_ash_rss <- function() {
    do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy,
                                       unmappable_effects = "ash",
                                       convergence_method = "pip")))
  }

  fitted_models <- list()
  if ("susie_inf_rss" %in% fit_methods || any_chained_init_rss) {
    inf_fit <- fit_one_susie_inf_rss()
    fitted_models[["susie_inf_rss"]] <- .set_finemapping_fit_class(inf_fit, "susie_inf_rss")
  }
  if ("susie_rss" %in% fit_methods ||
      identical(fit_methods, "single_effect") ||
      identical(fit_methods, "bayesian_conditional_regression")) {
    if (chain_inf_to_susie_rss) {
      chained_args <- prepare_susie_from_inf_args(
        list(L = L, L_greedy = L_greedy),
        fitted_models[["susie_inf_rss"]], refine_default = TRUE,
        unmappable_effects = "none"
      )
      rss_fit <- do.call(susie_rss, c(common, chained_args))
    } else {
      rss_fit <- fit_one_susie_rss()
    }
    # Label by analysis_method when in legacy single-method mode, else "susie_rss"
    rss_label <- if (is.null(methods)) analysis_method else "susie_rss"
    fitted_models[[rss_label]] <- .set_finemapping_fit_class(rss_fit, rss_label)
  }
  if ("susie_ash_rss" %in% fit_methods) {
    if (chain_inf_to_susie_ash_rss) {
      chained_args <- prepare_susie_from_inf_args(
        list(L = L, L_greedy = L_greedy),
        fitted_models[["susie_inf_rss"]], refine_default = NULL,
        unmappable_effects = "ash"
      )
      ash_fit <- do.call(susie_rss, c(common, chained_args))
    } else {
      ash_fit <- fit_one_susie_ash_rss()
    }
    fitted_models[["susie_ash_rss"]] <- .set_finemapping_fit_class(ash_fit, "susie_ash_rss")
  }

  # Drop SuSiE-inf-RSS from post-processing if it was only fit for init
  if (any_chained_init_rss && !("susie_inf_rss" %in% fit_methods)) {
    fitted_models[["susie_inf_rss"]] <- NULL
  }

  # For post-processing, pass genotype matrix X directly when available.
  if (!is.null(LD_mat)) {
    data_x <- LD_mat
    pp_cs_input <- "Xcorr"
  } else if (is.list(X_mat) && !is.matrix(X_mat)) {
    data_x <- do.call(rbind, X_mat)[, seq_along(z), drop = FALSE]
    pp_cs_input <- "X"
  } else {
    data_x <- X_mat[, seq_along(z), drop = FALSE]
    pp_cs_input <- "X"
  }

  post <- postprocess_finemapping_fits(
    fits = fitted_models,
    data_x = data_x,
    data_y = list(z = z),
    coverage = coverage,
    secondary_coverage = secondary_coverage,
    signal_cutoff = signal_cutoff,
    min_abs_corr = min_abs_corr,
    cs_input = pp_cs_input
  )
  # Primary method preference: "susie_rss" > other names > first fit
  primary <- if ("susie_rss" %in% names(fitted_models)) "susie_rss" else names(fitted_models)[1]
  format_finemapping_output(post, primary_method = primary)
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
