#' @title Heritability Estimation Entry Points and Converters
#' @description Top-level entry point for heritability estimation,
#'   LD score computation methods, H2Estimate accessors, and a converter
#'   to bridge H2Estimate into the sldsc_wrapper.R postprocessing pipeline.
#' @include AllGenerics.R
#' @importFrom stats median
NULL

# =============================================================================
# estimateH2 — main dispatch
# =============================================================================

#' @rdname estimateH2
#' @export
setMethod("estimateH2",
  signature(sumstats = "GWASSumStats", ld_ref = "LDStatistic"),
  function(sumstats, ld_ref, method = "lder", annotations = NULL,
           local = FALSE, ...) {

    method <- match.arg(method, c("lder", "gldsc", "hdl"))
    .validate_method_ref(method, ld_ref)

    z <- getZ(sumstats)
    n <- median(getN(sumstats))
    M <- nSnps(sumstats)

    # Apply var_y correction for case-control studies
    var_y <- getVarY(sumstats)
    if (!is.null(var_y)) {
      # Adjust effective N for liability scale: n_eff = n / var_y
      n <- n / var_y
    }

    # Dispatch to method-specific function
    result <- switch(method,
      "lder" = lder_univariate(z, n, ld_ref, annotations, local, ...),
      "gldsc" = gldsc_univariate(z, n, ld_ref, annotations, local, ...),
      "hdl" = hdl_univariate(z, n, ld_ref, annotations, local, ...)
    )

    # Wrap into H2Estimate S4 object
    new("H2Estimate",
      h2 = result$h2,
      h2_se = result$h2_se,
      intercept = result$intercept %||% NA_real_,
      intercept_se = result$intercept_se %||% NA_real_,
      local = result$local,
      enrichment = result$enrichment,
      tau_blocks = result$tau_blocks,
      score_stats = result$score_stats,
      method = method,
      n_snps = as.integer(M),
      trait_name = sumstats@trait_name
    )
  }
)

#' @keywords internal
.validate_method_ref <- function(method, ld_ref) {
  if (method %in% c("lder", "hdl") && !is(ld_ref, "LDEigen")) {
    stop("Method '", method, "' requires an LDEigen object, ",
         "got ", class(ld_ref))
  }
  if (method == "gldsc" && !is(ld_ref, "LDScore")) {
    stop("Method 'gldsc' requires an LDScore object, ",
         "got ", class(ld_ref))
  }
  invisible(TRUE)
}

# =============================================================================
# computeLdScores — LD score computation
# =============================================================================

#' @rdname computeLdScores
#' @export
setMethod("computeLdScores",
  signature(ld_ref = "LDEigen"),
  function(ld_ref, annotations = NULL, ...) {
    # Reconstruct LD scores from eigendecompositions
    # l2[j] = sum_k r^2_{jk} = sum_b sum_{eigenvalues in b} V[j,.]^2 * d
    n_snps <- nrow(ld_ref@snp_info)

    if (is.null(annotations)) {
      # Base LD scores only
      l2 <- numeric(n_snps)
      for (b in seq_along(ld_ref@eigen_list)) {
        block <- ld_ref@eigen_list[[b]]
        idx <- block$snp_idx
        V <- block$vectors
        d <- block$values
        # LD score for SNP j = sum_i V[j,i]^2 * d[i]^2
        # (since R = V D V', R^2_{jk} = sum_i V[j,i]^2 * d[i]^2 * V[k,i]^2)
        # Simplified: l2[j] = sum_i (V[j,i] * d[i])^2
        Vd <- sweep(V, 2, d, "*")
        l2[idx] <- rowSums(Vd^2)
      }
      return(matrix(l2, ncol = 1,
                     dimnames = list(NULL, "base_l2")))
    }

    # Annotation-stratified LD scores
    annot_mat <- annotations@annotations
    n_annot <- ncol(annot_mat)
    # Base + annotation-stratified columns
    l2_strat <- matrix(0, nrow = n_snps, ncol = 1 + n_annot)

    for (b in seq_along(ld_ref@eigen_list)) {
      block <- ld_ref@eigen_list[[b]]
      idx <- block$snp_idx
      V <- block$vectors
      d <- block$values

      # Base LD score
      Vd <- sweep(V, 2, d, "*")
      l2_strat[idx, 1] <- rowSums(Vd^2)

      # Annotation-stratified: l2_a[j] = sum_k r^2_{jk} * annot[k,a]
      # Using eigendecomposition: l2_a[j] = sum_i V[j,i]^2 * d[i]^2 *
      #   (sum_k V[k,i]^2 * annot[k,a])
      for (a in seq_len(n_annot)) {
        annot_col <- annot_mat[idx, a]
        # For each eigenvalue: weight = sum_k V[k,i]^2 * annot[k,a]
        annot_weights <- as.vector(crossprod(V^2, annot_col))
        l2_strat[idx, 1 + a] <- as.vector(Vd^2 %*% annot_weights)
      }
    }

    col_names <- c("base_l2", annotations@annotation_meta$name)
    colnames(l2_strat) <- col_names
    l2_strat
  }
)

#' @rdname computeLdScores
#' @export
setMethod("computeLdScores",
  signature(ld_ref = "LDScore"),
  function(ld_ref, annotations = NULL, ...) {
    if (is.null(annotations)) {
      return(ld_ref@ld_scores)
    }

    # Compute annotation-stratified LD scores using LD matrices
    if (length(ld_ref@ld_matrix_list) == 0) {
      stop("Annotation-stratified LD scores require ld_matrix_list in LDScore. ",
           "Recompute the LD reference with full LD matrices.")
    }

    n_snps <- nrow(ld_ref@snp_info)
    annot_mat <- annotations@annotations
    n_annot <- ncol(annot_mat)

    # Base L2 + annotation-stratified columns
    l2_strat <- matrix(0, nrow = n_snps, ncol = 1 + n_annot)
    l2_strat[, 1] <- ld_ref@ld_scores[, 1]

    for (b in seq_along(ld_ref@ld_matrix_list)) {
      block <- ld_ref@ld_matrix_list[[b]]
      R <- block$R
      idx <- block$snp_idx
      R2 <- R^2
      for (a in seq_len(n_annot)) {
        # l2_a[j] = sum_k R^2_{jk} * annot[k, a]
        l2_strat[idx, 1 + a] <- as.vector(R2 %*% annot_mat[idx, a])
      }
    }

    col_names <- c("base_l2", annotations@annotation_meta$name)
    colnames(l2_strat) <- col_names
    l2_strat
  }
)

# =============================================================================
# H2Estimate accessor methods
# =============================================================================

#' @rdname getLocal
#' @export
setMethod("getLocal", "H2Estimate", function(object) {
  object@local
})

#' @rdname getEnrichment
#' @export
setMethod("getEnrichment", "H2Estimate", function(object) {
  object@enrichment
})

#' @rdname getScoreStats
#' @export
setMethod("getScoreStats", "H2Estimate", function(object) {
  object@score_stats
})

# =============================================================================
# Converter: H2Estimate -> sldsc_wrapper list format
# =============================================================================

#' @title Convert H2Estimate to S-LDSC Trait Format
#' @description Convert an \code{H2Estimate} object into the list format
#'   expected by \code{\link{standardize_sldsc_trait}} and
#'   \code{\link{meta_sldsc_random}}. This bridges the h2 estimation
#'   methods (LDER, gLDSC, HDL) into the sldsc_wrapper.R postprocessing
#'   pipeline.
#' @param h2_est An \code{H2Estimate} object with enrichment and tau_blocks.
#' @return A named list matching the format of \code{\link{read_sldsc_trait}}:
#'   \describe{
#'     \item{categories}{Character vector of annotation names}
#'     \item{tau}{Named numeric vector of per-annotation coefficients}
#'     \item{tau_se}{Named numeric vector of tau standard errors}
#'     \item{enrichment}{Named numeric vector of enrichment ratios}
#'     \item{enrichment_se}{Named numeric vector of enrichment SEs}
#'     \item{enrichment_p}{Named numeric vector of enrichment p-values}
#'     \item{prop_h2}{Named numeric vector of proportion of h2}
#'     \item{prop_snps}{Named numeric vector of proportion of SNPs}
#'     \item{h2g}{Numeric scalar, global h2 estimate}
#'     \item{tau_blocks}{Matrix (n_blocks x n_categories) for jackknife}
#'     \item{n_blocks}{Integer, number of jackknife blocks}
#'   }
#' @export
h2estimate_to_sldsc_trait <- function(h2_est) {
  if (!is(h2_est, "H2Estimate")) {
    stop("h2_est must be an H2Estimate object")
  }

  enrich_df <- h2_est@enrichment
  if (is.null(enrich_df)) {
    stop("H2Estimate has no enrichment results. ",
         "Run estimateH2 with annotations to get enrichment estimates.")
  }

  cats <- as.character(enrich_df$annotation)
  n_cats <- length(cats)

  tau_blocks <- h2_est@tau_blocks
  if (is.null(tau_blocks)) {
    # Create a dummy single-block matrix from the point estimates
    tau_blocks <- matrix(enrich_df$tau, nrow = 1)
    colnames(tau_blocks) <- cats
    n_blocks <- 1L
  } else {
    n_blocks <- nrow(tau_blocks)
    if (is.null(colnames(tau_blocks))) {
      colnames(tau_blocks) <- cats
    }
  }

  list(
    categories    = cats,
    tau           = setNames(enrich_df$tau, cats),
    tau_se        = setNames(enrich_df$tau_se, cats),
    enrichment    = setNames(enrich_df$enrichment, cats),
    enrichment_se = setNames(enrich_df$enrichment_se, cats),
    enrichment_p  = setNames(enrich_df$enrichment_p, cats),
    prop_h2       = setNames(enrich_df$prop_h2, cats),
    prop_snps     = setNames(enrich_df$prop_snps, cats),
    h2g           = h2_est@h2,
    tau_blocks    = tau_blocks,
    n_blocks      = n_blocks
  )
}
