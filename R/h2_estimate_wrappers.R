#' @title Heritability Estimation Entry Points and Converters
#' @description Top-level entry point for heritability estimation,
#'   LD score computation methods, H2Estimate accessors, and a converter
#'   to bridge H2Estimate into the sldsc_wrapper.R postprocessing pipeline.
#' @name pecotmr-h2-wrappers
#' @keywords internal
#' @include AllGenerics.R
#' @importFrom stats median
NULL

# =============================================================================
# estimateH2 — main dispatch
# =============================================================================

#' @rdname estimateH2
#' @export
setMethod("estimateH2",
  signature(sumstats = "GwasSumStats", ldRef = "LdStatistic"),
  function(sumstats, ldRef, method = "lder", annotations = NULL,
           local = FALSE, ...) {

    method <- match.arg(method, c("lder", "gldsc", "hdl"))
    .validateMethodRef(method, ldRef)

    z <- getZ(sumstats)
    n <- median(getN(sumstats))
    M <- nSnps(sumstats)

    # Apply the legacy heritability-wrapper correction. This is separate from
    # the SuSiE RSS binaryTraitModel handling in the fine-mapping pipeline.
    varY <- getVarY(sumstats)
    if (!is.null(varY)) {
      n <- n / varY
    }

    # Dispatch to method-specific function
    result <- switch(method,
      "lder" = lderUnivariate(z, n, ldRef, annotations, local, ...),
      "gldsc" = gldscUnivariate(z, n, ldRef, annotations, local, ...),
      "hdl" = hdlUnivariate(z, n, ldRef, annotations, local, ...)
    )

    # Wrap into H2Estimate S4 object
    new("H2Estimate",
      h2 = result$h2,
      h2Se = result$h2_se,
      intercept = result$intercept %||% NA_real_,
      interceptSe = result$intercept_se %||% NA_real_,
      local = result$local,
      enrichment = result$enrichment,
      tauBlocks = result$tau_blocks,
      scoreStats = result$score_stats,
      method = method,
      nSnps = as.integer(M),
      traitName = sumstats@traitName
    )
  }
)

#' @keywords internal
.validateMethodRef <- function(method, ldRef) {
  if (method %in% c("lder", "hdl") && !is(ldRef, "LdEigen")) {
    stop("Method '", method, "' requires an LdEigen object, ",
         "got ", class(ldRef))
  }
  if (method == "gldsc" && !is(ldRef, "LdScore")) {
    stop("Method 'gldsc' requires an LdScore object, ",
         "got ", class(ldRef))
  }
  invisible(TRUE)
}

# =============================================================================
# computeLdScores — LD score computation
# =============================================================================

#' @rdname computeLdScores
#' @export
setMethod("computeLdScores",
  signature(ldRef = "LdEigen"),
  function(ldRef, annotations = NULL, ...) {
    # Reconstruct LD scores from eigendecompositions
    # l2[j] = sum_k r^2_{jk} = sum_b sum_{eigenvalues in b} V[j,.]^2 * d
    nSnps <- nrow(ldRef@snpInfo)

    if (is.null(annotations)) {
      # Base LD scores only
      l2 <- numeric(nSnps)
      for (b in seq_along(ldRef@eigenList)) {
        block <- ldRef@eigenList[[b]]
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
    annotMat <- annotations@annotations
    nAnnot <- ncol(annotMat)
    # Base + annotation-stratified columns
    l2Strat <- matrix(0, nrow = nSnps, ncol = 1 + nAnnot)

    for (b in seq_along(ldRef@eigenList)) {
      block <- ldRef@eigenList[[b]]
      idx <- block$snp_idx
      V <- block$vectors
      d <- block$values

      # Base LD score
      Vd <- sweep(V, 2, d, "*")
      l2Strat[idx, 1] <- rowSums(Vd^2)

      # Annotation-stratified: l2_a[j] = sum_k r^2_{jk} * annot[k,a]
      # Using eigendecomposition: l2_a[j] = sum_i V[j,i]^2 * d[i]^2 *
      #   (sum_k V[k,i]^2 * annot[k,a])
      for (a in seq_len(nAnnot)) {
        annotCol <- annotMat[idx, a]
        # For each eigenvalue: weight = sum_k V[k,i]^2 * annot[k,a]
        annotWeights <- as.vector(crossprod(V^2, annotCol))
        l2Strat[idx, 1 + a] <- as.vector(Vd^2 %*% annotWeights)
      }
    }

    colNames <- c("base_l2", annotations@annotationMeta$name)
    colnames(l2Strat) <- colNames
    l2Strat
  }
)

#' @rdname computeLdScores
#' @export
setMethod("computeLdScores",
  signature(ldRef = "LdScore"),
  function(ldRef, annotations = NULL, ...) {
    if (is.null(annotations)) {
      return(ldRef@ldScores)
    }

    # Compute annotation-stratified LD scores using LD matrices
    if (length(ldRef@ldMatrixList) == 0) {
      stop("Annotation-stratified LD scores require ldMatrixList in LdScore. ",
           "Recompute the LD reference with full LD matrices.")
    }

    nSnps <- nrow(ldRef@snpInfo)
    annotMat <- annotations@annotations
    nAnnot <- ncol(annotMat)

    # Base L2 + annotation-stratified columns
    l2Strat <- matrix(0, nrow = nSnps, ncol = 1 + nAnnot)
    l2Strat[, 1] <- ldRef@ldScores[, 1]

    for (b in seq_along(ldRef@ldMatrixList)) {
      block <- ldRef@ldMatrixList[[b]]
      R <- block$R
      idx <- block$snp_idx
      R2 <- R^2
      for (a in seq_len(nAnnot)) {
        # l2_a[j] = sum_k R^2_{jk} * annot[k, a]
        l2Strat[idx, 1 + a] <- as.vector(R2 %*% annotMat[idx, a])
      }
    }

    colNames <- c("base_l2", annotations@annotationMeta$name)
    colnames(l2Strat) <- colNames
    l2Strat
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
  object@scoreStats
})

# =============================================================================
# Converter: H2Estimate -> sldsc_wrapper list format
# =============================================================================

#' @title Convert H2Estimate to S-LDSC Trait Format
#' @description Convert an \code{H2Estimate} object into the list format
#'   expected by \code{\link{standardizeSldscTrait}} and
#'   \code{\link{metaSldscRandom}}. This bridges the h2 estimation
#'   methods (LDER, gLDSC, HDL) into the sldsc_wrapper.R postprocessing
#'   pipeline.
#' @param h2Est An \code{H2Estimate} object with enrichment and tauBlocks.
#' @return A named list matching the format of \code{\link{readSldscTrait}}:
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
#'     \item{tau_blocks}{Matrix (nBlocks x nCategories) for jackknife}
#'     \item{n_blocks}{Integer, number of jackknife blocks}
#'   }
#' @export
h2EstimateToSldscTrait <- function(h2Est) {
  if (!is(h2Est, "H2Estimate")) {
    stop("h2Est must be an H2Estimate object")
  }

  enrichDf <- h2Est@enrichment
  if (is.null(enrichDf)) {
    stop("H2Estimate has no enrichment results. ",
         "Run estimateH2 with annotations to get enrichment estimates.")
  }

  cats <- as.character(enrichDf$annotation)
  nCats <- length(cats)

  tauBlocks <- h2Est@tauBlocks
  if (is.null(tauBlocks)) {
    # Create a dummy single-block matrix from the point estimates
    tauBlocks <- matrix(enrichDf$tau, nrow = 1)
    colnames(tauBlocks) <- cats
    nBlocks <- 1L
  } else {
    nBlocks <- nrow(tauBlocks)
    if (is.null(colnames(tauBlocks))) {
      colnames(tauBlocks) <- cats
    }
  }

  list(
    categories    = cats,
    tau           = setNames(enrichDf$tau, cats),
    tau_se        = setNames(enrichDf$tau_se, cats),
    enrichment    = setNames(enrichDf$enrichment, cats),
    enrichment_se = setNames(enrichDf$enrichment_se, cats),
    enrichment_p  = setNames(enrichDf$enrichment_p, cats),
    prop_h2       = setNames(enrichDf$prop_h2, cats),
    prop_snps     = setNames(enrichDf$prop_snps, cats),
    h2g           = h2Est@h2,
    tau_blocks    = tauBlocks,
    n_blocks      = nBlocks
  )
}
