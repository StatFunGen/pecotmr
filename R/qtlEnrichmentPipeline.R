#' @title QTL Enrichment Pipeline (Genome-Wide)
#' @description Genome-wide pipeline that computes per-pair (GWAS study,
#'   QTL context) enrichment estimates by passing the GWAS PIP vector
#'   and the QTL credible-set posteriors to
#'   \code{\link{qtlEnrichment}}. The returned table feeds
#'   \code{\link{colocPipeline}} via its \code{enrichment} argument.
#'
#'   \strong{Not gene-parallelisable}: the enrichment estimator runs
#'   over the full genome of GWAS PIPs and the full collection of QTL
#'   fits at once.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{gwasFineMappingResult}: a genome-wide
#'     \code{\link{GwasFineMappingResult}} (one row per (study, LD
#'     block) tuple). Each entry's \code{FineMappingEntry$trimmedFit}
#'     must carry a \code{pip} vector.
#'   \item \code{qtlFineMappingResult}: the genome-wide
#'     \code{\link{QtlFineMappingResult}}. Each entry's
#'     \code{trimmedFit} must carry \code{alpha}, \code{pip}, and
#'     prior-variance fields (\code{V}).
#' }
#'
#' @section LD-sketch identity check:
#' The GWAS \code{\link{FineMappingResultBase}} must have a non-NULL
#' \code{ldSketch} (RSS-derived). If the QTL FMR also has a non-NULL
#' \code{ldSketch}, the two must match exactly. When the QTL FMR's
#' \code{ldSketch} is NULL (individual-level QTL fit), validation is
#' skipped on the QTL side.
#'
#' @param gwasFineMappingResult See above.
#' @param qtlFineMappingResult See above.
#' @param numGwas Number of GWAS variants used to estimate \code{piGwas}.
#'   When \code{NULL} (default) it is estimated from the data — bias
#'   warning applies if the input PIP vector is not genome-wide.
#' @param piQtl Per-variant prior of being a QTL causal variant.
#'   \code{NULL} (default) estimates from the data.
#' @param lambda Shrinkage parameter for the enrichment estimator.
#'   Default \code{1.0}.
#' @param impN Number of imputed samples used by the estimator. Default
#'   \code{25}.
#' @param numThreads Number of threads used by
#'   \code{qtlEnrichment}. Default \code{1}.
#' @param ... Additional arguments forwarded to
#'   \code{\link{qtlEnrichment}}.
#' @return A data frame with one row per (gwasStudy, qtlContext) pair
#'   and columns \code{gwasStudy}, \code{qtlContext},
#'   \code{enrichment}, \code{enrichmentSe}, \code{enrichmentLogOdds},
#'   plus any extras the underlying estimator emits. Suitable as the
#'   \code{enrichment} argument to \code{\link{colocPipeline}}.
#' @export
qtlEnrichmentPipeline <- function(gwasFineMappingResult,
                                  qtlFineMappingResult,
                                  numGwas = NULL,
                                  piQtl = NULL,
                                  lambda = 1.0,
                                  impN = 25,
                                  numThreads = 1L,
                                  ...) {
  if (!methods::is(gwasFineMappingResult, "GwasFineMappingResult"))
    stop("`gwasFineMappingResult` must be a GwasFineMappingResult.")
  if (!methods::is(qtlFineMappingResult, "QtlFineMappingResult"))
    stop("`qtlFineMappingResult` must be a QtlFineMappingResult.")

  gwasLd <- getLdSketch(gwasFineMappingResult)
  if (is.null(gwasLd))
    stop("qtlEnrichmentPipeline: the GWAS FineMappingResult must have a ",
         "non-NULL ldSketch (it should be RSS-derived).")
  qtlLd <- getLdSketch(qtlFineMappingResult)
  .colocRequireMatchingLdSketches(qtlLd, gwasLd)

  # Per-study genome-wide GWAS PIP vector (named by variant id).
  gwasStudies <- unique(as.character(gwasFineMappingResult$study))
  qtlContexts <- unique(as.character(qtlFineMappingResult$context))

  if (length(gwasStudies) == 0L || length(qtlContexts) == 0L)
    stop("qtlEnrichmentPipeline: no (gwasStudy, qtlContext) pairs ",
         "to compute (one of the inputs has zero rows).")

  results <- list()
  for (gStudy in gwasStudies) {
    gwasPip <- .enrBuildGwasPipVector(gwasFineMappingResult, gStudy)
    if (length(gwasPip) == 0L) {
      warning(sprintf(
        "qtlEnrichmentPipeline: no usable PIPs for gwasStudy='%s'; skipping.",
        gStudy))
      next
    }
    for (qContext in qtlContexts) {
      qtlRegions <- .enrBuildQtlRegionsList(qtlFineMappingResult, qContext)
      if (length(qtlRegions) == 0L) {
        warning(sprintf(
          "qtlEnrichmentPipeline: no usable QTL regions for qtlContext='%s'; skipping.",
          qContext))
        next
      }
      enr <- tryCatch(
        qtlEnrichment(
          gwasPip          = gwasPip,
          susieQtlRegions  = qtlRegions,
          numGwas          = numGwas,
          piQtl            = piQtl,
          lambda           = lambda,
          impN             = impN,
          numThreads       = numThreads,
          ...),
        error = function(e) {
          warning(sprintf(
            "qtlEnrichmentPipeline: qtlEnrichment failed for (gwasStudy='%s', qtlContext='%s'): %s",
            gStudy, qContext, conditionMessage(e)))
          NULL
        })
      if (is.null(enr)) next
      row <- .enrFlattenEnrichment(enr)
      row$gwasStudy  <- gStudy
      row$qtlContext <- qContext
      results[[length(results) + 1L]] <- row
    }
  }

  if (length(results) == 0L) {
    return(data.frame(
      gwasStudy  = character(0),
      qtlContext = character(0),
      enrichment = numeric(0),
      enrichmentSe = numeric(0),
      enrichmentLogOdds = numeric(0),
      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, lapply(results, as.data.frame,
                              stringsAsFactors = FALSE))
  rownames(out) <- NULL
  idCols <- c("gwasStudy", "qtlContext")
  other  <- setdiff(colnames(out), idCols)
  out[, c(idCols, other), drop = FALSE]
}

# =============================================================================
# Internal helpers
# =============================================================================

# Build a named GWAS PIP vector for one study. Walks every row of the
# GwasFineMappingResult tagged with that study, extracts the per-row
# pip from each FineMappingEntry, and concatenates with variant-id
# names. Errors if any single variant appears with conflicting PIP
# values across rows.
# @noRd
.enrBuildGwasPipVector <- function(gwasFmr, gStudy) {
  idx <- which(as.character(gwasFmr$study) == gStudy)
  if (length(idx) == 0L) return(numeric(0))
  pieces <- list()
  for (i in idx) {
    entry <- gwasFmr$entry[[i]]
    fit <- getTrimmedFit(entry)
    if (is.null(fit) || is.null(fit$pip)) next
    pip <- as.numeric(fit$pip)
    ids <- if (!is.null(names(fit$pip))) names(fit$pip)
           else getVariantIds(entry)
    if (length(ids) != length(pip)) next
    pieces[[length(pieces) + 1L]] <-
      stats::setNames(pip, as.character(ids))
  }
  if (length(pieces) == 0L) return(numeric(0))
  all <- unlist(pieces)
  uniqIds <- unique(names(all))
  if (length(uniqIds) != length(all)) {
    # Same variant id in multiple blocks. Verify the values agree.
    dedup <- tapply(all, names(all), function(v) {
      if (length(unique(round(v, 12))) > 1L) {
        stop("qtlEnrichmentPipeline: variant '", names(v)[[1L]],
             "' appears with conflicting PIPs across GWAS blocks for ",
             "study '", "?", "'; the GWAS fine-mapping must produce a ",
             "consistent PIP per variant.")
      }
      v[[1L]]
    })
    all <- as.numeric(dedup)
    names(all) <- names(dedup)
  }
  all
}

# Build the per-(qtl context) list of region fits in the shape that
# qtlEnrichment expects: list(d) where each d carries
# alpha, pip, prior_variance (V).
# @noRd
.enrBuildQtlRegionsList <- function(qtlFmr, qContext) {
  idx <- which(as.character(qtlFmr$context) == qContext)
  if (length(idx) == 0L) return(list())
  out <- list()
  for (i in idx) {
    entry <- qtlFmr$entry[[i]]
    fit <- getTrimmedFit(entry)
    if (is.null(fit) || is.null(fit$alpha) || is.null(fit$pip)) next
    pV <- if (!is.null(fit$V)) fit$V
          else if (!is.null(fit$prior_variance)) fit$prior_variance
          else NULL
    if (is.null(pV)) next
    if (is.null(names(fit$pip)))
      names(fit$pip) <- getVariantIds(entry)
    out[[length(out) + 1L]] <- list(
      alpha          = fit$alpha,
      pip            = fit$pip,
      prior_variance = pV)
  }
  out
}

# Coerce qtlEnrichment's variable-shape output into a single-row
# named list with the canonical columns the caller documents. The
# underlying estimator returns either a list with named numeric scalars
# (enrichment, enrichmentSe, enrichmentLogOdds, ...) or a matrix/df —
# this helper handles both.
# @noRd
.enrFlattenEnrichment <- function(enr) {
  if (is.list(enr) && is.null(dim(enr))) {
    pickScalar <- function(field) {
      v <- enr[[field]]
      if (is.null(v)) NA_real_
      else as.numeric(v[[1L]])
    }
    list(
      enrichment        = pickScalar("enrichment"),
      enrichmentSe      = pickScalar("enrichmentSe"),
      enrichmentLogOdds = pickScalar("enrichmentLogOdds"))
  } else if (is.matrix(enr) || is.data.frame(enr)) {
    df <- as.data.frame(enr, stringsAsFactors = FALSE)
    if (nrow(df) == 0L) {
      list(enrichment = NA_real_, enrichmentSe = NA_real_,
           enrichmentLogOdds = NA_real_)
    } else {
      list(
        enrichment        = .enrPickColumn(df, c("enrichment", "Enrichment")),
        enrichmentSe      = .enrPickColumn(df, c("enrichmentSe", "se", "stderr")),
        enrichmentLogOdds = .enrPickColumn(df, c("enrichmentLogOdds", "logOdds", "log_odds")))
    }
  } else {
    list(enrichment = NA_real_, enrichmentSe = NA_real_,
         enrichmentLogOdds = NA_real_)
  }
}

# @noRd
.enrPickColumn <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) return(NA_real_)
  as.numeric(df[[hit[[1L]]]][[1L]])
}


# =============================================================================
# qtlEnrichment: low-level enrichment estimation
# -----------------------------------------------------------------------------
# Per-(GWAS, QTL-region-list) enrichment estimator. Called per-(gwasStudy,
# qtlContext) pair by qtlEnrichmentPipeline above. Uses the fastenloc-style
# C++ kernel (qtlEnrichmentRcpp) under the hood.
# =============================================================================
#' @title Implementation of enrichment analysis described in https://doi.org/10.1371/journal.pgen.1006646
#'
#' @description Largely follows from fastenloc https://github.com/xqwen/fastenloc
#' but uses `susieR` fitted objects as input to estimate prior for use with `coloc` package (coloc v5, aka SuSiE-coloc).
#' The main differences are 1) now enrichment is based on all QTL variants whether or not they are inside signal clusters;
#' 2) Causal QTL are sampled from SuSiE single effects, not signal clusters;
#' 3) Allow a variant to be QTL for not only multiple conditions (eg cell types) but also multiple regions (eg genes).
#' Other minor improvements include 1) Make GSL RNG thread-safe; 2) Release memory from QTL binary annotation samples immediately after they are used.
#' @details Uses output of \code{\link[susieR]{susie}} from the
#'   \code{susieR} package.
#'
#' @param gwasPip This is a vector of GWAS PIP, genome-wide.
#' @param susieQtlRegions This is a list of SuSiE fitted objects per QTL unit analyzed
#' @param numGwas This parameter is highly important if GWAS input does not contain all SNPs interrogated (e.g., in some cases, only fine-mapped geomic regions are included).
#' Then users must pick a value of total_variants and estimate piGwas beforehand by: sum(gwasPip$pip)/numGwas. If numGwas is null, piGwas would be sum(gwasPip$pip)/total_variants.
#' @param piQtl This parameter can be safely left to default if your input QTL data has enough regions to estimate it.
#' @param lambda Similar to the shrinkage parameter used in ridge regression. It takes any non-negative value and shrinks the enrichment estimate towards 0.
#' When it is set to 0, no shrinkage will be applied. A large value indicates strong shrinkage. The default value is set to 1.0.
#' @param impN Rounds of multiple imputation to draw QTL from, default is 25.
#' @param numThreads Number of Simultaneous running CPU threads for multiple imputation, default is 1.
#' @return A list of enrichment parameter estimates
#'
#' @examples
#'
#' # Simulate fake data for gwasPip
#' nGwasPip <- 1000
#' gwasPip <- runif(nGwasPip)
#' names(gwasPip) <- paste0("snp", 1:nGwasPip)
#' gwasFit <- list(pip = gwasPip)
#' # Simulate fake data for a single SuSiEFit object
#' simulateSusiefit <- function(n, p) {
#'   pip <- runif(n)
#'   names(pip) <- paste0("snp", 1:n)
#'   alpha <- t(matrix(runif(n * p), nrow = n))
#'   alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
#'   list(
#'     pip = pip,
#'     alpha = alpha,
#'     prior_variance = runif(p)
#'   )
#' }
#' # Simulate multiple SuSiEFit objects
#' nSusieFits <- 2
#' susieFits <- replicate(nSusieFits, simulateSusiefit(nGwasPip, 10), simplify = FALSE)
#' # Add these fits to a list, providing names to each element
#' names(susieFits) <- paste0("fit", 1:length(susieFits))
#' # Set other parameters
#' impN <- 10
#' lambda <- 1
#' numThreads <- 1
#' library(pecotmr)
#' en <- qtlEnrichment(gwasFit, susieFits, lambda = lambda, impN = impN, numThreads = numThreads)
#'
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib pecotmr, .registration = TRUE
#' @export
#'
qtlEnrichment <- function(gwasPip, susieQtlRegions,
                                 numGwas = NULL, piQtl = NULL,
                                 lambda = 1.0, impN = 25,
                                 doubleShrinkage = FALSE,
                                 besselCorrection = TRUE,
                                 numThreads = 1, verbose = TRUE) {
  if (is.null(numGwas)) {
    warning("numGwas is not provided. Estimating piGwas from the data. Note that this estimate may be biased if the input gwasPip does not contain genome-wide variants.")
    piGwas <- sum(gwasPip) / length(gwasPip)
    if (verbose) {
      message(paste("Estimated piGwas: ", round(piGwas, 5), "\n"))
    }
  } else {
    piGwas <- sum(gwasPip) / numGwas
  }

  if (is.null(piQtl)) {
    warning("piQtl is not provided. Estimating piQtl from the data. Note that this estimate may be biased if either 1) the input susieQtlRegions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
    numSignal <- 0
    numTest <- 0
    for (d in susieQtlRegions) {
      numSignal <- numSignal + sum(d$pip)
      numTest <- numTest + length(d$pip)
    }
    piQtl <- numSignal / numTest
    if (verbose) {
      message(paste("Estimated piQtl: ", round(piQtl, 5), "\n"))
    }
  }

  if (piGwas == 0) stop("Cannot perform enrichment analysis. No association signal found in GWAS data.")
  if (piQtl == 0) stop("Cannot perform enrichment analysis. No QTL associated with the molecular phenotype.")

  # Check if names of gwasPip and susieQtlRegions$pip are both available
  if (is.null(names(gwasPip))) {
    stop("Variant names are missing in gwasPip. Please provide named gwasPip data.")
  }
  if (!all(sapply(susieQtlRegions, function(x) !is.null(names(x$pip))))) {
    stop("Variant names are missing in susieQtlRegions$pip. Please provide susieQtlRegions with named pip data.")
  }

  # Align the names of susieQtlRegions$pip to gwasPip names and document unmatched variants
  alignedSusieQtlRegions <- lapply(susieQtlRegions, function(x) {
    alignmentResult <- alignVariantNames(names(x$pip), names(gwasPip))
    names(x$pip) <- alignmentResult$alignedVariants
    if (length(alignmentResult$unmatchedIndices) > 0) {
      x$unmatched_variants <- names(x$pip)[alignmentResult$unmatchedIndices]
    }
    x
  })
  unmatchedVariants <- lapply(alignedSusieQtlRegions, function(x) x$unmatched_variants)

  # Update susieQtlRegions with the aligned variant names
  susieQtlRegions <- lapply(alignedSusieQtlRegions, function(x) {
    x$unmatched_variants <- NULL
    x
  })

  # cpp11 requires exact integer types for int parameters
  en <- qtlEnrichmentRcpp(
    rGwasPip = gwasPip,
    rQtlSusieFit = susieQtlRegions,
    piGwas = piGwas,
    piQtl = piQtl,
    ImpN = as.integer(impN),
    shrinkageLambda = lambda,
    doubleShrinkage = doubleShrinkage,
    besselCorrection = besselCorrection,
    numThreads = as.integer(numThreads)
  )

  # Add the unmatched variants to the output
  en <- list(en)
  en$unused_xqtl_variants <- unmatchedVariants

  return(en)
}

