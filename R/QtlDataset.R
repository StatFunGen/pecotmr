# =============================================================================
# QtlDataset S4 class
# -----------------------------------------------------------------------------
# Single-study individual-level QTL container: a GenotypeHandle, a named
# list of per-context phenotype SummarizedExperiments, optional genotype
# covariates, and constructor-level QC knobs (mafCutoff, macCutoff,
# xvarCutoff, imissCutoff, keepSamples, keepVariants). Backed by lazy
# getGenotypes / getResidualizedGenotypes accessors that apply QC at
# extraction time. The entry point for individual-level fine-mapping
# (fineMappingPipeline), TWAS weight learning (twasWeightsPipeline), and
# multi-study composition (MultiStudyQtlDataset).
# =============================================================================

#' @include allGenerics.R
NULL

setClass("QtlDataset",
  representation(
    study              = "character",
    genotypes          = "GenotypeHandle",
    phenotypes         = "list",
    genotypeCovariates = "matrix",
    scaleResiduals     = "logical",
    mafCutoff          = "numeric",
    macCutoff          = "numeric",
    xvarCutoff         = "numeric",
    imissCutoff        = "numeric",
    keepSamples        = "character",
    keepVariants       = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@study) != 1L || !nzchar(object@study))
      errors <- c(errors, "'study' must be a single non-empty character string")
    if (length(object@scaleResiduals) != 1L)
      errors <- c(errors, "'scaleResiduals' must be a single logical value")
    for (nm in c("mafCutoff", "macCutoff", "xvarCutoff", "imissCutoff")) {
      v <- methods::slot(object, nm)
      if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 0)
        errors <- c(errors, sprintf(
          "'%s' must be a single finite non-negative numeric", nm))
    }
    if (length(object@phenotypes) == 0L)
      errors <- c(errors, "'phenotypes' must not be empty")
    contextNames <- names(object@phenotypes)
    if (is.null(contextNames) || any(!nzchar(contextNames)) ||
        any(is.na(contextNames)))
      errors <- c(errors, "'phenotypes' must be a named list with non-empty names")
    else if (anyDuplicated(contextNames))
      errors <- c(errors, "context names in 'phenotypes' must be unique")
    for (ctx in seq_along(object@phenotypes)) {
      se <- object@phenotypes[[ctx]]
      if (!methods::is(se, "SummarizedExperiment")) {
        errors <- c(errors, sprintf(
          "phenotypes[[%d]] must be a SummarizedExperiment (got %s)",
          ctx, class(se)[[1L]]))
      }
    }
    # Trait-position consistency across shared traits in different contexts.
    if (length(object@phenotypes) > 1L &&
        all(vapply(object@phenotypes, methods::is, logical(1),
                   "SummarizedExperiment"))) {
      traitToRange <- list()
      for (ctx in seq_along(object@phenotypes)) {
        se <- object@phenotypes[[ctx]]
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- rownames(se)
        if (length(rr) != length(ids)) next
        for (i in seq_along(ids)) {
          tid <- ids[[i]]
          prev <- traitToRange[[tid]]
          this <- rr[i]
          if (is.null(prev)) {
            traitToRange[[tid]] <- this
          } else {
            if (!isTRUE(all.equal(
                  as.character(GenomicRanges::seqnames(prev)),
                  as.character(GenomicRanges::seqnames(this))
                )) ||
                GenomicRanges::start(prev) != GenomicRanges::start(this) ||
                GenomicRanges::end(prev) != GenomicRanges::end(this)) {
              errors <- c(errors, sprintf(
                "trait '%s' has inconsistent rowRanges across contexts", tid))
            }
          }
        }
      }
    }
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Per-entry payload classes for FineMappingResult and TwasWeights
# =============================================================================

#' @title Fine-Mapping Entry (per-tuple payload)
#' @description S4 container for a single fine-mapping fit attached to a
#'   \code{FineMappingResult} row. One entry corresponds to one
#'   \code{(study, context, trait, method)} tuple.
#'
#'   For joint fits (e.g., multi-trait mvSuSiE or fSuSiE), multiple
#'   \code{FineMappingEntry} objects in the same \code{FineMappingResult}
#'   collection may carry references to the same underlying R fit object
#'   (R's copy-on-modify semantics keep this memory-efficient).
#' @slot variantIds Character vector of variant IDs in the fit.
#' @slot trimmedFit The method-specific fit object (SuSiE list, mvSuSiE
#'   object, fSuSiE object, etc.). May be shared by reference across
#'   joint-fit entries.
#' @slot topLoci A long-format \code{data.frame} with at minimum
#'   \code{variant_id} and \code{pip} columns; optional \code{cs},
#'   \code{coverage}, \code{betahat}, \code{sd}, \code{csLog10bf},
#'   \code{z}.
#' @slot sumstats A list of summary statistics used in the fit, or
#'   \code{NULL}.
#' @export

# =============================================================================
# QtlDataset constructor and accessors
# =============================================================================

#' @title Create a QtlDataset Object
#' @description Construct a \code{QtlDataset} S4 object containing one
#'   study's individual-level QTL data: a genotype handle and a named list
#'   of \code{SummarizedExperiment} objects (one per QTL context), plus
#'   genotype-derived covariates and a residual-scaling flag.
#' @param study Character (length 1). Study identifier.
#' @param genotypes A \code{GenotypeHandle}.
#' @param phenotypes Named list of \code{SummarizedExperiment} objects,
#'   keyed by context. Each SE must have \code{rowRanges} carrying trait
#'   positions and \code{colData} carrying per-context phenotype covariates.
#' @param genotypeCovariates Numeric matrix of genotype-derived covariates
#'   (e.g., ancestry PCs); rows are samples.
#' @param scaleResiduals Logical (length 1). Default \code{TRUE}.
#' @return A \code{QtlDataset} object.
#' @export
QtlDataset <- function(study, genotypes, phenotypes,
                       genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0),
                       scaleResiduals = TRUE,
                       mafCutoff = 0,
                       macCutoff = 0,
                       xvarCutoff = 0,
                       imissCutoff = 0,
                       keepSamples = character(0),
                       keepVariants = character(0)) {
  obj <- new("QtlDataset",
             study              = as.character(study),
             genotypes          = genotypes,
             phenotypes         = phenotypes,
             genotypeCovariates = as.matrix(genotypeCovariates),
             scaleResiduals     = isTRUE(scaleResiduals),
             mafCutoff          = as.numeric(mafCutoff),
             macCutoff          = as.numeric(macCutoff),
             xvarCutoff         = as.numeric(xvarCutoff),
             imissCutoff        = as.numeric(imissCutoff),
             keepSamples        = as.character(keepSamples),
             keepVariants       = as.character(keepVariants))
  validObject(obj)
  obj
}

#' @rdname getStudy
#' @export
setMethod("getStudy", "QtlDataset", function(x) x@study)

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlDataset", function(x) names(x@phenotypes))

#' @rdname getGenotypeCovariates
#' @export
setMethod("getGenotypeCovariates", "QtlDataset",
          function(x) x@genotypeCovariates)

#' @rdname getScaleResiduals
#' @export
setMethod("getScaleResiduals", "QtlDataset", function(x) x@scaleResiduals)

# --- Internal: resolve the variant-selection region for the genotype handle.
# Returns a single GRanges. When `traitId` is supplied, expand each trait's
# rowRange by `cisWindow` bp and take the union span (per the multi-trait rule:
# `[min(start) - cisWindow, max(end) + cisWindow]`). When `region` is supplied,
# extend by `cisWindow` if given. Exactly one of (traitId, region) may be
# supplied; if neither is, return NULL meaning "all variants in handle".
.qtlResolveVariantRegion <- function(x, traitId = NULL, region = NULL,
                                     cisWindow = NULL) {
  if (!is.null(traitId) && !is.null(region)) {
    stop("Specify either `traitId` or `region`, not both.")
  }
  if (is.null(traitId) && is.null(region)) {
    return(NULL)
  }
  if (!is.null(traitId)) {
    if (is.null(cisWindow) || length(cisWindow) != 1L || cisWindow < 0) {
      stop("`cisWindow` is required (and must be non-negative) when ",
           "`traitId` is specified.")
    }
    # Build a GRanges from each trait's rowRanges across all contexts;
    # take the union span +/- cisWindow.
    perTraitRanges <- list()
    for (ctxIdx in seq_along(x@phenotypes)) {
      se <- x@phenotypes[[ctxIdx]]
      rr <- SummarizedExperiment::rowRanges(se)
      hits <- match(traitId, rownames(se))
      hits <- hits[!is.na(hits)]
      if (length(hits) > 0) {
        perTraitRanges[[length(perTraitRanges) + 1L]] <- rr[hits]
      }
    }
    if (length(perTraitRanges) == 0L) {
      stop("None of the requested traitId values were found in any context.")
    }
    allRanges <- do.call(c, perTraitRanges)
    chrs <- unique(as.character(GenomicRanges::seqnames(allRanges)))
    if (length(chrs) != 1L) {
      stop("Multi-trait variant extraction requires all selected traits to ",
           "share a chromosome (got: ",
           paste(chrs, collapse = ", "), ").")
    }
    spanStart <- max(1L, min(GenomicRanges::start(allRanges)) - cisWindow)
    spanEnd   <- max(GenomicRanges::end(allRanges)) + cisWindow
    return(GenomicRanges::GRanges(
      seqnames = chrs,
      ranges   = IRanges::IRanges(start = spanStart, end = spanEnd)
    ))
  }
  # region path
  if (!methods::is(region, "GRanges")) {
    stop("`region` must be a GRanges object.")
  }
  if (length(region) != 1L) {
    stop("`region` must be a single range.")
  }
  if (!is.null(cisWindow)) {
    if (length(cisWindow) != 1L || cisWindow < 0) {
      stop("`cisWindow` must be a single non-negative value.")
    }
    region <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(region),
      ranges   = IRanges::IRanges(
        start = max(1L, GenomicRanges::start(region) - cisWindow),
        end   = GenomicRanges::end(region) + cisWindow
      )
    )
  }
  region
}

# Internal: map a GRanges region into 1-based snpIdx into handle@snpInfo.
.qtlVariantIndices <- function(x, region = NULL) {
  handle <- x@genotypes
  if (is.null(region)) {
    return(seq_len(nrow(handle@snpInfo)))
  }
  chr <- as.character(GenomicRanges::seqnames(region))[[1L]]
  chrCanon <- sub("^chr", "", chr, ignore.case = TRUE)
  snpInfo <- handle@snpInfo
  siChr <- sub("^chr", "", as.character(snpInfo$CHR), ignore.case = TRUE)
  bp <- as.integer(snpInfo$BP)
  start <- GenomicRanges::start(region)
  end   <- GenomicRanges::end(region)
  which(siChr == chrCanon & bp >= start & bp <= end)
}

# Internal: extract the panel dosage block (samples x variants) for the
# requested region, narrow to the requested sample set, and apply lazy QC
# (per-sample imiss filter, then per-variant max(mafCutoff,
# macCutoff / (2 * n)) and xvarCutoff filters). Used by getGenotypes,
# getResidualizedGenotypes (via getGenotypes), and getMaf so all three
# share a single variant/sample selection result.
#
# Returns a list:
#   geno       : numeric matrix (kept samples x kept variants)
#   variantIds : character vector of kept variant IDs (= colnames(geno))
#   sampleIds  : character vector of kept sample IDs (= rownames(geno))
#   maf        : numeric vector of per-variant MAF for kept variants
.qtlExtractBlock <- function(x, traitId = NULL, region = NULL,
                             cisWindow = NULL, samples = NULL) {
  gr <- .qtlResolveVariantRegion(x, traitId = traitId, region = region,
                                 cisWindow = cisWindow)
  snpIdx <- .qtlVariantIndices(x, gr)
  if (length(snpIdx) == 0L) {
    return(list(
      geno       = matrix(numeric(0), nrow = x@genotypes@nSamples, ncol = 0L,
                          dimnames = list(x@genotypes@sampleIds, character(0))),
      variantIds = character(0),
      sampleIds  = x@genotypes@sampleIds,
      maf        = numeric(0)
    ))
  }

  # Apply keepVariants restriction before materialization so we do not
  # extract dosage we will immediately drop.
  if (length(x@keepVariants) > 0L) {
    snpAll <- as.character(x@genotypes@snpInfo$SNP[snpIdx])
    keepMask <- snpAll %in% x@keepVariants
    snpIdx <- snpIdx[keepMask]
    if (length(snpIdx) == 0L) {
      return(list(
        geno       = matrix(numeric(0), nrow = 0L, ncol = 0L,
                            dimnames = list(character(0), character(0))),
        variantIds = character(0),
        sampleIds  = character(0),
        maf        = numeric(0)
      ))
    }
  }

  block <- extractBlockGenotypes(x@genotypes, snpIdx, meanImpute = FALSE)
  # `block` is variants x samples (Bioc convention); transpose to
  # samples x variants for analysis-style operations.
  dosage <- t(SummarizedExperiment::assay(block, "dosage"))

  # Resolve the requested sample set: keepSamples (panel-level) intersected
  # with the per-call samples arg, then intersected with the panel sample IDs.
  panelSamples <- rownames(dosage)
  keep <- panelSamples
  if (length(x@keepSamples) > 0L) {
    keep <- intersect(keep, x@keepSamples)
  }
  if (!is.null(samples)) {
    keep <- intersect(keep, as.character(samples))
  }
  if (length(keep) == 0L) {
    return(list(
      geno       = dosage[integer(0), , drop = FALSE],
      variantIds = colnames(dosage),
      sampleIds  = character(0),
      maf        = rep(NA_real_, ncol(dosage))
    ))
  }
  dosage <- dosage[keep, , drop = FALSE]

  # Per-sample missingness filter.
  if (x@imissCutoff > 0 && nrow(dosage) > 0L && ncol(dosage) > 0L) {
    imiss <- rowMeans(is.na(dosage))
    keepSampleMask <- imiss <= x@imissCutoff
    dosage <- dosage[keepSampleMask, , drop = FALSE]
  }

  # Per-variant MAF / MAC / X-variance filters against the post-narrowing
  # sample count. We mean-impute internally for the variance / dosage
  # returned but compute MAF from the un-imputed values so missingness is
  # handled correctly.
  nSamp <- nrow(dosage)
  if (ncol(dosage) > 0L) {
    nObs <- colSums(!is.na(dosage))
    sumD <- colSums(dosage, na.rm = TRUE)
    p <- ifelse(nObs > 0L, sumD / (2 * nObs), NA_real_)
    mafVec <- pmin(p, 1 - p)
    effectiveMaf <- max(x@mafCutoff, if (nSamp > 0L)
      x@macCutoff / (2 * nSamp) else 0)
    keepVarMask <- !is.na(mafVec) & mafVec >= effectiveMaf
    if (x@xvarCutoff > 0 && nSamp > 1L) {
      # Compute variance with mean imputation per column so variance is
      # defined when missingness is present.
      mu <- ifelse(nObs > 0L, sumD / nObs, 0)
      centered <- sweep(dosage, 2L, mu, FUN = "-")
      centered[is.na(centered)] <- 0
      varVec <- colSums(centered * centered) / (nSamp - 1L)
      keepVarMask <- keepVarMask & varVec >= x@xvarCutoff
    }
    dosage <- dosage[, keepVarMask, drop = FALSE]
    mafVec <- mafVec[keepVarMask]
  } else {
    mafVec <- numeric(0)
  }

  # Mean-impute remaining missing dosage cells so downstream linear
  # algebra is well-defined; MAF was computed before imputation.
  if (anyNA(dosage)) {
    for (j in seq_len(ncol(dosage))) {
      col <- dosage[, j]
      na <- is.na(col)
      if (any(na)) {
        col[na] <- mean(col[!na])
        dosage[, j] <- col
      }
    }
  }

  list(
    geno       = dosage,
    variantIds = colnames(dosage),
    sampleIds  = rownames(dosage),
    maf        = mafVec
  )
}

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "QtlDataset",
  function(x, traitId = NULL, region = NULL, cisWindow = NULL,
           samples = NULL, ...) {
    .qtlExtractBlock(x, traitId = traitId, region = region,
                     cisWindow = cisWindow, samples = samples)$geno
  })

#' @rdname getMaf
#' @export
setMethod("getMaf", "QtlDataset",
  function(x, region = NULL, cisWindow = NULL, samples = NULL, ...) {
    block <- .qtlExtractBlock(x, traitId = NULL, region = region,
                              cisWindow = cisWindow, samples = samples)
    out <- block$maf
    names(out) <- block$variantIds
    out
  })

#' @rdname getPhenotypes
#' @export
setMethod("getPhenotypes", "QtlDataset",
  function(x, contexts, traitId = NULL, region = NULL,
           naAction = c("keep", "drop", "impute"),
           outlierAction = c("keep", "drop"),
           outlierPvalThreshold = 1e-3,
           ...) {
    naAction <- match.arg(naAction)
    outlierAction <- match.arg(outlierAction)
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required for getPhenotypes(QtlDataset). ",
           "Pass a character vector of one or more context names; ",
           "use getContexts(x) to list the available contexts.")
    }
    available <- names(x@phenotypes)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "),
           ". Available: ", paste(available, collapse = ", "))
    }
    out <- x@phenotypes[contexts]
    if (!is.null(traitId)) {
      out <- lapply(seq_along(out), function(i) {
        se <- out[[i]]
        ctx <- names(out)[[i]]
        present <- intersect(traitId, rownames(se))
        missing <- setdiff(traitId, rownames(se))
        if (length(missing) > 0L) {
          warning(sprintf("context '%s' is missing trait(s): %s",
                          ctx, paste(missing, collapse = ", ")))
        }
        se[present, , drop = FALSE]
      })
      names(out) <- contexts
    }
    if (!is.null(region)) {
      out <- lapply(out, function(se) {
        rr <- SummarizedExperiment::rowRanges(se)
        keep <- IRanges::overlapsAny(rr, region)
        se[keep, , drop = FALSE]
      })
      names(out) <- contexts
    }
    if (naAction != "keep") {
      out <- lapply(out, function(se) .qtlApplyPhenoNaAction(se, naAction))
      names(out) <- contexts
    }
    if (outlierAction != "keep") {
      out <- lapply(out, function(se)
        .qtlApplyPhenoOutliers(se, outlierAction, outlierPvalThreshold))
      names(out) <- contexts
    }
    if (length(contexts) == 1L) out[[1L]] else out
  })

# Internal: apply naAction to a SummarizedExperiment slice. SE assay rows
# are traits and columns are samples.
#   "drop"   -> drop samples (cols) where any selected trait is NA
#   "impute" -> mean-impute each trait (row) independently over its
#               non-NA sample values
# Operates jointly over the rows currently in `se` -- the caller is
# expected to have already subset the SE to the user's requested
# (traitId, region) subset.
.qtlApplyPhenoNaAction <- function(se, naAction) {
  assayName <- SummarizedExperiment::assayNames(se)[[1L]]
  Y <- SummarizedExperiment::assay(se, assayName)
  if (length(Y) == 0L) return(se)
  if (naAction == "drop") {
    keepSamp <- colSums(is.na(Y)) == 0L
    se <- se[, keepSamp, drop = FALSE]
  } else if (naAction == "impute") {
    if (anyNA(Y)) {
      for (j in seq_len(nrow(Y))) {
        row <- Y[j, ]
        na <- is.na(row)
        if (any(na)) {
          obs <- row[!na]
          row[na] <- if (length(obs) > 0L) mean(obs) else 0
          Y[j, ] <- row
        }
      }
      SummarizedExperiment::assay(se, assayName) <- Y
    }
  }
  se
}

# Multivariate-outlier keep mask via Mahalanobis distance against a
# (preferably robust) centre / covariance estimate. Returns a logical
# vector of length nrow(Y); TRUE = keep, FALSE = drop.
#
# When the `robustbase` package is installed, the centre and covariance
# come from `robustbase::covMcd` (minimum-covariance-determinant) so the
# detector itself is resistant to the outliers it's trying to find.
# Without robustbase we fall back to `colMeans` / `cov` with a one-shot
# message; the test then still works but its estimates are pulled by
# the very outliers it should be flagging.
#
# Significance: per-sample chi-squared(p) p-value with Bonferroni
# correction over the sample count. A sample is flagged when its
# corrected p-value falls below `pvalThreshold`. With single-trait Y
# (ncol == 1) this reduces to the standard z-test on (y - center)/sd.
#
# Returns all-TRUE (no-op) when there are too few samples to support
# a covariance estimate (n < p + 2).
.qtlOutlierKeepMask <- function(Y, pvalThreshold) {
  Y <- as.matrix(Y)
  n <- nrow(Y); p <- ncol(Y)
  if (n == 0L || p == 0L) return(rep(TRUE, n))
  if (n < p + 2L) {
    warning(sprintf(
      "outlier detection skipped: %d samples < %d traits + 2 needed for ",
      n, p), "a covariance estimate.")
    return(rep(TRUE, n))
  }
  if (requireNamespace("robustbase", quietly = TRUE)) {
    mcd <- tryCatch(robustbase::covMcd(Y), error = function(e) NULL)
    if (!is.null(mcd)) {
      ctr <- mcd$center
      covMat <- mcd$cov
    } else {
      ctr <- colMeans(Y); covMat <- stats::cov(Y)
    }
  } else {
    message("outlier detection: install 'robustbase' for an MCD-based ",
            "estimator; falling back to non-robust colMeans/cov.")
    ctr <- colMeans(Y); covMat <- stats::cov(Y)
  }
  invCov <- tryCatch(solve(covMat),
                     error = function(e) MASS::ginv(covMat))
  Yc <- sweep(Y, 2L, ctr)
  d2 <- rowSums((Yc %*% invCov) * Yc)
  raw <- stats::pchisq(d2, df = p, lower.tail = FALSE)
  raw >= (pvalThreshold / n)
}

# Wrapper: apply the keep-mask to a SummarizedExperiment slice. SE
# columns are samples; transpose the assay (traits x samples) before
# calling .qtlOutlierKeepMask which expects samples x traits.
.qtlApplyPhenoOutliers <- function(se, action, pvalThreshold) {
  if (action == "keep") return(se)
  assayName <- SummarizedExperiment::assayNames(se)[[1L]]
  Y <- t(SummarizedExperiment::assay(se, assayName))
  keep <- .qtlOutlierKeepMask(Y, pvalThreshold)
  if (all(keep)) return(se)
  se[, keep, drop = FALSE]
}

#' @rdname getPhenotypeCovariates
#' @export
setMethod("getPhenotypeCovariates", "QtlDataset",
  function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required.")
    }
    available <- names(x@phenotypes)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "))
    }
    out <- lapply(contexts, function(ctx) {
      se <- x@phenotypes[[ctx]]
      cd <- SummarizedExperiment::colData(se)
      as.matrix(as.data.frame(cd))
    })
    names(out) <- contexts
    out
  })

# Internal: residualize a numeric matrix Y (n x k) against a covariate
# matrix C (n x p) via pivoted QR decomposition. Adds an intercept column
# to C. When C is rank-deficient (e.g., union of all contexts' phenotype
# covariates includes collinear / duplicate columns), the pivoted QR drops
# the redundant columns automatically. Optionally rescales each residual
# column to unit standard deviation; constant-valued columns are left
# unchanged.
.qtlResidualizeQR <- function(Y, C, scaleResiduals = TRUE) {
  X <- if (is.null(C) || ncol(C) == 0L) {
    matrix(1, nrow = nrow(Y), ncol = 1L,
           dimnames = list(rownames(Y), "intercept"))
  } else {
    cbind(intercept = 1, C)
  }
  # `qr.resid` does not support LAPACK pivoted QR, so use `lm.fit`. It
  # handles rank-deficient designs gracefully via base-R's pivoted QR
  # internally — same effect the LAPACK path was meant to deliver.
  res <- stats::lm.fit(x = X, y = Y)$residuals
  res <- as.matrix(res)
  rownames(res) <- rownames(Y)
  colnames(res) <- colnames(Y)
  if (isTRUE(scaleResiduals)) {
    sds <- apply(res, 2L, function(v) stats::sd(v, na.rm = TRUE))
    # `sds == 0` exact-zero test is unreliable for residuals coming out of
    # lm.fit on a constant Y: roundoff gives sd ~ 1e-16 instead of 0, and
    # dividing the (also-tiny) residuals by it amplifies floating-point
    # noise to unit-scale. Treat anything below sqrt(.Machine$double.eps)
    # as effectively zero (column is constant) and skip rescaling.
    nearZero <- !is.finite(sds) | sds < sqrt(.Machine$double.eps)
    sds[nearZero] <- 1
    res[, nearZero] <- 0
    res <- sweep(res, 2L, sds, FUN = "/")
  }
  res
}

# Internal: validate and resolve the `*ToResidualize` argument against a
# set of contexts and the covariates actually present in those contexts'
# colData. Accepts either NULL (use all), a character vector (apply to all
# listed contexts), or a named list keyed by context. Returns a named list
# keyed by context giving the actual character vector of covariate names
# to use for that context (or character(0) if none). Errors when:
#   - a named-list key is not in `contexts`
#   - `contexts` contains entries missing from a supplied named-list
#     (per the rule: named-list keys must equal `contexts`)
#   - an explicitly requested name matches no actual covariate
.qtlResolvePhenoSelection <- function(x, contexts, toResidualize) {
  resolveOne <- function(ctx, requested) {
    se <- x@phenotypes[[ctx]]
    avail <- colnames(SummarizedExperiment::colData(se))
    if (is.null(requested)) return(avail)
    keep <- intersect(requested, avail)
    if (length(keep) != length(requested)) {
      missingNames <- setdiff(requested, avail)
      stop(sprintf(
        "phenotypeCovariatesToResidualize: context '%s' has no covariate(s) named: %s",
        ctx, paste(missingNames, collapse = ", ")))
    }
    keep
  }
  if (is.null(toResidualize)) {
    out <- lapply(contexts, resolveOne, requested = NULL)
    names(out) <- contexts
    return(out)
  }
  if (is.list(toResidualize)) {
    if (is.null(names(toResidualize)) ||
        any(!nzchar(names(toResidualize)))) {
      stop("phenotypeCovariatesToResidualize: when supplied as a list, ",
           "it must be named with context names.")
    }
    badKeys <- setdiff(names(toResidualize), contexts)
    if (length(badKeys) > 0L) {
      stop("phenotypeCovariatesToResidualize: list key(s) not in `contexts`: ",
           paste(badKeys, collapse = ", "))
    }
    missingKeys <- setdiff(contexts, names(toResidualize))
    if (length(missingKeys) > 0L) {
      stop("phenotypeCovariatesToResidualize: list does not cover all ",
           "`contexts`. Per-context lists must have exactly the same ",
           "context set as `contexts`. Missing keys: ",
           paste(missingKeys, collapse = ", "))
    }
    out <- lapply(contexts, function(ctx) resolveOne(ctx, toResidualize[[ctx]]))
    names(out) <- contexts
    return(out)
  }
  if (is.character(toResidualize)) {
    out <- lapply(contexts, resolveOne, requested = toResidualize)
    names(out) <- contexts
    return(out)
  }
  stop("phenotypeCovariatesToResidualize must be NULL, a character vector, ",
       "or a named list keyed by context.")
}

# Internal: validate the genotype-covariate selection vector. Returns
# character(0) when nothing selected, the resolved set otherwise.
.qtlResolveGenoSelection <- function(x, toResidualize) {
  avail <- colnames(x@genotypeCovariates)
  if (is.null(avail)) avail <- character(0)
  if (is.null(toResidualize)) return(avail)
  keep <- intersect(toResidualize, avail)
  if (length(keep) != length(toResidualize)) {
    missingNames <- setdiff(toResidualize, avail)
    stop("genotypeCovariatesToResidualize: no covariate(s) named: ",
         paste(missingNames, collapse = ", "))
  }
  keep
}

# Internal: build the covariate matrix used for residualization, given a
# set of contexts, the resolved per-context phenotype selections, and the
# resolved genotype-covariate selection. Honors the inclusion flags. For
# `length(contexts) == 1` (per-context mode), the per-context phenotype
# covariates and the genotype covariates are taken with no cross-context
# alignment. For `length(contexts) >= 2` (joint mode), per-context
# phenotype covariates from all listed contexts are concatenated
# (prefixed with "{context}." to keep same-named columns distinct) and
# the sample set is intersected across all contributing matrices.
# Returns a single matrix (with rownames = sample IDs) or NULL.
.qtlBuildResidualizationDesign <- function(x, contexts,
                                           phenoSelection,
                                           genoSelection,
                                           includePheno, includeGeno) {
  perContext <- list()
  if (includePheno) {
    for (ctx in contexts) {
      keep <- phenoSelection[[ctx]]
      if (length(keep) == 0L) next
      se <- x@phenotypes[[ctx]]
      cd <- as.matrix(as.data.frame(SummarizedExperiment::colData(se)))
      cdMat <- cd[, keep, drop = FALSE]
      colnames(cdMat) <- paste0(ctx, ".", colnames(cdMat))
      if (is.null(rownames(cdMat))) {
        rownames(cdMat) <- as.character(
          rownames(SummarizedExperiment::colData(se)))
      }
      perContext[[ctx]] <- cdMat
    }
  }
  gCov <- if (includeGeno && length(genoSelection) > 0L) {
    x@genotypeCovariates[, genoSelection, drop = FALSE]
  } else {
    matrix(numeric(0), nrow = 0, ncol = 0)
  }
  haveAny <- length(perContext) > 0L ||
    (!is.null(gCov) && ncol(gCov) > 0L)
  if (!haveAny) return(NULL)

  sampleSets <- list()
  for (mat in perContext) {
    if (!is.null(rownames(mat))) {
      sampleSets[[length(sampleSets) + 1L]] <- rownames(mat)
    }
  }
  if (!is.null(gCov) && ncol(gCov) > 0L && !is.null(rownames(gCov))) {
    sampleSets[[length(sampleSets) + 1L]] <- rownames(gCov)
  }
  common <- if (length(sampleSets) == 0L) character(0)
            else Reduce(intersect, sampleSets)
  if (length(common) == 0L) return(NULL)

  blocks <- list()
  for (mat in perContext) {
    blocks[[length(blocks) + 1L]] <- mat[common, , drop = FALSE]
  }
  if (!is.null(gCov) && ncol(gCov) > 0L) {
    blocks[[length(blocks) + 1L]] <- gCov[common, , drop = FALSE]
  }
  do.call(cbind, blocks)
}

# Internal: resolve a (convenience, precise) flag pair to a single boolean.
# `missing*` arguments are passed as the result of `missing()` evaluated in
# the calling method to detect whether the user explicitly set the value.
# Rules:
#   - both missing: returns TRUE (the documented default)
#   - only convenience set: returns convenience
#   - only precise set: returns precise
#   - both set: must agree, else error
.qtlResolveResidualizationFlag <- function(conveniencePassed, convenienceMissing,
                                           precisePassed, preciseMissing,
                                           convenienceName, preciseName) {
  if (preciseMissing && convenienceMissing) return(TRUE)
  if (preciseMissing) return(isTRUE(conveniencePassed))
  if (convenienceMissing) return(isTRUE(precisePassed))
  if (isTRUE(conveniencePassed) != isTRUE(precisePassed)) {
    stop(sprintf(
      "Conflicting values: `%s` = %s and `%s` = %s. Set only one, or ",
      convenienceName, conveniencePassed,
      preciseName, precisePassed),
      "pass consistent values.")
  }
  isTRUE(precisePassed)
}

#' @rdname getResidualizedGenotypes
#' @export
setMethod("getResidualizedGenotypes", "QtlDataset",
  function(x, contexts, traitId = NULL, region = NULL, cisWindow = NULL,
           samples = NULL,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize = NULL,
           residualizePhenotypeCovariates = TRUE,
           residualizeGenotypeCovariates  = TRUE,
           residualizePhenotypeCovariatesFromGenotypes = NULL,
           residualizeGenotypeCovariatesFromGenotypes  = NULL,
           ...) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required for getResidualizedGenotypes(QtlDataset). ",
           "Use getContexts(x) to list the available contexts. ",
           "Pass a single context for per-context mode or multiple ",
           "contexts for joint mode (sample intersection).")
    }
    bad <- setdiff(contexts, names(x@phenotypes))
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "))
    }

    # Resolve inclusion flags (convenience vs precise).
    convPhenoMissing <- missing(residualizePhenotypeCovariates)
    convGenoMissing  <- missing(residualizeGenotypeCovariates)
    precPhenoMissing <- missing(residualizePhenotypeCovariatesFromGenotypes) ||
                       is.null(residualizePhenotypeCovariatesFromGenotypes)
    precGenoMissing  <- missing(residualizeGenotypeCovariatesFromGenotypes) ||
                       is.null(residualizeGenotypeCovariatesFromGenotypes)
    includePheno <- .qtlResolveResidualizationFlag(
      residualizePhenotypeCovariates, convPhenoMissing,
      residualizePhenotypeCovariatesFromGenotypes, precPhenoMissing,
      "residualizePhenotypeCovariates",
      "residualizePhenotypeCovariatesFromGenotypes")
    includeGeno <- .qtlResolveResidualizationFlag(
      residualizeGenotypeCovariates, convGenoMissing,
      residualizeGenotypeCovariatesFromGenotypes, precGenoMissing,
      "residualizeGenotypeCovariates",
      "residualizeGenotypeCovariatesFromGenotypes")

    # Resolve the covariate selections.
    phenoSel <- .qtlResolvePhenoSelection(x, contexts,
                                          phenotypeCovariatesToResidualize)
    genoSel  <- .qtlResolveGenoSelection(x, genotypeCovariatesToResidualize)

    G <- getGenotypes(x, traitId = traitId, region = region,
                      cisWindow = cisWindow, samples = samples)
    if (ncol(G) == 0L) return(G)

    C <- .qtlBuildResidualizationDesign(
      x, contexts = contexts,
      phenoSelection = phenoSel,
      genoSelection  = genoSel,
      includePheno   = includePheno,
      includeGeno    = includeGeno)
    if (!is.null(C)) {
      common <- intersect(rownames(G), rownames(C))
      if (length(common) == 0L) {
        stop("No samples in common between the genotype matrix and the ",
             "covariate matrix for contexts: ",
             paste(contexts, collapse = ", "))
      }
      G <- G[common, , drop = FALSE]
      C <- C[common, , drop = FALSE]
    }
    .qtlResidualizeQR(G, C, scaleResiduals = x@scaleResiduals)
  })

#' @rdname getResidualizedPhenotypes
#' @export
setMethod("getResidualizedPhenotypes", "QtlDataset",
  function(x, contexts, traitId = NULL, region = NULL,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize = NULL,
           residualizePhenotypeCovariates = TRUE,
           residualizeGenotypeCovariates  = TRUE,
           residualizePhenotypeCovariatesFromPhenotypes = NULL,
           residualizeGenotypeCovariatesFromPhenotypes  = NULL,
           naAction = c("keep", "drop", "impute"),
           outlierAction = c("keep", "drop"),
           outlierPvalThreshold = 1e-3,
           ...) {
    naAction <- match.arg(naAction)
    outlierAction <- match.arg(outlierAction)
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required for getResidualizedPhenotypes().")
    }
    bad <- setdiff(contexts, names(x@phenotypes))
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "))
    }

    convPhenoMissing <- missing(residualizePhenotypeCovariates)
    convGenoMissing  <- missing(residualizeGenotypeCovariates)
    precPhenoMissing <- missing(residualizePhenotypeCovariatesFromPhenotypes) ||
                       is.null(residualizePhenotypeCovariatesFromPhenotypes)
    precGenoMissing  <- missing(residualizeGenotypeCovariatesFromPhenotypes) ||
                       is.null(residualizeGenotypeCovariatesFromPhenotypes)
    includePheno <- .qtlResolveResidualizationFlag(
      residualizePhenotypeCovariates, convPhenoMissing,
      residualizePhenotypeCovariatesFromPhenotypes, precPhenoMissing,
      "residualizePhenotypeCovariates",
      "residualizePhenotypeCovariatesFromPhenotypes")
    includeGeno <- .qtlResolveResidualizationFlag(
      residualizeGenotypeCovariates, convGenoMissing,
      residualizeGenotypeCovariatesFromPhenotypes, precGenoMissing,
      "residualizeGenotypeCovariates",
      "residualizeGenotypeCovariatesFromPhenotypes")

    phenoSel <- .qtlResolvePhenoSelection(x, contexts,
                                          phenotypeCovariatesToResidualize)
    genoSel  <- .qtlResolveGenoSelection(x, genotypeCovariatesToResidualize)

    # NA-handle Y on the raw phenotype side, before residualization.
    Yraw <- getPhenotypes(x, contexts = contexts, traitId = traitId,
                          region = region, naAction = naAction)
    # getPhenotypes auto-unwraps single context; re-wrap so the lapply
    # below sees a consistent list shape.
    if (length(contexts) == 1L) Yraw <- setNames(list(Yraw), contexts)
    C <- .qtlBuildResidualizationDesign(
      x, contexts = contexts,
      phenoSelection = phenoSel,
      genoSelection  = genoSel,
      includePheno   = includePheno,
      includeGeno    = includeGeno)

    out <- lapply(contexts, function(ctx) {
      se <- Yraw[[ctx]]
      Y <- t(SummarizedExperiment::assay(se))  # samples x traits
      if (!is.null(C)) {
        common <- intersect(rownames(Y), rownames(C))
        if (length(common) == 0L) {
          stop(sprintf(
            "context '%s': no samples shared between phenotype data and ",
            "the resolved covariate matrix.", ctx))
        }
        Y <- Y[common, , drop = FALSE]
        Cctx <- C[common, , drop = FALSE]
      } else {
        Cctx <- NULL
      }
      Yres <- .qtlResidualizeQR(Y, Cctx, scaleResiduals = x@scaleResiduals)
      # Outlier detection on the residualized scale (samples whose
      # residualized phenotype is unusual *given* their covariates).
      if (outlierAction != "keep") {
        keep <- .qtlOutlierKeepMask(Yres, outlierPvalThreshold)
        if (!all(keep)) Yres <- Yres[keep, , drop = FALSE]
      }
      Yres
    })
    names(out) <- contexts
    if (length(contexts) == 1L) out[[1L]] else out
  })


#' @export
setMethod("show", "QtlDataset", function(object) {
  nCtx <- length(object@phenotypes)
  ctxNames <- names(object@phenotypes)
  totalTraits <- length(unique(unlist(
    lapply(object@phenotypes, rownames), use.names = FALSE)))
  cat(sprintf("QtlDataset for study '%s'\n", object@study))
  cat(sprintf("  %d context(s): %s\n", nCtx,
              paste(ctxNames, collapse = ", ")))
  cat(sprintf("  %d unique traits across contexts\n", totalTraits))
  cat(sprintf("  Genotypes: %s\n",
              paste0(object@genotypes@format, " @ ", object@genotypes@path)))
  cat(sprintf("  Genotype covariates: %d cols\n",
              ncol(object@genotypeCovariates)))
  cat(sprintf("  Scale residuals: %s\n", object@scaleResiduals))
})
