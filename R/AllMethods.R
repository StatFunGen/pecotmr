#' @title S4 Method Implementations
#' @description Constructors and accessor method implementations for all
#'   S4 classes: LdData, RegionalData, FineMappingResult, TwasWeights.
#' @name pecotmr-methods
#' @keywords internal
#' @include AllGenerics.R
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom GenomicRanges seqnames GRanges
NULL

# =============================================================================
# LdData constructor and accessors
# =============================================================================

#' @title Create an LdData Object
#' @description Construct an \code{LdData} from a correlation matrix and/or
#'   genotype handle, plus variant metadata as a GRanges.
#' @param correlation A correlation matrix, list of matrices, or NULL.
#' @param genotypeHandle A GenotypeHandle, list of GenotypeHandles, or NULL.
#' @param snpIdx Integer vector of SNP indices, or NULL.
#' @param variants A GRanges with variant metadata (must have variant_id in
#'   mcols, plus A1, A2).
#' @param blockMetadata LdBlocks or data.frame with block info.
#' @param nRef Integer, reference panel sample size.
#' @return An \code{LdData} object.
#' @export
LdData <- function(correlation = NULL, genotypeHandle = NULL,
                   snpIdx = NULL, variants, blockMetadata,
                   nRef = 0L) {
  obj <- new("LdData",
    correlation = correlation,
    genotypeHandle = genotypeHandle,
    snpIdx = snpIdx,
    variants = variants,
    blockMetadata = blockMetadata,
    nRef = as.integer(nRef)
  )
  validObject(obj)
  obj
}

#' @rdname getCorrelation
#' @export
setMethod("getCorrelation", "LdData", function(x) {
  if (!is.null(x@correlation)) return(x@correlation)
  if (is.null(x@genotypeHandle)) {
    stop("No correlation matrix or genotype handle available")
  }
  if (is.list(x@genotypeHandle)) {
    stop("Cannot compute single correlation matrix from mixture panels. ",
         "Use getGenotypes() and compute LD per-panel, or pass X directly ",
         "to susie_rss().")
  }
  geno <- extractBlockGenotypes(x@genotypeHandle, x@snpIdx)
  X <- t(assay(geno, "dosage"))
  computeLd(X, method = "sample")
})

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "LdData", function(x) {
  if (is.null(x@genotypeHandle)) return(NULL)
  # Plain matrix stored directly (e.g. from loadLdSketch after filtering)
  if (is.matrix(x@genotypeHandle)) return(x@genotypeHandle)
  if (is.list(x@genotypeHandle)) {
    lapply(x@genotypeHandle, function(h) {
      geno <- extractBlockGenotypes(h, x@snpIdx)
      t(assay(geno, "dosage"))
    })
  } else {
    geno <- extractBlockGenotypes(x@genotypeHandle, x@snpIdx)
    t(assay(geno, "dosage"))
  }
})

#' @rdname hasGenotypes
#' @export
setMethod("hasGenotypes", "LdData", function(x) {
  !is.null(x@genotypeHandle)
})

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "LdData", function(x) {
  mcols(x@variants)$variant_id
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "LdData", function(x) {
  x@variants
})

#' @rdname getBlockMetadata
#' @export
setMethod("getBlockMetadata", "LdData", function(x) {
  x@blockMetadata
})

#' @rdname getRefPanel
#' @export
setMethod("getRefPanel", "LdData", function(x) {
  mc <- as.data.frame(mcols(x@variants))
  mc$chrom <- as.character(seqnames(x@variants))
  mc$pos <- start(x@variants)
  mc
})

# =============================================================================
# Helper: build variant GRanges from refPanel data.frame
# =============================================================================

#' @title Build Variant GRanges
#' @description Convert a refPanel data.frame to a GRanges object suitable
#'   for the LdData variants slot.
#' @param refPanel data.frame with columns: chrom, pos, A1, A2, variant_id.
#'   May also have allele_freq, variance, n_nomiss.
#' @return A GRanges object.
#' @keywords internal
#' @noRd
.refPanelToGranges <- function(refPanel) {
  chr <- as.character(refPanel$chrom)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- paste0("chr", chr)
  pos <- as.integer(refPanel$pos)

  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = pos, width = 1L)
  )

  mcolsData <- DataFrame(
    variant_id = refPanel$variant_id,
    A1 = refPanel$A1,
    A2 = refPanel$A2
  )

  optional <- c("allele_freq", "variance", "n_nomiss")
  for (col in optional) {
    if (col %in% names(refPanel)) {
      mcolsData[[col]] <- refPanel[[col]]
    }
  }
  mcols(gr) <- mcolsData
  gr
}

# =============================================================================
# RegionalData constructor and accessors
# =============================================================================

#' @title Create a RegionalData Object
#' @description Construct a \code{RegionalData} from genotype matrix,
#'   phenotype and covariate lists.
#' @param genotypeMatrix Numeric matrix (samples x variants).
#' @param phenotypes Named list of phenotype matrices.
#' @param covariates Named list of covariate matrices.
#' @param scaleResiduals Logical.
#' @param maf Named list of MAF vectors.
#' @param region GRanges or NULL.
#' @param droppedSamples List.
#' @param coordinates data.frame or NULL.
#' @return A \code{RegionalData} object.
#' @export
RegionalData <- function(genotypeMatrix, phenotypes, covariates,
                         scaleResiduals = FALSE, maf = list(),
                         region = NULL, droppedSamples = list(),
                         coordinates = NULL) {
  obj <- new("RegionalData",
    genotypeMatrix = genotypeMatrix,
    phenotypes = phenotypes,
    covariates = covariates,
    scaleResiduals = scaleResiduals,
    maf = maf,
    region = region,
    droppedSamples = droppedSamples,
    coordinates = coordinates
  )
  validObject(obj)
  obj
}

#' @rdname getResidualX
#' @export
setMethod("getResidualX", "RegionalData", function(x, condition = 1L) {
  condition <- as.integer(condition)
  X <- x@genotypeMatrix
  covar <- x@covariates[[condition]]
  # Subset X to samples present in this condition's covariate
  common <- intersect(rownames(X), rownames(covar))
  XSub <- X[common, , drop = FALSE]
  covarSub <- covar[common, , drop = FALSE]
  res <- .lm.fit(x = cbind(1, covarSub), y = XSub)$residuals
  res <- as.matrix(res)
  if (x@scaleResiduals) res <- scale(res)
  colnames(res) <- colnames(XSub)
  rownames(res) <- common
  res
})

#' @rdname getResidualY
#' @export
setMethod("getResidualY", "RegionalData", function(x, condition = 1L) {
  condition <- as.integer(condition)
  Y <- x@phenotypes[[condition]]
  covar <- x@covariates[[condition]]
  common <- intersect(rownames(Y), rownames(covar))
  YSub <- as.matrix(Y[common, , drop = FALSE])
  covarSub <- covar[common, , drop = FALSE]
  res <- .lm.fit(x = cbind(1, covarSub), y = YSub)$residuals
  res <- as.matrix(res)
  colnames(res) <- colnames(YSub)
  if (x@scaleResiduals) res <- scale(res)
  rownames(res) <- common
  res
})

#' @rdname getResidualXScalar
#' @export
setMethod("getResidualXScalar", "RegionalData", function(x, condition = 1L) {
  if (!x@scaleResiduals) return(rep(1, ncol(x@genotypeMatrix)))
  res <- getResidualX(x, condition)
  apply(res, 2, sd)
})

#' @rdname getResidualYScalar
#' @export
setMethod("getResidualYScalar", "RegionalData", function(x, condition = 1L) {
  if (!x@scaleResiduals) return(1)
  res <- getResidualY(x, condition)
  apply(res, 2, sd)
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "RegionalData", function(x) {
  colnames(x@genotypeMatrix)
})

#' @rdname getPhenotypes
#' @export
setMethod("getPhenotypes", "RegionalData", function(x) x@phenotypes)

#' @rdname getCovariates
#' @export
setMethod("getCovariates", "RegionalData", function(x) x@covariates)

#' @rdname getGenotypeMatrix
#' @export
setMethod("getGenotypeMatrix", "RegionalData", function(x) x@genotypeMatrix)

#' @rdname getGenotypeMatrix
#' @export
setMethod("getGenotypeMatrix", "MultivariateRegionalData", function(x) x@genotypeMatrix)

# ----- MultivariateRegionalData constructor and accessors -----

#' @title Construct a MultivariateRegionalData object
#' @description Build a \code{MultivariateRegionalData} S4 object capturing
#'   regional association data prepared for multivariate modeling (single
#'   joint Y matrix across conditions).
#' @param genotypeMatrix Numeric matrix (samples x variants).
#' @param Y Numeric matrix (samples x conditions).
#' @param scaling Numeric vector of per-condition scaling factors.
#' @param droppedSamples Character vector or list of dropped sample IDs.
#' @param region A \code{GRanges} or NULL.
#' @param coordinates A data.frame of phenotype coordinates, or NULL.
#' @return A \code{MultivariateRegionalData} object.
#' @export
MultivariateRegionalData <- function(genotypeMatrix, Y, scaling,
                                     droppedSamples = NULL,
                                     region = NULL,
                                     coordinates = NULL) {
  obj <- new("MultivariateRegionalData",
             genotypeMatrix = genotypeMatrix,
             Y = Y,
             scaling = as.numeric(scaling),
             droppedSamples = droppedSamples,
             region = region,
             coordinates = coordinates)
  validObject(obj)
  obj
}

#' @rdname getY
#' @export
setMethod("getY", "MultivariateRegionalData", function(x) x@Y)

#' @rdname getScaling
#' @export
setMethod("getScaling", "MultivariateRegionalData", function(x) x@scaling)

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "MultivariateRegionalData", function(x) {
  colnames(x@genotypeMatrix)
})

#' @rdname getChrom
#' @export
setMethod("getChrom", "MultivariateRegionalData", function(x) {
  if (is.null(x@region)) return(NULL)
  as.character(GenomicRanges::seqnames(x@region))[1]
})

#' @rdname getGrange
#' @export
setMethod("getGrange", "MultivariateRegionalData", function(x) {
  if (is.null(x@region)) return(NULL)
  as.character(c(GenomicRanges::start(x@region),
                 GenomicRanges::end(x@region)))
})

#' @rdname getMaf
#' @export
setMethod("getMaf", "MultivariateRegionalData", function(x) {
  apply(x@genotypeMatrix, 2, computeMaf)
})

#' @rdname getXVariance
#' @export
setMethod("getXVariance", "MultivariateRegionalData", function(x, condition = 1L) {
  matrixStats::colVars(x@genotypeMatrix)
})

#' @rdname getXVariance
#' @export
setMethod("getXVariance", "RegionalData", function(x, condition = 1L) {
  res <- getResidualX(x, condition)
  matrixStats::colVars(res)
})

#' @rdname getChrom
#' @export
setMethod("getChrom", "RegionalData", function(x) {
  if (is.null(x@region)) return(NULL)
  as.character(GenomicRanges::seqnames(x@region))[1]
})

#' @rdname getGrange
#' @export
setMethod("getGrange", "RegionalData", function(x) {
  if (is.null(x@region)) return(NULL)
  as.character(c(GenomicRanges::start(x@region),
                 GenomicRanges::end(x@region)))
})

#' @title Combine Two RegionalData Objects
#' @description Concatenate two \code{RegionalData} objects by appending
#'   their per-condition slots (phenotypes, covariates, maf, droppedSamples).
#'   Used by multi-panel pipelines that load per-LD-panel data and aggregate
#'   them. The \code{genotypeMatrix} of \code{x} is retained as the
#'   canonical genotype reference; the \code{region} is taken from \code{y}
#'   (mirrors prior list-merge behavior).
#' @param x First \code{RegionalData} object.
#' @param y Second \code{RegionalData} object.
#' @return A merged \code{RegionalData}.
#' @export
setMethod("c", "RegionalData", function(x, ...) {
  others <- list(...)
  if (length(others) == 0L) return(x)
  result <- x
  for (y in others) {
    if (!is(y, "RegionalData")) stop("All arguments to c() must be RegionalData")
    result <- RegionalData(
      genotypeMatrix = result@genotypeMatrix,
      phenotypes = c(result@phenotypes, y@phenotypes),
      covariates = c(result@covariates, y@covariates),
      scaleResiduals = result@scaleResiduals,
      maf = c(result@maf, y@maf),
      region = y@region,
      droppedSamples = list(
        X = c(result@droppedSamples$X, y@droppedSamples$X),
        Y = c(result@droppedSamples$Y, y@droppedSamples$Y),
        covar = c(result@droppedSamples$covar, y@droppedSamples$covar)
      ),
      coordinates = result@coordinates
    )
  }
  result
})

# =============================================================================
# FineMappingResult constructor and accessors
# =============================================================================

#' @title Create a FineMappingResult Object
#' @description Construct a \code{FineMappingResult} from fine-mapping output.
#' @param variantNames Character vector.
#' @param trimmedFit List.
#' @param topLoci data.frame in long format.
#' @param method Character.
#' @param sumstats List or NULL.
#' @return A \code{FineMappingResult} object.
#' @export
FineMappingResult <- function(variantNames, trimmedFit, topLoci,
                              method, sumstats = NULL) {
  obj <- new("FineMappingResult",
    variantNames = variantNames,
    trimmedFit = trimmedFit,
    topLoci = topLoci,
    method = method,
    sumstats = sumstats
  )
  validObject(obj)
  obj
}

#' @rdname getPip
#' @export
setMethod("getPip", "FineMappingResult", function(x) {
  tl <- x@topLoci
  if (nrow(tl) == 0 || !"pip" %in% names(tl)) return(numeric(0))
  setNames(tl$pip, tl$variant_id)
})

#' @rdname getCs
#' @export
setMethod("getCs", "FineMappingResult", function(x, coverage = 0.95) {
  tl <- x@topLoci
  if (nrow(tl) == 0) return(data.frame())
  csCol <- grep(paste0("^cs.*", coverage * 100), names(tl), value = TRUE)
  if (length(csCol) == 0 && "cs" %in% names(tl)) csCol <- "cs"
  if (length(csCol) == 0) return(data.frame())
  tl[tl[[csCol[1]]] > 0, , drop = FALSE]
})

#' @rdname getLbf
#' @export
setMethod("getLbf", "FineMappingResult", function(x) {
  fit <- x@trimmedFit
  if (is.null(fit)) return(data.frame())

  # Extract lbf_variable matrix (effects x variants)
  lbf <- NULL
  if (!is.null(fit$lbf_variable)) {
    lbf <- fit$lbf_variable
  }
  if (is.null(lbf)) return(data.frame())

  # Build data.frame: variant_id + one column per effect
  variantIds <- x@variantNames
  if (length(variantIds) == 0 && !is.null(names(fit$pip)))
    variantIds <- names(fit$pip)

  df <- data.frame(variant_id = variantIds, stringsAsFactors = FALSE)
  lbfT <- t(lbf)  # transpose to variants x effects
  colnames(lbfT) <- paste0("L", seq_len(ncol(lbfT)))
  cbind(df, as.data.frame(lbfT))
})

#' @rdname getEffects
#' @export
setMethod("getEffects", "FineMappingResult", function(x) {
  fit <- x@trimmedFit
  if (is.null(fit)) return(data.frame())

  variantIds <- x@variantNames
  if (length(variantIds) == 0 && !is.null(names(fit$pip)))
    variantIds <- names(fit$pip)

  # Number of effects
  nEffects <- if (!is.null(fit$V)) length(fit$V) else if (!is.null(fit$alpha)) nrow(fit$alpha) else 0L
  if (nEffects == 0) return(data.frame())

  effectIds <- paste0("L", seq_len(nEffects))

  # Prior variance
  V <- if (!is.null(fit$V)) fit$V else rep(NA_real_, nEffects)

  # Per-effect log BF
  csLog10bf <- if (!is.null(fit$lbf)) fit$lbf else rep(NA_real_, nEffects)

  # Credible set info
  csSets <- fit$sets$cs
  csCoverage <- fit$sets$coverage
  csPurity <- fit$sets$purity

  csVariants <- character(nEffects)
  coverage <- numeric(nEffects)
  csMinR2 <- numeric(nEffects)
  csAvgR2 <- numeric(nEffects)

  for (i in seq_len(nEffects)) {
    eid <- effectIds[i]
    if (!is.null(csSets) && eid %in% names(csSets)) {
      idx <- csSets[[eid]]
      csVariants[i] <- paste(variantIds[idx], collapse = ";")
      csIdx <- which(names(csSets) == eid)
      if (!is.null(csCoverage))
        coverage[i] <- csCoverage[csIdx]
      if (!is.null(csPurity) && is.matrix(csPurity)) {
        csMinR2[i] <- csPurity[eid, 1]
        csAvgR2[i] <- csPurity[eid, 2]
      }
    } else {
      csVariants[i] <- "None"
    }
  }

  data.frame(
    effect_id = effectIds,
    V = V,
    cs_log10bf = csLog10bf,
    cs_min_r2 = csMinR2,
    cs_avg_r2 = csAvgR2,
    coverage = coverage,
    cs = csVariants,
    stringsAsFactors = FALSE)
})

# =============================================================================
# TwasWeights constructor and accessors
# =============================================================================

#' @title Create a TwasWeights Object
#' @description Construct a \code{TwasWeights} from weight matrices.
#' @param weights Named list of matrices.
#' @param variantIds Character vector.
#' @param fits Named list or NULL.
#' @param cvPerformance Named list or NULL.
#' @param standardized Logical. If TRUE, weights are on the standardized
#'   (correlation) scale and do not need variance scaling in harmonizeTwas.
#'   Defaults to FALSE (individual-level / raw genotype scale).
#' @return A \code{TwasWeights} object.
#' @export
TwasWeights <- function(weights, variantIds, fits = NULL,
                        cvPerformance = NULL, standardized = FALSE,
                        molecularId = character(0), dataType = NULL) {
  obj <- new("TwasWeights",
    weights = weights,
    variantIds = variantIds,
    methods = names(weights),
    fits = fits,
    cvPerformance = cvPerformance,
    standardized = standardized,
    molecularId = molecularId,
    dataType = dataType
  )
  validObject(obj)
  obj
}

#' @rdname getWeights
#' @export
setMethod("getWeights", "TwasWeights", function(x, method = NULL) {
  if (is.null(method)) return(x@weights)
  if (!method %in% x@methods)
    stop("Method '", method, "' not found. Available: ",
         paste(x@methods, collapse = ", "))
  x@weights[[method]]
})

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TwasWeights", function(x) x@standardized)

#' @rdname getCvPerformance
#' @export
setMethod("getCvPerformance", "TwasWeights", function(x, method = NULL) {
  if (is.null(method)) return(x@cvPerformance)
  x@cvPerformance[[method]]
})

#' @rdname getFits
#' @export
setMethod("getFits", "TwasWeights", function(x, method = NULL) {
  if (is.null(method)) return(x@fits)
  x@fits[[method]]
})

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "TwasWeights", function(x) x@methods)

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TwasWeights", function(x) x@variantIds)

#' @rdname getMolecularId
#' @export
setMethod("getMolecularId", "TwasWeights", function(x) x@molecularId)

#' @rdname getDataType
#' @export
setMethod("getDataType", "TwasWeights", function(x) x@dataType)

# =============================================================================
# FineMappingResult additional accessors
# =============================================================================

#' @rdname getTrimmedFit
#' @export
setMethod("getTrimmedFit", "FineMappingResult", function(x) x@trimmedFit)

#' @rdname getVariantNames
#' @export
setMethod("getVariantNames", "FineMappingResult", function(x) x@variantNames)

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "FineMappingResult", function(x) x@topLoci)

# =============================================================================
# AlleleQcResult constructor and accessors
# =============================================================================

#' @title Construct an AlleleQcResult object
#' @description Build an \code{AlleleQcResult} S4 object wrapping the post-QC
#'   harmonized variants and the full per-variant QC diagnostics.
#' @param harmonizedData Data frame of variants retained after allele QC.
#' @param qcSummary Data frame of per-variant diagnostic columns.
#' @return An \code{AlleleQcResult} object.
#' @export
AlleleQcResult <- function(harmonizedData, qcSummary) {
  obj <- new("AlleleQcResult",
      harmonizedData = as.data.frame(harmonizedData),
      qcSummary = as.data.frame(qcSummary))
  validObject(obj)
  obj
}

#' @rdname getHarmonizedData
#' @export
setMethod("getHarmonizedData", "AlleleQcResult", function(x) x@harmonizedData)

#' @rdname getQcSummary
#' @export
setMethod("getQcSummary", "AlleleQcResult", function(x) x@qcSummary)

# =============================================================================
# QcResult constructor and accessors
# =============================================================================

#' @title Construct a QcResult object
#' @description Build a \code{QcResult} S4 object capturing the output of
#'   summary-statistic QC. Validates that \code{ldData} is an \code{LdData}
#'   or NULL.
#' @param ldData An \code{LdData} or NULL.
#' @param rssInput List with \code{sumstats}, \code{n}, \code{varY}.
#' @param preprocess List with \code{sumstats} and \code{ldData}.
#' @param outlierNumber Integer count of LD-mismatch outliers removed.
#' @param skipped Single logical indicating a short-circuit.
#' @param skipReason Character explanation; defaults to empty.
#' @return A \code{QcResult} object.
#' @export
QcResult <- function(ldData = NULL,
                     rssInput = list(),
                     preprocess = list(),
                     outlierNumber = 0L,
                     skipped = FALSE,
                     skipReason = "") {
  reason <- if (length(skipReason) == 0L) "" else as.character(skipReason)[[1]]
  obj <- new("QcResult",
      ldData = ldData,
      rssInput = rssInput,
      preprocess = preprocess,
      outlierNumber = as.integer(outlierNumber),
      skipped = isTRUE(skipped),
      skipReason = reason)
  validObject(obj)
  obj
}

#' @rdname getLdData
#' @export
setMethod("getLdData", "QcResult", function(x) x@ldData)

#' @rdname getRssInput
#' @export
setMethod("getRssInput", "QcResult", function(x) x@rssInput)

#' @rdname getPreprocess
#' @export
setMethod("getPreprocess", "QcResult", function(x) x@preprocess)

#' @rdname getOutlierNumber
#' @export
setMethod("getOutlierNumber", "QcResult", function(x) x@outlierNumber)

#' @rdname isSkipped
#' @export
setMethod("isSkipped", "QcResult", function(x) x@skipped)

#' @rdname getSkipReason
#' @export
setMethod("getSkipReason", "QcResult", function(x) x@skipReason)

# =============================================================================
# topLoci GRanges conversion
# =============================================================================

#' @title Convert topLoci to GRanges
#' @description Convert a long-format topLoci data.frame (with variant_id
#'   in chr:pos:A2:A1 format) to a GRanges object with all metadata columns.
#' @param topLoci A data.frame with a variant_id column encoding
#'   chr:pos:A2:A1.
#' @return A GRanges object with all original columns as metadata.
#' @export
topLociToGranges <- function(topLoci) {
  if (is.null(topLoci) || nrow(topLoci) == 0) {
    return(GRanges())
  }
  parsed <- parseVariantId(topLoci$variant_id)
  chr <- paste0("chr", parsed$chrom)
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = parsed$pos, width = 1L)
  )
  mcols(gr) <- DataFrame(topLoci)
  gr
}
