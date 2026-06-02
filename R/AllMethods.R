#' @title S4 Method Implementations
#' @description Constructors and accessor method implementations for all
#'   S4 classes: LDData, RegionalData, FineMappingResult, TWASWeights.
#' @include AllGenerics.R
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom GenomicRanges seqnames GRanges
NULL

# =============================================================================
# LDData constructor and accessors
# =============================================================================

#' @title Create an LDData Object
#' @description Construct an \code{LDData} from a correlation matrix and/or
#'   genotype handle, plus variant metadata as a GRanges.
#' @param correlation A correlation matrix, list of matrices, or NULL.
#' @param genotype_handle A GenotypeHandle, list of GenotypeHandles, or NULL.
#' @param snp_idx Integer vector of SNP indices, or NULL.
#' @param variants A GRanges with variant metadata (must have variant_id in
#'   mcols, plus A1, A2).
#' @param block_metadata LDBlocks or data.frame with block info.
#' @param n_ref Integer, reference panel sample size.
#' @return An \code{LDData} object.
#' @export
LDData <- function(correlation = NULL, genotype_handle = NULL,
                   snp_idx = NULL, variants, block_metadata,
                   n_ref = 0L) {
  obj <- new("LDData",
    correlation = correlation,
    genotype_handle = genotype_handle,
    snp_idx = snp_idx,
    variants = variants,
    block_metadata = block_metadata,
    n_ref = as.integer(n_ref)
  )
  validObject(obj)
  obj
}

#' @rdname getCorrelation
#' @export
setMethod("getCorrelation", "LDData", function(x) {
  if (!is.null(x@correlation)) return(x@correlation)
  if (is.null(x@genotype_handle)) {
    stop("No correlation matrix or genotype handle available")
  }
  if (is.list(x@genotype_handle)) {
    stop("Cannot compute single correlation matrix from mixture panels. ",
         "Use getGenotypes() and compute LD per-panel, or pass X directly ",
         "to susie_rss().")
  }
  geno <- extractBlockGenotypes(x@genotype_handle, x@snp_idx)
  X <- t(assay(geno, "dosage"))
  compute_LD(X, method = "sample")
})

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "LDData", function(x) {
  if (is.null(x@genotype_handle)) return(NULL)
  # Plain matrix stored directly (e.g. from load_ld_sketch after filtering)
  if (is.matrix(x@genotype_handle)) return(x@genotype_handle)
  if (is.list(x@genotype_handle)) {
    lapply(x@genotype_handle, function(h) {
      geno <- extractBlockGenotypes(h, x@snp_idx)
      t(assay(geno, "dosage"))
    })
  } else {
    geno <- extractBlockGenotypes(x@genotype_handle, x@snp_idx)
    t(assay(geno, "dosage"))
  }
})

#' @rdname hasGenotypes
#' @export
setMethod("hasGenotypes", "LDData", function(x) {
  !is.null(x@genotype_handle)
})

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "LDData", function(x) {
  mcols(x@variants)$variant_id
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "LDData", function(x) {
  x@variants
})

#' @rdname getBlockMetadata
#' @export
setMethod("getBlockMetadata", "LDData", function(x) {
  x@block_metadata
})

#' @rdname getRefPanel
#' @export
setMethod("getRefPanel", "LDData", function(x) {
  mc <- as.data.frame(mcols(x@variants))
  mc$chrom <- as.character(seqnames(x@variants))
  mc$pos <- start(x@variants)
  mc
})

# =============================================================================
# Helper: build variant GRanges from ref_panel data.frame
# =============================================================================

#' @title Build Variant GRanges
#' @description Convert a ref_panel data.frame to a GRanges object suitable
#'   for the LDData variants slot.
#' @param ref_panel data.frame with columns: chrom, pos, A1, A2, variant_id.
#'   May also have allele_freq, variance, n_nomiss.
#' @return A GRanges object.
#' @keywords internal
#' @noRd
.ref_panel_to_granges <- function(ref_panel) {
  chr <- as.character(ref_panel$chrom)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- paste0("chr", chr)
  pos <- as.integer(ref_panel$pos)

  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = pos, width = 1L)
  )

  mcols_data <- DataFrame(
    variant_id = ref_panel$variant_id,
    A1 = ref_panel$A1,
    A2 = ref_panel$A2
  )

  optional <- c("allele_freq", "variance", "n_nomiss")
  for (col in optional) {
    if (col %in% names(ref_panel)) {
      mcols_data[[col]] <- ref_panel[[col]]
    }
  }
  mcols(gr) <- mcols_data
  gr
}

# =============================================================================
# RegionalData constructor and accessors
# =============================================================================

#' @title Create a RegionalData Object
#' @description Construct a \code{RegionalData} from genotype matrix,
#'   phenotype and covariate lists.
#' @param genotype_matrix Numeric matrix (samples x variants).
#' @param phenotypes Named list of phenotype matrices.
#' @param covariates Named list of covariate matrices.
#' @param scale_residuals Logical.
#' @param maf Named list of MAF vectors.
#' @param region GRanges or NULL.
#' @param dropped_samples List.
#' @param Y_coordinates data.frame or NULL.
#' @return A \code{RegionalData} object.
#' @export
RegionalData <- function(genotype_matrix, phenotypes, covariates,
                         scale_residuals = FALSE, maf = list(),
                         region = NULL, dropped_samples = list(),
                         Y_coordinates = NULL) {
  obj <- new("RegionalData",
    genotype_matrix = genotype_matrix,
    phenotypes = phenotypes,
    covariates = covariates,
    scale_residuals = scale_residuals,
    maf = maf,
    region = region,
    dropped_samples = dropped_samples,
    Y_coordinates = Y_coordinates
  )
  validObject(obj)
  obj
}

#' @rdname getResidualX
#' @export
setMethod("getResidualX", "RegionalData", function(x, condition = 1L) {
  condition <- as.integer(condition)
  X <- x@genotype_matrix
  covar <- x@covariates[[condition]]
  # Subset X to samples present in this condition's covariate
  common <- intersect(rownames(X), rownames(covar))
  X_sub <- X[common, , drop = FALSE]
  covar_sub <- covar[common, , drop = FALSE]
  res <- .lm.fit(x = cbind(1, covar_sub), y = X_sub)$residuals
  res <- as.matrix(res)
  if (x@scale_residuals) res <- scale(res)
  colnames(res) <- colnames(X_sub)
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
  Y_sub <- as.matrix(Y[common, , drop = FALSE])
  covar_sub <- covar[common, , drop = FALSE]
  res <- .lm.fit(x = cbind(1, covar_sub), y = Y_sub)$residuals
  res <- as.matrix(res)
  colnames(res) <- colnames(Y_sub)
  if (x@scale_residuals) res <- scale(res)
  rownames(res) <- common
  res
})

#' @rdname getResidualXScalar
#' @export
setMethod("getResidualXScalar", "RegionalData", function(x, condition = 1L) {
  if (!x@scale_residuals) return(rep(1, ncol(x@genotype_matrix)))
  res <- getResidualX(x, condition)
  apply(res, 2, sd)
})

#' @rdname getResidualYScalar
#' @export
setMethod("getResidualYScalar", "RegionalData", function(x, condition = 1L) {
  if (!x@scale_residuals) return(1)
  res <- getResidualY(x, condition)
  apply(res, 2, sd)
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "RegionalData", function(x) {
  colnames(x@genotype_matrix)
})

#' @rdname getPhenotypes
#' @export
setMethod("getPhenotypes", "RegionalData", function(x) x@phenotypes)

#' @rdname getCovariates
#' @export
setMethod("getCovariates", "RegionalData", function(x) x@covariates)

#' @rdname getGenotypeMatrix
#' @export
setMethod("getGenotypeMatrix", "RegionalData", function(x) x@genotype_matrix)

#' @rdname getGenotypeMatrix
#' @export
setMethod("getGenotypeMatrix", "MultivariateRegionalData", function(x) x@genotype_matrix)

# ----- MultivariateRegionalData constructor and accessors -----

#' @title Construct a MultivariateRegionalData object
#' @description Build a \code{MultivariateRegionalData} S4 object capturing
#'   regional association data prepared for multivariate modeling (single
#'   joint Y matrix across conditions).
#' @param genotype_matrix Numeric matrix (samples x variants).
#' @param Y_matrix Numeric matrix (samples x conditions).
#' @param Y_scalar Numeric vector of per-condition scaling factors.
#' @param dropped_samples Character vector or list of dropped sample IDs.
#' @param region A \code{GRanges} or NULL.
#' @param Y_coordinates A data.frame of phenotype coordinates, or NULL.
#' @return A \code{MultivariateRegionalData} object.
#' @export
MultivariateRegionalData <- function(genotype_matrix, Y_matrix, Y_scalar,
                                     dropped_samples = NULL,
                                     region = NULL,
                                     Y_coordinates = NULL) {
  obj <- new("MultivariateRegionalData",
             genotype_matrix = genotype_matrix,
             Y_matrix = Y_matrix,
             Y_scalar = as.numeric(Y_scalar),
             dropped_samples = dropped_samples,
             region = region,
             Y_coordinates = Y_coordinates)
  validObject(obj)
  obj
}

#' @rdname getYMatrix
#' @export
setMethod("getYMatrix", "MultivariateRegionalData", function(x) x@Y_matrix)

#' @rdname getYScalar
#' @export
setMethod("getYScalar", "MultivariateRegionalData", function(x) x@Y_scalar)

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "MultivariateRegionalData", function(x) {
  colnames(x@genotype_matrix)
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
  apply(x@genotype_matrix, 2, compute_maf)
})

#' @rdname getXVariance
#' @export
setMethod("getXVariance", "MultivariateRegionalData", function(x, condition = 1L) {
  matrixStats::colVars(x@genotype_matrix)
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
#'   their per-condition slots (phenotypes, covariates, maf, dropped_samples).
#'   Used by multi-panel pipelines that load per-LD-panel data and aggregate
#'   them. The \code{genotype_matrix} of \code{x} is retained as the
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
      genotype_matrix = result@genotype_matrix,
      phenotypes = c(result@phenotypes, y@phenotypes),
      covariates = c(result@covariates, y@covariates),
      scale_residuals = result@scale_residuals,
      maf = c(result@maf, y@maf),
      region = y@region,
      dropped_samples = list(
        X = c(result@dropped_samples$X, y@dropped_samples$X),
        Y = c(result@dropped_samples$Y, y@dropped_samples$Y),
        covar = c(result@dropped_samples$covar, y@dropped_samples$covar)
      ),
      Y_coordinates = result@Y_coordinates
    )
  }
  result
})

# =============================================================================
# FineMappingResult constructor and accessors
# =============================================================================

#' @title Create a FineMappingResult Object
#' @description Construct a \code{FineMappingResult} from fine-mapping output.
#' @param variant_names Character vector.
#' @param trimmed_fit List.
#' @param top_loci data.frame in long format.
#' @param method Character.
#' @param sumstats List or NULL.
#' @return A \code{FineMappingResult} object.
#' @export
FineMappingResult <- function(variant_names, trimmed_fit, top_loci,
                              method, sumstats = NULL) {
  obj <- new("FineMappingResult",
    variant_names = variant_names,
    trimmed_fit = trimmed_fit,
    top_loci = top_loci,
    method = method,
    sumstats = sumstats
  )
  validObject(obj)
  obj
}

#' @rdname getPIP
#' @export
setMethod("getPIP", "FineMappingResult", function(x) {
  tl <- x@top_loci
  if (nrow(tl) == 0 || !"pip" %in% names(tl)) return(numeric(0))
  setNames(tl$pip, tl$variant_id)
})

#' @rdname getCS
#' @export
setMethod("getCS", "FineMappingResult", function(x, coverage = 0.95) {
  tl <- x@top_loci
  if (nrow(tl) == 0) return(data.frame())
  cs_col <- grep(paste0("^cs.*", coverage * 100), names(tl), value = TRUE)
  if (length(cs_col) == 0 && "cs" %in% names(tl)) cs_col <- "cs"
  if (length(cs_col) == 0) return(data.frame())
  tl[tl[[cs_col[1]]] > 0, , drop = FALSE]
})

#' @rdname getLBF
#' @export
setMethod("getLBF", "FineMappingResult", function(x) {
  fit <- x@trimmed_fit
  if (is.null(fit)) return(data.frame())

  # Extract lbf_variable matrix (effects x variants)
  lbf <- NULL
  if (!is.null(fit$lbf_variable)) {
    lbf <- fit$lbf_variable
  }
  if (is.null(lbf)) return(data.frame())

  # Build data.frame: variant_id + one column per effect
  variant_ids <- x@variant_names
  if (length(variant_ids) == 0 && !is.null(names(fit$pip)))
    variant_ids <- names(fit$pip)

  df <- data.frame(variant_id = variant_ids, stringsAsFactors = FALSE)
  lbf_t <- t(lbf)  # transpose to variants x effects
  colnames(lbf_t) <- paste0("L", seq_len(ncol(lbf_t)))
  cbind(df, as.data.frame(lbf_t))
})

#' @rdname getEffects
#' @export
setMethod("getEffects", "FineMappingResult", function(x) {
  fit <- x@trimmed_fit
  if (is.null(fit)) return(data.frame())

  variant_ids <- x@variant_names
  if (length(variant_ids) == 0 && !is.null(names(fit$pip)))
    variant_ids <- names(fit$pip)

  # Number of effects
  n_effects <- if (!is.null(fit$V)) length(fit$V) else if (!is.null(fit$alpha)) nrow(fit$alpha) else 0L
  if (n_effects == 0) return(data.frame())

  effect_ids <- paste0("L", seq_len(n_effects))

  # Prior variance
  V <- if (!is.null(fit$V)) fit$V else rep(NA_real_, n_effects)

  # Per-effect log BF
  cs_log10bf <- if (!is.null(fit$lbf)) fit$lbf else rep(NA_real_, n_effects)

  # Credible set info
  cs_sets <- fit$sets$cs
  cs_coverage <- fit$sets$coverage
  cs_purity <- fit$sets$purity

  cs_variants <- character(n_effects)
  coverage <- numeric(n_effects)
  cs_min_r2 <- numeric(n_effects)
  cs_avg_r2 <- numeric(n_effects)

  for (i in seq_len(n_effects)) {
    eid <- effect_ids[i]
    if (!is.null(cs_sets) && eid %in% names(cs_sets)) {
      idx <- cs_sets[[eid]]
      cs_variants[i] <- paste(variant_ids[idx], collapse = ";")
      cs_idx <- which(names(cs_sets) == eid)
      if (!is.null(cs_coverage))
        coverage[i] <- cs_coverage[cs_idx]
      if (!is.null(cs_purity) && is.matrix(cs_purity)) {
        cs_min_r2[i] <- cs_purity[eid, 1]
        cs_avg_r2[i] <- cs_purity[eid, 2]
      }
    } else {
      cs_variants[i] <- "None"
    }
  }

  data.frame(
    effect_id = effect_ids,
    V = V,
    cs_log10bf = cs_log10bf,
    cs_min_r2 = cs_min_r2,
    cs_avg_r2 = cs_avg_r2,
    coverage = coverage,
    cs = cs_variants,
    stringsAsFactors = FALSE)
})

# =============================================================================
# TWASWeights constructor and accessors
# =============================================================================

#' @title Create a TWASWeights Object
#' @description Construct a \code{TWASWeights} from weight matrices.
#' @param weights Named list of matrices.
#' @param variant_ids Character vector.
#' @param fits Named list or NULL.
#' @param cv_performance Named list or NULL.
#' @param standardized Logical. If TRUE, weights are on the standardized
#'   (correlation) scale and do not need variance scaling in harmonize_twas.
#'   Defaults to FALSE (individual-level / raw genotype scale).
#' @return A \code{TWASWeights} object.
#' @export
TWASWeights <- function(weights, variant_ids, fits = NULL,
                        cv_performance = NULL, standardized = FALSE,
                        molecular_id = character(0), data_type = NULL) {
  obj <- new("TWASWeights",
    weights = weights,
    variant_ids = variant_ids,
    methods = names(weights),
    fits = fits,
    cv_performance = cv_performance,
    standardized = standardized,
    molecular_id = molecular_id,
    data_type = data_type
  )
  validObject(obj)
  obj
}

#' @rdname getWeights
#' @export
setMethod("getWeights", "TWASWeights", function(x, method = NULL) {
  if (is.null(method)) return(x@weights)
  if (!method %in% x@methods)
    stop("Method '", method, "' not found. Available: ",
         paste(x@methods, collapse = ", "))
  x@weights[[method]]
})

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TWASWeights", function(x) x@standardized)

#' @rdname getCVPerformance
#' @export
setMethod("getCVPerformance", "TWASWeights", function(x, method = NULL) {
  if (is.null(method)) return(x@cv_performance)
  x@cv_performance[[method]]
})

#' @rdname getFits
#' @export
setMethod("getFits", "TWASWeights", function(x, method = NULL) {
  if (is.null(method)) return(x@fits)
  x@fits[[method]]
})

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "TWASWeights", function(x) x@methods)

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TWASWeights", function(x) x@variant_ids)

#' @rdname getMolecularId
#' @export
setMethod("getMolecularId", "TWASWeights", function(x) x@molecular_id)

#' @rdname getDataType
#' @export
setMethod("getDataType", "TWASWeights", function(x) x@data_type)

# =============================================================================
# FineMappingResult additional accessors
# =============================================================================

#' @rdname getTrimmedFit
#' @export
setMethod("getTrimmedFit", "FineMappingResult", function(x) x@trimmed_fit)

#' @rdname getVariantNames
#' @export
setMethod("getVariantNames", "FineMappingResult", function(x) x@variant_names)

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "FineMappingResult", function(x) x@top_loci)

# =============================================================================
# AlleleQCResult constructor and accessors
# =============================================================================

#' @title Construct an AlleleQCResult object
#' @description Build an \code{AlleleQCResult} S4 object wrapping the post-QC
#'   harmonized variants and the full per-variant QC diagnostics.
#' @param harmonized_data Data frame of variants retained after allele QC.
#' @param qc_summary Data frame of per-variant diagnostic columns.
#' @return An \code{AlleleQCResult} object.
#' @export
AlleleQCResult <- function(harmonized_data, qc_summary) {
  obj <- new("AlleleQCResult",
      harmonized_data = as.data.frame(harmonized_data),
      qc_summary = as.data.frame(qc_summary))
  validObject(obj)
  obj
}

#' @rdname getHarmonizedData
#' @export
setMethod("getHarmonizedData", "AlleleQCResult", function(x) x@harmonized_data)

#' @rdname getQCSummary
#' @export
setMethod("getQCSummary", "AlleleQCResult", function(x) x@qc_summary)

# =============================================================================
# QCResult constructor and accessors
# =============================================================================

#' @title Construct a QCResult object
#' @description Build a \code{QCResult} S4 object capturing the output of
#'   summary-statistic QC. Validates that \code{ld_data} is an \code{LDData}
#'   or NULL.
#' @param ld_data An \code{LDData} or NULL.
#' @param rss_input List with \code{sumstats}, \code{n}, \code{var_y}.
#' @param preprocess List with \code{sumstats} and \code{ld_data}.
#' @param outlier_number Integer count of LD-mismatch outliers removed.
#' @param skipped Single logical indicating a short-circuit.
#' @param skip_reason Character explanation; defaults to empty.
#' @return A \code{QCResult} object.
#' @export
QCResult <- function(ld_data = NULL,
                     rss_input = list(),
                     preprocess = list(),
                     outlier_number = 0L,
                     skipped = FALSE,
                     skip_reason = "") {
  reason <- if (length(skip_reason) == 0L) "" else as.character(skip_reason)[[1]]
  obj <- new("QCResult",
      ld_data = ld_data,
      rss_input = rss_input,
      preprocess = preprocess,
      outlier_number = as.integer(outlier_number),
      skipped = isTRUE(skipped),
      skip_reason = reason)
  validObject(obj)
  obj
}

#' @rdname getLDData
#' @export
setMethod("getLDData", "QCResult", function(x) x@ld_data)

#' @rdname getRSSInput
#' @export
setMethod("getRSSInput", "QCResult", function(x) x@rss_input)

#' @rdname getPreprocess
#' @export
setMethod("getPreprocess", "QCResult", function(x) x@preprocess)

#' @rdname getOutlierNumber
#' @export
setMethod("getOutlierNumber", "QCResult", function(x) x@outlier_number)

#' @rdname isSkipped
#' @export
setMethod("isSkipped", "QCResult", function(x) x@skipped)

#' @rdname getSkipReason
#' @export
setMethod("getSkipReason", "QCResult", function(x) x@skip_reason)

# =============================================================================
# top_loci GRanges conversion
# =============================================================================

#' @title Convert top_loci to GRanges
#' @description Convert a long-format top_loci data.frame (with variant_id
#'   in chr:pos:A2:A1 format) to a GRanges object with all metadata columns.
#' @param top_loci A data.frame with a variant_id column encoding
#'   chr:pos:A2:A1.
#' @return A GRanges object with all original columns as metadata.
#' @export
top_loci_to_granges <- function(top_loci) {
  if (is.null(top_loci) || nrow(top_loci) == 0) {
    return(GRanges())
  }
  parsed <- parse_variant_id(top_loci$variant_id)
  chr <- paste0("chr", parsed$chrom)
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = parsed$pos, width = 1L)
  )
  mcols(gr) <- DataFrame(top_loci)
  gr
}
