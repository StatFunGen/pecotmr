#' xQTL GWAS Enrichment Analysis
#'
#' This function processes GWAS and xQTL finemapped data files and then computes QTL enrichment.
#' For details on the parameters `piGwas`, `piQtl`, `lambda`, `impN`, and `numThreads`,
#' refer to the documentation of the `computeQtlEnrichment` function.
#'
#' @param xqtlFiles Vector of xQTL RDS file paths.
#' @param gwasFiles Vector of GWAS RDS file paths.
#' @param xqtlFinemappingObj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwasFinemappingObj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtlVarnameObj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwasVarnameObj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param piGwas Optional parameter for GWAS enrichment estimation (see `computeQtlEnrichment`).
#' @param piQtl Optional parameter for xQTL enrichment estimation (see `computeQtlEnrichment`).
#' @param lambda Shrinkage parameter for enrichment computation (see `computeQtlEnrichment`).
#' @param impN Importance parameter for enrichment computation (see `computeQtlEnrichment`).
#' @param numThreads Number of threads for parallel processing (see `computeQtlEnrichment`).
#' @return The output from the computeQtlEnrichment function.
#' @examples
#' gwasFiles <- c("gwas_file1.rds", "gwas_file2.rds")
#' xqtlFiles <- c("xqtl_file1.rds", "xqtl_file2.rds")
#' result <- xqtlEnrichmentWrapper(gwasFiles, xqtlFiles)
#' @export
xqtlEnrichmentWrapper <- function(xqtlFiles, gwasFiles,
                                  xqtlFinemappingObj = NULL, gwasFinemappingObj = NULL,
                                  xqtlVarnameObj = NULL, gwasVarnameObj = NULL,
                                  numGwas = NULL, piQtl = NULL,
                                  lambda = 1.0, impN = 25,
                                  doubleShrinkage = FALSE,
                                  besselCorrection = TRUE,
                                  numThreads = 1) {
  processFinemappedData <- function(xqtlFiles, gwasFiles,
                                    xqtlFinemappingObj = NULL, gwasFinemappingObj = NULL,
                                    xqtlVarnameObj = NULL, gwasVarnameObj = NULL) {
    # Load and process GWAS data
    gwasPipList <- map(gwasFiles, function(file) {
      rawData <- readRDS(file)[[1]]
      gwasData <- if (!is.null(gwasFinemappingObj)) getNestedElement(rawData, gwasFinemappingObj) else rawData
      pip <- gwasData$pip
      if (!is.null(gwasVarnameObj)) names(pip) <- getNestedElement(rawData, gwasVarnameObj)
      pip
    })

    # Check for unique variant names in GWAS pip vectors
    allVariantNames <- unique(unlist(map(gwasPipList, names)))
    if (length(unique(allVariantNames)) != length(allVariantNames)) {
      stop("Non-unique variant names found in GWAS data with different pip values.")
    }
    gwasPip <- unlist(gwasPipList)

    # Process xQTL data
    xqtlData <- lapply(xqtlFiles, function(file) {
      rawData <- readRDS(file)[[1]]
      xqtlData <- tryCatch(
        {
          if (!is.null(xqtlFinemappingObj)) getNestedElement(rawData, xqtlFinemappingObj) else rawData
        },
        error = function(e) {
          return(NULL)
        }
      )
      if (!is.null(xqtlData)) {
        list(
          alpha = xqtlData$alpha,
          pip = setNames(xqtlData$pip, getNestedElement(rawData, xqtlVarnameObj)),
          prior_variance = xqtlData$V
        )
      } else {
        NULL
      }
    })

    # Return results as a list
    return(list(gwas_pip = gwasPip, xqtl_data = xqtlData))
  }

  # Load data
  dat <- processFinemappedData(xqtlFiles, gwasFiles, xqtlFinemappingObj, gwasFinemappingObj, xqtlVarnameObj, gwasVarnameObj)
  # Compute QTL enrichment
  return(computeQtlEnrichment(
    gwasPip = dat$gwas_pip, susieQtlRegions = dat$xqtl_data,
    numGwas = numGwas, piQtl = piQtl,
    lambda = lambda, impN = impN,
    doubleShrinkage = doubleShrinkage,
    besselCorrection = besselCorrection,
    numThreads = numThreads
  ))
}

#' Function to filter and order colocalization results
#' @noRd
filterAndOrderColocResults <- function(colocResultsFil) {
  # Ensure the input has more than one column
  if (ncol(colocResultsFil) <= 1) {
    stop("Insufficient number of columns in colocalization results")
  }

  csNum <- ncol(colocResultsFil) - 1
  map(seq_len(csNum), function(n) {
    colocResultsFil[, c(1, n + 1)] %>% .[order(.[, 2], decreasing = TRUE), ]
  })
}

#' Function to calculate cumulative sum
#' @noRd
calculateCumsum <- function(colocResults) {
  cumsum(colocResults[, 2])
}

#' Load LD matrix for a set of variants, narrowing the region and aligning names.
#' @importFrom stringr str_split
#' @noRd
extractLdForVariants <- function(ldMetaFilePath, analysisRegion, variants) {
  varPos <- as.numeric(str_split(variants, ":", simplify = TRUE)[, 2])
  chr <- str_split(analysisRegion, ":", simplify = TRUE)[, 1]
  regionNarrow <- paste0(chr, ":", min(varPos), "-", max(varPos))
  ldData <- loadLdMatrix(ldMetaFilePath, region = regionNarrow,
                         returnGenotype = "auto")
  if (!is(ldData, "LdData")) {
    stop("loadLdMatrix must return an LdData object")
  }
  ldVariants <- getVariantIds(ldData)
  hasGeno <- hasGenotypes(ldData)
  aligned <- alignVariantNames(ldVariants, variants)
  # When genotypes available, compute R only for the needed variant subset
  if (hasGeno) {
    X <- getGenotypes(ldData)
    colnames(X) <- aligned$aligned_variants
    xSub <- X[, variants, drop = FALSE]
    ldMatrix <- computeLd(xSub, method = "sample")
  } else {
    ldMatrix <- getCorrelation(ldData)
    colnames(ldMatrix) <- rownames(ldMatrix) <- aligned$aligned_variants
    ldMatrix <- ldMatrix[variants, variants]
  }
  ldMatrix
}

#' Function to calculate purity
#' @noRd
calculatePurity <- function(variants, extLd, squared = FALSE) {
  # This is a placeholder for calculating purity, adjust as per your actual function
  purity <- matrix(susieR:::get_purity(variants, Xcorr = extLd, squared), 1, 3)
  purity
}

#' Main processing function
#' This function is designed to summarize coloc results based on the following criteria:
#' 1. Among the colocalized variant pairs, PPH4 has the highest value compared to PPH0-PPH3.
#' 2. PPH4 exceeds threshold, default as 0 since we advocate not using PPH4 concept but rather use CoS
#' 3. We aggregate variants and cumulatively sum their PPH4 values to form a credible set until the threshold, default as 0.95.
#' 4. The cs's purity is computed with the `get_purity` function from the `gaow/susieR` package, and the same purity criteria are employed to filter the credibility set.
#' @noRd
processColocResults <- function(colocResult, ldMetaFilePath, analysisRegion, pph4Thres = 0, coverage = 0.95, minAbsCorr = 0.8, nullIndex = 0, colocIndex = "PP.H4.abf") {
  # Extract PIP values from colocResult summary
  colocSummary <- as.data.frame(colocResult$summary)
  colocPip <- colocSummary[, grepl("PP", colnames(colocSummary))]

  # Filter and extract relevant columns from colocResult results
  # PP.H4 is highest and > 0.8
  colocResultsDf <- as.data.frame(colocResult$results)
  colocFilter <- apply(colocPip, 1, function(row) {
    maxIndex <- which.max(row)
    maxValue <- row[maxIndex]
    return(maxValue > pph4Thres && colnames(colocPip)[maxIndex] == colocIndex)
  })

  colocRes <- list()

  if (sum(colocFilter) > 0) {
    colocResultsFil <- colocResultsDf[, c(1, which(colocFilter) + 1), drop = FALSE]
    colocSummaryFil <- colocSummary[which(colocFilter), , drop = FALSE]

    # prepare to calculate purity
    orderedResults <- filterAndOrderColocResults(colocResultsFil)
    cs <- map(orderedResults, function(res) {
      csm <- calculateCumsum(res)
      res[, 1][1:min(which(csm > coverage))]
    })

    purity <- map_dfr(seq_along(cs), function(n) {
      variants <- normalizeVariantId(cs[[n]])
      if (nullIndex > 0 && nullIndex %in% variants) {
        data.frame(min.abs.corr = -9, mean.abs.corr = -9, median.abs.corr = -9)
      } else {
        extLd <- extractLdForVariants(ldMetaFilePath, analysisRegion, variants)
        p <- calculatePurity(variants, extLd)
        data.frame(min.abs.corr = p[1, 1], mean.abs.corr = p[1, 2], median.abs.corr = p[1, 3])
      }
    })
    isPure <- which(purity[, 1] >= minAbsCorr)

    # Finalize the result
    if (length(isPure) > 0) {
      cs <- cs[isPure]
      purity <- purity[isPure, ]
      trueSummary <- colocSummaryFil[isPure, ]
      colocRes$sets <- list(cs = cs, purity = purity, true_summary = trueSummary)
    }
  } else {
    message("Coloc results did not find any variants that satisfy the condition of PP.H4 being the highest value and > ", pph4Thres)
    colocRes$sets <- list(cs = NULL)
  }

  return(colocRes)
}

# Extract and filter an LBF matrix from a finemapped data object.
# @noRd
.extractLbfMatrix <- function(rawData, finemappingObj, varnameObj,
                              filterLbfCs, filterLbfCsSecondary, priorTol) {
  fmData <- if (!is.null(finemappingObj)) {
    tryCatch(getNestedElement(rawData, finemappingObj),
      error = function(e) {
        message(paste("no", finemappingObj[2], "in", finemappingObj[1]))
        NULL
      }
    )
  } else {
    rawData
  }
  if (is.null(fmData)) return(NULL)

  lbfMatrix <- as.data.frame(fmData$lbf_variable)
  # fSuSiE has a different structure
  if (is.null(lbfMatrix) || nrow(lbfMatrix) == 0) {
    fsusieLbf <- NULL
    if (is.list(rawData) && length(rawData) >= 1 && is.list(rawData[[1]])) {
      fsusieLbf <- rawData[[1]]$fsusie_result$lBF
    }
    if (is.list(fsusieLbf) && length(fsusieLbf) > 0) {
      lbfMatrix <- do.call(rbind, fsusieLbf) %>% as.data.frame()
      if (nrow(lbfMatrix) > 0) message("This is a fSuSiE case")
    }
  }

  # Filter rows
  if (filterLbfCs && is.null(filterLbfCsSecondary)) {
    lbfMatrix <- lbfMatrix[fmData$sets$cs_index, , drop = FALSE]
  } else if (!is.null(filterLbfCsSecondary)) {
    lbfMatrix <- lbfMatrix[getFilterLbfIndex(fmData, coverage = filterLbfCsSecondary), , drop = FALSE]
  } else {
    if ("V" %in% names(fmData)) {
      lbfMatrix <- lbfMatrix[fmData$V > priorTol, , drop = FALSE]
    } else {
      message("No V found in original data.")
    }
  }

  # Set variant names and remove NA columns
  if (!is.null(varnameObj)) colnames(lbfMatrix) <- getNestedElement(rawData, varnameObj)
  lbfMatrix <- lbfMatrix[, !is.na(colnames(lbfMatrix))]

  list(lbf_matrix = lbfMatrix, fm_data = fmData)
}

# Extract LBF matrix from an rssAnalysisPipeline result object.
# Unlike .extractLbfMatrix which navigates RDS-loaded nested lists,
# this works directly with the in-memory pipeline output structure.
# @noRd
.extractLbfFromPipelineResult <- function(pipelineResult,
                                          filterLbfCs, filterLbfCsSecondary,
                                          priorTol) {
  methodNames <- setdiff(names(pipelineResult), "rss_data_analyzed")
  if (length(methodNames) == 0) return(NULL)

  methodResult <- pipelineResult[[methodNames[1]]]
  fmResult <- methodResult$finemapping_result
  if (is.null(fmResult) || !is(fmResult, "FineMappingResult")) return(NULL)
  fmData <- getTrimmedFit(fmResult)
  variantNames <- getVariantNames(fmResult)
  if (is.null(fmData) || is.null(fmData$lbf_variable)) return(NULL)

  lbfMatrix <- as.data.frame(fmData$lbf_variable)

  # Row filtering — same logic as .extractLbfMatrix
  if (filterLbfCs && is.null(filterLbfCsSecondary)) {
    lbfMatrix <- lbfMatrix[fmData$sets$cs_index, , drop = FALSE]
  } else if (!is.null(filterLbfCsSecondary)) {
    lbfMatrix <- lbfMatrix[getFilterLbfIndex(fmData, coverage = filterLbfCsSecondary), , drop = FALSE]
  } else if ("V" %in% names(fmData)) {
    lbfMatrix <- lbfMatrix[fmData$V > priorTol, , drop = FALSE]
  }

  if (!is.null(variantNames) && length(variantNames) == ncol(lbfMatrix)) {
    colnames(lbfMatrix) <- variantNames
  }
  lbfMatrix <- lbfMatrix[, !is.na(colnames(lbfMatrix))]
  list(lbf_matrix = lbfMatrix, fm_data = fmData)
}

# Save inline fine-mapping result to disk in a format compatible with the
# file-based reading path (readRDS(file)[[1]] + gwasFinemappingObj/gwasVarnameObj).
# @noRd
.saveFinemappingResult <- function(pipelineResult, savePath) {
  if (is.null(savePath) || is.null(pipelineResult)) return(invisible(NULL))
  methodNames <- setdiff(names(pipelineResult), "rss_data_analyzed")
  if (length(methodNames) == 0) return(invisible(NULL))
  methodResult <- pipelineResult[[methodNames[1]]]
  fmResult <- methodResult$finemapping_result
  if (is.null(fmResult) || !is(fmResult, "FineMappingResult")) return(invisible(NULL))
  saveData <- list(
    susie_fit = getTrimmedFit(fmResult),
    variant_names = getVariantNames(fmResult)
  )
  saveRDS(list(saveData), savePath)
  message("Fine-mapping result saved to: ", savePath,
          "\n  Reuse with: gwasFiles = '", savePath,
          "', gwasFinemappingObj = 'susie_fit', gwasVarnameObj = 'variant_names'")
  invisible(savePath)
}

#' Colocalization Analysis Wrapper
#'
#' Processes xQTL and GWAS finemapped data for colocalization analysis.
#' GWAS data can come from pre-computed RDS files or from inline fine-mapping
#' via \code{\link{rssAnalysisPipeline}}.
#'
#' @param xqtlFile Path to the xQTL RDS file.
#' @param gwasFiles Vector of paths to GWAS RDS files. Required when
#'   \code{runFinemapping = FALSE}. Ignored when \code{runFinemapping = TRUE}.
#' @param xqtlFinemappingObj Optional path in xQTL RDS to the finemapping object.
#' @param gwasFinemappingObj Optional path in GWAS RDS to the finemapping object.
#' @param xqtlVarnameObj Optional path in xQTL RDS to variant names.
#' @param gwasVarnameObj Optional path in GWAS RDS to variant names.
#' @param xqtlRegionObj Optional path in xQTL RDS to region info.
#' @param gwasRegionObj Optional path in GWAS RDS to region info.
#' @param filterLbfCs Logical. Filter LBF rows by credible set index.
#' @param filterLbfCsSecondary Coverage for secondary LBF filtering.
#' @param priorTol Minimum prior variance to retain an effect (default 1e-9).
#' @param p1 Prior probability a SNP is associated with trait 1 (default 1e-4).
#' @param p2 Prior probability a SNP is associated with trait 2 (default 1e-4).
#' @param p12 Prior probability a SNP is associated with both traits (default 5e-6).
#' @param runFinemapping Logical. If TRUE, run GWAS fine-mapping inline via
#'   \code{\link{rssAnalysisPipeline}}. Default FALSE.
#' @param sumstatPath Path to GWAS summary statistics file. Required when
#'   \code{runFinemapping = TRUE}.
#' @param columnFilePath Path to column mapping file for summary statistics.
#' @param ldData LD reference data (LdData object or list). Required when
#'   \code{runFinemapping = TRUE}.
#' @param nSample Sample size for GWAS.
#' @param nCase Number of cases for binary traits.
#' @param nControl Number of controls for binary traits.
#' @param region Genomic region string (e.g., "chr1:1000-2000").
#' @param zMismatchQc Z-score / LD-mismatch QC selector forwarded to
#'   \code{\link{rssAnalysisPipeline}}: "slalom", "dentist", or "none".
#'   Default "slalom". (Hard rename of the former \code{zMismatchQc}; no alias.)
#' @param finemappingMethod Fine-mapping method. Default "susie_rss".
#' @param finemappingOpts List of fine-mapping options passed to
#'   \code{\link{rssAnalysisPipeline}}.
#' @param impute Logical. Run RAISS imputation. Default TRUE.
#' @param imputeOpts List of imputation options.
#' @param saveFinemappingPath Path to save fine-mapping result as RDS. The
#'   saved file can be reused via \code{gwasFiles} with
#'   \code{gwasFinemappingObj = "susie_fit"} and
#'   \code{gwasVarnameObj = "variant_names"}.
#' @param returnFinemapping Logical. If TRUE and \code{runFinemapping = TRUE},
#'   include full fine-mapping result under \code{$gwas_finemapping}.
#' @param ... Additional arguments (currently unused).
#' @return A list containing the coloc results and the summarized sets.
#' @seealso \code{\link{rssAnalysisPipeline}}, \code{\link{colocPostProcessor}}
#' @importFrom dplyr bind_rows mutate across
#' @importFrom tidyr replace_na
#' @importFrom coloc coloc.bf_bf
#' @importFrom purrr map map_dfr
#' @export
colocWrapper <- function(xqtlFile, gwasFiles = NULL,
                         xqtlFinemappingObj = NULL, xqtlVarnameObj = NULL, xqtlRegionObj = NULL,
                         gwasFinemappingObj = NULL, gwasVarnameObj = NULL, gwasRegionObj = NULL,
                         filterLbfCs = FALSE, filterLbfCsSecondary = NULL,
                         priorTol = 1e-9, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                         runFinemapping = FALSE,
                         sumstatPath = NULL, columnFilePath = NULL,
                         ldData = NULL,
                         nSample = 0, nCase = 0, nControl = 0,
                         region = NULL,
                         zMismatchQc = "slalom",
                         finemappingMethod = "susie_rss",
                         finemappingOpts = list(
                           L = 20, L_greedy = 5,
                           coverage = c(0.95, 0.7, 0.5),
                           signal_cutoff = 0.025,
                           min_abs_corr = 0.8
                         ),
                         impute = TRUE,
                         imputeOpts = list(rcond = 0.01, R2_threshold = 0.6,
                                           minimum_ld = 5, lamb = 0.01),
                         saveFinemappingPath = NULL,
                         returnFinemapping = FALSE,
                         ...) {
  # --- Input validation ---
  if (!runFinemapping && is.null(gwasFiles)) {
    stop("Either set runFinemapping = TRUE with GWAS sumstat inputs, or provide gwasFiles paths to pre-computed results.")
  }
  if (runFinemapping && !is.null(gwasFiles)) {
    warning("Both runFinemapping = TRUE and gwasFiles provided. Inline fine-mapping will be used; gwasFiles ignored.")
    gwasFiles <- NULL
  }
  if (runFinemapping) {
    if (is.null(sumstatPath)) stop("sumstatPath is required when runFinemapping = TRUE.")
    if (is.null(ldData)) stop("ldData is required when runFinemapping = TRUE.")
  }

  gwasPipelineResult <- NULL

  if (runFinemapping) {
    # --- Inline fine-mapping path: QC runs inside rssAnalysisPipeline ---
    gwasPipelineResult <- rssAnalysisPipeline(
      sumstatPath = sumstatPath, columnFilePath = columnFilePath,
      ldData = ldData,
      nSample = nSample, nCase = nCase, nControl = nControl,
      region = region,
      zMismatchQc = zMismatchQc, finemappingMethod = finemappingMethod,
      finemappingOpts = finemappingOpts,
      impute = impute, imputeOpts = imputeOpts
    )

    # Save to disk before extraction (useful even if extraction fails)
    .saveFinemappingResult(gwasPipelineResult, saveFinemappingPath)

    gwasExtracted <- .extractLbfFromPipelineResult(
      gwasPipelineResult, filterLbfCs, filterLbfCsSecondary, priorTol
    )
    if (is.null(gwasExtracted)) {
      colocRes <- list("No GWAS fine-mapping results produced by inline pipeline.")
      result <- c(colocRes, analysis_region = region)
      if (returnFinemapping) result$gwas_finemapping <- gwasPipelineResult
      return(result)
    }
    combinedGwasLbfMatrix <- gwasExtracted$lbf_matrix %>%
      as.data.frame() %>% mutate(across(everything(), ~ replace_na(., 0)))
  } else {
    # --- File-based path (unchanged) ---
    gwasLbfMatrices <- map(gwasFiles, function(file) {
      rawData <- readRDS(file)[[1]]
      .extractLbfMatrix(rawData, gwasFinemappingObj, gwasVarnameObj,
                        filterLbfCs, filterLbfCsSecondary, priorTol)$lbf_matrix
    })
    combinedGwasLbfMatrix <- bind_rows(gwasLbfMatrices) %>%
      mutate(across(everything(), ~ replace_na(., 0)))
  }

  # Process xQTL data
  xqtlRawData <- readRDS(xqtlFile)[[1]]
  xqtlExtracted <- .extractLbfMatrix(xqtlRawData, xqtlFinemappingObj, xqtlVarnameObj,
                                     filterLbfCs, filterLbfCsSecondary, priorTol)

  if (!is.null(xqtlExtracted)) {
    xqtlLbfMatrix <- xqtlExtracted$lbf_matrix
    if (nrow(combinedGwasLbfMatrix) > 0 && nrow(xqtlLbfMatrix) > 0) {
      colnames(xqtlLbfMatrix) <- alignVariantNames(colnames(xqtlLbfMatrix), colnames(combinedGwasLbfMatrix))$aligned_variants
      commonColnames <- intersect(colnames(xqtlLbfMatrix), colnames(combinedGwasLbfMatrix))

      numDroppedCols <- ncol(xqtlLbfMatrix) - length(commonColnames)
      if (numDroppedCols > 0) {
        message("Number of columns dropped from xQTL matrix: ", numDroppedCols)
      }

      xqtlLbfMatrix <- xqtlLbfMatrix[, commonColnames, drop = FALSE] %>% as.matrix()
      combinedGwasLbfMatrix <- combinedGwasLbfMatrix[, commonColnames, drop = FALSE] %>% as.matrix()

      convertToString <- function(df) paste0("chr", df$chrom, ":", df$start, "-", df$end)
      analysisRegionOut <- if (!is.null(xqtlRegionObj)) {
        getNestedElement(xqtlRawData, xqtlRegionObj) %>% convertToString()
      } else {
        region
      }

      colocRes <- coloc.bf_bf(xqtlLbfMatrix, combinedGwasLbfMatrix, p1 = p1, p2 = p2, p12 = p12)
    } else {
      colocRes <- list("No coloc results due to the absence of a GWAS log Bayes factor matrix filtered by prior tolerance.")
      analysisRegionOut <- region
    }
  } else {
    colocRes <- list(paste("no", xqtlFinemappingObj[2], "in", xqtlFinemappingObj[1]))
    analysisRegionOut <- region
  }

  result <- c(colocRes, analysis_region = analysisRegionOut)
  if (returnFinemapping && !is.null(gwasPipelineResult)) {
    result$gwas_finemapping <- gwasPipelineResult
  }
  return(result)
}

#' colocPostProcessor function
#' @param colocRes coloc results from coloc.susie.
#' @param ldMetaFilePath Path to the metadata of LD reference.
#' @param analysisRegion Path to the analysis region of coloc result.
#' @return A list containing the coloc results and post processed coloc sets.
#' @export
colocPostProcessor <- function(colocRes, ldMetaFilePath = NULL, analysisRegion = NULL, ...) {
  if (!is.null(ldMetaFilePath)) {
    if (is.null(analysisRegion)) {
      stop("ldMetaFilePath is provided but analysisRegion is not provided. Please provide analysisRegion for purity filter.")
    }
    # Perform purity filter using ldMetaFilePath and analysisRegion
    colocRes <- c(colocRes, processColocResults(colocRes, ldMetaFilePath, analysisRegion = analysisRegion))
  } else {
    if (!is.null(analysisRegion)) {
      warning("analysisRegion is provided but will not be used as ldMetaFilePath is not provided.")
    }
    warning("ldMetaFilePath not provided. Purity filter cannot be applied.")
  }
  return(colocRes)
}

#' @export

# In practice, analysis will contain two lines:
# res <- colocWrapper(...)
# postProcessedRes <- colocPostProcessor
