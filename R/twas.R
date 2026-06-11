# Internal: RAISS-impute GWAS z-scores for LD-sketch variants missing from
# the harmonized GWAS subset. Returns the (possibly widened) sumstats data
# frame. Imputed rows fill `z` from RAISS; `beta` becomes the imputed z and
# `se` becomes 1 when those columns are present in the input. Other columns
# are filled with NA. Imputed variants with R^2 below the threshold are
# dropped by RAISS's internal filter.
imputeMissingGwasForSketch <- function(gwasDataSumstats, sketchRefPanel,
                                       sketchX, imputeOpts, contextLabel = "") {
  missingIds <- setdiff(sketchRefPanel$variant_id, gwasDataSumstats$variant_id)
  if (length(missingIds) == 0) return(gwasDataSumstats)

  refCols <- c("chrom", "pos", "variant_id", "A1", "A2")
  if (!all(refCols %in% colnames(sketchRefPanel))) {
    warning("imputeMissingGwasForSketch: sketch refPanel missing required columns; skipping imputation.")
    return(gwasDataSumstats)
  }
  if (!all(refCols %in% colnames(gwasDataSumstats)) || !"z" %in% colnames(gwasDataSumstats)) {
    warning("imputeMissingGwasForSketch: gwas sumstats missing required columns; skipping imputation.")
    return(gwasDataSumstats)
  }

  # RAISS requires inputs sorted by position (within each chromosome)
  refSorted <- sketchRefPanel[order(sketchRefPanel$chrom, sketchRefPanel$pos), refCols, drop = FALSE]
  knownSorted <- gwasDataSumstats[order(gwasDataSumstats$chrom, gwasDataSumstats$pos), c(refCols, "z"), drop = FALSE]
  # Reorder genotype matrix columns to match the sorted refPanel
  vidOrder <- match(refSorted$variant_id, colnames(sketchX))
  vidOrder <- vidOrder[!is.na(vidOrder)]
  sketchXSorted <- sketchX[, vidOrder, drop = FALSE]
  # Translate snake_case imputeOpts keys to camelCase raiss() arguments.
  imputeOptsRenamed <- imputeOpts
  if ("R2_threshold" %in% names(imputeOptsRenamed)) {
    imputeOptsRenamed$r2Threshold <- imputeOptsRenamed$R2_threshold
    imputeOptsRenamed$R2_threshold <- NULL
  }
  if ("minimum_ld" %in% names(imputeOptsRenamed)) {
    imputeOptsRenamed$minimumLd <- imputeOptsRenamed$minimum_ld
    imputeOptsRenamed$minimum_ld <- NULL
  }
  raissArgs <- c(list(
    refPanel = refSorted,
    knownZscores = knownSorted,
    genotypeMatrix = sketchXSorted,
    verbose = FALSE
  ), imputeOptsRenamed)
  raissOut <- tryCatch(do.call(raiss, raissArgs),
                       error = function(e) {
                         warning(sprintf("RAISS missing-variant imputation failed (%s): %s",
                                         contextLabel, e$message))
                         NULL
                       })
  if (is.null(raissOut) || is.null(raissOut$result_filter)) return(gwasDataSumstats)

  imputedDf <- raissOut$result_filter
  newRows <- imputedDf[!imputedDf$variant_id %in% gwasDataSumstats$variant_id, , drop = FALSE]
  if (nrow(newRows) == 0) return(gwasDataSumstats)

  added <- newRows[, c("variant_id", "chrom", "pos", "A1", "A2", "z"), drop = FALSE]
  if ("beta" %in% colnames(gwasDataSumstats)) added$beta <- newRows$z
  if ("se"   %in% colnames(gwasDataSumstats)) added$se   <- 1
  for (col in setdiff(colnames(gwasDataSumstats), colnames(added))) {
    added[[col]] <- NA
  }
  added <- added[, colnames(gwasDataSumstats), drop = FALSE]
  message(sprintf("RAISS imputed %d missing GWAS variants (%s).", nrow(added), contextLabel))
  rbind(gwasDataSumstats, added)
}

#' Function to perform allele flip QC and harmonization on the weights and GWAS against LD for a region.
#' FIXME: GWAS loading function from Haochen for both tabix & column-mapping yml application
#'
#' Function Conditions:
#' - processes data in the format of either the output from load_twas_weights/generate_twas_db or
#'   refined_twas_weights_data from twas pipeline.
#' - For the first format, we expect there is only one gene/events's information, that can be accessed through `region_info_obj`
#'   and refined_twas_weights_data contains per region multiple gene/event's refined weights data.
#'
#' Main Steps:
#' 1. allele QC for TWAS weights against the LD meta
#' 2. allele QC for GWA summary stats against the LD meta
#' 3. adjust susie/mvsusie weights based on the overlap variants
#'
#' @param twasWeightsData List of list of twas weights output from from generate_twas_db function.
#' @param gwasMetaFile A file path for a dataframe table with column of "study_id", "chrom" (integer), "file_path",
#' "column_mapping_file". Each file in "file_path" column is tab-delimited dataframe of GWAS summary statistics with column name
#' "chrom" (or #chrom" if tabix-indexed), "pos", "A2", "A1".
#' @param ldMetaFilePath Path to LD reference: either a PLINK2/PLINK1 prefix, or a tab-delimited
#'   metadata file with columns "#chrom", "start", "end", "path" (auto-detected).
#' @param ldReferenceSampleSize Sample size of the LD reference panel (integer). Required.
#'   Used to compute per-variant variance as 2*p*(1-p)*n/(n-1). For ADSP R4, use 17000.
#' @param imputeMissing Logical. When \code{TRUE}, RAISS imputes GWAS z-scores
#'   for variants that are present in the LD sketch but missing from the GWAS
#'   summary statistics. This widens GWAS coverage so weight variants with LD
#'   neighbors but no GWAS hit are no longer silently dropped at the
#'   weight-vs-GWAS intersection. Default \code{FALSE}.
#' @param imputeOpts Named list of RAISS imputation parameters. Used when
#'   \code{imputeMissing = TRUE}. Defaults:
#'   \code{list(rcond = 0.01, r2Threshold = 0.6, minimumLd = 5, lamb = 0.01)}.
#'   Imputed variants with \code{R2 < r2Threshold} are dropped.
#' @return A list of list for harmonized weights and dataframe of gwas summary statistics that is add to the original input of
#' twasWeightsData under each context.
#' @importFrom vroom vroom
#' @importFrom readr parse_number
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges findOverlaps start end reduce
#' @export
harmonizeTwas <- function(twasWeightsData, ldMetaFilePath, gwasMetaFile,
                          ldReferenceSampleSize, columnFilePath = NULL, commentString = "#",
                          imputeMissing = FALSE,
                          imputeOpts = list(rcond = 0.01, r2Threshold = 0.6,
                                             minimumLd = 5, lamb = 0.01)) {
  # Step 1: Normalize twasWeightsData -- accept bare TwasWeights or wrapper lists
  molecularIds <- names(twasWeightsData)
  for (molId in molecularIds) {
    entry <- twasWeightsData[[molId]]
    if (is(entry, "TwasWeights")) {
      # Already a bare TwasWeights, use directly
    } else if (is.list(entry) && is(entry$twas_weights, "TwasWeights")) {
      # Wrapper list with snake-case key -- extract the TwasWeights
      twasWeightsData[[molId]] <- entry$twas_weights
    } else if (is.list(entry) && is(entry$twasWeights, "TwasWeights")) {
      # Wrapper list with camelCase key -- extract the TwasWeights
      twasWeightsData[[molId]] <- entry$twasWeights
    } else {
      stop("Each element of twasWeightsData must be a TwasWeights S4 object ",
           "or a list with a $twas_weights TwasWeights element")
    }
  }
  firstTw <- twasWeightsData[[1]]
  chrom <- as.integer(parse_number(gsub(":.*$", "", getVariantIds(firstTw)[1])))
  gwasMetaDf <- as.data.frame(vroom(gwasMetaFile))
  gwasFiles <- unique(gwasMetaDf$file_path[gwasMetaDf$chrom == chrom])
  names(gwasFiles) <- unique(gwasMetaDf$study_id[gwasMetaDf$chrom == chrom])
  results <- list()

  # Per-gene loop: each gene loads its own LD sketch independently
  for (molecularId in molecularIds) {
    tw <- twasWeightsData[[molecularId]]
    molRes <- list(chrom = chrom, variant_names = list())
    molRes[["data_type"]] <- getDataType(tw)
    contexts <- getMethodNames(tw)

    # Step 2: Build gene window from all contexts' variant positions
    allWeightVariants <- getVariantIds(tw)
    variantPositions <- parseVariantId(allWeightVariants)$pos
    geneRegion <- paste0(chrom, ":", min(variantPositions), "-", max(variantPositions))

    # Step 3: Load LD sketch for this gene's window and compute SVD
    sketch <- loadLdSketch(ldMetaFilePath, geneRegion, nSample = ldReferenceSampleSize)
    sketchX <- getGenotypes(sketch)
    sketchRefPanel <- getRefPanel(sketch)
    sketchVariantIds <- getVariantIds(sketch)
    sketchN <- nrow(sketchX)
    xStd <- standardizeGenotypeHwe(sketchX, sketchRefPanel$allele_freq)
    svdResult <- safeSvd(xStd, tol = 0)

    # Warn when weight variants have no LD-reference counterpart at all
    # (cannot be imputed by RAISS; will be dropped at the weights-vs-sketch step).
    weightNoLd <- setdiff(allWeightVariants, sketchVariantIds)
    if (length(weightNoLd) > 0) {
      warning(sprintf(
        "harmonizeTwas: %d of %d weight variants for %s have no LD-reference counterpart and will be dropped.",
        length(weightNoLd), length(allWeightVariants), molecularId
      ))
    }

    # Step 4: Harmonize GWAS and weights against sketch variants
    for (study in names(gwasFiles)) {
      gwasFile <- gwasFiles[study]
      gwasDataSumstats <- harmonizeGwas(gwasFile, queryRegion = geneRegion,
                                        sketchVariantIds, c("beta", "z"),
                                        matchMinProp = 0, columnFilePath = columnFilePath,
                                        commentString = commentString)
      if (is.null(gwasDataSumstats)) next

      # Optional RAISS imputation: fill GWAS z-scores for sketch variants
      # absent from the harmonized GWAS. Widens GWAS so the downstream
      # weight-vs-GWAS intersection no longer drops weight variants that
      # have LD neighbors but no observed GWAS hit.
      if (isTRUE(imputeMissing)) {
        gwasDataSumstats <- imputeMissingGwasForSketch(
          gwasDataSumstats = gwasDataSumstats,
          sketchRefPanel = sketchRefPanel,
          sketchX = sketchX,
          imputeOpts = imputeOpts,
          contextLabel = sprintf("study=%s, gene=%s", study, molecularId)
        )
      }

      for (context in contexts) {
        weightsMatrix <- getWeights(tw, context)
        originalWeightVariants <- rownames(weightsMatrix)

        # Harmonize weights against sketch reference
        weightsMatrix <- cbind(variantIdToDf(rownames(weightsMatrix)), weightsMatrix)
        weightsMatrixQced <- matchRefPanel(weightsMatrix, sketchVariantIds,
          colnames(weightsMatrix)[!colnames(weightsMatrix) %in% c("chrom", "pos", "A2", "A1")],
          matchMinProp = 0
        )
        qcedData <- getHarmonizedData(weightsMatrixQced)
        weightsMatrixSubset <- as.matrix(qcedData[, !colnames(qcedData) %in% c(
          "chrom", "pos", "A2", "A1", "variant_id", "variants_id_original"
        ), drop = FALSE])
        rownames(weightsMatrixSubset) <- qcedData$variant_id

        # Ensure consistent chr prefix convention before intersecting
        chrMatched <- ensureChrMatch(gwasDataSumstats$variant_id, sketchVariantIds)
        gwasDataSumstats$variant_id <- chrMatched$ids_a
        rownames(weightsMatrixSubset) <- ensureChrMatch(rownames(weightsMatrixSubset), gwasDataSumstats$variant_id)$ids_a
        weightsMatrixSubset <- weightsMatrixSubset[rownames(weightsMatrixSubset) %in% gwasDataSumstats$variant_id, , drop = FALSE]
        if (nrow(weightsMatrixSubset) == 0) next
        postqcWeightVariants <- rownames(weightsMatrixSubset)

        # Step 5: adjust SuSiE weights based on available variants
        twWeightsCtx <- getWeights(tw, context)
        if ("susie_weights" %in% colnames(twWeightsCtx)) {
          # For adjustSusieWeights, wrap TwasWeights in the list format it expects
          molDataForAdjust <- list(
            susie_results = getFits(tw),
            weights = getWeights(tw),
            variant_names = lapply(getWeights(tw), function(w) if (is.matrix(w)) rownames(w) else names(w))
          )
          adjustedSusieWeights <- adjustSusieWeights(molDataForAdjust,
            keepVariants = postqcWeightVariants, runAlleleQc = TRUE,
            variableNameObj = c("variant_names", context),
            susieObj = c("susie_results", context),
            twasWeightsTable = c("weights", context), postqcWeightVariants, matchMinProp = 0
          )
          weightsMatrixSubset <- cbind(
            susie_weights = setNames(adjustedSusieWeights$adjusted_susie_weights, adjustedSusieWeights$remained_variants_ids),
            weightsMatrixSubset[adjustedSusieWeights$remained_variants_ids, !colnames(weightsMatrixSubset) %in% "susie_weights", drop = FALSE]
          )
          susieResults <- getFits(tw, context)
          susieIntermediate <- susieResults[c("pip", "cs_variants", "cs_purity")]
          names(susieIntermediate[["pip"]]) <- originalWeightVariants # original variants not yet qced
          pip <- susieIntermediate[["pip"]]
          pipQced <- matchRefPanel(cbind(parseVariantId(names(pip)), pip), sketchVariantIds, "pip", matchMinProp = 0)
          pipQcedDf <- getHarmonizedData(pipQced)
          susieIntermediate[["pip"]] <- abs(pipQcedDf$pip)
          names(susieIntermediate[["pip"]]) <- pipQcedDf$variant_id
          susieIntermediate[["cs_variants"]] <- lapply(susieIntermediate[["cs_variants"]], function(x) {
            variantQc <- matchRefPanel(x, sketchVariantIds, matchMinProp = 0)
            variantQcDf <- getHarmonizedData(variantQc)
            variantQcDf$variant_id[variantQcDf$variant_id %in% postqcWeightVariants]
          })
          molRes[["susie_weights_intermediate_qced"]][[context]] <- susieIntermediate
        }
        rm(weightsMatrix)

        if (nrow(weightsMatrixSubset) == 0) {
          warning("weightsMatrixSubset is empty. Skipping this context.")
          next
        }
        molRes[["variant_names"]][[context]][[study]] <- rownames(weightsMatrixSubset)

        # Step 6: scale weights by variance (from sketch ref_panel)
        # RSS/standardized weights are already on the correlation scale and
        # do not need sqrt(variance) scaling.
        isStandardized <- isTRUE(getStandardized(tw))
        if (isStandardized) {
          scaled <- weightsMatrixSubset
        } else {
          variance <- sketchRefPanel$variance[match(rownames(weightsMatrixSubset), sketchRefPanel$variant_id)]
          scaled <- weightsMatrixSubset * sqrt(variance)
        }
        molRes[["weights_qced"]][[context]][[study]] <- list(scaled_weights = scaled, weights = weightsMatrixSubset)
      }
      # Combine GWAS sumstats for this study (filter to variants used by any context)
      usedVariants <- unique(findData(molRes[["variant_names"]], c(2, study)))
      if (!is.null(usedVariants)) {
        gwasSubset <- gwasDataSumstats[gwasDataSumstats$variant_id %in% usedVariants, , drop = FALSE]
        molRes[["gwas_qced"]][[study]] <- rbind(molRes[["gwas_qced"]][[study]], gwasSubset)
        gwasQced <- molRes[["gwas_qced"]][[study]]
        molRes[["gwas_qced"]][[study]] <- gwasQced[!duplicated(gwasQced[, c("variant_id", "z")]), ]
      }
    }

    twasWeightsData[[molecularId]] <- NULL
    # Store SVD components for this gene
    if (is.null(molRes[["gwas_qced"]]) || length(molRes[["gwas_qced"]]) == 0) {
      results[[molecularId]] <- NULL
    } else {
      molRes[["svd_V"]] <- svdResult$v
      molRes[["svd_D"]] <- svdResult$d
      molRes[["n_sketch"]] <- sketchN
      molRes[["ld_variant_ids"]] <- sketchVariantIds
      results[[molecularId]] <- molRes
    }
  }
  return(list(twas_data_qced = results, ref_panel = sketchRefPanel))
}

#' Harmonize GWAS Summary Statistics
#' perform harmonization on gwas summary statistics for a chromosome data or specific queried region
#' @param gwasFile A string for the file path of gwas summary statistics file that is already tabix indexed
#' @param queryRegion A string for region of query for tabix-indexed gwas summary statistics file in the format of chr:start-end
#' @noRd
#' @export
harmonizeGwas <- function(gwasFile, queryRegion, ldVariants, colToFlip=NULL, matchMinProp=0, columnFilePath=NULL, commentString="#"){
    if(is.null(gwasFile)| is.na(gwasFile)) stop("No GWAS file path provided. ")
    if (!is.null(columnFilePath)) {
      rssResult <- loadRssData(
        sumstatPath = gwasFile,
        columnFilePath = columnFilePath,
        region = queryRegion,
        commentString = commentString
      )
      gwasDataSumstats <- rssResult$sumstats
    } else {
      gwasDataSumstats <- as.data.frame(tabixRegion(gwasFile, queryRegion))
      if (nrow(gwasDataSumstats) > 0) {
        gwasDataSumstats <- standardiseSumstatsColumns(gwasDataSumstats)
      }
    }
    if (nrow(gwasDataSumstats) == 0) {
        if (length(names(gwasFile))==0) names(gwasFile) <- gwasFile
        warning(paste0("No GWAS summary statistics found for the region of ", queryRegion, " in ", names(gwasFile), ". "))
        return(NULL)
    }
    # Check if sumstats has z-scores or (beta and se)
    if (!is.null(gwasDataSumstats$z)) {
      # z-scores already present, nothing to do
    } else if (!is.null(gwasDataSumstats$beta) && !is.null(gwasDataSumstats$se)) {
      gwasDataSumstats$z <- gwasDataSumstats$beta / gwasDataSumstats$se
    } else {
      stop("gwasDataSumstats should have 'z' or ('beta' and 'se') columns")
    }
    # check for overlapping variants
    if (!any(gwasDataSumstats$pos %in% gsub("\\:.*$", "", sub("^.*?\\:", "", ldVariants)))) return(NULL)
    gwasAlleleFlip <- matchRefPanel(gwasDataSumstats, ldVariants, colToFlip=colToFlip, matchMinProp = matchMinProp)
    gwasDataSumstats <- getHarmonizedData(gwasAlleleFlip) # post-qc gwas data that is flipped and corrected - gwas study level
    gwasDataSumstats <- gwasDataSumstats[!is.na(gwasDataSumstats$z) & !is.infinite(gwasDataSumstats$z), ]
    return(gwasDataSumstats)
}

#' Function to perform TWAS analysis for across multiple contexts.
#' This function peforms TWAS analysis for multiple contexts for imputable genes within an LD region and summarize the twas results.
#' @param twasWeightsData List of list of twas weights output from generate_twas_db function.
#' @param regionBlock A string with LD region informaiton of chromosome number, star and end position of LD block conneced with "_".
#' @param imputeMissing Logical. Passed to \code{\link{harmonizeTwas}}. When
#'   \code{TRUE}, RAISS imputes GWAS z-scores for variants present in the LD
#'   sketch but missing from the GWAS summary statistics, so weight variants
#'   with LD neighbors but no observed GWAS hit are not silently dropped at
#'   the weight-vs-GWAS intersection. Default \code{FALSE}.
#' @param imputeOpts Named list of RAISS imputation parameters used when
#'   \code{imputeMissing = TRUE}. Defaults to
#'   \code{list(rcond = 0.01, r2Threshold = 0.6, minimumLd = 5, lamb = 0.01)};
#'   imputed variants with \code{R2 < r2Threshold} are dropped.
#' @return A list of list containing twas result table and formatted TWAS data compatible with ctwas_sumstats() function.
#' \itemize{
#'   \item{twas_table}{ A dataframe of twas results summary is generated for each gene-contexts-method pair of all methods for imputable genes.}
#'   \item{twas_data}{ A list of list containing formatted TWAS data.}
#' }
# Shared shape for twasAnalysis() result rows. Internal.
buildTwasScoreRow <- function(twasRs, weightDb, context, study) {
  if (is.null(twasRs)) return(data.frame())
  # Strip trailing "_<suffix>" (snake_case) or "Weights" (camelCase) from
  # method keys to produce a short method name (e.g. enetWeights -> enet,
  # enet_weights -> enet).
  methodLabels <- sub("(_[^_]+|Weights)$", "", names(twasRs))
  data.frame(
    gwas_study   = study,
    method       = methodLabels,
    twas_z       = findData(twasRs, c(2, "z")),
    twas_pval    = findData(twasRs, c(2, "pval")),
    context      = context,
    molecular_id = weightDb
  )
}

# Internal: for each gene-context-study group, if the selected method produced
# NA/Inf TWAS z-scores, fall back to the next best method by rsq_cv.
applyMethodFallback <- function(df) {
  if (nrow(df) == 0 || !all(c("molecular_id", "context", "gwas_study", "is_selected_method", "twas_z", "rsq_cv", "is_imputable") %in% names(df))) {
    return(df)
  }
  groups <- split(seq_len(nrow(df)), list(df$molecular_id, df$context, df$gwas_study), drop = TRUE)
  for (idxs in groups) {
    selIdx <- idxs[df$is_selected_method[idxs]]
    if (length(selIdx) != 1) next
    zVal <- df$twas_z[selIdx]
    if (!is.na(zVal) && is.finite(zVal)) next
    # Selected method has invalid z — try fallback
    otherIdxs <- setdiff(idxs, selIdx)
    validMask <- !is.na(df$twas_z[otherIdxs]) & is.finite(df$twas_z[otherIdxs])
    if (any(validMask)) {
      candidates <- otherIdxs[validMask]
      best <- candidates[which.max(df$rsq_cv[candidates])]
      df$is_selected_method[selIdx] <- FALSE
      df$is_selected_method[best] <- TRUE
      message(paste0("TWAS method fallback for ", df$molecular_id[selIdx],
                     " / ", df$context[selIdx], " / ", df$gwas_study[selIdx],
                     ": ", df$method[selIdx], " -> ", df$method[best]))
    } else {
      # No method has valid z — mark group as non-imputable
      df$is_imputable[idxs] <- FALSE
    }
  }
  df
}

#' @importFrom stringr str_remove
#' @importFrom purrr list_flatten
#' @export
twasPipeline <- function(twasWeightsData,
                         ldMetaFilePath,
                         gwasMetaFile,
                         regionBlock,
                         ldReferenceSampleSize,
                         rsqCutoff = 0.01,
                         rsqPvalCutoff = 0.05,
                         rsqOption = c("rsq", "adj_rsq"),
                         rsqPvalOption = c("pval", "adj_rsq_pval"),
                         mrPvalCutoff = 0.05,
                         mrCoverageColumn = NULL,
                         mrMethod = "susie",
                         mrCoverage = 0.95,
                         outputTwasData = FALSE,
                         eventFilters=NULL,
                         columnFilePath = NULL,
                         commentString="#",
                         imputeMissing = FALSE,
                         imputeOpts = list(rcond = 0.01, r2Threshold = 0.6,
                                            minimumLd = 5, lamb = 0.01)) {
  # internal function to format TWAS output
  formatTwasData <- function(postQcTwasData, twasTable) {
    weightsList <- map(names(postQcTwasData), function(molecularId) {
      mol <- postQcTwasData[[molecularId]]
      contexts <- names(mol[["weights_qced"]])
      molChrom <- mol[["chrom"]]
      modelSel <- mol[["model_selection"]]

      map(contexts, function(context) {
        dataType <- mol[["data_type"]][[context]]
        if (!is.null(modelSel) && is.list(modelSel) && length(modelSel) > 0) {
          isImputable <- modelSel[[context]]$is_imputable
          modelSelected <- if (isTRUE(isImputable)) modelSel[[context]]$selected_model else NA
        } else {
          modelSelected <- NA
          isImputable <- NA
        }
        if (is.null(modelSelected) || !isTRUE(isImputable)) return(NULL)

        gwasStudies <- names(mol[["weights_qced"]][[context]])
        weightKey <- paste0(molecularId, "|", dataType, "_", context)
        studyEntries <- map(gwasStudies, function(study) {
          ctxWeights <- mol[["weights_qced"]][[context]][[study]]
          scaledWgt <- ctxWeights[["scaled_weights"]][, paste0(modelSelected, "_weights"), drop = FALSE]
          colnames(scaledWgt) <- "weight"
          contextVariants <- rownames(ctxWeights[["scaled_weights"]])
          contextRange <- parseVariantId(contextVariants)$pos
          entry <- list(list(
            chrom = molChrom, p0 = min(contextRange), p1 = max(contextRange),
            wgt = scaledWgt, molecular_id = molecularId,
            weight_name = paste0(dataType, "_", context), type = dataType,
            context = context, n_wgt = length(contextVariants)
          ))
          names(entry) <- study
          result <- list(entry)
          names(result) <- weightKey
          result
        }) %>% list_flatten()
        studyEntries
      }) %>% compact() %>% list_flatten()
    }) %>% list_flatten()
    weights <- compact(weightsList)
    # Optional susie_weights_intermediate_qced processing
    if ("susie_weights_intermediate_qced" %in% names(postQcTwasData[[1]])) {
      susieWeightsIntermediateQced <- setNames(lapply(
        names(postQcTwasData),
        function(x) postQcTwasData[[x]]$susie_weights_intermediate_qced
      ), names(postQcTwasData))
    } else {
      susieWeightsIntermediateQced <- NULL
    }

    # gene_z table
    if ("is_selected_method" %in% colnames(twasTable)) {
      twasTable <- twasTable[na.omit(twasTable$is_selected_method), , drop = FALSE]
    }
    if (nrow(twasTable) > 0) {
      twasTable$id <- paste0(twasTable$molecular_id, "|", twasTable$type, "_", twasTable$context)
      twasTable$group <- paste0(twasTable$context, "|", twasTable$type)

      twasTable$z <- twasTable$twas_z

      outputColumns <- c("id", "z", "type", "context", "group", "gwas_study")
      twasTable <- twasTable[, intersect(outputColumns, colnames(twasTable)), drop = FALSE]
      studies <- unique(twasTable$gwas_study)
      zGeneList <- list()
      zSnp <- list()
      for (study in studies) {
        zGeneList[[study]] <- twasTable[twasTable$gwas_study == study, , drop = FALSE]
      }
      result <- list(weights = weights, z_gene = zGeneList)
      if (!is.null(susieWeightsIntermediateQced)) {
        result$susie_weights_intermediate_qced <- susieWeightsIntermediateQced
      }
      return(result)
    } else {
      return(NULL)
    }
  }
  pickBestModel <- function(tw, molecularId, rsqCutoff, rsqPvalCutoff, rsqOption, rsqPvalOption) {
    bestRsq <- rsqCutoff
    cvPerf <- getCvPerformance(tw)
    methodNames <- getMethodNames(tw)
    # SS-TWAS path: no CV performance, all methods are valid
    if (is.null(cvPerf) || length(cvPerf) == 0) {
      modelSelection <- lapply(methodNames, function(context) {
        list(selected_model = NA, is_imputable = TRUE, all_methods = TRUE)
      })
      names(modelSelection) <- methodNames
      return(modelSelection)
    }
    # Determine if a gene/region is imputable and select the best model
    modelSelection <- lapply(methodNames, function(context) {
      selectedModel <- NULL
      availableModels <- do.call(c, lapply(names(cvPerf[[context]]), function(model) {
        if (!is.na(cvPerf[[context]][[model]][, rsqOption])) {
          return(model)
        }
      }))
      if (length(availableModels) <= 0) {
        message(paste0("No model provided TWAS cross validation performance metrics information at context ", context, ". "))
        return(NULL)
      }
      for (model in availableModels) {
        modelData <- cvPerf[[context]][[model]]
        if (modelData[, rsqOption] >= bestRsq & modelData[, colnames(modelData)[which(colnames(modelData) %in% rsqPvalOption)]] < rsqPvalCutoff) {
          bestRsq <- modelData[, rsqOption]
          selectedModel <- model
        }
      }
      if (is.null(selectedModel)) {
        message(paste0(
          "No model has p-value < ", rsqPvalCutoff, " and r2 >= ", rsqCutoff, ", skipping context ", context,
          " at region ", molecularId, ". "
        ))
        return(list(selected_model = c("context_non_imputable"), is_imputable = FALSE)) # No significant model found
      } else {
        selectedModel <- unlist(strsplit(selectedModel, "_performance"))
        message(paste0("The selected best performing model for context ", context, " at region ", molecularId, " is ", selectedModel, ". "))
        return(list(selected_model = selectedModel, is_imputable = TRUE))
      }
    })
    names(modelSelection) <- methodNames
    return(modelSelection)
  }

  # Step 1: TWAS and MR analysis for all methods for imputable gene
  rsqOption <- match.arg(rsqOption)

  # Normalize twasWeightsData entries to TwasWeights S4
  for (wdb in names(twasWeightsData)) {
    entry <- twasWeightsData[[wdb]]
    if (is(entry, "TwasWeights")) next
    if (is.list(entry) && is(entry[["twas_weights"]], "TwasWeights")) {
      # Wrapper list with $twas_weights — unwrap but merge metadata into S4
      twInner <- entry[["twas_weights"]]
      twasWeightsData[[wdb]] <- TwasWeights(
        weights = getWeights(twInner),
        variantIds = getVariantIds(twInner),
        fits = getFits(twInner),
        cvPerformance = getCvPerformance(twInner),
        standardized = getStandardized(twInner),
        molecularId = if (!is.null(entry[["molecular_id"]])) entry[["molecular_id"]] else getMolecularId(twInner),
        dataType = if (!is.null(entry[["data_type"]])) entry[["data_type"]] else getDataType(twInner)
      )
    } else if (is.list(entry) && !is.null(entry[["weights"]])) {
      # Legacy list from load_twas_weights or test fixtures
      wts <- entry[["weights"]]
      vid <- if (!is.null(names(wts)) && length(wts) > 0 && !is.null(rownames(wts[[1]]))) {
        Reduce(union, lapply(wts, rownames))
      } else character(0)
      twasWeightsData[[wdb]] <- TwasWeights(
        weights = wts,
        variantIds = vid,
        fits = entry[["susie_results"]],
        cvPerformance = entry[["twas_cv_performance"]],
        molecularId = if (!is.null(entry[["molecular_id"]])) entry[["molecular_id"]] else character(0),
        dataType = entry[["data_type"]]
      )
    }
  }

  # filter events
  if (!is.null(eventFilters)) {
    for (weightDb in names(twasWeightsData)) {
      tw <- twasWeightsData[[weightDb]]
      contexts <- getMethodNames(tw)
      filteredEvents <- filterMolecularEvents(contexts, eventFilters, removeAllGroup = TRUE)
      if (length(filteredEvents) != 0) {
        # Rebuild TwasWeights with only the filtered contexts
        twasWeightsData[[weightDb]] <- TwasWeights(
          weights = getWeights(tw)[filteredEvents],
          variantIds = getVariantIds(tw),
          fits = if (!is.null(getFits(tw))) getFits(tw)[intersect(filteredEvents, names(getFits(tw)))] else NULL,
          cvPerformance = if (!is.null(getCvPerformance(tw))) getCvPerformance(tw)[intersect(filteredEvents, names(getCvPerformance(tw)))] else NULL,
          standardized = getStandardized(tw),
          molecularId = getMolecularId(tw),
          dataType = getDataType(tw)
        )
      } else {
        twasWeightsData[[weightDb]] <- NULL
      }
    }
  }
  if (length(twasWeightsData)==0) {
    return(list(NULL))
  }

  # harmonize twas weights and gwas sumstats against LD
  twasDataQcedResult <- harmonizeTwas(twasWeightsData, ldMetaFilePath, gwasMetaFile,
                                      ldReferenceSampleSize = ldReferenceSampleSize,
                                      columnFilePath = columnFilePath, commentString = commentString,
                                      imputeMissing = imputeMissing,
                                      imputeOpts = imputeOpts)
  twasResultsDb <- lapply(names(twasWeightsData), function(weightDb) {
    tw <- twasWeightsData[[weightDb]]
    twMethods <- getMethodNames(tw)
    twCv <- getCvPerformance(tw)
    twFits <- getFits(tw)
    twasDataQced <- twasDataQcedResult$twas_data_qced
    if (length(twasDataQced[[weightDb]]) == 0 | is.null(twasDataQced[[weightDb]])) {
      warning(paste0("No data harmonized for ", weightDb, ". Returning NULL for TWAS result for this region."))
      return(NULL)
    }
    if (rsqCutoff > 0) {
      message("Selecting the best model based on criteria...")
      bestModelSelection <- pickBestModel(
        tw, molecularId = weightDb,
        rsqCutoff = rsqCutoff,
        rsqPvalCutoff = rsqPvalCutoff,
        rsqOption = rsqOption,
        rsqPvalOption = rsqPvalOption
      )
      twasDataQced[[weightDb]][["model_selection"]] <- setNames(bestModelSelection, twMethods)
    } else {
      message("Skipping best model selection. Assigning NA of model_selection to all weights.")
      twasDataQced[[weightDb]][["model_selection"]] <- setNames(
        rep(NA, length(twMethods)), twMethods
      )
    }
    dt <- getDataType(tw)
    if (is.null(dt)) {
      twasDataQced[[weightDb]][["data_type"]] <- setNames(
        rep(list(NA), length(twMethods)), twMethods
      )
    }
    if (length(weightDb) < 1) stop(paste0("No data harmonized for ", weightDb, ". "))
    contexts <- names(twasDataQced[[weightDb]][["weights_qced"]])
    gwasStudies <- names(twasDataQced[[weightDb]][["gwas_qced"]])

    # Combined loop for TWAS and MR analysis
    mrCols <- c("gene_name", "num_CS", "num_IV", "cpip", "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval", "I2")

    # Nested lapply for contexts and gwas studies
    twasGeneResults <- lapply(contexts, function(context) {
      studyResults <- lapply(gwasStudies, function(study) {
        twasVariants <- Reduce(intersect, list(rownames(twasDataQced[[weightDb]][["weights_qced"]][[context]][[study]][["weights"]]),
          twasDataQced[[weightDb]][["variant_names"]][[context]][[study]],
          twasDataQced[[weightDb]][["gwas_qced"]][[study]]$variant_id)
        )
        if (length(twasVariants) == 0) {
          return(list(twas_rs_df = data.frame(), mr_rs_df = data.frame()))
        }
        # twas analysis -- enable omnibus when no CV performance available
        hasCv <- !is.null(twCv) && length(twCv) > 0
        twasRs <- twasAnalysis(
          twasDataQced[[weightDb]][["weights_qced"]][[context]][[study]][["weights"]],
          twasDataQced[[weightDb]][["gwas_qced"]][[study]],
          extractVariantsObjs = twasVariants,
          V = twasDataQced[[weightDb]][["svd_V"]],
          D = twasDataQced[[weightDb]][["svd_D"]],
          nSketch = twasDataQced[[weightDb]][["n_sketch"]],
          ldVariantIds = twasDataQced[[weightDb]][["ld_variant_ids"]],
          combineIfNoCv = !hasCv
        )
        if (is.null(twasRs)) {
          return(list(twas_rs_df = data.frame(), mr_rs_df = data.frame()))
        }
        twasRsDf <- buildTwasScoreRow(twasRs, weightDb, context, study)
        # MR analysis
        if (!is.null(twFits) &&
          any(na.omit(twasRsDf$twas_pval) < mrPvalCutoff) &&
          !is.null(twFits[[context]]) && "top_loci" %in% names(twFits[[context]])) {
          if (!"effect_allele_frequency" %in% colnames(twasDataQced[[weightDb]][["gwas_qced"]][[study]])) {
            warning(paste0("skip MR for ", weightDb, " for ", study, ", the effect_allele_frequency information is not available."))
            return(list(twas_rs_df = twasRsDf, mr_rs_df = data.frame()))
          }
          combinedLdMetaDf <- twasDataQcedResult$ref_panel
          # mrFormat expects a nested list with $molecular_id and $susie_results
          mrInput <- list(molecular_id = weightDb, susie_results = twFits)
          mrFormattedInput <- mrFormat(mrInput, context, twasDataQced[[weightDb]][["gwas_qced"]][[study]],
            coverage = mrCoverageColumn, runAlleleQc = TRUE, method = mrMethod,
            coverageLevel = mrCoverage, molecularNameObj = c("molecular_id"),
            ldMetaDf = combinedLdMetaDf
          )
          if (all(is.na(mrFormattedInput$bhat_y))) {
            # FIXME: after updating gwas beta and se NA problem, mr analysis will be restored
            mrRsDf <- as.data.frame(matrix(rep(NA, length(mrCols)), nrow = 1))
            colnames(mrRsDf) <- mrCols
          } else {
            mrRsDf <- as.data.frame(mrAnalysis(mrFormattedInput, cpipCutoff = 0.1))
          }
        } else {
          mrRsDf <- as.data.frame(matrix(rep(NA, length(mrCols)), nrow = 1))
          colnames(mrRsDf) <- mrCols
        }
        mrRsDf$context <- context
        mrRsDf$gwas_study <- study
        mrRsDf$gene_name <- weightDb
        return(list(twas_rs_df = twasRsDf, mr_rs_df = mrRsDf))
      })
      twasContextTable <- do.call(rbind, lapply(studyResults, function(x) x$twas_rs_df))
      mrContextTable <- do.call(rbind, lapply(studyResults, function(x) x$mr_rs_df))
      return(list(twas_context_table = twasContextTable, mr_context_table = mrContextTable))
    })
    twasDataQced[[weightDb]][["svd_V"]] <- NULL
    twasDataQced[[weightDb]][["svd_D"]] <- NULL
    twasDataQced[[weightDb]][["n_sketch"]] <- NULL
    twasDataQced[[weightDb]][["ld_variant_ids"]] <- NULL
    twasWeightsData[[weightDb]] <- NULL
    twasGeneTable <- do.call(rbind, lapply(twasGeneResults, function(x) x$twas_context_table))
    mrGeneTable <- do.call(rbind, lapply(twasGeneResults, function(x) x$mr_context_table))
    return(list(twas_table = twasGeneTable, twas_data_qced = twasDataQced[weightDb], mr_result = mrGeneTable))
  })
  rm(twasDataQcedResult)
  gc()
  twasResultsDb <- twasResultsDb[!sapply(twasResultsDb, function(x) is.null(x) || (is.list(x) && all(sapply(x, is.null))))]
  if (length(twasResultsDb) == 0) {
    return(list(NULL))
  }
  twasResultsTable <- do.call(rbind, lapply(twasResultsDb, function(x) x$twas_table))
  mrResults <- do.call(rbind, lapply(twasResultsDb, function(x) x$mr_result))
  twasData <- do.call(c, lapply(twasResultsDb, function(x) x$twas_data_qced))
  # snp_info <- do.call(c, lapply(twasResultsDb, function(x) x$snp_info))
  rm(twasResultsDb)
  gc()

  # Step 2: Summarize and merge twas cv results and region information for all methods for all contexts for imputable genes.
  twasTable <- do.call(rbind, lapply(names(twasData), function(molecularId) {
    twMol <- twasWeightsData[[molecularId]]
    contexts <- getMethodNames(twMol)
    twMolCv <- getCvPerformance(twMol)
    twMolDt <- getDataType(twMol)
    # merge twas_cv information for same gene across all weight db files, loop through each context for all methods
    geneTable <- do.call(rbind, lapply(contexts, function(context) {
      cvPerf <- if (!is.null(twMolCv)) twMolCv[[context]] else NULL
      modelSel <- twasData[[molecularId]][["model_selection"]][[context]]
      isImputable <- if (!is.null(modelSel)) modelSel$is_imputable else TRUE

      if (is.null(cvPerf) || length(cvPerf) == 0) {
        # SS-TWAS path: no CV, derive methods from weight matrix columns
        wtMat <- getWeights(twMol, context)
        methods <- if (is.matrix(wtMat)) colnames(wtMat) else names(wtMat)
        if (is.null(methods)) methods <- "unknown"
        dtVal <- if (!is.null(twMolDt)) twMolDt[[context]] else NA
        contextTable <- data.frame(
          context = context, method = methods,
          is_imputable = isImputable,
          is_selected_method = FALSE,
          rsq_cv = NA_real_, pval_cv = NA_real_,
          type = dtVal
        )
      } else {
        methods <- sub("_[^_]+$", "", names(cvPerf))
        selectedMethod <- if (!is.null(modelSel)) modelSel$selected_model else NA
        if (is.null(selectedMethod)) selectedMethod <- NA
        isSelectedMethod <- ifelse(methods == selectedMethod, TRUE, FALSE)

        cvRsqs <- sapply(cvPerf, function(x) x[, rsqOption])
        cvPvals <- sapply(cvPerf, function(x) x[, colnames(x)[which(colnames(x) %in% rsqPvalOption)]])

        dtVal <- if (!is.null(twMolDt)) twMolDt[[context]] else NA
        contextTable <- data.frame(
          context = context, method = methods,
          is_imputable = isImputable,
          is_selected_method = isSelectedMethod,
          rsq_cv = cvRsqs, pval_cv = cvPvals,
          type = dtVal
        )
      }
      return(contextTable)
    }))
    geneTable$molecular_id <- molecularId
    return(geneTable)
  }))
  twasTable$chr <- as.integer(stripChrPrefix(gsub("\\_.*", "", regionBlock)))
  twasTable$block <- regionBlock

  # Step 3. merge twas result table and twas input into twasData to output
  colnameOrdered <- c("chr", "molecular_id", "context", "gwas_study", "method", "is_imputable", "is_selected_method", "rsq_cv", "pval_cv", "twas_z", "twas_pval", "type", "block")
  if (nrow(twasResultsTable) == 0) {
    return(list(twas_result = NULL, twas_data = NULL, mr_result = NULL))
  }
  twasTable <- merge(twasTable, twasResultsTable, by = c("molecular_id", "context", "method"))
  twasTable <- applyMethodFallback(twasTable)
  twasTable <- twasTable[twasTable$is_imputable, , drop = FALSE]
  if (outputTwasData & nrow(twasTable) > 0) {
    twasDataSubset <- formatTwasData(twasData, twasTable)
    # if (!is.null(twasDataSubset)) twasDataSubset$snp_info <- snp_info
  } else {
    twasDataSubset <- NULL
  }
  return(list(twas_result = twasTable[, colnameOrdered], twas_data = twasDataSubset, mr_result = mrResults))
}

#' Calculate TWAS z-score and p-value
#'
#' This function calculates the TWAS z-score and p-value given the weights, z-scores,
#' and optionally the correlation matrix (R) or the genotype matrix (X).
#'
#' @param weights A numeric vector of weights.
#' @param z A numeric vector of z-scores.
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item z: The TWAS z-score.
#'   \item pval: The corresponding p-value.
#' }
#'
#' @importFrom stats cor pchisq
#'
#' @export
twasZ <- function(weights, z, R = NULL, X = NULL, V = NULL, D = NULL, nSketch = NULL) {
  # Check that weights and z-scores have the same length
  if (length(weights) != length(z)) {
    stop("Weights and z-scores must have the same length.")
  }

  stat <- t(weights) %*% z

  if (!is.null(V) && !is.null(D) && !is.null(nSketch)) {
    # SVD path: denom = wᵀRw = sum(Lambda * (Vᵀw)²) where Lambda = D²/(nSketch-1)
    Lambda <- D^2 / (nSketch - 1)
    Vw <- crossprod(V, weights)
    denom <- sum(Lambda * Vw^2)
  } else {
    if (is.null(R)) R <- computeLd(X)
    denom <- t(weights) %*% R %*% weights
  }

  zscore <- stat / sqrt(denom)
  pval <- pchisq(zscore * zscore, 1, lower.tail = FALSE)

  return(list(z = zscore, pval = pval))
}

#' Multi-condition TWAS joint test
#'
#' Computes per-condition TWAS z-scores from a variants x conditions weight
#' matrix and an LD sketch (eigenvalues / eigenvectors), and combines them
#' into a joint p-value across conditions using any of the p-value
#' combination methods supported elsewhere in the package.
#'
#' Per-condition test statistics use the cross-condition correlation
#' matrix induced by the weights and LD sketch; methods that need the
#' correlation (\code{"fisher"}, \code{"stouffer"}, \code{"invchisq"},
#' \code{"gbj"}, \code{"aspu"}, \code{"gates"}) consume it directly.
#'
#' @param weights A matrix of weights, one column per condition.
#' @param z A numeric vector of GWAS z-scores aligned to the rows of
#'   \code{weights}.
#' @param V SVD right-singular vectors (variants x components) of the
#'   LD sketch.
#' @param dSvd SVD singular values (vector) of the LD sketch.
#' @param nSketch Sample size of the LD sketch.
#' @param combineMethod Cross-condition p-value combination method. One
#'   of \code{"acat"} (default), \code{"hmp"}, \code{"fisher"},
#'   \code{"stouffer"}, \code{"invchisq"}, \code{"gbj"}, \code{"aspu"},
#'   or \code{"gates"}.
#' @param R,X Legacy alternatives to the LD sketch SVD path; supplying
#'   either still works but is no longer recommended. Documented
#'   workflows use \code{V}, \code{dSvd}, \code{nSketch}.
#'
#' @return A list with:
#' \describe{
#'   \item{Z}{Per-condition Z-score and p-value matrix
#'     (one row per condition).}
#'   \item{combined}{List with \code{method} (the requested
#'     \code{combineMethod}) and \code{pval} (the joint p-value).}
#' }
#'
#' @importFrom stats pnorm
#' @export
twasJointZ <- function(weights, z, R = NULL, X = NULL,
                       V = NULL, dSvd = NULL, nSketch = NULL,
                       combineMethod = c("acat", "hmp", "fisher",
                                          "stouffer", "invchisq",
                                          "gbj", "aspu", "gates")) {
  combineMethod <- match.arg(combineMethod)
  if (nrow(weights) != length(z)) {
    stop("Number of rows in weights must match the length of z-scores.")
  }

  useSvd <- !is.null(V) && !is.null(dSvd) && !is.null(nSketch)

  if (useSvd) {
    # Eigendecomposition path: R = V diag(Lambda) V' with
    # Lambda_i = dSvd_i^2 / (nSketch - 1). Avoid ever forming R.
    Lambda <- dSvd^2 / (nSketch - 1)
    idx <- which(rownames(V) %in% rownames(weights))
    vSub <- V[idx, , drop = FALSE]
    VtW <- crossprod(vSub, weights)  # r x k
    covY <- crossprod(VtW * sqrt(Lambda))  # k x k
  } else {
    # Legacy R / X path (kept for backwards compatibility).
    if (is.null(R)) R <- computeLd(X)
    idx <- which(rownames(R) %in% rownames(weights))
    rSub <- R[idx, idx]
    covY <- crossprod(weights, rSub) %*% weights
  }

  ySd <- sqrt(diag(covY))
  xSd <- rep(1, nrow(weights))  # standardized genotype scale

  # Gamma scaling per condition: gamma_k = diag(xSd / ySd[k])
  g <- setNames(lapply(colnames(weights), function(cond) {
    diag(xSd / ySd[cond], length(xSd), length(xSd))
  }), colnames(weights))

  # Per-condition Z-score and two-sided p-value
  zMatrix <- do.call(rbind, lapply(colnames(weights), function(cond) {
    Zi <- crossprod(weights[, cond], g[[cond]]) %*% as.numeric(z)
    pval <- 2 * pnorm(abs(Zi), lower.tail = FALSE)
    setNames(c(Zi, pval), c("Z", "pval"))
  }))
  rownames(zMatrix) <- colnames(weights)

  # Cross-condition correlation sig[i,j] from weighted LD sketch.
  lam <- matrix(NA_real_, nrow = ncol(weights), ncol = nrow(weights),
                dimnames = list(colnames(weights), NULL))
  for (cond in colnames(weights)) {
    lam[cond, ] <- as.numeric(weights[, cond] %*% g[[cond]])
  }
  if (useSvd) {
    LV <- lam %*% vSub                               # k x r
    sig <- tcrossprod(sweep(LV, 2, Lambda, "*"), LV) # k x k
  } else {
    sig <- tcrossprod((lam %*% rSub), lam)
  }

  # Dispatch to the requested combination method. Methods reuse the same
  # helpers as twasAnalysis's cross-method omnibus.
  zscores <- as.numeric(zMatrix[, "Z"])
  pvals   <- as.numeric(zMatrix[, "pval"])
  valid <- is.finite(pvals) & pvals > 0 & pvals < 1
  combinedPval <- if (sum(valid) < 2L) {
    NA_real_
  } else {
    sigSub <- sig[valid, valid, drop = FALSE]
    tryCatch(
      switch(combineMethod,
        acat     = pvalAcat(pvals[valid]),
        hmp      = pvalHmp(pvals[valid]),
        fisher   = ,
        stouffer = ,
        invchisq = pvalPoolr(pvals[valid], method = combineMethod, R = sigSub),
        gbj      = pvalGbj(zscores[valid], R = sigSub, method = combineMethod),
        aspu     = ,
        gates    = pvalAspu(zscores[valid], pvals[valid],
                              R = sigSub, method = combineMethod)
      ),
      error = function(e) {
        warning(sprintf("twasJointZ combineMethod = '%s' failed: %s",
                        combineMethod, e$message))
        NA_real_
      }
    )
  }

  list(Z = zMatrix,
       combined = list(method = combineMethod, pval = combinedPval))
}

#' TWAS Analysis
#'
#' Performs TWAS analysis using the provided weights matrix, GWAS summary statistics database,
#' and LD matrix. It extracts the necessary GWAS summary statistics and LD matrix based on the
#' specified variants and computes the z-score and p-value for each gene.
#'
#' When \code{combineIfNoCv = TRUE} and there are at least two methods with
#' valid p-values, an omnibus p-value is computed via the method specified in
#' \code{combineMethod} and appended as an \code{"omnibus"} entry. This is
#' intended for summary-statistics TWAS where cross-validation performance is
#' not available for model selection.
#'
#' @param weightsMatrix A matrix containing weights for all methods.
#' @param gwasSumstatsDb A data frame containing the GWAS summary statistics.
#' @param ldMatrix A matrix representing linkage disequilibrium between variants.
#' @param extractVariantsObjs A vector of variant identifiers to extract from the GWAS and LD matrix.
#' @param V SVD right-singular vectors from LD sketch (optional).
#' @param D SVD singular values from LD sketch (optional).
#' @param nSketch Sample size of LD sketch (optional).
#' @param ldVariantIds Variant IDs in the LD sketch (optional).
#' @param combineMethod P-value combination method: \code{"acat"} (default),
#'   \code{"hmp"}, \code{"fisher"}, \code{"stouffer"}, \code{"invchisq"},
#'   \code{"gbj"}, \code{"aspu"}, or \code{"gates"}.
#' @param combineIfNoCv Logical. If TRUE and no CV performance is available,
#'   combine per-method p-values into an omnibus result.
#'
#' @return A list with TWAS z-scores and p-values across methods for each gene.
#'   When omnibus combination is enabled, includes an additional \code{"omnibus"}
#'   entry.
#' @export
twasAnalysis <- function(weightsMatrix, gwasSumstatsDb, ldMatrix = NULL,
                         extractVariantsObjs, V = NULL, D = NULL,
                         nSketch = NULL, ldVariantIds = NULL,
                         combineMethod = "acat",
                         combineIfNoCv = FALSE) {
  # Extract gwas_sumstats
  gwasSumstatsSubset <- gwasSumstatsDb[match(extractVariantsObjs, gwasSumstatsDb$variant_id), ]
  # Validate that the GWAS subset is not empty
  if (nrow(gwasSumstatsSubset) == 0 | all(is.na(gwasSumstatsSubset))) {
    warning("No GWAS summary statistics found for the specified variants.")
    return(NULL)
  }

  # SVD path
  if (!is.null(V) && !is.null(D) && !is.null(nSketch) && !is.null(ldVariantIds)) {
    validIndices <- extractVariantsObjs %in% ldVariantIds
    if (!any(validIndices)) {
      warning("None of the specified variants are present in the LD sketch. Skipping this context.")
      return(NULL)
    }
    validVariantsObjs <- extractVariantsObjs[validIndices]
    # Subset V rows to match the valid variants
    vRowIdx <- match(validVariantsObjs, ldVariantIds)
    vSubset <- V[vRowIdx, , drop = FALSE]
    weightsMatrix <- weightsMatrix[validVariantsObjs, , drop = FALSE]
    gwasSumstatsSubset <- gwasSumstatsDb[match(validVariantsObjs, gwasSumstatsDb$variant_id), ]
    twasZPval <- apply(
      as.matrix(weightsMatrix), 2,
      function(x) twasZ(x, gwasSumstatsSubset$z, V = vSubset, D = D, nSketch = nSketch)
    )
    return(.maybeAddOmnibus(twasZPval, weightsMatrix, ldMatrix,
                            combineMethod, combineIfNoCv))
  }

  # LD matrix path
  validIndices <- extractVariantsObjs %in% rownames(ldMatrix)
  if (!any(validIndices)) {
    warning("None of the specified variants are present in the LD matrix. Skipping this context.")
    return(NULL)
  }
  validVariantsObjs <- extractVariantsObjs[validIndices]
  ldMatrixSubset <- ldMatrix[validVariantsObjs, validVariantsObjs]
  weightsMatrix <- weightsMatrix[validVariantsObjs, , drop = FALSE]
  gwasSumstatsSubset <- gwasSumstatsDb[match(validVariantsObjs, gwasSumstatsDb$variant_id), ]
  twasZPval <- apply(
    as.matrix(weightsMatrix), 2,
    function(x) twasZ(x, gwasSumstatsSubset$z, R = ldMatrixSubset)
  )
  return(.maybeAddOmnibus(twasZPval, weightsMatrix, ldMatrixSubset,
                          combineMethod, combineIfNoCv))
}

#' Add omnibus p-value combination to TWAS results
#' @noRd
.maybeAddOmnibus <- function(twasZPval, weightsMatrix, ldMatrix,
                             combineMethod, combineIfNoCv) {
  if (!isTRUE(combineIfNoCv) || length(twasZPval) < 2) {
    return(twasZPval)
  }

  pvals <- vapply(twasZPval, function(x) as.numeric(x$pval), numeric(1))
  zscores <- vapply(twasZPval, function(x) as.numeric(x$z), numeric(1))
  valid <- !is.na(pvals) & is.finite(pvals) & pvals > 0 & pvals < 1

  if (sum(valid) < 2) return(twasZPval)

  combinedPval <- tryCatch({
    switch(combineMethod,
      acat = pvalAcat(pvals[valid]),
      hmp = pvalHmp(pvals[valid]),
      fisher = , stouffer = , invchisq = {
        methodCor <- twasMethodCor(
          lapply(which(valid), function(i) weightsMatrix[, i]),
          ldMatrix)
        pvalPoolr(pvals[valid], method = combineMethod, R = methodCor)
      },
      gbj = {
        methodCor <- twasMethodCor(
          lapply(which(valid), function(i) weightsMatrix[, i]),
          ldMatrix)
        pvalGbj(zscores[valid], R = methodCor, method = combineMethod)
      },
      aspu = , gates = {
        methodCor <- twasMethodCor(
          lapply(which(valid), function(i) weightsMatrix[, i]),
          ldMatrix)
        pvalAspu(zscores[valid], pvals[valid], R = methodCor, method = combineMethod)
      },
      pvalAcat(pvals[valid])  # fallback
    )
  }, error = function(e) {
    warning(sprintf("Omnibus combination (%s) failed: %s", combineMethod, e$message))
    NA_real_
  })

  twasZPval[["omnibus"]] <- list(z = NA_real_, pval = combinedPval)
  twasZPval
}
