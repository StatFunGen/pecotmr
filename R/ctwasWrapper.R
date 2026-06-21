#' Function to select variants for ctwas weights input
#' @param regionData A list of list containing weights list and snp_info list data for multiple genes/events within a single LD block region.
#' @param exportTwasWeightDb A list of list of fine-mapping result data formatted by generate_twas_db function.
#' @param regionBlock A string for region information for region_weights, consisted of chromosome number, star and end position of LD block conneced with "_".
#' @export
trimCtwasVariants <- function(regionData, twasWeightCutoff = 1e-5, csMinCor = 0.8,
                              minPipCutoff = 0.1, maxNumVariants = 1000) {
  # internal functions to select variants for a gene-context pair weight list
  selectVariants <- function(groupName, studyName, regionData, csMinCor, minPipCutoff, maxNumVariants) {
    weightList <- regionData$weights[[groupName]][[studyName]]
    context <- weightList$context
    selectedVariantsByContext <- c()
    molecularId <- gsub("\\|.*", "", groupName)

    if ("csVariants" %in% names(regionData$susieWeightsIntermediate[[molecularId]][[context]]) & length(regionData$susieWeightsIntermediate[[molecularId]][[context]][["csVariants"]]) != 0) {
      csMinAbsCor <- regionData$susieWeightsIntermediate[[molecularId]][[context]]$csPurity$minAbsCorr
      for (L in seq_along(regionData$susieWeightsIntermediate[[molecularId]][[context]]$csVariants)) {
        # we includ all variants in $cs_variant if min_abs_corr > csMinCor for the set
        if (csMinAbsCor[L] >= csMinCor) {
          csVariants <- regionData$susieWeightsIntermediate[[molecularId]][[context]]$csVariants[[L]]
          selectedVariantsByContext <- csVariants[csVariants %in% rownames(weightList$wgt)]
        }
      }
    }
    contextPip <- regionData$susieWeightsIntermediate[[molecularId]][[context]]$pip
    # variant IDs are in canonical chr-prefix format from allele_qc pipeline
    highPipVariants <- names(contextPip[contextPip > minPipCutoff])[names(contextPip[contextPip > minPipCutoff]) %in% rownames(weightList$wgt)]
    selectedVariantsByContext <- unique(c(selectedVariantsByContext, highPipVariants))

    # prioritize SNPs based on PIP if maxNumVariants different from Inf
    availableVariants <- intersect(rownames(weightList$wgt), names(contextPip))
    prioritized <- unique(c(selectedVariantsByContext, setdiff(availableVariants, selectedVariantsByContext)))
    prioritized <- prioritized[order(-contextPip[prioritized])]
    selectedVariantsByContext <- head(prioritized, maxNumVariants)
    weightList$wgt <- weightList$wgt[selectedVariantsByContext, , drop = FALSE]
    return(weightList)
  }
  mergeByStudy <- function(weights) {
    weightList <- list()
    for (group in names(weights)) {
      for (study in names(weights[[group]])) {
        weightList[[study]][[group]] <- weights[[group]][[study]]
      }
    }
    return(weightList)
  }

  weights <- setNames(lapply(names(regionData$weights), function(group) {
    for (study in names(regionData$weights[[group]])) {
      regionData$weights[[group]][[study]]$wgt <- regionData$weights[[group]][[study]]$wgt[abs(regionData$weights[[group]][[study]]$wgt[, 1]) >= twasWeightCutoff, , drop = FALSE]
      if (nrow(regionData$weights[[group]][[study]]$wgt) < 1) {
        regionData$weights[[group]][[study]] <- NULL
        next
      }
      if (all(is.na(regionData$weights[[group]][[study]]$wgt[, 1])) || all(is.nan(regionData$weights[[group]][[study]]$wgt[, 1]))) {
        regionData$weights[[group]][[study]] <- NULL
        next
      }
      if (nrow(regionData$weights[[group]][[study]]$wgt) < maxNumVariants) {
        regionData$weights[[group]][[study]]$nWgt <- nrow(regionData$weights[[group]][[study]]$wgt)
      } else {
        regionData$weights[[group]][[study]] <- selectVariants(group, study, regionData, csMinCor = csMinCor, minPipCutoff = minPipCutoff, maxNumVariants = maxNumVariants)
        regionData$weights[[group]][[study]]$nWgt <- nrow(regionData$weights[[group]][[study]]$wgt)
      }
      regionData$weights[[group]] <- Filter(Negate(is.null), regionData$weights[[group]])
      contextRange <- as.integer(sapply(rownames(regionData$weights[[group]][[study]]$wgt), function(variant) strsplit(variant, "\\:")[[1]][2]))
      if(twasWeightCutoff!=0 | csMinCor!=0 | minPipCutoff!=0 | maxNumVariants!=Inf){
        regionData$weights[[group]][[study]][["p0"]] = min(contextRange)# update min max position
        regionData$weights[[group]][[study]][["p1"]] = max(contextRange)
      }
    }
    return(regionData$weights[[group]])
  }), names(regionData$weights))
  weights <- Filter(Negate(is.null), weights)
  weights <- mergeByStudy(weights)
  return(weights)
}

