#' Load a PLINK .bim file for cTWAS
#'
#' @description
#' \strong{Deprecated.} Use [readBim()] via the standard I/O path
#' instead. This wrapper remains for backwards compatibility and calls
#' [readBim()] internally, mapping its output to the legacy column names.
#'
#' @param bimFilePath Path to a PLINK \code{.bim} file (or a \code{.bed}
#'   file - the \code{.bim} extension is resolved automatically).
#'
#' @return A data.frame with columns \code{chrom}, \code{id}, \code{GD},
#'   \code{pos}, \code{A1}, \code{A2}. Variant IDs are normalised via
#'   [normalizeVariantId()].
#'
#' @export
ctwasBimfileLoader <- function(bimFilePath) {
  .Deprecated("readBim", package = "pecotmr",
              msg = "ctwasBimfileLoader() is deprecated. Use readBim() instead.")
  # readBim() expects a .bed path and derives .bim from it.
  # Accept either .bim or .bed and normalise to .bed.
  bedPath <- sub("\\.bim$", ".bed", bimFilePath)
  bim <- readBim(bedPath)
  # Map new column names back to legacy names
  snpInfo <- data.frame(
    chrom = bim$chrom,
    id    = normalizeVariantId(bim$id),
    GD    = bim$gpos,
    pos   = bim$pos,
    A1    = bim$a1,
    A2    = bim$a0,
    stringsAsFactors = FALSE
  )
  return(snpInfo)
}

#' Load cTWAS LD meta-data
#'
#' @description
#' \strong{Deprecated.} Use [ldLoader()] with its \code{ldInfo}
#' argument instead. This wrapper remains for backwards compatibility and
#' produces the same \code{list(LD_info, region_info)} output as the original.
#'
#' @param ldMetaDataFile Path to the LD meta-data TSV file.
#' @param subsetRegionIds Optional character vector of region IDs
#'   (\code{"chrom_start_end"}) to subset to.
#'
#' @return A list with components:
#' \describe{
#'   \item{LD_info}{Data.frame with columns \code{region_id}, \code{LD_file},
#'     \code{SNP_file}.}
#'   \item{region_info}{Data.frame with columns \code{chrom}, \code{start},
#'     \code{stop}, \code{region_id}.}
#' }
#'
#' @importFrom vroom vroom
#' @export
getCtwasMetaData <- function(ldMetaDataFile, subsetRegionIds = NULL) {
  .Deprecated("ldLoader", package = "pecotmr",
              msg = "getCtwasMetaData() is deprecated. Use ldLoader() with ldInfo instead.")
  LD_info <- as.data.frame(vroom(ldMetaDataFile))
  colnames(LD_info)[1] <- "chrom"
  LD_info$region_id <- paste(as.integer(stripChrPrefix(LD_info$chrom)),
                             LD_info$start, LD_info$end, sep = "_")
  LD_info$LD_file <- paste0(dirname(ldMetaDataFile), "/",
                            gsub(",.*$", "", LD_info$path))
  LD_info$SNP_file <- paste0(LD_info$LD_file, ".bim")
  LD_info <- LD_info[, c("region_id", "LD_file", "SNP_file")]
  region_info <- LD_info[, "region_id", drop = FALSE]
  region_info$chrom <- as.integer(gsub("\\_.*$", "", region_info$region_id))
  region_info$start <- as.integer(gsub("\\_.*$", "",
                                       sub("^.*?\\_", "", region_info$region_id)))
  region_info$stop <- as.integer(sub("^.*?\\_", "",
                                      sub("^.*?\\_", "", region_info$region_id)))
  region_info$region_id <- paste0(region_info$chrom, "_",
                                   region_info$start, "_",
                                   region_info$stop)
  region_info <- region_info[, c("chrom", "start", "stop", "region_id")]
  if (!is.null(subsetRegionIds)) {
    region_info <- region_info[region_info$region_id %in% subsetRegionIds, ]
  }
  return(list(LD_info = LD_info, region_info = region_info))
}

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

    if ("cs_variants" %in% names(regionData$susie_weights_intermediate[[molecularId]][[context]]) & length(regionData$susie_weights_intermediate[[molecularId]][[context]][["cs_variants"]]) != 0) {
      csMinAbsCor <- regionData$susie_weights_intermediate[[molecularId]][[context]]$cs_purity$min.abs.corr
      for (L in seq_along(regionData$susie_weights_intermediate[[molecularId]][[context]]$cs_variants)) {
        # we includ all variants in $cs_variant if min_abs_corr > csMinCor for the set
        if (csMinAbsCor[L] >= csMinCor) {
          csVariants <- regionData$susie_weights_intermediate[[molecularId]][[context]]$cs_variants[[L]]
          selectedVariantsByContext <- csVariants[csVariants %in% rownames(weightList$wgt)]
        }
      }
    }
    contextPip <- regionData$susie_weights_intermediate[[molecularId]][[context]]$pip
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
        regionData$weights[[group]][[study]]$n_wgt <- nrow(regionData$weights[[group]][[study]]$wgt)
      } else {
        regionData$weights[[group]][[study]] <- selectVariants(group, study, regionData, csMinCor = csMinCor, minPipCutoff = minPipCutoff, maxNumVariants = maxNumVariants)
        regionData$weights[[group]][[study]]$n_wgt <- nrow(regionData$weights[[group]][[study]]$wgt)
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

