#' Deduplicate and sort genomic regions by chromosome and start position.
#' @importFrom dplyr distinct arrange
#' @importFrom magrittr %>%
#' @noRd
orderDedupRegions <- function(df) {
  df$chrom <- as.integer(stripChrPrefix(df$chrom))
  df <- distinct(df, chrom, start, .keep_all = TRUE) %>%
    arrange(chrom, start)
  df
}

#' Find the first and last rows of genomicData that overlap a query region.
#' Clamps the query to the available data range before searching.
#' @importFrom dplyr filter arrange slice
#' @noRd
findIntersectionRows <- function(genomicData, regionChrom, regionStart, regionEnd) {
  chromData <- genomicData %>% filter(chrom == regionChrom)
  if (nrow(chromData) == 0) stop("No data for chromosome ", regionChrom)

  # Clamp query to available range
  regionStart <- max(regionStart, min(chromData$start))
  regionEnd   <- min(regionEnd,   max(chromData$end))

  startRow <- genomicData %>%
    filter(chrom == regionChrom, start <= regionStart, end > regionStart) %>%
    slice(1)
  endRow <- genomicData %>%
    filter(chrom == regionChrom, start < regionEnd, end >= regionEnd) %>%
    arrange(desc(end)) %>%
    slice(1)

  if (nrow(startRow) == 0 || nrow(endRow) == 0) {
    stop("Region ", regionChrom, ":", regionStart, "-", regionEnd,
         " is not covered by any rows in the LD metadata.")
  }
  list(start_row = startRow, end_row = endRow)
}

#' Validate that startRow..endRow fully covers [regionStart, regionEnd].
#' @noRd
validateSelectedRegion <- function(startRow, endRow, regionStart, regionEnd) {
  if (startRow$start > regionStart || endRow$end < regionEnd) {
    stop("Region ", regionStart, "-", regionEnd, " is not fully covered by the LD metadata ",
         "(available: ", startRow$start, "-", endRow$end, ").")
  }
}

#' Extract values of a column for rows spanning the intersection range.
#' @noRd
extractFilePaths <- function(genomicData, intersectionRows, columnToExtract) {
  if (!columnToExtract %in% names(genomicData)) {
    stop("Column '", columnToExtract, "' not found in genomic data.")
  }
  idx <- which(genomicData$chrom == intersectionRows$start_row$chrom &
               genomicData$start >= intersectionRows$start_row$start &
               genomicData$start <= intersectionRows$end_row$start)
  genomicData[[columnToExtract]][idx]
}

#' Find LD blocks overlapping a query region from a metadata TSV file.
#'
#' @param ldReferenceMetaFile TSV with columns chrom, start, end, path.
#'   The path column may be comma-separated: "ld_file,bim_file".
#' @param region "chr:start-end" string or data.frame with chrom/start/end.
#' @param completeCoverageRequired If TRUE, error when the region extends
#'   beyond available LD blocks.
#' @return A list with: intersections (LD_file_paths, bim_file_paths),
#'   ld_meta_data, and parsed region.
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom vroom vroom
#' @noRd
getRegionalLdMeta <- function(ldReferenceMetaFile, region, completeCoverageRequired = FALSE) {
  genomicData <- vroom(ldReferenceMetaFile)
  region <- parseRegion(region)
  # Set column names
  names(genomicData) <- c("chrom", "start", "end", "path")
  names(region) <- c("chrom", "start", "end")

  # Treat start=0, end=0 as "covers all regions" (used for whole-chromosome PLINK files)
  wholeChrom <- genomicData$start == 0 & genomicData$end == 0
  if (any(wholeChrom)) genomicData$end[wholeChrom] <- Inf

  # Order and deduplicate regions
  genomicData <- orderDedupRegions(genomicData)
  region <- orderDedupRegions(region)

  # Process file paths
  filePath <- genomicData$path %>%
    str_split(",", simplify = TRUE) %>%
    data.frame() %>%
    `colnames<-`(if (ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))

  genomicData <- cbind(genomicData, filePath) %>% select(-path)

  # Find intersection rows
  intersectionRows <- findIntersectionRows(genomicData, region$chrom, region$start, region$end)

  # Validate region
  if (completeCoverageRequired) {
    validateSelectedRegion(intersectionRows$start_row, intersectionRows$end_row, region$start, region$end)
  }

  # Extract file paths
  ldPaths <- findValidFilePaths(ldReferenceMetaFile, extractFilePaths(genomicData, intersectionRows, "LD_file_path"))
  bimPaths <- if ("bim_file_path" %in% names(genomicData)) {
    findValidFilePaths(ldReferenceMetaFile, extractFilePaths(genomicData, intersectionRows, "bim_file_path"))
  } else {
    NULL
  }

  return(list(
    intersections = list(
      start_index = intersectionRows$start_row,
      end_index = intersectionRows$end_row,
      LD_file_paths = ldPaths,
      bim_file_paths = bimPaths
    ),
    ld_meta_data = genomicData,
    region = region
  ))
}

#' Read a pre-computed LD matrix (.cor.xz) and its bim file, returning a
#' symmetric matrix with variants ordered by position.
#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom stats setNames
#' @noRd
processLdMatrix <- function(ldFilePath, snpFilePath = NULL) {
  # Read .cor.xz matrix
  ldFileCon <- xzfile(ldFilePath)
  ldMatrix <- scan(ldFileCon, quiet = TRUE)
  close(ldFileCon)
  ldMatrix <- matrix(ldMatrix, ncol = sqrt(length(ldMatrix)), byrow = TRUE)

  # Auto-detect variant metadata file: .bim (PLINK1) or .pvar/.pvar.zst (PLINK2)
  if (is.null(snpFilePath)) {
    candidates <- paste0(ldFilePath, c(".bim", ".pvar", ".pvar.zst"))
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) stop("No variant file found for: ", ldFilePath,
                                  " (tried .bim, .pvar, .pvar.zst)")
    snpFilePath <- found[1]
  }

  ldVariants <- readVariantMetadata(snpFilePath)
  isPvar <- !("gpos" %in% names(ldVariants))
  ldVariants <- ldVariants %>%
    mutate(chrom = as.character(as.integer(stripChrPrefix(chrom))),
           variants = normalizeVariantId(id))
  if (isPvar) {
    ldVariants <- rename(ldVariants, GD = pos)
    ldVariants$GD <- ldVariants$pos <- as.integer(
      sapply(ldVariants$variants, function(v) strsplit(v, ":")[[1]][2]))
  } else {
    ldVariants <- rename(ldVariants, GD = gpos)
  }

  # Label and symmetrize the matrix
  colnames(ldMatrix) <- rownames(ldMatrix) <- ldVariants$variants
  if (all(ldMatrix[lower.tri(ldMatrix)] == 0)) {
    ldMatrix[lower.tri(ldMatrix)] <- t(ldMatrix)[lower.tri(ldMatrix)]
  } else {
    ldMatrix[upper.tri(ldMatrix)] <- t(ldMatrix)[upper.tri(ldMatrix)]
  }

  # Order variants by genomic position
  posOrder <- order(sapply(ldVariants$variants, function(v) as.integer(strsplit(v, ":")[[1]][2])))
  ldVariants <- ldVariants[posOrder, ]
  ldMatrix <- ldMatrix[ldVariants$variants, ldVariants$variants]

  list(LD_matrix = ldMatrix, LD_variants = ldVariants)
}

#' Subset an LD matrix and variant info to a genomic region, optionally
#' further restricted to specific coordinates.
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @noRd
extractLdForRegion <- function(ldMatrix, variants, region, extractCoordinates) {
  extracted <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)

  if (!is.null(extractCoordinates)) {
    extractCoordinates <- extractCoordinates %>%
      mutate(chrom = as.integer(stripChrPrefix(chrom))) %>%
      select(chrom, pos)
    extracted <- extracted %>%
      mutate(chrom = as.integer(stripChrPrefix(chrom))) %>%
      merge(extractCoordinates, by = c("chrom", "pos"))
    keepCols <- intersect(c("chrom", "variants", "pos", "GD", "A1", "A2",
                             "variance", "allele_freq", "n_nomiss"), names(extracted))
    extracted <- select(extracted, all_of(keepCols))
  }

  mat <- ldMatrix[extracted$variants, extracted$variants, drop = FALSE]
  list(extracted_LD_matrix = mat, extracted_LD_variants = extracted)
}

#' Combine multiple block-level LD matrices into one, handling boundary overlaps.
#' @importFrom utils tail
#' @noRd
createLdMatrix <- function(ldMatrices, variants) {
  # Merge variant lists, deduplicating boundary overlaps
  mergeVariants <- function(variantList) {
    merged <- character(0)
    for (v in variantList) {
      ids <- if (is.list(v) && !is.null(v$variants)) v$variants else v
      if (length(ids) == 0) next
      if (length(merged) > 0 && tail(merged, 1) == ids[1]) ids <- ids[-1]
      merged <- c(merged, ids)
    }
    merged
  }

  allVariants <- mergeVariants(variants)
  combined <- matrix(0, nrow = length(allVariants), ncol = length(allVariants),
                     dimnames = list(allVariants, allVariants))

  # Place each block into the combined matrix
  for (i in seq_along(ldMatrices)) {
    v <- rownames(ldMatrices[[i]])
    idx <- match(v, allVariants)
    combined[idx, idx] <- ldMatrices[[i]]
  }
  combined
}

#' Load and Process Linkage Disequilibrium (LD) Matrix
#'
#' Unified entry point for loading LD data from a metadata TSV file.
#'
#' The metadata TSV must have columns: chrom, start, end, path. Two formats:
#' \itemize{
#'   \item Pre-computed LD blocks: many rows per chromosome with block boundaries
#'     in start/end and path pointing to .cor.xz files (optionally comma-separated
#'     with a .bim path).
#'   \item PLINK genotype files: one row per chromosome with start=0, end=0, and
#'     path pointing to a per-chromosome PLINK prefix (.pgen/.pvar[.zst]/.psam or
#'     .bed/.bim/.fam). LD is computed on the fly via \code{computeLd()}.
#' }
#'
#' @param ldMetaFilePath Path to the LD metadata TSV file.
#' @param region Region of interest: "chr:start-end" string or data.frame with chrom/start/end.
#' @param extractCoordinates Optional data.frame with columns "chrom" and "pos" for
#'   specific coordinates extraction (only for pre-computed LD blocks).
#' @param returnGenotype Controls what LD_matrix contains in the return value.
#'   FALSE (default): always return correlation matrix R.
#'   TRUE: return genotype matrix X (only valid for PLINK sources).
#'   "auto": return X for PLINK sources, R for pre-computed sources.
#' @param nSample Optional sample size for computing variance (= 2*p*(1-p)*n/(n-1)).
#'   If NULL, ref_panel will not include variance or n_nomiss columns.
#'   Only used for PLINK genotype sources.
#'
#' @return A list with:
#' \describe{
#'   \item{LD_variants}{Character vector of variant IDs (canonical format).}
#'   \item{LD_matrix}{LD correlation matrix R (or genotype matrix X when returnGenotype is TRUE or "auto" with PLINK source).}
#'   \item{ref_panel}{Data.frame with variant metadata (chrom, pos, A2, A1, variant_id,
#'     and optionally allele_freq, variance, n_nomiss).}
#'   \item{is_genotype}{Logical: TRUE if LD_matrix contains genotype X, FALSE if correlation R.}
#'   \item{block_metadata}{Data.frame with region/block info. For pre-computed LD: one row per block.
#'     For PLINK: a single row spanning the loaded region.}
#' }
#' @export
loadLdMatrix <- function(ldMetaFilePath, region, extractCoordinates = NULL,
                           returnGenotype = FALSE, nSample = NULL) {
  source <- resolveLdSource(ldMetaFilePath)
  isGeno <- source$type %in% c("plink2", "plink1", "vcf", "gds")

  # "auto": return X for genotype sources, R for pre-computed
  if (identical(returnGenotype, "auto")) returnGenotype <- isGeno

  if (isGeno) {
    genoPath <- resolveGenotypePathForRegion(source$meta_path, region)
    result <- loadLdFromGenotype(genoPath, region,
                                    returnGenotype = returnGenotype,
                                    nSample = nSample)
  } else {
    # Pre-computed LD blocks (.cor.xz)
    if (returnGenotype) {
      stop("returnGenotype=TRUE requires genotype files, not pre-computed LD matrices.")
    }
    result <- loadLdFromBlocks(source$meta_path, region, extractCoordinates, nSample = nSample)
  }

  # Remove any duplicate variant IDs (safety net for boundary overlaps)
  variantIds <- getVariantIds(result)
  if (!is.null(variantIds)) {
    dupIdx <- which(duplicated(variantIds))
    if (length(dupIdx) > 0) {
      variantIdsClean <- variantIds[-dupIdx]
      corr <- getCorrelation(result)
      if (!is.null(corr)) {
        corr <- corr[-dupIdx, -dupIdx, drop = FALSE]
      }
      variantsGr <- result@variants[-dupIdx]
      result <- LdData(
        correlation = corr,
        genotypeHandle = result@genotypeHandle,
        snpIdx = result@snpIdx,
        variants = variantsGr,
        blockMetadata = result@blockMetadata,
        nRef = result@nRef
      )
    }
  }

  result
}

# ---------- Internal: resolve LD source type ----------

#' @noRd
hasPlink2Files <- function(prefix) {
  file.exists(paste0(prefix, ".pgen")) &&
    (file.exists(paste0(prefix, ".pvar")) || file.exists(paste0(prefix, ".pvar.zst"))) &&
    file.exists(paste0(prefix, ".psam"))
}

#' @noRd
hasPlink1Files <- function(prefix) {
  file.exists(paste0(prefix, ".bed")) &&
    file.exists(paste0(prefix, ".bim")) &&
    file.exists(paste0(prefix, ".fam"))
}

#' @noRd
isVcfPath <- function(path) {
  grepl("\\.(vcf|vcf\\.gz|bcf)$", path) && file.exists(path)
}

#' @noRd
isGdsPath <- function(path) {
  grepl("\\.gds$", path) && file.exists(path)
}

#' Check whether a path points to a genotype source (PLINK, VCF, or GDS).
#' @noRd
isGenotypeSource <- function(path) {
  hasPlink2Files(path) || hasPlink1Files(path) || isVcfPath(path) || isGdsPath(path)
}

#' Resolve an LD source metadata TSV to its actual data type.
#'
#' The metadata TSV has columns: chrom, start, end, path. Three categories are
#' supported:
#' \itemize{
#'   \item Pre-computed LD blocks (.cor.xz): many rows per chromosome, each with
#'     specific start/end block boundaries and path pointing to .cor.xz files.
#'   \item Genotype files (PLINK2, PLINK1, VCF, or GDS): one row per chromosome
#'     with start=0, end=0, and path pointing to a per-chromosome genotype file
#'     or prefix. The actual region filter is applied by the genotype loader.
#' }
#'
#' This function peeks at the first row to determine the data type.
#' The actual per-chromosome path is resolved later by
#' \code{resolveGenotypePathForRegion()} at load time.
#'
#' @param path Path to a metadata TSV file with columns chrom, start, end, path.
#' @return A list with:
#'   \item{type}{"plink2", "plink1", "vcf", "gds", or "precomputed"}
#'   \item{data_path}{Genotype path from first row (for type detection only; actual
#'     per-chromosome path is resolved at load time)}
#'   \item{meta_path}{The metadata TSV path (always set)}
#' @importFrom vroom vroom
#' @noRd
resolveLdSource <- function(path) {
  if (!file.exists(path)) {
    stop("LD metadata file not found: ", path,
         "\n  Expected: a TSV file with columns chrom, start, end, path.")
  }

  # Peek at first row to determine underlying data type
  meta <- as.data.frame(vroom(path, show_col_types = FALSE, n_max = 1))
  if (ncol(meta) < 4) stop("LD metadata file must have at least 4 columns (chrom, start, end, path): ", path)
  colnames(meta)[1:4] <- c("chrom", "start", "end", "path")
  rawPath <- gsub(",.*$", "", meta$path[1])  # strip comma-separated bim path
  resolved <- file.path(dirname(path), rawPath)

  if (hasPlink2Files(resolved)) return(list(type = "plink2", data_path = resolved, meta_path = path))
  if (hasPlink1Files(resolved)) return(list(type = "plink1", data_path = resolved, meta_path = path))
  if (isVcfPath(resolved)) return(list(type = "vcf", data_path = resolved, meta_path = path))
  if (isGdsPath(resolved)) return(list(type = "gds", data_path = resolved, meta_path = path))

  # Pre-computed .cor.xz blocks - verify not using 0:0 sentinel
  if (!is.na(meta$start) && !is.na(meta$end) && meta$start == 0 && meta$end == 0) {
    stop("Metadata has start=0, end=0 but path does not resolve to genotype files: ", resolved,
         "\n  The 0:0 sentinel is only valid for whole-chromosome genotype files.")
  }

  list(type = "precomputed", meta_path = path)
}

#' Resolve the correct genotype path for a given region from a metadata TSV.
#' Reads the TSV, finds the row matching the query region's chromosome,
#' and returns the resolved genotype file path or prefix.
#' @importFrom vroom vroom
#' @noRd
resolveGenotypePathForRegion <- function(metaPath, region) {
  parsed <- parseRegion(region)
  meta <- as.data.frame(vroom(metaPath, show_col_types = FALSE))
  colnames(meta) <- c("chrom", "start", "end", "path")
  meta$chrom <- as.integer(stripChrPrefix(meta$chrom))
  queryChrom <- as.integer(stripChrPrefix(parsed$chrom))

  matching <- meta[meta$chrom == queryChrom, , drop = FALSE]
  if (nrow(matching) == 0) {
    stop("No entry for chromosome ", queryChrom, " in metadata file: ", metaPath)
  }
  rawPath <- gsub(",.*$", "", matching$path[1])
  file.path(dirname(metaPath), rawPath)
}

# ---------- Internal: load LD from genotype files ----------

#' Load genotype data and compute LD or return genotype matrix.
#' @noRd
loadLdFromGenotype <- function(genotypePath, region,
                                  returnGenotype = FALSE, nSample = NULL) {
  # Load genotype matrix and variant info via the unified loader
  result <- loadGenotypeRegion(genotypePath, region = region,
                                 returnVariantInfo = TRUE)
  X <- result$X
  variantInfo <- result$variant_info

  # Normalize variant IDs to canonical format (chr:pos:A2:A1)
  variantIds <- normalizeVariantId(
    formatVariantId(variantInfo$chrom, variantInfo$pos, variantInfo$A2, variantInfo$A1)
  )
  colnames(X) <- variantIds

  # Build ref_panel
  refPanel <- parseVariantId(variantIds)
  refPanel$variant_id <- variantIds

  # Load allele frequency from .afreq file if available, otherwise compute from genotypes
  afreq <- readAfreq(genotypePath)
  if (!is.null(afreq)) {
    freqMatch <- match(variantInfo$id, afreq$id)
    nUnmatched <- sum(is.na(freqMatch))
    if (nUnmatched > 0) {
      warning(nUnmatched, " out of ", length(freqMatch),
              " variants have no allele frequency in .afreq file.")
    }
    refPanel$allele_freq <- afreq$alt_freq[freqMatch]
  } else {
    # Compute ALT allele frequency directly from the dosage matrix
    refPanel$allele_freq <- colMeans(X, na.rm = TRUE) / 2
  }

  # Compute variance if sample size provided
  if (!is.null(nSample)) {
    p <- refPanel$allele_freq
    refPanel$variance <- 2 * p * (1 - p) * nSample / (nSample - 1)
    refPanel$n_nomiss <- nSample
  }

  # Block metadata (single block spanning the loaded region)
  positions <- variantInfo$pos
  blockMetadata <- data.frame(
    block_id = 1L,
    chrom = as.character(variantInfo$chrom[1]),
    block_start = min(positions),
    block_end = max(positions),
    size = length(variantIds),
    start_idx = 1L,
    end_idx = length(variantIds),
    stringsAsFactors = FALSE
  )

  # Build variant GRanges for LdData
  variantsGr <- .refPanelToGranges(refPanel)

  if (returnGenotype) {
    # Store genotype handle + snp_idx for lazy access
    handle <- readGenotypes(genotypePath)
    snpIdx <- .regionToSnpIdx(handle@snpInfo, region)
    return(LdData(
      correlation = NULL,
      genotypeHandle = handle,
      snpIdx = snpIdx,
      variants = variantsGr,
      blockMetadata = blockMetadata,
      nRef = as.integer(nrow(X))
    ))
  }

  R <- computeLd(X, method = "sample")

  LdData(
    correlation = R,
    genotypeHandle = NULL,
    snpIdx = NULL,
    variants = variantsGr,
    blockMetadata = blockMetadata,
    nRef = as.integer(nrow(X))
  )
}

# ---------- LD sketch: genotype loading ----------

#' HWE-based standardization of a genotype matrix
#'
#' Centers by 2*alleleFreq, scales by sqrt(2*alleleFreq*(1-alleleFreq)).
#' Assumes monomorphic variants have already been removed.
#'
#' @param X Numeric genotype matrix (n x p).
#' @param alleleFreq Numeric vector of allele frequencies (length p).
#' @return Standardized matrix (n x p).
#' @noRd
standardizeGenotypeHwe <- function(X, alleleFreq) {
  Xstd <- sweep(X, 2, 2 * alleleFreq)
  sweep(Xstd, 2, sqrt(2 * alleleFreq * (1 - alleleFreq)), "/")
}

#' Load LD sketch genotypes for a region
#'
#' Loads genotype data for a region via \code{loadLdMatrix(returnGenotype=TRUE)}
#' and removes monomorphic variants. Returns the raw genotype matrix and metadata,
#' which callers can use to derive either a correlation matrix R (for summary-based
#' weight training or fine-mapping) or an SVD (for TWAS z-score computation).
#'
#' @param ldMetaFilePath Path to the LD metadata TSV file.
#' @param region Region of interest: "chr:start-end" string or data.frame with chrom/start/end.
#' @param nSample Optional original panel sample size for computing variance
#'   (= 2*p*(1-p)*n/(n-1)). Passed through to \code{loadLdMatrix()}.
#'
#' @return An \code{LdData} S4 object with monomorphic variants removed.
#'   Consumers should use S4 accessors: \code{getGenotypes()}, \code{getRefPanel()},
#'   \code{getVariantIds()}. The number of sketch samples is
#'   \code{nrow(getGenotypes(result))}.
#' @export
loadLdSketch <- function(ldMetaFilePath, region, nSample = NULL) {
  result <- loadLdMatrix(ldMetaFilePath, region, returnGenotype = TRUE, nSample = nSample)
  if (!is(result, "LdData")) {
    stop("loadLdMatrix must return an LdData object")
  }
  X <- getGenotypes(result)
  refPanel <- getRefPanel(result)

  # Remove monomorphic variants (zero variance under HWE)
  p <- refPanel$allele_freq
  polymorphic <- p > 0 & p < 1
  if (!all(polymorphic)) {
    X <- X[, polymorphic, drop = FALSE]
    refPanel <- refPanel[polymorphic, , drop = FALSE]
  }

  # Rebuild LdData with the extracted (and filtered) genotype matrix stored
  # directly in genotypeHandle so getGenotypes() returns it without needing
  # the original file handle.
  variantsGr <- .refPanelToGranges(refPanel)
  LdData(
    correlation = NULL,
    genotypeHandle = X,
    snpIdx = NULL,
    variants = variantsGr,
    blockMetadata = getBlockMetadata(result),
    nRef = result@nRef
  )
}

# ---------- Internal: load LD from pre-computed blocks ----------

#' Load pre-computed LD from block-based metadata files.
#' @noRd
loadLdFromBlocks <- function(ldMetaFilePath, region, extractCoordinates = NULL, nSample = NULL) {
  # Intersect LD metadata with specified regions
  intersectedLdFiles <- getRegionalLdMeta(ldMetaFilePath, region)

  ldFilePaths <- intersectedLdFiles$intersections$LD_file_paths
  bimFilePaths <- intersectedLdFiles$intersections$bim_file_paths

  extractedLdMatricesList <- list()
  extractedLdVariantsList <- list()
  blockChroms <- character(length(ldFilePaths))

  for (j in seq_along(ldFilePaths)) {
    ldMatrixProcessed <- processLdMatrix(ldFilePaths[j], bimFilePaths[j])
    extractedLdList <- extractLdForRegion(
      ldMatrix = ldMatrixProcessed$LD_matrix,
      variants = ldMatrixProcessed$LD_variants,
      region = intersectedLdFiles$region,
      extractCoordinates = extractCoordinates
    )
    extractedLdMatricesList[[j]] <- extractedLdList$extracted_LD_matrix
    extractedLdVariantsList[[j]] <- extractedLdList$extracted_LD_variants
    if (nrow(extractedLdVariantsList[[j]]) > 0) {
      blockChroms[j] <- as.character(extractedLdVariantsList[[j]]$chrom[1])
    } else {
      blockChroms[j] <- as.character(intersectedLdFiles$region$chrom)
    }
  }

  # Filter out empty blocks before combining
  nonEmpty <- sapply(extractedLdVariantsList, function(v) nrow(v) > 0)
  if (!any(nonEmpty)) {
    stop("No variants found in any LD block for the specified region.")
  }
  if (any(!nonEmpty)) {
    message(paste(
      "Removing", sum(!nonEmpty), "empty LD block(s) with no variants in the region."
    ))
    extractedLdMatricesList <- extractedLdMatricesList[nonEmpty]
    extractedLdVariantsList <- extractedLdVariantsList[nonEmpty]
    blockChroms <- blockChroms[nonEmpty]
    ldFilePaths <- ldFilePaths[nonEmpty]
  }

  ldMatrix <- createLdMatrix(
    ldMatrices = extractedLdMatricesList,
    variants = extractedLdVariantsList
  )
  ldVariants <- rownames(ldMatrix)

  blockVariants <- lapply(extractedLdVariantsList, function(v) v$variants)
  blockPositions <- lapply(extractedLdVariantsList, function(v) v$pos)
  blockMetadata <- data.frame(
    block_id = seq_along(ldFilePaths),
    chrom = blockChroms,
    block_start = sapply(blockPositions, min),
    block_end = sapply(blockPositions, max),
    size = sapply(blockVariants, length),
    start_idx = sapply(blockVariants, function(v) min(match(v, ldVariants))),
    end_idx = sapply(blockVariants, function(v) max(match(v, ldVariants))),
    stringsAsFactors = FALSE
  )

  rm(extractedLdMatricesList)

  refPanel <- parseVariantId(rownames(ldMatrix))
  mergedVariantList <- do.call(rbind, extractedLdVariantsList)
  refPanel$variant_id <- rownames(ldMatrix)

  if ("allele_freq" %in% colnames(mergedVariantList)) {
    refPanel$allele_freq <- mergedVariantList$allele_freq[match(rownames(ldMatrix), mergedVariantList$variants)]
  }
  if ("variance" %in% colnames(mergedVariantList)) {
    refPanel$variance <- mergedVariantList$variance[match(rownames(ldMatrix), mergedVariantList$variants)]
  }
  if ("n_nomiss" %in% colnames(mergedVariantList)) {
    refPanel$n_nomiss <- mergedVariantList$n_nomiss[match(rownames(ldMatrix), mergedVariantList$variants)]
  }

  # Compute variance from nSample + allele_freq if not already present
  if (!is.null(nSample) && (!"variance" %in% colnames(refPanel) || all(is.na(refPanel$variance)))) {
    if ("allele_freq" %in% colnames(refPanel)) {
      p <- refPanel$allele_freq
      refPanel$variance <- 2 * p * (1 - p) * nSample / (nSample - 1)
      refPanel$n_nomiss <- nSample
    }
  }

  variantsGr <- .refPanelToGranges(refPanel)

  LdData(
    correlation = ldMatrix,
    genotypeHandle = NULL,
    snpIdx = NULL,
    variants = variantsGr,
    blockMetadata = blockMetadata,
    nRef = if (is.null(nSample)) 0L else as.integer(nSample)
  )
}

#' Filter variants by LD Reference
#'
#' Filters a vector of variant IDs to those present in the LD reference panel.
#' Auto-detects the reference type (PLINK2, PLINK1, or pre-computed LD metadata).
#'
#' @param variantIds variant names in the format chr:pos:ref:alt.
#' @param ldReferenceMetaFile Path to LD metadata file or PLINK prefix.
#' @param keepIndel Whether to keep indel variants. Default TRUE.
#' @return A list with:
#'   \item{data}{Character vector of filtered variant IDs.}
#'   \item{idx}{Integer vector of indices into the original variantIds.}
#' @importFrom dplyr group_by summarise
#' @importFrom vroom vroom
#' @importFrom magrittr %>%
#' @export
filterVariantsByLdReference <- function(variantIds, ldReferenceMetaFile, keepIndel = TRUE) {
  variantsDf <- parseVariantId(variantIds)

  # Derive region to scope the reference lookup
  regionDf <- variantsDf %>%
    group_by(chrom) %>%
    summarise(start = min(pos), end = max(pos))

  # Use shared helper -- no genotype loading
  refInfo <- getRefVariantInfo(ldReferenceMetaFile, regionDf)
  refChrom <- as.integer(stripChrPrefix(refInfo$chrom))
  refKey <- paste0(refChrom, ":", refInfo$pos)

  variantKey <- paste0(variantsDf$chrom, ":", variantsDf$pos)
  keepIndices <- which(variantKey %in% refKey)

  if (!keepIndel) {
    snpIdx <- which(isSnpAlleles(variantsDf$A1, variantsDf$A2))
    keepIndices <- intersect(keepIndices, snpIdx)
  }

  nDropped <- length(variantIds) - length(keepIndices)
  if (nDropped > 0) {
    message(nDropped, " out of ", length(variantIds),
            " total variants dropped due to absence on the reference LD panel.")
  }

  list(data = variantIds[keepIndices], idx = keepIndices)
}

#' Partition LD Matrix into Block-Specific Matrices
#'
#' This function takes the output from loadLdMatrix and partitions the combined LD matrix
#' into a list of smaller matrices based on the block_indices, making it easier to work with
#' large LD matrices that span multiple blocks.
#'
#' @param ldData An \code{LdData} S4 object as returned by \code{loadLdMatrix()}.
#' @param mergeSmallBlocks Logical, whether to merge blocks smaller than minMergedBlockSize (default: TRUE).
#' @param minMergedBlockSize Integer, minimum number of variants for a block after merging (default: 500).
#' @param maxMergedBlockSize Integer, maximum number of variants in a block after merging (default: 10000).
#'
#' @return returns a list containing:
#' \describe{
#' \item{ld_matrices}{A list of matrices, each representing LD for a specific block.}
#' \item{variant_indices}{A data frame that maps variant IDs to their corresponding block.}
#' \item{block_metadata}{Information about each block including size, chromosome, start and end positions.}
#' }
#' @noRd
partitionLdMatrix <- function(ldData, mergeSmallBlocks = TRUE,
                                minMergedBlockSize = 500, maxMergedBlockSize = 10000) {
  if (!is(ldData, "LdData")) {
    stop("ldData must be an LdData object")
  }
  combinedMatrix <- getCorrelation(ldData)
  blockMetadata <- ldData@blockMetadata
  if (is(blockMetadata, "LdBlocks")) {
    blockMetadata <- as.data.frame(blockMetadata@blocks)
  }
  variantIds <- getVariantIds(ldData)

  # Error if matrix is empty
  if (is.null(combinedMatrix) || nrow(combinedMatrix) == 0 || ncol(combinedMatrix) == 0) {
    stop("Empty or NULL LD matrix provided.")
  }

  # Ensure the row and column names of the matrix match the variantIds
  if (is.null(rownames(combinedMatrix)) || is.null(colnames(combinedMatrix)) ||
    !identical(rownames(combinedMatrix), variantIds) || !identical(colnames(combinedMatrix), variantIds)) {
    rownames(combinedMatrix) <- variantIds
    colnames(combinedMatrix) <- variantIds
  }

  # Filter out blocks with invalid indices (empty blocks, out-of-range, NA, Inf)
  nVariants <- length(variantIds)
  validBlocks <- sapply(seq_len(nrow(blockMetadata)), function(i) {
    s <- blockMetadata$start_idx[i]
    e <- blockMetadata$end_idx[i]
    sz <- blockMetadata$size[i]
    # Block is valid if: size > 0, indices are finite integers, and within range
    !is.na(s) && !is.na(e) && is.finite(s) && is.finite(e) &&
      sz > 0 && s >= 1 && e >= s && e <= nVariants
  })

  if (!any(validBlocks)) {
    stop("No valid LD blocks found. All block indices are out of range or empty.")
  }

  if (any(!validBlocks)) {
    message(paste(
      "Removing", sum(!validBlocks),
      "LD block(s) with invalid or out-of-range indices."
    ))
    blockMetadata <- blockMetadata[validBlocks, , drop = FALSE]
    blockMetadata$block_id <- seq_len(nrow(blockMetadata))
  }

  # Validate the block structure of the matrix (skip if only one block)
  if (nrow(blockMetadata) > 1) {
    validateBlockStructure(combinedMatrix, blockMetadata, variantIds)
  }

  # Optionally merge small blocks
  if (mergeSmallBlocks && any(blockMetadata$size < minMergedBlockSize) && nrow(blockMetadata) > 1) {
    blockMetadata <- mergeBlocks(blockMetadata, minMergedBlockSize, maxMergedBlockSize)
  }

  # Partition the matrix based on block metadata
  result <- extractBlockMatrices(combinedMatrix, blockMetadata, variantIds)
  return(result)
}

#' Validate that cross-block entries are zero (excluding boundary variants).
#' @noRd
validateBlockStructure <- function(matrix, blockMetadata, variantIds) {
  msgs <- character(0)
  n <- length(variantIds)

  for (i in 1:(nrow(blockMetadata) - 1)) {
    for (j in (i + 1):nrow(blockMetadata)) {
      si <- blockMetadata$start_idx[i]; ei <- blockMetadata$end_idx[i]
      sj <- blockMetadata$start_idx[j]; ej <- blockMetadata$end_idx[j]
      if (si > n || ei > n || sj > n || ej > n) {
        msgs <- c(msgs, paste("Block indices out of range for blocks", i, "and", j))
        next
      }
      # Exclude boundary variants (potential overlaps)
      vi <- variantIds[si:(ei - 1)]
      vj <- variantIds[(sj + 1):ej]
      if (length(vi) > 0 && length(vj) > 0) {
        maxVal <- max(abs(matrix[vi, vj, drop = FALSE]))
        if (maxVal > 1e-10) {
          msgs <- c(msgs, paste("Non-zero correlation between blocks", i, "and", j,
                                "- max:", maxVal))
        }
      }
    }
  }
  if (length(msgs) > 0) stop("Matrix lacks expected block structure:\n", paste(msgs, collapse = "\n"))
}

#' @noRd
canMerge <- function(block1, block2, maxSize) {
  block1$chrom == block2$chrom && (block1$size + block2$size) <= maxSize
}

#' @noRd
mergeTwoBlocks <- function(blockMetadata, idx1, idx2) {
  if (idx1 > idx2) { tmp <- idx1; idx1 <- idx2; idx2 <- tmp }
  result <- blockMetadata
  result$end_idx[idx1] <- blockMetadata$end_idx[idx2]
  result$size[idx1] <- blockMetadata$size[idx1] + blockMetadata$size[idx2]
  result <- result[-idx2, ]
  result$block_id <- seq_len(nrow(result))
  result
}

#' Find blocks below minSize and identify the best neighbor to merge with.
#' @noRd
findMergeCandidates <- function(blockMetadata, minSize, maxSize) {
  candidates <- data.frame(block_idx = integer(), merge_with = integer(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(blockMetadata))) {
    if (blockMetadata$size[i] >= minSize) next
    prevOk <- i > 1 && canMerge(blockMetadata[i, ], blockMetadata[i - 1, ], maxSize)
    nextOk <- i < nrow(blockMetadata) && canMerge(blockMetadata[i, ], blockMetadata[i + 1, ], maxSize)
    mergeWith <- if (prevOk && nextOk) {
      if (blockMetadata$size[i - 1] <= blockMetadata$size[i + 1]) i - 1 else i + 1
    } else if (prevOk) i - 1
      else if (nextOk) i + 1
      else next
    candidates <- rbind(candidates, data.frame(block_idx = i, merge_with = mergeWith))
  }
  candidates
}

#' Iteratively merge blocks below minSize with their smallest neighbor.
#' @noRd
mergeBlocks <- function(blockMetadata, minSize, maxSize) {
  if (nrow(blockMetadata) <= 1) return(blockMetadata)
  repeat {
    candidates <- findMergeCandidates(blockMetadata, minSize, maxSize)
    if (nrow(candidates) == 0) break
    blockMetadata <- mergeTwoBlocks(blockMetadata, candidates$block_idx[1], candidates$merge_with[1])
  }
  blockMetadata
}

# Helper function to extract block matrices
extractBlockMatrices <- function(matrix, blockMetadata, variantIds) {
  ldMatrices <- list()
  variantMapping <- data.frame(
    variant_id = character(),
    block_id = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(blockMetadata))) {
    startIdx <- blockMetadata$start_idx[i]
    endIdx <- blockMetadata$end_idx[i]

    # Skip empty blocks
    if (endIdx < startIdx) next

    # Ensure indices are within bounds
    if (startIdx > length(variantIds) || endIdx > length(variantIds)) {
      warning(paste("Block", i, "has indices outside the range of variantIds. Skipping."))
      next
    }

    # Extract variant IDs for this block
    blockVariants <- variantIds[startIdx:endIdx]

    # Extract submatrix for this block - use named indexing
    blockMatrix <- matrix[blockVariants, blockVariants, drop = FALSE]

    # Store in list
    ldMatrices[[i]] <- blockMatrix

    # Update variant mapping
    blockMapping <- data.frame(
      variant_id = blockVariants,
      block_id = i,
      stringsAsFactors = FALSE
    )
    variantMapping <- rbind(variantMapping, blockMapping)

  }

  return(list(
    ld_matrices = ldMatrices,
    variant_indices = variantMapping,
    block_metadata = blockMetadata
  ))
}


#' Check and optionally repair LD matrix quality
#'
#' Diagnoses positive-definiteness of an LD correlation matrix and optionally
#' repairs it. Downstream methods like PRS-CS require positive-definite LD
#' (Cholesky decomposition), while others (lassosum, SDPR) handle non-PD
#' matrices internally via their own regularization.
#'
#' Three modes are available:
#' \describe{
#'   \item{\code{"check"}}{Diagnostic only - returns eigenvalue statistics
#'     without modifying the matrix.}
#'   \item{\code{"shrink"}}{Apply shrinkage toward identity:
#'     \code{R_s = (1 - shrinkage) * R + shrinkage * I}. Simple and fast;
#'     always produces a positive-definite matrix when \code{shrinkage > 0}.}
#'   \item{\code{"eigenfix"}}{Set negative eigenvalues to zero and
#'     reconstruct the matrix. Matches the approach used in susieR's
#'     \code{rss_lambda_constructor} and is the closest positive
#'     semidefinite matrix in the Frobenius norm. Does not inflate the
#'     diagonal like shrinkage does.}
#' }
#'
#' @param R Symmetric correlation matrix.
#' @param method One of \code{"check"}, \code{"shrink"}, or \code{"eigenfix"}.
#' @param rTol Eigenvalue tolerance. Eigenvalues with absolute value below
#'   \code{rTol} are treated as zero. Default: \code{1e-8}.
#' @param shrinkage Shrinkage parameter for \code{method = "shrink"}.
#'   Default: \code{0.01}.
#'
#' @return A list with components:
#' \describe{
#'   \item{R}{The (possibly repaired) LD matrix.}
#'   \item{is_pd}{Logical: is the matrix positive definite?}
#'   \item{is_psd}{Logical: is the matrix positive semidefinite (within rTol)?}
#'   \item{min_eigenvalue}{Smallest eigenvalue of the original matrix.}
#'   \item{n_negative}{Number of negative eigenvalues (below -rTol).}
#'   \item{condition_number}{Ratio of largest to smallest positive eigenvalue
#'     (\code{Inf} if any eigenvalue is zero).}
#'   \item{method_applied}{Character: \code{"none"}, \code{"shrink"}, or
#'     \code{"eigenfix"}.}
#' }
#'
#' @examples
#' # A well-conditioned matrix
#' R_good <- diag(5)
#' checkLd(R_good)$is_pd  # TRUE
#'
#' # A matrix with negative eigenvalues
#' R_bad <- matrix(0.9, 3, 3); diag(R_bad) <- 1; R_bad[1,3] <- R_bad[3,1] <- -0.5
#' checkLd(R_bad)$is_psd  # FALSE
#' R_fixed <- checkLd(R_bad, method = "eigenfix")$R
#' checkLd(R_fixed)$is_psd  # TRUE
#'
#' @export
checkLd <- function(R,
                     method = c("check", "shrink", "eigenfix"),
                     rTol = 1e-8,
                     shrinkage = 0.01) {
  method <- match.arg(method)
  p <- nrow(R)

  # Eigen decomposition (symmetric)
  eig <- eigen(R, symmetric = TRUE)
  vals <- eig$values

  # Diagnostics
  minEval <- min(vals)
  nNeg <- sum(vals < -rTol)
  posVals <- vals[vals > rTol]
  condNum <- if (length(posVals) > 0) max(posVals) / min(posVals) else Inf
  isPsd <- !any(vals < -rTol)
  isPd <- all(vals > rTol)

  methodApplied <- "none"
  Rout <- R

  if (method == "shrink" && !isPd) {
    Rout <- (1 - shrinkage) * R + shrinkage * diag(p)
    methodApplied <- "shrink"
  } else if (method == "eigenfix" && !isPd) {
    # Set negative eigenvalues to a small positive value and reconstruct.
    # Using rTol (not zero) ensures the result is strictly positive
    # definite, which is required by methods that use Cholesky decomposition
    # (PRS-CS, SDPR). Setting to exactly zero would produce PSD but not PD.
    valsFixed <- pmax(vals, rTol)
    Rout <- eig$vectors %*% diag(valsFixed) %*% t(eig$vectors)
    # Restore exact symmetry and unit diagonal
    Rout <- (Rout + t(Rout)) / 2
    diag(Rout) <- 1
    methodApplied <- "eigenfix"
  }

  list(
    R = Rout,
    is_pd = isPd,
    is_psd = isPsd,
    min_eigenvalue = minEval,
    n_negative = nNeg,
    condition_number = condNum,
    method_applied = methodApplied
  )
}

#' Prune columns by pairwise correlation (LD-style prune)
#'
#' Performs LD pruning using one of two backends. The default \code{"hclust"}
#' backend computes the full correlation matrix, builds a single-linkage
#' hierarchical clustering on the distance (1 - |cor|), and keeps one
#' representative column per cluster. The \code{"snprelate"} backend delegates
#' to \code{SNPRelate::snpgdsLDpruning}, which performs a sliding-window
#' greedy prune directly on a temporary GDS file.
#'
#' @param X Numeric matrix. Columns are the variables to prune (typically SNP
#'   genotype dosages); rows are observations.
#' @param corThres Numeric in (0, 1). Absolute correlation threshold.
#'   Columns whose pairwise |cor| exceeds this are grouped; one survivor is
#'   kept per group. Default 0.8.
#' @param backend Character, one of \code{"hclust"} (default) or
#'   \code{"snprelate"}. Controls the pruning algorithm:
#'   \describe{
#'     \item{\code{"hclust"}}{Uses the internal hierarchical-clustering approach
#'       with \code{Rfast::cora} (if available) or base \code{cor()}.}
#'     \item{\code{"snprelate"}}{Requires \pkg{SNPRelate} and \pkg{gdsfmt}.
#'       Creates a temporary GDS file and runs
#'       \code{SNPRelate::snpgdsLDpruning(method = "corr")}.}
#'   }
#' @param verbose Logical. If TRUE, print progress messages. Default FALSE.
#'
#' @return A list with:
#'   \describe{
#'     \item{X.new}{Matrix containing the retained columns of \code{X}.}
#'     \item{filter.id}{Integer vector of the column indices of \code{X} that
#'       were retained (in original order).}
#'   }
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' X[, 2] <- X[, 1] + rnorm(100, sd = 0.01)   # near-duplicate of col 1
#' res <- ldPruneByCorrelation(X, corThres = 0.9)
#' ncol(res$X.new)
#'
#' @importFrom stats as.dist hclust cutree cor
#' @export
ldPruneByCorrelation <- function(X, corThres = 0.8,
                                    backend = c("hclust", "snprelate"),
                                    verbose = FALSE) {
  backend <- match.arg(backend)
  p <- ncol(X)

  if (backend == "snprelate") {
    return(.ldPruneSnprelate(X, corThres = corThres, verbose = verbose))
  }

  # ---- hclust backend (default) ----
  if (requireNamespace("Rfast", quietly = TRUE)) {
    cor.X <- Rfast::cora(X, large = TRUE)
  } else {
    cor.X <- cor(X)
  }
  Sigma.distance <- as.dist(1 - abs(cor.X))
  fit <- hclust(Sigma.distance, method = "single")
  clusters <- cutree(fit, h = 1 - corThres)
  groups <- unique(clusters)
  ind.delete <- NULL
  X.new <- X
  filter.id <- seq_len(p)
  for (ig in seq_along(groups)) {
    temp.group <- which(clusters == groups[ig])
    if (length(temp.group) > 1) {
      ind.delete <- c(ind.delete, temp.group[-1])
    }
  }
  ind.delete <- unique(ind.delete)
  if (length(ind.delete) > 0) {
    X.new <- as.matrix(X[, -ind.delete])
    filter.id <- filter.id[-ind.delete]
    if (verbose) {
      message("ldPruneByCorrelation: pruned ", length(ind.delete),
              " of ", p, " columns at |cor| > ", corThres)
    }
  } else if (verbose) {
    message("ldPruneByCorrelation: no columns pruned at |cor| > ", corThres)
  }

  if (ncol(X.new) == 1) {
    colnames(X.new) <- colnames(X)[-ind.delete]
  }

  list(X.new = X.new, filter.id = filter.id)
}

#' SNPRelate-based LD pruning helper
#' @noRd
.ldPruneSnprelate <- function(X, corThres, verbose) {
  if (!requireNamespace("SNPRelate", quietly = TRUE) ||
      !requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Packages 'SNPRelate' and 'gdsfmt' are required for backend='snprelate'.")
  }
  p <- ncol(X)
  snpNames <- colnames(X) %||% paste0("snp", seq_len(p))

  # Round dosages to integer genotype codes for GDS
  genoInt <- round(X)
  storage.mode(genoInt) <- "integer"

  tmpGds <- tempfile(fileext = ".gds")
  on.exit(unlink(tmpGds), add = TRUE)

  SNPRelate::snpgdsCreateGeno(
    gds.fn = tmpGds,
    genmat = t(genoInt),
    sample.id = seq_len(nrow(X)),
    snp.id = seq_len(p),
    snp.rs.id = snpNames,
    snp.chromosome = rep(1L, p),
    snp.position = seq_len(p),
    snpfirstdim = TRUE
  )

  gds <- SNPRelate::snpgdsOpen(tmpGds, allow.duplicate = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  keepList <- SNPRelate::snpgdsLDpruning(
    gds,
    method = "corr",
    ld.threshold = corThres,
    verbose = verbose
  )

  keepIds <- sort(unlist(keepList, use.names = FALSE))
  filter.id <- keepIds
  X.new <- X[, keepIds, drop = FALSE]

  if (verbose) {
    message("ldPruneByCorrelation (snprelate): kept ", length(keepIds),
            " of ", p, " columns at |cor| > ", corThres)
  }

  list(X.new = X.new, filter.id = filter.id)
}

#' Drop collinear columns from a design matrix by a chosen strategy
#'
#' Given a numeric matrix \code{X} and a set of column names known to be
#' involved in linear dependencies, remove one column using one of three
#' strategies. Designed to be called iteratively by
#' \code{\link{enforceDesignFullRank}}, but can be used standalone.
#'
#' @param X Numeric matrix. Must have column names covering
#'   \code{problematicCols}.
#' @param problematicCols Character vector of column names in \code{X} that
#'   are candidates for removal. If empty, \code{X} is returned unchanged.
#' @param strategy One of \code{"correlation"} (remove the column with the
#'   largest sum of absolute pairwise correlations among the candidates;
#'   when only two candidates, one is picked at random), \code{"variance"}
#'   (remove the lowest-variance candidate), or \code{"response_correlation"}
#'   (remove the candidate whose correlation with \code{response} has the
#'   smallest magnitude).
#' @param response Numeric vector required when \code{strategy =
#'   "response_correlation"}; the outcome to correlate against.
#' @param verbose Logical. If TRUE, print which column was removed. Default
#'   FALSE.
#'
#' @return \code{X} with exactly one column removed (or unchanged if
#'   \code{problematicCols} is empty).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 3), 100, 3)
#' X[, 3] <- X[, 1] + X[, 2]
#' colnames(X) <- c("a", "b", "c")
#' dropCollinearColumns(X, problematicCols = c("a", "b", "c"),
#'                        strategy = "variance")
#'
#' @importFrom stats var cor
#' @keywords internal
#' @noRd
dropCollinearColumns <- function(X, problematicCols,
                                   strategy = c("correlation", "variance", "response_correlation"),
                                   response = NULL, verbose = FALSE) {
  strategy <- match.arg(strategy)

  if (length(problematicCols) == 0) {
    return(X)
  }

  if (length(problematicCols) == 1) {
    colToRemove <- problematicCols[1]
    if (verbose) message("dropCollinearColumns: removing single column ", colToRemove)
    X <- X[, !(colnames(X) %in% colToRemove), drop = FALSE]
    return(X)
  }

  if (strategy == "variance") {
    variances <- apply(X[, problematicCols, drop = FALSE], 2, var)
    colToRemove <- problematicCols[which.min(variances)]
    if (verbose) message("dropCollinearColumns: smallest variance -> removing ", colToRemove)
  } else if (strategy == "correlation") {
    corMatrix <- abs(cor(X[, problematicCols, drop = FALSE]))
    diag(corMatrix) <- 0

    if (length(problematicCols) == 2) {
      colToRemove <- sample(problematicCols, 1)
      if (verbose) message("dropCollinearColumns: two candidates, randomly removing ", colToRemove)
    } else {
      corSums <- colSums(corMatrix)
      colToRemove <- problematicCols[which.max(corSums)]
      if (verbose) message("dropCollinearColumns: highest sum |cor| -> removing ", colToRemove)
    }
  } else if (strategy == "response_correlation") {
    if (is.null(response)) {
      stop("response must be supplied for strategy = 'response_correlation'")
    }
    corWithResponse <- apply(X[, problematicCols, drop = FALSE], 2,
                               function(col) cor(col, response))
    colToRemove <- problematicCols[which.min(abs(corWithResponse))]
    if (verbose) message("dropCollinearColumns: smallest |cor| with response -> removing ", colToRemove)
  }

  X[, !(colnames(X) %in% colToRemove), drop = FALSE]
}

#' Iteratively enforce full column rank on a design matrix
#'
#' Given a candidate predictor matrix \code{X} and an optional unnamed
#' covariate matrix \code{C}, builds the design \code{[1, X, C]} and removes
#' rank-deficient columns from \code{X} until the design has full column rank.
#' Rank-deficient columns are identified via the pivot of
#' \code{qr([1, X, C])}. On each iteration, one problematic column is dropped
#' using \code{\link{dropCollinearColumns}}. If iterative pruning does not
#' achieve full rank, falls back to \code{\link{ldPruneByCorrelation}} at a
#' descending sequence of correlation thresholds.
#'
#' @param X Numeric matrix with column names (the predictors subject to
#'   pruning).
#' @param C Numeric matrix of covariates (can be unnamed) that will be kept.
#'   Pass \code{NULL} or a zero-column matrix when there are no covariates.
#' @param strategy Passed through to \code{\link{dropCollinearColumns}}.
#' @param response Passed through to \code{\link{dropCollinearColumns}}
#'   when \code{strategy = "response_correlation"}.
#' @param maxIterations Integer. Hard cap on the iterative-prune loop.
#'   Default 300.
#' @param corrThresholds Numeric vector of |cor| thresholds used for the
#'   \code{\link{ldPruneByCorrelation}} fallback, tried in order.
#'   Default \code{seq(0.75, 0.5, by = -0.05)}.
#' @param verbose Logical. If TRUE, print per-iteration progress. Default
#'   FALSE.
#'
#' @return The pruned predictor matrix \code{X} (covariates \code{C} are not
#'   modified).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 4), 100, 4)
#' X[, 4] <- X[, 1] + X[, 2]          # rank-deficient
#' colnames(X) <- c("a", "b", "c", "d")
#' C <- matrix(rnorm(100), 100, 1)
#' X2 <- enforceDesignFullRank(X, C, strategy = "variance")
#' qr(cbind(1, X2, C))$rank == ncol(cbind(1, X2, C))
#'
#' @export
enforceDesignFullRank <- function(X, C,
                                     strategy = c("correlation", "variance", "response_correlation"),
                                     response = NULL,
                                     maxIterations = 300L,
                                     corrThresholds = seq(0.75, 0.5, by = -0.05),
                                     verbose = FALSE) {
  strategy <- match.arg(strategy)
  originalColnames <- colnames(X)
  initialNcol <- ncol(X)
  iteration <- 0L

  buildDesign <- function(X) {
    XD <- cbind(1, X, C)
    colnames(XD)[seq_len(ncol(X) + 1L)] <- c("Intercept", colnames(X))
    XD
  }

  Xdesign <- buildDesign(X)
  matrixRank <- qr(Xdesign)$rank
  if (verbose) {
    message("enforceDesignFullRank: initial rank ", matrixRank,
            " / ", ncol(Xdesign))
  }

  skipIterative <- FALSE

  # Fast path: try removing all QR-pivot-flagged columns at once.
  if (matrixRank < ncol(Xdesign)) {
    qrd <- qr(Xdesign)
    problematicCols <- qrd$pivot[(qrd$rank + 1L):ncol(Xdesign)]
    problematicColnames <- colnames(Xdesign)[problematicCols]
    problematicColnames <- problematicColnames[problematicColnames %in% colnames(X)]

    if (length(problematicColnames) > 0) {
      Xtemp <- X[, !(colnames(X) %in% problematicColnames), drop = FALSE]
      if (qr(buildDesign(Xtemp))$rank == ncol(buildDesign(Xtemp))) {
        if (verbose) {
          message("enforceDesignFullRank: full rank after batch-removing ",
                  length(problematicColnames), " column(s)")
        }
      } else {
        skipIterative <- TRUE
        if (verbose) {
          message("enforceDesignFullRank: batch removal insufficient, ",
                  "skipping to correlation-pruning fallback")
        }
      }
    }
  }

  # Iterative path.
  if (!skipIterative) {
    while (matrixRank < ncol(Xdesign) && iteration < maxIterations) {
      qrd <- qr(Xdesign)
      problematicCols <- qrd$pivot[(qrd$rank + 1L):ncol(Xdesign)]
      problematicColnames <- colnames(Xdesign)[problematicCols]
      problematicColnames <- problematicColnames[problematicColnames %in% colnames(X)]

      if (length(problematicColnames) == 0) break

      X <- dropCollinearColumns(X, problematicColnames, strategy = strategy,
                                  response = response, verbose = verbose)

      Xdesign <- buildDesign(X)
      matrixRank <- qr(Xdesign)$rank
      iteration <- iteration + 1L
      if (verbose) {
        message("enforceDesignFullRank: iter ", iteration,
                " rank ", matrixRank, " / ", ncol(Xdesign))
      }
    }

    if (iteration == maxIterations) {
      warning("enforceDesignFullRank: maxIterations reached; design may still be rank-deficient")
    }
  }

  # Correlation-threshold fallback.
  Xdesign <- buildDesign(X)
  matrixRank <- qr(Xdesign)$rank
  if (matrixRank < ncol(Xdesign)) {
    if (verbose) {
      message("enforceDesignFullRank: applying ldPruneByCorrelation fallback")
    }
    for (threshold in corrThresholds) {
      filterResult <- ldPruneByCorrelation(X, corThres = threshold,
                                               verbose = verbose)
      X <- filterResult$X.new
      Xdesign <- buildDesign(X)
      matrixRank <- qr(Xdesign)$rank
      if (verbose) {
        message("enforceDesignFullRank: threshold ", threshold,
                " -> rank ", matrixRank, " / ", ncol(Xdesign))
      }
      if (matrixRank == ncol(Xdesign)) break
    }
  }

  if (ncol(X) == 1L && initialNcol == 1L) {
    colnames(X) <- originalColnames
  }
  X
}

#' LD clumping by a per-variant score using bigsnpr
#'
#' Wraps \code{bigsnpr::snp_clumping} with the boilerplate of wrapping a
#' numeric dosage matrix into a \code{bigstatsr::FBM.code256} object and of
#' handling the common pitfall of a single-variant input.
#'
#' @param X Numeric matrix of 0/1/2 allele dosages, n rows by p variants.
#'   Column names are expected to be variant IDs but are not required.
#' @param score Numeric vector of length \code{ncol(X)}. Higher values favour
#'   retention during clumping (e.g. -log10 p, |Z|, MAF). May be \code{NULL},
#'   in which case bigsnpr falls back to minor allele frequency computed from
#'   \code{X}.
#' @param chr Integer or character vector of length \code{ncol(X)} giving the
#'   chromosome for each variant.
#' @param pos Integer vector of length \code{ncol(X)} giving the base-pair
#'   position for each variant.
#' @param r2 Numeric in (0, 1]. r-squared threshold for clumping (variants
#'   within \code{windowKb} whose r2 exceeds \code{r2} and have lower
#'   \code{score} are removed). Default 0.2.
#' @param windowKb Numeric. Window size in kilobases. Default is
#'   \code{100 / r2}, matching the common "ld-clump size = 100/r2" heuristic
#'   used in many GWAS pipelines.
#' @param verbose Logical. If TRUE, print the number of retained variants.
#'   Default FALSE.
#'
#' @return An integer vector of indices (into \code{X} columns) kept after
#'   clumping. For a single-column \code{X}, returns \code{1L}.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 500; p <- 20
#'   X <- matrix(rbinom(n * p, 2, 0.3), n, p)
#'   colnames(X) <- paste0("chr1:", seq_len(p) * 1000, ":A:G")
#'   s <- runif(p)
#'   chr <- rep(1L, p); pos <- seq_len(p) * 1000L
#'   keep <- ldClumpByScore(X, score = s, chr = chr, pos = pos, r2 = 0.2)
#' }
#'
#' @export
ldClumpByScore <- function(X, score, chr, pos, r2 = 0.2,
                              windowKb = 100 / r2, verbose = FALSE) {
  if (!requireNamespace("bigsnpr", quietly = TRUE)) {
    stop("Package 'bigsnpr' is required. Install from CRAN: install.packages('bigsnpr')")
  }
  if (!requireNamespace("bigstatsr", quietly = TRUE)) {
    stop("Package 'bigstatsr' is required. Install from CRAN: install.packages('bigstatsr')")
  }

  if (ncol(X) < 1L) stop("ldClumpByScore: X must have at least one column")
  if (!is.null(score) && length(score) != ncol(X)) {
    stop("ldClumpByScore: length(score) must equal ncol(X)")
  }
  if (length(chr) != ncol(X) || length(pos) != ncol(X)) {
    stop("ldClumpByScore: chr and pos must have length equal to ncol(X)")
  }

  if (ncol(X) == 1L) {
    if (verbose) message("ldClumpByScore: single variant, skipping clumping")
    return(1L)
  }

  if (inherits(X, "FBM")) {
    G <- X
  } else {
    codeVec <- c(0, 1, 2, rep(NA, 256L - 3L))
    G <- bigstatsr::FBM.code256(
      nrow = nrow(X), ncol = ncol(X),
      init = X, code = codeVec
    )
  }

  keep <- bigsnpr::snp_clumping(
    G = G,
    infos.chr = as.integer(chr),
    infos.pos = as.integer(pos),
    S = score,
    thr.r2 = r2,
    size = windowKb
  )

  if (verbose) {
    message("ldClumpByScore: ", length(keep), " / ", ncol(X),
            " variants retained at r2 <= ", r2)
  }
  keep
}
