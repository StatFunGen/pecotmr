#' Deduplicate and sort genomic regions by chromosome and start position.
#' @importFrom dplyr distinct arrange
#' @importFrom magrittr %>%
#' @noRd
orderDedupRegions <- function(df) {
    df$chrom <- canonChrom(df$chrom)
    df <- distinct(df, chrom, start, .keep_all = TRUE) %>%
        arrange(chromOrder(chrom), start)
    df
}

#' Find the first and last rows of genomicData that overlap a query region.
#' Clamps the query to the available data range before searching.
#' @importFrom dplyr filter arrange slice desc
#' @noRd
findIntersectionRows <- function(
    genomicData,
    regionChrom,
    regionStart,
    regionEnd
) {
    chromData <- genomicData %>% filter(chrom == regionChrom)
    if (nrow(chromData) == 0) {
        stop("No data for chromosome ", regionChrom)
    }

    # Clamp query to available range
    regionStart <- max(regionStart, min(chromData$start))
    regionEnd <- min(regionEnd, max(chromData$end))

    startRow <- genomicData %>%
        filter(
            chrom == regionChrom,
            start <= regionStart,
            end > regionStart
        ) %>%
        slice(1)
    endRow <- genomicData %>%
        filter(chrom == regionChrom, start < regionEnd, end >= regionEnd) %>%
        arrange(desc(end)) %>%
        slice(1)

    if (nrow(startRow) == 0 || nrow(endRow) == 0) {
        stop(
            "Region ",
            regionChrom,
            ":",
            regionStart,
            "-",
            regionEnd,
            " is not covered by any rows in the LD metadata."
        )
    }
    list(startRow = startRow, endRow = endRow)
}

#' Validate that startRow..endRow fully covers [regionStart, regionEnd].
#' @noRd
validateSelectedRegion <- function(startRow, endRow, regionStart, regionEnd) {
    if (startRow$start > regionStart || endRow$end < regionEnd) {
        stop(
            "Region ",
            regionStart,
            "-",
            regionEnd,
            " is not fully covered by the LD metadata ",
            "(available: ",
            startRow$start,
            "-",
            endRow$end,
            ")."
        )
    }
}

#' Extract values of a column for rows spanning the intersection range.
#' @noRd
extractFilePaths <- function(genomicData, intersectionRows, columnToExtract) {
    if (!columnToExtract %in% names(genomicData)) {
        stop("Column '", columnToExtract, "' not found in genomic data.")
    }
    idx <- which(
        genomicData$chrom == intersectionRows$startRow$chrom &
            genomicData$start >= intersectionRows$startRow$start &
            genomicData$start <= intersectionRows$endRow$start
    )
    genomicData[[columnToExtract]][idx]
}

# Internal: resolve a sidecar file path declared in an LD-meta TSV.
# Paths in the TSV's `path` column are conventionally written relative
# to the TSV's own directory, not the analysis CWD; try the path as
# given first, then `dirname(ldReferenceMetaFile)/<path>`, then the
# manifest itself as a fallback.
# @noRd
.findValidFilePath <- function(referenceFilePath, targetFilePath) {
    if (file.exists(targetFilePath)) {
        return(targetFilePath)
    }
    targetFullPath <- file.path(dirname(referenceFilePath), targetFilePath)
    if (file.exists(targetFullPath)) {
        return(targetFullPath)
    }
    if (file.exists(referenceFilePath)) {
        return(referenceFilePath)
    }
    stop(sprintf(
        paste0(
            "Both reference and target file paths do not work. Tried ",
            "paths: '%s' and '%s'"
        ),
        referenceFilePath,
        targetFullPath
    ))
}

# Vectorised .findValidFilePath over a vector of target paths.
# @noRd
.findValidFilePaths <- function(referenceFilePath, targetFilePaths) {
    map_chr(targetFilePaths, function(x) {
        .findValidFilePath(referenceFilePath, x)
    })
}

#' Find LD blocks overlapping a query region from a metadata TSV file.
#'
#' @param ldReferenceMetaFile TSV with columns chrom, start, end, path. The path
#'   column may be comma-separated: "ld_file,bim_file".
#' @param region "chr:start-end" string or data.frame with chrom/start/end.
#' @param completeCoverageRequired If TRUE, error when the region extends beyond
#'   available LD blocks.
#' @return A list with: intersections (LD_file_paths, bimFilePaths), ldMetaData,
#'   and parsed region.
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom vroom vroom
#' @noRd
# Split the comma-joined path column into LD (+ optional bim) path columns.
.regionalLdParsePaths <- function(genomicData) {
    filePath <- genomicData$path %>%
        str_split(",", simplify = TRUE) %>%
        data.frame() %>%
        `colnames<-`(
            if (ncol(.) == 2) {
                c("LD_file_path", "bim_file_path")
            } else {
                c("LD_file_path")
            }
        )
    cbind(genomicData, filePath) %>% select(-path)
}

# Resolve the LD (and optional bim) file paths for the intersected rows.
.regionalLdExtractPaths <- function(
    ldReferenceMetaFile,
    genomicData,
    intersectionRows
) {
    ldPaths <- .findValidFilePaths(
        ldReferenceMetaFile,
        extractFilePaths(genomicData, intersectionRows, "LD_file_path")
    )
    bimPaths <- if ("bim_file_path" %in% names(genomicData)) {
        .findValidFilePaths(
            ldReferenceMetaFile,
            extractFilePaths(genomicData, intersectionRows, "bim_file_path")
        )
    } else {
        NULL
    }
    list(ldPaths = ldPaths, bimPaths = bimPaths)
}

getRegionalLdMeta <- function(
    ldReferenceMetaFile,
    region,
    completeCoverageRequired = FALSE
) {
    genomicData <- vroom(ldReferenceMetaFile)
    region <- parseRegion(region)
    names(genomicData) <- c("chrom", "start", "end", "path")
    names(region) <- c("chrom", "start", "end")
    # Treat start=0, end=0 as "covers all regions" (whole-chromosome files).
    wholeChrom <- genomicData$start == 0 & genomicData$end == 0
    if (any(wholeChrom)) {
        genomicData$end[wholeChrom] <- Inf
    }
    genomicData <- orderDedupRegions(genomicData)
    region <- orderDedupRegions(region)
    genomicData <- .regionalLdParsePaths(genomicData)
    intersectionRows <- findIntersectionRows(
        genomicData,
        region$chrom,
        region$start,
        region$end
    )
    if (completeCoverageRequired) {
        validateSelectedRegion(
            intersectionRows$startRow,
            intersectionRows$endRow,
            region$start,
            region$end
        )
    }
    paths <- .regionalLdExtractPaths(
        ldReferenceMetaFile,
        genomicData,
        intersectionRows
    )
    list(
        intersections = list(
            startIndex = intersectionRows$startRow,
            endIndex = intersectionRows$endRow,
            LD_file_paths = paths$ldPaths,
            bimFilePaths = paths$bimPaths
        ),
        ldMetaData = genomicData,
        region = region
    )
}

#' Read a pre-computed LD matrix (.cor.xz) and its bim file, returning a
#' symmetric matrix with variants ordered by position.
#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom stats setNames
#' @noRd
# Auto-detect the variant-metadata file (.bim / .pvar / .pvar.zst).
.processLdSnpFile <- function(ldFilePath, snpFilePath) {
    if (!is.null(snpFilePath)) {
        return(snpFilePath)
    }
    candidates <- paste0(ldFilePath, c(".bim", ".pvar", ".pvar.zst"))
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) {
        stop(
            "No variant file found for: ",
            ldFilePath,
            " (tried .bim, .pvar, .pvar.zst)"
        )
    }
    found[1]
}

# Read + normalise the LD variant metadata (canonical chrom / variant id / GD).
.processLdVariants <- function(snpFilePath) {
    ldVariants <- readVariantMetadata(snpFilePath)
    isPvar <- !("gpos" %in% names(ldVariants))
    ldVariants <- ldVariants %>%
        mutate(chrom = canonChrom(chrom), variants = normalizeVariantId(id))
    if (isPvar) {
        ldVariants <- rename(ldVariants, GD = pos)
        ldVariants$GD <- ldVariants$pos <- as.integer(
            map_chr(ldVariants$variants, function(v) strsplit(v, ":")[[1]][2])
        )
    } else {
        ldVariants <- rename(ldVariants, GD = gpos)
    }
    ldVariants
}

processLdMatrix <- function(ldFilePath, snpFilePath = NULL) {
    ldFileCon <- xzfile(ldFilePath)
    ldMatrix <- scan(ldFileCon, quiet = TRUE)
    close(ldFileCon)
    ldMatrix <- matrix(ldMatrix, ncol = sqrt(length(ldMatrix)), byrow = TRUE)
    snpFilePath <- .processLdSnpFile(ldFilePath, snpFilePath)
    ldVariants <- .processLdVariants(snpFilePath)
    # Label and symmetrize the matrix.
    colnames(ldMatrix) <- rownames(ldMatrix) <- ldVariants$variants
    if (all(ldMatrix[lower.tri(ldMatrix)] == 0)) {
        ldMatrix[lower.tri(ldMatrix)] <- t(ldMatrix)[lower.tri(ldMatrix)]
    } else {
        ldMatrix[upper.tri(ldMatrix)] <- t(ldMatrix)[upper.tri(ldMatrix)]
    }
    # Order variants by genomic position.
    posOrder <- order(map_int(ldVariants$variants, function(v) {
        as.integer(strsplit(v, ":")[[1]][2])
    }))
    ldVariants <- ldVariants[posOrder, ]
    ldMatrix <- ldMatrix[ldVariants$variants, ldVariants$variants]
    list(ldMatrix = ldMatrix, ldVariants = ldVariants)
}

#' Subset an LD matrix and variant info to a genomic region, optionally further
#' restricted to specific coordinates.
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @noRd
extractLdForRegion <- function(ldMatrix, variants, region, extractCoordinates) {
    extracted <- subset(
        variants,
        chrom == region$chrom & pos >= region$start & pos <= region$end
    )

    if (!is.null(extractCoordinates)) {
        extractCoordinates <- extractCoordinates %>%
            mutate(chrom = canonChrom(chrom)) %>%
            select(chrom, pos)
        extracted <- extracted %>%
            mutate(chrom = canonChrom(chrom)) %>%
            merge(extractCoordinates, by = c("chrom", "pos"))
        keepCols <- intersect(
            c(
                "chrom",
                "variants",
                "pos",
                "GD",
                "A1",
                "A2",
                "variance",
                "allele_freq",
                "n_nomiss"
            ),
            names(extracted)
        )
        extracted <- select(extracted, all_of(keepCols))
    }

    mat <- ldMatrix[extracted$variants, extracted$variants, drop = FALSE]
    list(extractedLdMatrix = mat, extractedLdVariants = extracted)
}

# Concatenate per-block variant-id lists into one deduplicated vector, dropping
# a repeated boundary variant shared between adjacent blocks.
# @noRd
.ldMergeVariants <- function(variantList) {
    merged <- character(0)
    for (v in variantList) {
        ids <- if (is.list(v) && !is.null(v$variants)) v$variants else v
        if (length(ids) == 0) {
            next
        }
        if (length(merged) > 0 && tail(merged, 1) == ids[1]) {
            ids <- ids[-1]
        }
        merged <- c(merged, ids)
    }
    merged
}

#' Combine multiple block-level LD matrices into one, handling boundary
#' overlaps.
#' @importFrom utils tail
#' @noRd
createLdMatrix <- function(ldMatrices, variants) {
    allVariants <- .ldMergeVariants(variants)
    combined <- matrix(
        0,
        nrow = length(allVariants),
        ncol = length(allVariants),
        dimnames = list(allVariants, allVariants)
    )

    # Place each block into the combined matrix
    for (i in seq_along(ldMatrices)) {
        v <- rownames(ldMatrices[[i]])
        idx <- match(v, allVariants)
        combined[idx, idx] <- ldMatrices[[i]]
    }
    combined
}

# Dispatch to the genotype- or pre-computed-block LD loader for the source.
.loadLdDispatch <- function(
    source,
    isGeno,
    region,
    extractCoordinates,
    returnGenotype,
    nSample
) {
    if (isGeno) {
        genoPath <- resolveGenotypePathForRegion(source$metaPath, region)
        return(loadLdFromGenotype(
            genoPath,
            region,
            returnGenotype = returnGenotype,
            nSample = nSample
        ))
    }
    if (returnGenotype) {
        stop(
            "returnGenotype=TRUE requires genotype files, not ",
            "pre-computed LD matrices."
        )
    }
    loadLdFromBlocks(
        source$metaPath,
        region,
        extractCoordinates,
        nSample = nSample
    )
}

# Drop duplicate variant ids (boundary-overlap safety net).
.loadLdDedup <- function(result) {
    variantIds <- getVariantIds(result)
    if (is.null(variantIds)) {
        return(result)
    }
    dupIdx <- which(duplicated(variantIds))
    if (length(dupIdx) == 0) {
        return(result)
    }
    corr <- getCorrelation(result)
    if (!is.null(corr)) {
        corr <- corr[-dupIdx, -dupIdx, drop = FALSE]
    }
    LdData(
        correlation = corr,
        genotypeHandle = getGenotypeHandle(result),
        snpIdx = getSnpIdx(result),
        variants = getVariantInfo(result)[-dupIdx],
        blockMetadata = getBlockMetadata(result),
        nRef = getNRef(result)
    )
}

#' Load and Process Linkage Disequilibrium (LD) Matrix
#'
#' Unified entry point for loading LD data from a metadata TSV file.
#'
#' The metadata TSV must have columns: chrom, start, end, path. Two formats:
#' \itemize{
#'   \item Pre-computed LD blocks: many rows per chromosome with block
#'   boundaries
#'     in start/end and path pointing to .cor.xz files (optionally
#'     comma-separated
#'     with a .bim path).
#'   \item PLINK genotype files: one row per chromosome with start=0, end=0, and
#'     path pointing to a per-chromosome PLINK prefix (.pgen/.pvar[.zst]/.psam
#'     or
#'     .bed/.bim/.fam). LD is computed on the fly via \code{computeLd()}.
#' }
#'
#' @param ldMetaFilePath Path to the LD metadata TSV file.
#' @param region Region of interest: "chr:start-end" string or data.frame with
#'   chrom/start/end.
#' @param extractCoordinates Optional data.frame with columns "chrom" and "pos"
#'   for specific coordinates extraction (only for pre-computed LD blocks).
#' @param returnGenotype Controls what ldMatrix contains in the return value.
#'   FALSE (default): always return correlation matrix R. TRUE: return genotype
#'   matrix X (only valid for PLINK sources). "auto": return X for PLINK
#'   sources, R for pre-computed sources.
#' @param nSample Optional sample size for computing variance (=
#'   2*p*(1-p)*n/(n-1)). If NULL, ref_panel will not include variance or
#'   n_nomiss columns. Only used for PLINK genotype sources.
#'
#' @return A list with:
#' \describe{
#'   \item{ldVariants}{Character vector of variant IDs (canonical format).}
#'   \item{ldMatrix}{LD correlation matrix R (or genotype matrix X when
#'   returnGenotype is TRUE or "auto" with PLINK source).}
#'   \item{ref_panel}{Data.frame with variant metadata (chrom, pos, A2, A1,
#'   variant_id,
#'     and optionally allele_freq, variance, n_nomiss).}
#'   \item{is_genotype}{Logical: TRUE if ldMatrix contains genotype X, FALSE if
#'   correlation R.}
#'   \item{blockMetadata}{Data.frame with region/block info. For pre-computed
#'   LD: one row per block.
#'     For PLINK: a single row spanning the loaded region.}
#' }
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' loadLdMatrix(ldMetaFilePath = meta, region = "chr22:16000000-18000000")
#' @export
loadLdMatrix <- function(
    ldMetaFilePath,
    region,
    extractCoordinates = NULL,
    returnGenotype = FALSE,
    nSample = NULL
) {
    source <- resolveLdSource(ldMetaFilePath)
    isGeno <- source$type %in% c("plink2", "plink1", "vcf", "gds")
    # "auto": return X for genotype sources, R for pre-computed.
    if (identical(returnGenotype, "auto")) {
        returnGenotype <- isGeno
    }
    result <- .loadLdDispatch(
        source,
        isGeno,
        region,
        extractCoordinates,
        returnGenotype,
        nSample
    )
    .loadLdDedup(result)
}

# ---------- Internal: resolve LD source type ----------

#' @noRd
hasPlink2Files <- function(prefix) {
    file.exists(paste0(prefix, ".pgen")) &&
        (file.exists(paste0(prefix, ".pvar")) ||
            file.exists(paste0(prefix, ".pvar.zst"))) &&
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
    hasPlink2Files(path) ||
        hasPlink1Files(path) ||
        isVcfPath(path) ||
        isGdsPath(path)
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
#' This function peeks at the first row to determine the data type. The actual
#' per-chromosome path is resolved later by
#' \code{resolveGenotypePathForRegion()} at load time.
#'
#' @param path Path to a metadata TSV file with columns chrom, start, end, path.
#' @return A list with:
#'   \item{type}{"plink2", "plink1", "vcf", "gds", or "precomputed"}
#'   \item{dataPath}{Genotype path from first row (for type detection only;
#'   actual
#'     per-chromosome path is resolved at load time)}
#'   \item{metaPath}{The metadata TSV path (always set)}
#' @importFrom vroom vroom
#' @noRd
# Read + validate the first row of an LD metadata TSV (>=4 columns).
.resolveLdReadMeta <- function(path) {
    if (!file.exists(path)) {
        stop(
            "LD metadata file not found: ",
            path,
            "\n  Expected: a TSV file with columns chrom, start, end, path."
        )
    }
    meta <- as.data.frame(vroom(path, show_col_types = FALSE, n_max = 1))
    if (ncol(meta) < 4) {
        stop(
            "LD metadata file must have at least 4 columns (chrom, ",
            "start, end, path): ",
            path
        )
    }
    colnames(meta)[1:4] <- c("chrom", "start", "end", "path")
    meta
}

# Genotype source descriptor for the resolved path, or NULL if pre-computed.
.resolveLdGenotypeType <- function(resolved, path) {
    if (hasPlink2Files(resolved)) {
        return(list(type = "plink2", dataPath = resolved, metaPath = path))
    }
    if (hasPlink1Files(resolved)) {
        return(list(type = "plink1", dataPath = resolved, metaPath = path))
    }
    if (isVcfPath(resolved)) {
        return(list(type = "vcf", dataPath = resolved, metaPath = path))
    }
    if (isGdsPath(resolved)) {
        return(list(type = "gds", dataPath = resolved, metaPath = path))
    }
    NULL
}

resolveLdSource <- function(path) {
    meta <- .resolveLdReadMeta(path)
    # Strip the comma-separated bim path, then resolve relative to the meta dir.
    rawPath <- gsub(",.*$", "", meta$path[1])
    resolved <- file.path(dirname(path), rawPath)
    genoType <- .resolveLdGenotypeType(resolved, path)
    if (!is.null(genoType)) {
        return(genoType)
    }
    if (
        !is.na(meta$start) &&
            !is.na(meta$end) &&
            meta$start == 0 &&
            meta$end == 0
    ) {
        stop(
            "Metadata has start=0, end=0 but path does not resolve to ",
            "genotype files: ",
            resolved,
            "\n  The 0:0 sentinel is only valid for whole-chromosome ",
            "genotype files."
        )
    }
    list(type = "precomputed", metaPath = path)
}

#' Resolve the correct genotype path for a given region from a metadata TSV.
#' Reads the TSV, finds the row matching the query region's chromosome, and
#' returns the resolved genotype file path or prefix.
#' @importFrom vroom vroom
#' @noRd
resolveGenotypePathForRegion <- function(metaPath, region) {
    parsed <- parseRegion(region)
    meta <- as.data.frame(vroom(metaPath, show_col_types = FALSE))
    colnames(meta) <- c("chrom", "start", "end", "path")
    meta$chrom <- canonChrom(meta$chrom)
    queryChrom <- canonChrom(parsed$chrom)

    matching <- meta[meta$chrom == queryChrom, , drop = FALSE]
    if (nrow(matching) == 0) {
        stop(
            "No entry for chromosome ",
            queryChrom,
            " in metadata file: ",
            metaPath
        )
    }
    rawPath <- gsub(",.*$", "", matching$path[1])
    file.path(dirname(metaPath), rawPath)
}

# ---------- Internal: load LD from genotype files ----------

#' Load genotype data and compute LD or return genotype matrix.
#' @noRd
# --- loadLdFromGenotype helpers ---------------------------------------------

# Reference panel from variant ids + allele frequency (.afreq or dosage-derived)
# + variance when a sample size is supplied.
.loadLdGtRefPanel <- function(
    X,
    variantInfo,
    variantIds,
    genotypePath,
    nSample
) {
    refPanel <- parseVariantId(variantIds)
    refPanel$variant_id <- variantIds
    afreq <- readAfreq(genotypePath)
    if (!is.null(afreq)) {
        freqMatch <- match(variantInfo$id, afreq$id)
        nUnmatched <- sum(is.na(freqMatch))
        if (nUnmatched > 0) {
            warning(
                nUnmatched,
                " out of ",
                length(freqMatch),
                " variants have no allele frequency in .afreq file."
            )
        }
        refPanel$allele_freq <- afreq$alt_freq[freqMatch]
    } else {
        refPanel$allele_freq <- colMeans(X, na.rm = TRUE) / 2
    }
    if (!is.null(nSample)) {
        p <- refPanel$allele_freq
        refPanel$variance <- 2 * p * (1 - p) * nSample / (nSample - 1)
        refPanel$n_nomiss <- nSample
    }
    refPanel
}

# Single-block metadata spanning the loaded region.
.loadLdGtBlockMeta <- function(variantInfo, variantIds) {
    positions <- variantInfo$pos
    data.frame(
        blockId = 1L,
        chrom = as.character(variantInfo$chrom[1]),
        blockStart = min(positions),
        blockEnd = max(positions),
        size = length(variantIds),
        startIdx = 1L,
        endIdx = length(variantIds),
        stringsAsFactors = FALSE
    )
}

# Lazy-genotype LdData result (handle + region snp index, no correlation).
.loadLdGtGenotypeResult <- function(
    genotypePath,
    region,
    variantsGr,
    blockMetadata,
    X
) {
    handle <- readGenotypes(genotypePath)
    snpIdx <- .regionToSnpIdx(getSnpInfo(handle), region)
    LdData(
        correlation = NULL,
        genotypeHandle = handle,
        snpIdx = snpIdx,
        variants = variantsGr,
        blockMetadata = blockMetadata,
        nRef = as.integer(nrow(X))
    )
}

loadLdFromGenotype <- function(
    genotypePath,
    region,
    returnGenotype = FALSE,
    nSample = NULL
) {
    result <- loadGenotypeRegion(
        genotypePath,
        region = region,
        returnVariantInfo = TRUE
    )
    X <- result$X
    variantInfo <- result$variant_info
    variantIds <- normalizeVariantId(formatVariantId(
        variantInfo$chrom,
        variantInfo$pos,
        variantInfo$A2,
        variantInfo$A1
    ))
    colnames(X) <- variantIds
    refPanel <- .loadLdGtRefPanel(
        X,
        variantInfo,
        variantIds,
        genotypePath,
        nSample
    )
    blockMetadata <- .loadLdGtBlockMeta(variantInfo, variantIds)
    variantsGr <- .refPanelToGranges(refPanel)
    if (returnGenotype) {
        return(.loadLdGtGenotypeResult(
            genotypePath,
            region,
            variantsGr,
            blockMetadata,
            X
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

# ---------- LD sketch: per-variant LD matrix ----------

# Internal: build a sample-correlation LD matrix for a specified variant
# subset of an `ldSketch` `GenotypeHandle`. Shared by twasWeightsPipeline,
# fineMappingPipeline, causalInferencePipeline, and colocboostPipeline. The
# four sites differed only in their error message prefix and in whether
# variants absent from the panel raise an error or get silently dropped.
#
# Arguments:
#   ldSketch    A GenotypeHandle.
#   variantIds  Character vector of SNP IDs to extract.
#   label       Error-message prefix, e.g. ".twasLdFromSketch".
#   onMissing   "error" (default) -> any unmatched id stops the call;
#               "drop"           -> unmatched ids are silently filtered.
#               When "drop" leaves no surviving variants the function
#               returns NULL.
#
# Returns:
#   A `length(variantIds) x length(variantIds)` symmetric LD matrix with
#   rows/cols named by the (possibly filtered) `variantIds`.
#   With `onMissing = "drop"` the returned matrix carries an attribute
#   `"keptVariantIds"` so callers can recover which ids survived.
# Require a non-NULL GenotypeHandle ldSketch.
.ldFromSketchValidate <- function(ldSketch, label) {
    if (is.null(ldSketch)) {
        stop(sprintf(
            paste0(
                "%s: the SumStats/collection carries no ldSketch ",
                "(ldSketch = NULL); this step needs an LD reference."
            ),
            label
        ))
    }
    if (!methods::is(ldSketch, "GenotypeHandle")) {
        stop(sprintf("%s: ldSketch must be a GenotypeHandle.", label))
    }
}

# Match requested ids to the panel; NULL if none match. Returns kept ids/order.
.ldFromSketchMatch <- function(ldSketch, variantIds, label, onMissing) {
    snpInfo <- getSnpInfo(ldSketch)
    # Match by (chrom, pos, allele) tuple with an exact id-string fallback for
    # rsID panels; the caller's original ids and order are preserved.
    m <- matchVariants(
        variantIds,
        as.character(snpInfo$SNP),
        removeStrandAmbiguous = FALSE
    )
    nMissing <- length(variantIds) - length(m$idxA)
    if (nMissing > 0L && onMissing == "error") {
        stop(sprintf(
            "%s: %d variant id(s) not present in the LD sketch panel.",
            label,
            nMissing
        ))
    }
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    o <- order(m$idxA) # restore the caller's requested order
    list(keptIds = variantIds[m$idxA[o]], idx = m$idxB[o])
}

.ldFromSketch <- function(
    ldSketch,
    variantIds,
    label = ".ldFromSketch",
    onMissing = c("error", "drop")
) {
    .ldFromSketchValidate(ldSketch, label)
    onMissing <- match.arg(onMissing)
    matched <- .ldFromSketchMatch(ldSketch, variantIds, label, onMissing)
    if (is.null(matched)) {
        return(NULL)
    }
    geno <- .dosageMatrix(ldSketch, matched$idx, meanImpute = TRUE)
    colnames(geno) <- matched$keptIds
    ldMat <- computeLd(geno, method = "sample")
    dimnames(ldMat) <- list(matched$keptIds, matched$keptIds)
    if (onMissing == "drop") {
        attr(ldMat, "keptVariantIds") <- matched$keptIds
    }
    ldMat
}

# ---------- LD sketch: cross-pipeline LD-panel equality check ----------

# Internal: assert that two `GenotypeHandle` LD sketches describe the same
# reference panel: same variant identity (chr-agnostic CHR via canonChrom,
# exact BP/A1/A2, in the same order) and the same sampleIds. The SNP label is
# not compared, so a pure chr-prefix difference does not fail; an allele swap
# (different A1/A2) still does, since it means a different LD coding. Shared by
# causalInferencePipeline, colocPipeline,
# qtlEnrichmentPipeline, ctwasPipeline, and
# colocboostPipeline.
#
# NULL handling:
#   nullPolicy = "qtl-required" (default): a NULL qtlLd skips the check; a
#     non-NULL qtlLd with a NULL gwasLd is an error. Used by cip / coloc /
#     ctwas / enloc / qtlEnrichment.
#   nullPolicy = "lenient": a NULL on either side skips the check. Used by
#     colocboostPipeline.
#
# `label` is the human-readable name of the QTL-side input (e.g.
# "twasWeights" or "fineMappingResult"); it is woven into the error
# messages when provided. `pipelineName` prefixes every error so the
# failure source remains discoverable.
# --- .requireMatchingLdSketches helpers -------------------------------------

# Null handling: TRUE if the caller should return early (qtl NULL, or gwas NULL
# under a lenient policy); errors when a non-NULL qtl faces a NULL gwas sketch.
.ldSketchNullGuard <- function(qtlLd, gwasLd, pipelineName, label, nullPolicy) {
    if (is.null(qtlLd)) {
        return(TRUE)
    }
    if (is.null(gwasLd)) {
        if (nullPolicy == "lenient") {
            return(TRUE)
        }
        stop(
            pipelineName,
            ": ",
            if (!is.null(label)) {
                sprintf("ldSketch on `%s` is non-NULL ", label)
            } else {
                "qtl ldSketch is non-NULL "
            },
            "but the GWAS ldSketch is NULL."
        )
    }
    FALSE
}

# Both must be GenotypeHandles with matching panel size.
.ldSketchCheckShape <- function(qtlLd, gwasLd, pipelineName, between) {
    if (
        !methods::is(qtlLd, "GenotypeHandle") ||
            !methods::is(gwasLd, "GenotypeHandle")
    ) {
        stop(
            pipelineName,
            ": ldSketch slots",
            between,
            " must both be GenotypeHandle objects ",
            "for the cross-pipeline LD reference check."
        )
    }
    qSnp <- getSnpInfo(qtlLd)
    gSnp <- getSnpInfo(gwasLd)
    if (nrow(qSnp) != nrow(gSnp)) {
        stop(
            pipelineName,
            ": ldSketch panels differ in size (",
            nrow(qSnp),
            " vs ",
            nrow(gSnp),
            " variants)",
            between,
            "; the two ldSketch GenotypeHandles must match exactly."
        )
    }
}

# Panels must agree on CHR/BP/A1/A2 columns and on the sample set.
.ldSketchCheckContent <- function(qtlLd, gwasLd, pipelineName, between) {
    qSnp <- getSnpInfo(qtlLd)
    gSnp <- getSnpInfo(gwasLd)
    if (!identical(canonChrom(qSnp$CHR), canonChrom(gSnp$CHR))) {
        stop(
            pipelineName,
            ": ldSketch panels differ in column CHR",
            between,
            "; use the same ldSketch on both."
        )
    }
    for (col in c("BP", "A1", "A2")) {
        if (!identical(as.character(qSnp[[col]]), as.character(gSnp[[col]]))) {
            stop(
                pipelineName,
                ": ldSketch panels differ in column ",
                col,
                between,
                "; use the same ldSketch on both."
            )
        }
    }
    if (!identical(getSampleIds(qtlLd), getSampleIds(gwasLd))) {
        stop(
            pipelineName,
            ": ldSketch panels have different sample sets",
            between,
            "; use the same ldSketch on both."
        )
    }
}

.requireMatchingLdSketches <- function(
    qtlLd,
    gwasLd,
    pipelineName,
    label = NULL,
    nullPolicy = c("qtl-required", "lenient")
) {
    nullPolicy <- match.arg(nullPolicy)
    if (.ldSketchNullGuard(qtlLd, gwasLd, pipelineName, label, nullPolicy)) {
        return(invisible(NULL))
    }
    between <- if (!is.null(label)) {
        sprintf(" between `%s` and gwas inputs", label)
    } else {
        ""
    }
    .ldSketchCheckShape(qtlLd, gwasLd, pipelineName, between)
    .ldSketchCheckContent(qtlLd, gwasLd, pipelineName, between)
    invisible(NULL)
}

# ---------- LD sketch: genotype loading ----------

#' HWE-based standardization of a genotype matrix
#'
#' Centers by 2*alleleFreq, scales by sqrt(2*alleleFreq*(1-alleleFreq)). Assumes
#' monomorphic variants have already been removed.
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
#' and removes monomorphic variants. Returns the raw genotype matrix and
#' metadata, which callers can use to derive either a correlation matrix R (for
#' summary-based weight training or fine-mapping) or an SVD (for TWAS z-score
#' computation).
#'
#' @param ldMetaFilePath Path to the LD metadata TSV file.
#' @param region Region of interest: "chr:start-end" string or data.frame with
#'   chrom/start/end.
#' @param nSample Optional original panel sample size for computing variance (=
#'   2*p*(1-p)*n/(n-1)). Passed through to \code{loadLdMatrix()}.
#'
#' @return An \code{LdData} S4 object with monomorphic variants removed.
#'   Consumers should use S4 accessors: \code{getGenotypes()},
#'   \code{getRefPanel()}, \code{getVariantIds()}. The number of sketch samples
#'   is \code{nrow(getGenotypes(result))}.
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' loadLdSketch(ldMetaFilePath = meta, region = "chr22:16000000-18000000")
#' @export
loadLdSketch <- function(ldMetaFilePath, region, nSample = NULL) {
    result <- loadLdMatrix(
        ldMetaFilePath,
        region,
        returnGenotype = TRUE,
        nSample = nSample
    )
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
        nRef = getNRef(result)
    )
}

# ---------- Internal: load LD from pre-computed blocks ----------

#' Load pre-computed LD from block-based metadata files.
#' @importFrom purrr map_chr map_int map_dbl
#' @noRd
# --- loadLdFromBlocks helpers -----------------------------------------------

# Load + region-extract each LD block; track each block's chromosome.
.loadLdBlocksLoop <- function(
    ldFilePaths,
    bimFilePaths,
    intersectedLdFiles,
    extractCoordinates
) {
    matrices <- list()
    variants <- list()
    blockChroms <- character(length(ldFilePaths))
    for (j in seq_along(ldFilePaths)) {
        proc <- processLdMatrix(ldFilePaths[j], bimFilePaths[j])
        extracted <- extractLdForRegion(
            ldMatrix = proc$ldMatrix,
            variants = proc$ldVariants,
            region = intersectedLdFiles$region,
            extractCoordinates = extractCoordinates
        )
        matrices[[j]] <- extracted$extractedLdMatrix
        variants[[j]] <- extracted$extractedLdVariants
        blockChroms[j] <- if (nrow(variants[[j]]) > 0) {
            as.character(variants[[j]]$chrom[1])
        } else {
            as.character(intersectedLdFiles$region$chrom)
        }
    }
    list(matrices = matrices, variants = variants, blockChroms = blockChroms)
}

# Drop blocks with no variants in the region (error if none remain).
.loadLdFilterEmpty <- function(blocks, ldFilePaths) {
    nonEmpty <- map_lgl(blocks$variants, function(v) nrow(v) > 0)
    if (!any(nonEmpty)) {
        stop("No variants found in any LD block for the specified region.")
    }
    if (any(!nonEmpty)) {
        message(
            "Removing ",
            sum(!nonEmpty),
            " empty LD block(s) with no variants in the region."
        )
    }
    list(
        matrices = blocks$matrices[nonEmpty],
        variants = blocks$variants[nonEmpty],
        blockChroms = blocks$blockChroms[nonEmpty],
        ldFilePaths = ldFilePaths[nonEmpty]
    )
}

# Per-block metadata (id, chrom, span, size, index range in the merged matrix).
.loadLdBlockMetadata <- function(
    variants,
    ldFilePaths,
    blockChroms,
    ldVariants
) {
    blockVariants <- lapply(variants, function(v) v$variants)
    blockPositions <- lapply(variants, function(v) v$pos)
    data.frame(
        blockId = seq_along(ldFilePaths),
        chrom = blockChroms,
        blockStart = map_dbl(blockPositions, min),
        blockEnd = map_dbl(blockPositions, max),
        size = map_int(blockVariants, length),
        startIdx = map_int(blockVariants, function(v) {
            min(match(v, ldVariants))
        }),
        endIdx = map_int(blockVariants, function(v) max(match(v, ldVariants))),
        stringsAsFactors = FALSE
    )
}

# Build the reference panel: variant ids + merged per-variant annotations,
# deriving variance from nSample + allele_freq when it is otherwise absent.
.loadLdRefPanel <- function(ldMatrix, extractedLdVariantsList, nSample) {
    refPanel <- parseVariantId(rownames(ldMatrix))
    mergedVariantList <- do.call(rbind, extractedLdVariantsList)
    ids <- rownames(ldMatrix)
    refPanel$variant_id <- ids
    for (col in c("allele_freq", "variance", "n_nomiss")) {
        if (col %in% colnames(mergedVariantList)) {
            refPanel[[col]] <- mergedVariantList[[col]][
                match(ids, mergedVariantList$variants)
            ]
        }
    }
    needVar <- !"variance" %in% colnames(refPanel) ||
        all(is.na(refPanel$variance))
    if (!is.null(nSample) && needVar && "allele_freq" %in% colnames(refPanel)) {
        p <- refPanel$allele_freq
        refPanel$variance <- 2 * p * (1 - p) * nSample / (nSample - 1)
        refPanel$n_nomiss <- nSample
    }
    refPanel
}

loadLdFromBlocks <- function(
    ldMetaFilePath,
    region,
    extractCoordinates = NULL,
    nSample = NULL
) {
    intersectedLdFiles <- getRegionalLdMeta(ldMetaFilePath, region)
    ldFilePaths <- intersectedLdFiles$intersections$LD_file_paths
    bimFilePaths <- intersectedLdFiles$intersections$bimFilePaths
    blocks <- .loadLdBlocksLoop(
        ldFilePaths,
        bimFilePaths,
        intersectedLdFiles,
        extractCoordinates
    )
    filtered <- .loadLdFilterEmpty(blocks, ldFilePaths)
    ldMatrix <- createLdMatrix(
        ldMatrices = filtered$matrices,
        variants = filtered$variants
    )
    ldVariants <- rownames(ldMatrix)
    blockMetadata <- .loadLdBlockMetadata(
        filtered$variants,
        filtered$ldFilePaths,
        filtered$blockChroms,
        ldVariants
    )
    refPanel <- .loadLdRefPanel(ldMatrix, filtered$variants, nSample)
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
#' Auto-detects the reference type (PLINK2, PLINK1, or pre-computed LD
#' metadata).
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
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' filterVariantsByLdReference(
#'   variantIds = c("chr22:16050000:A:G", "chr22:17000000:C:T"),
#'   ldReferenceMetaFile = meta)
#' @export
filterVariantsByLdReference <- function(
    variantIds,
    ldReferenceMetaFile,
    keepIndel = TRUE
) {
    variantsDf <- parseVariantId(variantIds)

    # Derive region to scope the reference lookup
    regionDf <- variantsDf %>%
        group_by(chrom) %>%
        summarise(start = min(pos), end = max(pos))

    # Use shared helper -- no genotype loading
    refInfo <- getRefVariantInfo(ldReferenceMetaFile, regionDf)
    refChrom <- canonChrom(refInfo$chrom)
    refKey <- paste0(refChrom, ":", refInfo$pos)

    variantKey <- paste0(variantsDf$chrom, ":", variantsDf$pos)
    keepIndices <- which(variantKey %in% refKey)

    if (!keepIndel) {
        snpIdx <- which(isSnpAlleles(variantsDf$A1, variantsDf$A2))
        keepIndices <- intersect(keepIndices, snpIdx)
    }

    nDropped <- length(variantIds) - length(keepIndices)
    if (nDropped > 0) {
        message(
            nDropped,
            " out of ",
            length(variantIds),
            " total variants dropped due to absence on the reference LD panel."
        )
    }

    list(data = variantIds[keepIndices], idx = keepIndices)
}

#' Partition LD Matrix into Block-Specific Matrices
#'
#' This function takes the output from loadLdMatrix and partitions the combined
#' LD matrix into a list of smaller matrices based on the block_indices, making
#' it easier to work with large LD matrices that span multiple blocks.
#'
#' @param ldData An \code{LdData} S4 object as returned by
#'   \code{loadLdMatrix()}.
#' @param mergeSmallBlocks Logical, whether to merge blocks smaller than
#'   minMergedBlockSize (default: TRUE).
#' @param minMergedBlockSize Integer, minimum number of variants for a block
#'   after merging (default: 500).
#' @param maxMergedBlockSize Integer, maximum number of variants in a block
#'   after merging (default: 10000).
#'
#' @return returns a list containing:
#' \describe{
#' \item{ldMatrices}{A list of matrices, each representing LD for a specific
#' block.}
#' \item{variantIndices}{A data frame that maps variant IDs to their
#' corresponding block.}
#' \item{blockMetadata}{Information about each block including size,
#' chromosome, start and end positions.}
#' }
#' @noRd
# --- partitionLdMatrix helpers ----------------------------------------------

# Reject empty/NULL matrices; align row/col names to the variant ids.
.partitionValidateMatrix <- function(combinedMatrix, variantIds) {
    if (
        is.null(combinedMatrix) ||
            nrow(combinedMatrix) == 0 ||
            ncol(combinedMatrix) == 0
    ) {
        stop("Empty or NULL LD matrix provided.")
    }
    if (
        is.null(rownames(combinedMatrix)) ||
            is.null(colnames(combinedMatrix)) ||
            !identical(rownames(combinedMatrix), variantIds) ||
            !identical(colnames(combinedMatrix), variantIds)
    ) {
        rownames(combinedMatrix) <- variantIds
        colnames(combinedMatrix) <- variantIds
    }
    combinedMatrix
}

# Drop blocks with invalid/out-of-range indices; renumber the survivors.
.partitionFilterBlocks <- function(blockMetadata, nVariants) {
    validBlocks <- map_lgl(seq_len(nrow(blockMetadata)), function(i) {
        s <- blockMetadata$startIdx[i]
        e <- blockMetadata$endIdx[i]
        sz <- blockMetadata$size[i]
        !is.na(s) &&
            !is.na(e) &&
            is.finite(s) &&
            is.finite(e) &&
            sz > 0 &&
            s >= 1 &&
            e >= s &&
            e <= nVariants
    })
    if (!any(validBlocks)) {
        stop(
            "No valid LD blocks found. All block indices are out of ",
            "range or empty."
        )
    }
    if (any(!validBlocks)) {
        message(
            "Removing ",
            sum(!validBlocks),
            " LD block(s) with invalid or out-of-range indices."
        )
        blockMetadata <- blockMetadata[validBlocks, , drop = FALSE]
        blockMetadata$blockId <- seq_len(nrow(blockMetadata))
    }
    blockMetadata
}

partitionLdMatrix <- function(
    ldData,
    mergeSmallBlocks = TRUE,
    minMergedBlockSize = 500,
    maxMergedBlockSize = 10000
) {
    if (!is(ldData, "LdData")) {
        stop("ldData must be an LdData object")
    }
    combinedMatrix <- getCorrelation(ldData)
    blockMetadata <- getBlockMetadata(ldData)
    if (is(blockMetadata, "LdBlocks")) {
        blockMetadata <- as.data.frame(getBlocks(blockMetadata))
    }
    variantIds <- getVariantIds(ldData)
    combinedMatrix <- .partitionValidateMatrix(combinedMatrix, variantIds)
    blockMetadata <- .partitionFilterBlocks(blockMetadata, length(variantIds))
    # Validate the block structure of the matrix (skip if only one block).
    if (nrow(blockMetadata) > 1) {
        validateBlockStructure(combinedMatrix, blockMetadata, variantIds)
    }
    if (
        mergeSmallBlocks &&
            any(blockMetadata$size < minMergedBlockSize) &&
            nrow(blockMetadata) > 1
    ) {
        blockMetadata <- mergeBlocks(
            blockMetadata,
            minMergedBlockSize,
            maxMergedBlockSize
        )
    }
    extractBlockMatrices(combinedMatrix, blockMetadata, variantIds)
}

#' Validate that cross-block entries are zero (excluding boundary variants).
#' @noRd
validateBlockStructure <- function(matrix, blockMetadata, variantIds) {
    msgs <- character(0)
    n <- length(variantIds)

    for (i in 1:(nrow(blockMetadata) - 1)) {
        for (j in (i + 1):nrow(blockMetadata)) {
            si <- blockMetadata$startIdx[i]
            ei <- blockMetadata$endIdx[i]
            sj <- blockMetadata$startIdx[j]
            ej <- blockMetadata$endIdx[j]
            if (si > n || ei > n || sj > n || ej > n) {
                msgs <- c(
                    msgs,
                    paste("Block indices out of range for blocks", i, "and", j)
                )
                next
            }
            # Exclude boundary variants (potential overlaps)
            vi <- variantIds[si:(ei - 1)]
            vj <- variantIds[(sj + 1):ej]
            if (length(vi) > 0 && length(vj) > 0) {
                maxVal <- max(abs(matrix[vi, vj, drop = FALSE]))
                if (maxVal > 1e-10) {
                    msgs <- c(
                        msgs,
                        paste(
                            "Non-zero correlation between blocks",
                            i,
                            "and",
                            j,
                            "- max:",
                            maxVal
                        )
                    )
                }
            }
        }
    }
    if (length(msgs) > 0) {
        stop(
            "Matrix lacks expected block structure:\n",
            paste(msgs, collapse = "\n")
        )
    }
}

#' @noRd
canMerge <- function(block1, block2, maxSize) {
    block1$chrom == block2$chrom && (block1$size + block2$size) <= maxSize
}

#' @noRd
mergeTwoBlocks <- function(blockMetadata, idx1, idx2) {
    if (idx1 > idx2) {
        tmp <- idx1
        idx1 <- idx2
        idx2 <- tmp
    }
    result <- blockMetadata
    result$endIdx[idx1] <- blockMetadata$endIdx[idx2]
    result$size[idx1] <- blockMetadata$size[idx1] + blockMetadata$size[idx2]
    result <- result[-idx2, ]
    result$blockId <- seq_len(nrow(result))
    result
}

#' Find blocks below minSize and identify the best neighbor to merge with.
#' @noRd
findMergeCandidates <- function(blockMetadata, minSize, maxSize) {
    candidates <- data.frame(
        block_idx = integer(),
        merge_with = integer(),
        stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(blockMetadata))) {
        if (blockMetadata$size[i] >= minSize) {
            next
        }
        prevOk <- i > 1 &&
            canMerge(blockMetadata[i, ], blockMetadata[i - 1, ], maxSize)
        nextOk <- i < nrow(blockMetadata) &&
            canMerge(blockMetadata[i, ], blockMetadata[i + 1, ], maxSize)
        mergeWith <- if (prevOk && nextOk) {
            if (blockMetadata$size[i - 1] <= blockMetadata$size[i + 1]) {
                i - 1
            } else {
                i + 1
            }
        } else if (prevOk) {
            i - 1
        } else if (nextOk) {
            i + 1
        } else {
            next
        }
        candidates <- rbind(
            candidates,
            data.frame(block_idx = i, merge_with = mergeWith)
        )
    }
    candidates
}

#' Iteratively merge blocks below minSize with their smallest neighbor.
#' @noRd
mergeBlocks <- function(blockMetadata, minSize, maxSize) {
    if (nrow(blockMetadata) <= 1) {
        return(blockMetadata)
    }
    repeat {
        candidates <- findMergeCandidates(blockMetadata, minSize, maxSize)
        if (nrow(candidates) == 0) {
            break
        }
        blockMetadata <- mergeTwoBlocks(
            blockMetadata,
            candidates$block_idx[1],
            candidates$merge_with[1]
        )
    }
    blockMetadata
}

# Helper function to extract block matrices
# Extract one block's submatrix + mapping (NULL to skip empty/OOB blocks).
.extractOneBlock <- function(matrix, variantIds, startIdx, endIdx, i) {
    if (endIdx < startIdx) {
        return(NULL)
    }
    if (startIdx > length(variantIds) || endIdx > length(variantIds)) {
        warning(
            "Block ",
            i,
            " has indices outside the range of variantIds. Skipping."
        )
        return(NULL)
    }
    blockVariants <- variantIds[startIdx:endIdx]
    list(
        matrix = matrix[blockVariants, blockVariants, drop = FALSE],
        mapping = data.frame(
            variant_id = blockVariants,
            blockId = i,
            stringsAsFactors = FALSE
        )
    )
}

extractBlockMatrices <- function(matrix, blockMetadata, variantIds) {
    ldMatrices <- list()
    variantMapping <- data.frame(
        variant_id = character(),
        blockId = integer(),
        stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(blockMetadata))) {
        block <- .extractOneBlock(
            matrix,
            variantIds,
            blockMetadata$startIdx[i],
            blockMetadata$endIdx[i],
            i
        )
        if (is.null(block)) {
            next
        }
        ldMatrices[[i]] <- block$matrix
        variantMapping <- rbind(variantMapping, block$mapping)
    }
    list(
        ldMatrices = ldMatrices,
        variantIndices = variantMapping,
        blockMetadata = blockMetadata
    )
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
#' @param shrinkage Shrinkage parameter for \code{method = "shrink"}. Default:
#'   \code{0.01}.
#'
#' @return A list with components:
#' \describe{
#'   \item{R}{The (possibly repaired) LD matrix.}
#'   \item{isPd}{Logical: is the matrix positive definite?}
#'   \item{isPsd}{Logical: is the matrix positive semidefinite (within rTol)?}
#'   \item{minEigenvalue}{Smallest eigenvalue of the original matrix.}
#'   \item{nNegative}{Number of negative eigenvalues (below -rTol).}
#'   \item{conditionNumber}{Ratio of largest to smallest positive eigenvalue
#'     (\code{Inf} if any eigenvalue is zero).}
#'   \item{methodApplied}{Character: \code{"none"}, \code{"shrink"}, or
#'     \code{"eigenfix"}.}
#' }
#'
#' @examples
#' # A well-conditioned matrix
#' R_good <- diag(5)
#' checkLd(R_good)$isPd  # TRUE
#'
#' # A matrix with negative eigenvalues
#' R_bad <- matrix(0.9, 3, 3); diag(R_bad) <- 1
#' R_bad[1, 3] <- R_bad[3, 1] <- -0.5
#' checkLd(R_bad)$isPsd  # FALSE
#' R_fixed <- checkLd(R_bad, method = "eigenfix")$R
#' checkLd(R_fixed)$isPsd  # TRUE
#'
#' @export
checkLd <- function(
    R,
    method = c("check", "shrink", "eigenfix"),
    rTol = 1e-8,
    shrinkage = 0.01
) {
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
        isPd = isPd,
        isPsd = isPsd,
        minEigenvalue = minEval,
        nNegative = nNeg,
        conditionNumber = condNum,
        methodApplied = methodApplied
    )
}

# hclust (single-linkage) LD pruning: keep one representative per |cor| cluster.
.ldPruneHclust <- function(X, corThres, verbose) {
    p <- ncol(X)
    if (requireNamespace("Rfast", quietly = TRUE)) {
        cor.X <- Rfast::cora(X, large = TRUE)
    } else {
        cor.X <- cor(X)
    }
    Sigma.distance <- as.dist(1 - abs(cor.X))
    fit <- hclust(Sigma.distance, method = "single")
    clusters <- cutree(fit, h = 1 - corThres)
    ind.delete <- NULL
    for (ig in unique(clusters)) {
        temp.group <- which(clusters == ig)
        if (length(temp.group) > 1) {
            ind.delete <- c(ind.delete, temp.group[-1])
        }
    }
    ind.delete <- unique(ind.delete)
    X.new <- X
    filter.id <- seq_len(p)
    if (length(ind.delete) > 0) {
        X.new <- as.matrix(X[, -ind.delete])
        filter.id <- filter.id[-ind.delete]
        if (verbose) {
            message(
                "ldPruneByCorrelation: pruned ",
                length(ind.delete),
                " of ",
                p,
                " columns at |cor| > ",
                corThres
            )
        }
    } else if (verbose) {
        message(
            "ldPruneByCorrelation: no columns pruned at |cor| > ",
            corThres
        )
    }
    if (ncol(X.new) == 1) {
        colnames(X.new) <- colnames(X)[-ind.delete]
    }
    list(X.new = X.new, filter.id = filter.id)
}

#' Prune columns by pairwise correlation (LD-style prune)
#'
#' Performs LD pruning using one of two backends. The default \code{"hclust"}
#' backend computes the full correlation matrix, builds a single-linkage
#' hierarchical clustering on the distance (1 - |cor|), and keeps one
#' representative column per cluster. The \code{"snprelate"} backend delegates
#' to \code{SNPRelate::snpgdsLDpruning}, which performs a sliding-window greedy
#' prune directly on a temporary GDS file.
#'
#' @param X Numeric matrix. Columns are the variables to prune (typically SNP
#'   genotype dosages); rows are observations.
#' @param corThres Numeric in (0, 1). Absolute correlation threshold. Columns
#'   whose pairwise |cor| exceeds this are grouped; one survivor is kept per
#'   group. Default 0.8.
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
ldPruneByCorrelation <- function(
    X,
    corThres = 0.8,
    backend = c("hclust", "snprelate"),
    verbose = FALSE
) {
    backend <- match.arg(backend)
    if (backend == "snprelate") {
        return(.ldPruneSnprelate(X, corThres = corThres, verbose = verbose))
    }
    .ldPruneHclust(X, corThres, verbose)
}

#' SNPRelate-based LD pruning helper
#' @noRd
#' SNPRelate-based LD pruning helper
#' @noRd
.ldPruneSnprelateDeps <- function() {
    if (
        !requireNamespace("SNPRelate", quietly = TRUE) ||
            !requireNamespace("gdsfmt", quietly = TRUE)
    ) {
        # nocov start
        stop(
            "Packages 'SNPRelate' and 'gdsfmt' are required for ",
            "backend='snprelate'."
        )
        # nocov end
    }
}

# Write X (rounded to integer genotype codes) to a temporary GDS for SNPRelate.
.ldPruneSnprelateCreateGds <- function(tmpGds, X, snpNames, p) {
    genoInt <- round(X)
    storage.mode(genoInt) <- "integer"
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
}

.ldPruneSnprelate <- function(X, corThres, verbose) {
    .ldPruneSnprelateDeps()
    p <- ncol(X)
    snpNames <- colnames(X) %||% paste0("snp", seq_len(p))
    tmpGds <- tempfile(fileext = ".gds")
    on.exit(unlink(tmpGds), add = TRUE)
    .ldPruneSnprelateCreateGds(tmpGds, X, snpNames, p)
    gds <- SNPRelate::snpgdsOpen(tmpGds, allow.duplicate = TRUE)
    on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)
    keepList <- SNPRelate::snpgdsLDpruning(
        gds,
        method = "corr",
        ld.threshold = corThres,
        verbose = verbose
    )
    keepIds <- sort(unlist(keepList, use.names = FALSE))
    X.new <- X[, keepIds, drop = FALSE]
    if (verbose) {
        message(
            "ldPruneByCorrelation (snprelate): kept ",
            length(keepIds),
            " of ",
            p,
            " columns at |cor| > ",
            corThres
        )
    }
    list(X.new = X.new, filter.id = keepIds)
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
#' @param problematicCols Character vector of column names in \code{X} that are
#'   candidates for removal. If empty, \code{X} is returned unchanged.
#' @param strategy One of \code{"correlation"} (remove the column with the
#'   largest sum of absolute pairwise correlations among the candidates; when
#'   only two candidates, one is picked at random), \code{"variance"} (remove
#'   the lowest-variance candidate), or \code{"response_correlation"} (remove
#'   the candidate whose correlation with \code{response} has the smallest
#'   magnitude).
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
# Correlation strategy: drop the most-connected column (random tie-break at 2).
.dropCollinearPickCor <- function(X, problematicCols, verbose) {
    corMatrix <- abs(cor(X[, problematicCols, drop = FALSE]))
    diag(corMatrix) <- 0
    if (length(problematicCols) == 2) {
        colToRemove <- sample(problematicCols, 1)
        if (verbose) {
            message(
                "dropCollinearColumns: two candidates, randomly removing ",
                colToRemove
            )
        }
        return(colToRemove)
    }
    colToRemove <- problematicCols[which.max(colSums(corMatrix))]
    if (verbose) {
        message(
            "dropCollinearColumns: highest sum |cor| -> removing ",
            colToRemove
        )
    }
    colToRemove
}

# Choose which of >=2 collinear columns to drop, per the requested strategy.
.dropCollinearPick <- function(
    X,
    problematicCols,
    strategy,
    response,
    verbose
) {
    if (strategy == "variance") {
        variances <- apply(X[, problematicCols, drop = FALSE], 2, var)
        colToRemove <- problematicCols[which.min(variances)]
        if (verbose) {
            message(
                "dropCollinearColumns: smallest variance -> removing ",
                colToRemove
            )
        }
        return(colToRemove)
    }
    if (strategy == "correlation") {
        return(.dropCollinearPickCor(X, problematicCols, verbose))
    }
    if (is.null(response)) {
        stop(
            "response must be supplied for strategy = ",
            "'response_correlation'"
        )
    }
    corWithResponse <- apply(
        X[, problematicCols, drop = FALSE],
        2,
        function(col) cor(col, response)
    )
    colToRemove <- problematicCols[which.min(abs(corWithResponse))]
    if (verbose) {
        message(
            "dropCollinearColumns: smallest |cor| with response -> removing ",
            colToRemove
        )
    }
    colToRemove
}

dropCollinearColumns <- function(
    X,
    problematicCols,
    strategy = c("correlation", "variance", "response_correlation"),
    response = NULL,
    verbose = FALSE
) {
    strategy <- match.arg(strategy)
    if (length(problematicCols) == 0) {
        return(X)
    }
    if (length(problematicCols) == 1) {
        colToRemove <- problematicCols[1]
        if (verbose) {
            message(
                "dropCollinearColumns: removing single column ",
                colToRemove
            )
        }
        return(X[, !(colnames(X) %in% colToRemove), drop = FALSE])
    }
    colToRemove <- .dropCollinearPick(
        X,
        problematicCols,
        strategy,
        response,
        verbose
    )
    X[, !(colnames(X) %in% colToRemove), drop = FALSE]
}

# Design matrix [1 | X | C] with the intercept + X columns named.
# @noRd
.ldBuildDesign <- function(X, C) {
    XD <- cbind(1, X, C)
    colnames(XD)[seq_len(ncol(X) + 1L)] <- c("Intercept", colnames(X))
    XD
}

# --- enforceDesignFullRank helpers ------------------------------------------

# QR-pivot columns of the design that are collinear (and present in X).
.edfrProblematicColnames <- function(Xdesign, X) {
    qrd <- qr(Xdesign)
    if (qrd$rank >= ncol(Xdesign)) {
        return(character(0))
    }
    cols <- qrd$pivot[(qrd$rank + 1L):ncol(Xdesign)]
    nms <- colnames(Xdesign)[cols]
    nms[nms %in% colnames(X)]
}

# Fast pre-check: would batch-removing the flagged columns restore full rank?
# Returns TRUE to skip the (slow) iterative path in favour of the fallback.
.edfrCheckBatch <- function(X, C, Xdesign, matrixRank, verbose) {
    if (matrixRank >= ncol(Xdesign)) {
        return(FALSE)
    }
    problematicColnames <- .edfrProblematicColnames(Xdesign, X)
    if (length(problematicColnames) == 0) {
        return(FALSE)
    }
    Xtemp <- X[, !(colnames(X) %in% problematicColnames), drop = FALSE]
    tempDesign <- .ldBuildDesign(Xtemp, C)
    if (qr(tempDesign)$rank == ncol(tempDesign)) {
        if (verbose) {
            message(
                "enforceDesignFullRank: full rank after batch-removing ",
                length(problematicColnames),
                " column(s)"
            )
        }
        return(FALSE)
    }
    if (verbose) {
        message(
            "enforceDesignFullRank: batch removal insufficient, ",
            "skipping to correlation-pruning fallback"
        )
    }
    TRUE
}

# Iteratively drop collinear columns until the design is full rank.
.edfrIterate <- function(X, C, strategy, response, maxIterations, verbose) {
    iteration <- 0L
    Xdesign <- .ldBuildDesign(X, C)
    matrixRank <- qr(Xdesign)$rank
    while (matrixRank < ncol(Xdesign) && iteration < maxIterations) {
        problematicColnames <- .edfrProblematicColnames(Xdesign, X)
        if (length(problematicColnames) == 0) {
            break
        }
        X <- dropCollinearColumns(
            X,
            problematicColnames,
            strategy = strategy,
            response = response,
            verbose = verbose
        )
        Xdesign <- .ldBuildDesign(X, C)
        matrixRank <- qr(Xdesign)$rank
        iteration <- iteration + 1L
        if (verbose) {
            message(
                "enforceDesignFullRank: iter ",
                iteration,
                " rank ",
                matrixRank,
                " / ",
                ncol(Xdesign)
            )
        }
    }
    if (iteration == maxIterations) {
        warning(
            "enforceDesignFullRank: maxIterations reached; design may ",
            "still be rank-deficient"
        )
    }
    X
}

# Correlation-threshold pruning fallback when the design is still deficient.
.edfrCorrelationFallback <- function(X, C, corrThresholds, verbose) {
    Xdesign <- .ldBuildDesign(X, C)
    matrixRank <- qr(Xdesign)$rank
    if (matrixRank >= ncol(Xdesign)) {
        return(X)
    }
    if (verbose) {
        message("enforceDesignFullRank: applying ldPruneByCorrelation fallback")
    }
    for (threshold in corrThresholds) {
        filterResult <- ldPruneByCorrelation(
            X,
            corThres = threshold,
            verbose = verbose
        )
        X <- filterResult$X.new
        Xdesign <- .ldBuildDesign(X, C)
        matrixRank <- qr(Xdesign)$rank
        if (verbose) {
            message(
                "enforceDesignFullRank: threshold ",
                threshold,
                " -> rank ",
                matrixRank,
                " / ",
                ncol(Xdesign)
            )
        }
        if (matrixRank == ncol(Xdesign)) break
    }
    X
}

#' Iteratively enforce full column rank on a design matrix
#'
#' Given a candidate predictor matrix \code{X} and an optional unnamed covariate
#' matrix \code{C}, builds the design \code{[1, X, C]} and removes
#' rank-deficient columns from \code{X} until the design has full column rank.
#' Rank-deficient columns are identified via the pivot of \code{qr([1, X, C])}.
#' On each iteration, one problematic column is dropped using
#' \code{dropCollinearColumns}. If iterative pruning does not achieve full rank,
#' falls back to \code{\link{ldPruneByCorrelation}} at a descending sequence of
#' correlation thresholds.
#'
#' @param X Numeric matrix with column names (the predictors subject to
#'   pruning).
#' @param C Numeric matrix of covariates (can be unnamed) that will be kept.
#'   Pass \code{NULL} or a zero-column matrix when there are no covariates.
#' @param strategy Passed through to \code{dropCollinearColumns}.
#' @param response Passed through to \code{dropCollinearColumns} when
#'   \code{strategy = "response_correlation"}.
#' @param maxIterations Integer. Hard cap on the iterative-prune loop. Default
#'   300.
#' @param corrThresholds Numeric vector of |cor| thresholds used for the
#'   \code{\link{ldPruneByCorrelation}} fallback, tried in order. Default
#'   \code{seq(0.75, 0.5, by = -0.05)}.
#' @param verbose Logical. If TRUE, print per-iteration progress. Default FALSE.
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
enforceDesignFullRank <- function(
    X,
    C,
    strategy = c("correlation", "variance", "response_correlation"),
    response = NULL,
    maxIterations = 300L,
    corrThresholds = seq(0.75, 0.5, by = -0.05),
    verbose = FALSE
) {
    strategy <- match.arg(strategy)
    originalColnames <- colnames(X)
    initialNcol <- ncol(X)
    Xdesign <- .ldBuildDesign(X, C)
    matrixRank <- qr(Xdesign)$rank
    if (verbose) {
        message(
            "enforceDesignFullRank: initial rank ",
            matrixRank,
            " / ",
            ncol(Xdesign)
        )
    }
    skipIterative <- .edfrCheckBatch(X, C, Xdesign, matrixRank, verbose)
    if (!skipIterative) {
        X <- .edfrIterate(X, C, strategy, response, maxIterations, verbose)
    }
    X <- .edfrCorrelationFallback(X, C, corrThresholds, verbose)
    if (ncol(X) == 1L && initialNcol == 1L) {
        colnames(X) <- originalColnames
    }
    X
}

# Require the bigsnpr/bigstatsr packages used for score-based LD clumping.
.ldClumpCheckDeps <- function() {
    if (!requireNamespace("bigsnpr", quietly = TRUE)) {
        # nocov start
        stop(
            "Package 'bigsnpr' is required. Install from CRAN: ",
            "install.packages('bigsnpr')"
        )
        # nocov end
    }
    if (!requireNamespace("bigstatsr", quietly = TRUE)) {
        # nocov start
        stop(
            "Package 'bigstatsr' is required. Install from CRAN: ",
            "install.packages('bigstatsr')"
        )
        # nocov end
    }
}

# Validate the clumping inputs (dimensions of score/chr/pos vs X).
.ldClumpValidate <- function(X, score, chr, pos) {
    if (ncol(X) < 1L) {
        stop("ldClumpByScore: X must have at least one column")
    }
    if (!is.null(score) && length(score) != ncol(X)) {
        stop("ldClumpByScore: length(score) must equal ncol(X)")
    }
    if (length(chr) != ncol(X) || length(pos) != ncol(X)) {
        stop("ldClumpByScore: chr and pos must have length equal to ncol(X)")
    }
}

# Wrap X as a bigstatsr FBM (pass through if already one).
.ldClumpFbm <- function(X) {
    if (inherits(X, "FBM")) {
        return(X)
    }
    codeVec <- c(0, 1, 2, rep(NA, 256L - 3L))
    bigstatsr::FBM.code256(
        nrow = nrow(X),
        ncol = ncol(X),
        init = X,
        code = codeVec
    )
}

#' LD clumping by a per-variant score using bigsnpr
#'
#' Wraps \code{bigsnpr::snp_clumping} with the boilerplate of wrapping a numeric
#' dosage matrix into a \code{bigstatsr::FBM.code256} object and of handling the
#' common pitfall of a single-variant input.
#'
#' @param X Numeric matrix of 0/1/2 allele dosages, n rows by p variants. Column
#'   names are expected to be variant IDs but are not required.
#' @param score Numeric vector of length \code{ncol(X)}. Higher values favour
#'   retention during clumping (e.g. -log10 p, |Z|, MAF). May be \code{NULL}, in
#'   which case bigsnpr falls back to minor allele frequency computed from
#'   \code{X}.
#' @param chr Integer or character vector of length \code{ncol(X)} giving the
#'   chromosome for each variant.
#' @param pos Integer vector of length \code{ncol(X)} giving the base-pair
#'   position for each variant.
#' @param r2 Numeric in (0, 1]. r-squared threshold for clumping (variants
#'   within \code{windowKb} whose r2 exceeds \code{r2} and have lower
#'   \code{score} are removed). Default 0.2.
#' @param windowKb Numeric. Window size in kilobases. Default is \code{100 /
#'   r2}, matching the common "ld-clump size = 100/r2" heuristic used in many
#'   GWAS pipelines.
#' @param verbose Logical. If TRUE, print the number of retained variants.
#'   Default FALSE.
#'
#' @return An integer vector of indices (into \code{X} columns) kept after
#'   clumping. For a single-column \code{X}, returns \code{1L}.
#'
#' @examples
#'   set.seed(1)
#'   n <- 500; p <- 20
#'   X <- matrix(rbinom(n * p, 2, 0.3), n, p)
#'   colnames(X) <- paste0("chr1:", seq_len(p) * 1000, ":A:G")
#'   s <- runif(p)
#'   chr <- rep(1L, p); pos <- seq_len(p) * 1000L
#'   keep <- ldClumpByScore(X, score = s, chr = chr, pos = pos, r2 = 0.2)
#'
#' @export
ldClumpByScore <- function(
    X,
    score,
    chr,
    pos,
    r2 = 0.2,
    windowKb = 100 / r2,
    verbose = FALSE
) {
    .ldClumpCheckDeps()
    .ldClumpValidate(X, score, chr, pos)
    if (ncol(X) == 1L) {
        if (verbose) {
            message("ldClumpByScore: single variant, skipping clumping")
        }
        return(1L)
    }
    G <- .ldClumpFbm(X)
    keep <- bigsnpr::snp_clumping(
        G = G,
        infos.chr = as.integer(chr),
        infos.pos = as.integer(pos),
        S = score,
        thr.r2 = r2,
        size = windowKb
    )
    if (verbose) {
        message(
            "ldClumpByScore: ",
            length(keep),
            " / ",
            ncol(X),
            " variants retained at r2 <= ",
            r2
        )
    }
    keep
}


# =============================================================================
# Block-wise LD loaders
# -----------------------------------------------------------------------------
# High-level helpers that sit on top of `loadLdMatrix` / `processLdMatrix`
# to retrieve per-block LD or genotype matrices on demand. Used by
# downstream pipelines (cTWAS, etc.) that need to walk many LD blocks
# without materializing them all in memory at once.
# =============================================================================

#' Extract the LD or genotype matrix from an LdData S4 object.
#' @param ld An LdData object.
#' @param wantGenotype Logical; if TRUE, extract the genotype matrix (via
#'   \code{getGenotypes()}).
#' @return A matrix.
#' @noRd
extractLdMatrix <- function(ld, wantGenotype = FALSE) {
    if (!is(ld, "LdData")) {
        stop("ld must be an LdData object")
    }
    if (wantGenotype && hasGenotypes(ld)) {
        return(getGenotypes(ld))
    }
    getCorrelation(ld)
}

# Validate an ldLoader spec: exactly one source + per-mode requirements.
.ldLoaderValidate <- function(rList, xList, ldMetaPath, regions, ldInfo) {
    nSources <- sum(
        !is.null(rList),
        !is.null(xList),
        !is.null(ldMetaPath),
        !is.null(ldInfo)
    )
    if (nSources != 1) {
        stop("Provide exactly one of rList, xList, ldMetaPath, or ldInfo.")
    }
    if (!is.null(ldMetaPath) && is.null(regions)) {
        stop("'regions' is required when using ldMetaPath.")
    }
    if (
        !is.null(ldInfo) &&
            (!is.data.frame(ldInfo) || !"LD_file" %in% colnames(ldInfo))
    ) {
        stop("ldInfo must be a data.frame with column 'LD_file'.")
    }
}

#' Create an LD loader for on-demand block-wise LD retrieval
#'
#' Constructs a loader function that retrieves per-block LD matrices on demand.
#' This avoids loading all blocks into memory simultaneously, which is critical
#' for genome-wide analyses with hundreds of blocks.
#'
#' Four modes are supported:
#'
#' \describe{
#'   \item{list mode (R)}{Pre-loaded list of LD correlation matrices.
#'     Simple but uses more memory. Set \code{R_list}.}
#'   \item{list mode (X)}{Pre-loaded list of genotype matrices (n x p_g).
#'     Set \code{X_list}.}
#'   \item{region mode}{Loads LD from a pecotmr metadata TSV file on the fly
#'     via \code{\link{loadLdMatrix}}. Memory-efficient for large datasets.
#'     Set \code{ld_meta_path} and \code{regions}.}
#'   \item{ldInfo mode}{Loads pre-computed LD blocks from \code{.cor.xz}
#'     files listed in an \code{ldInfo} data.frame (as returned by
#'     cTWAS meta-data utilities). Set \code{ldInfo}.}
#' }
#'
#' @param rList List of G precomputed LD correlation matrices (p_g x p_g).
#' @param xList List of G genotype matrices (n x p_g).
#' @param ldMetaPath Path to a pecotmr LD metadata TSV file (as used by
#'   \code{\link{loadLdMatrix}}).
#' @param regions Character vector of G region strings (e.g.,
#'   \code{"chr22:17238266-19744294"}). Required when \code{ldMetaPath} is used.
#' @param ldInfo A data.frame with column \code{LD_file} (paths to genotype
#'   files or \code{.cor.xz} LD matrix files) and optionally \code{SNP_file}
#'   (paths to companion \code{.bim} files for pre-computed blocks; defaults to
#'   \code{paste0(LD_file, ".bim")} if absent). Genotype paths can be PLINK2
#'   prefixes, PLINK1 prefixes, VCF files, or GDS files. As returned by cTWAS
#'   meta-data utilities.
#' @param returnGenotype Logical. When using region mode, return the genotype
#'   matrix X (\code{TRUE}) or LD correlation R (\code{FALSE}, default).
#' @param maxVariants Integer or \code{NULL}. If set, randomly subsample blocks
#'   larger than this to control memory usage.
#'
#' @return An \code{ldLoaderSpec} object (an opaque list describing the source).
#'   Pass it with a block index to \code{\link{loadLdBlock}} to load one block.
#'
#' @seealso \code{\link{loadLdBlock}}
#' @examples
#' # List mode with pre-computed LD
#' R1 <- diag(10)
#' R2 <- diag(15)
#' spec <- ldLoader(rList = list(R1, R2))
#' loadLdBlock(spec, 1)  # returns R1
#' loadLdBlock(spec, 2)  # returns R2
#'
#' @export
ldLoader <- function(
    rList = NULL,
    xList = NULL,
    ldMetaPath = NULL,
    regions = NULL,
    ldInfo = NULL,
    returnGenotype = FALSE,
    maxVariants = NULL
) {
    .ldLoaderValidate(rList, xList, ldMetaPath, regions, ldInfo)
    mode <- if (!is.null(rList)) {
        "rList"
    } else if (!is.null(xList)) {
        "xList"
    } else if (!is.null(ldMetaPath)) {
        "meta"
    } else {
        "info"
    }
    structure(
        list(
            mode = mode,
            rList = rList,
            xList = xList,
            ldMetaPath = ldMetaPath,
            regions = regions,
            ldInfo = ldInfo,
            returnGenotype = returnGenotype,
            maxVariants = maxVariants
        ),
        class = "ldLoaderSpec"
    )
}


# ---- Per-block loaders (one per ldLoader source mode) -----------------------
# Top-level workers dispatched by .ldLoadBlock() on the spec's $mode; formerly
# branch closures inside ldLoader().

# @noRd
.ldLoadRList <- function(g, rList, maxVariants) {
    R <- rList[[g]]
    if (!is.null(maxVariants) && ncol(R) > maxVariants) {
        keep <- sort(sample(ncol(R), maxVariants))
        R <- R[keep, keep]
    }
    R
}

# @noRd
.ldLoadXList <- function(g, xList, maxVariants) {
    X <- xList[[g]]
    if (!is.null(maxVariants) && ncol(X) > maxVariants) {
        keep <- sort(sample(ncol(X), maxVariants))
        X <- X[, keep]
    }
    X
}

# @noRd
.ldLoadRegionMeta <- function(
    g,
    ldMetaPath,
    regions,
    returnGenotype,
    maxVariants
) {
    ld <- loadLdMatrix(
        ldMetaPath,
        region = regions[g],
        returnGenotype = returnGenotype
    )
    mat <- extractLdMatrix(ld, wantGenotype = returnGenotype)
    if (!is.null(maxVariants) && ncol(mat) > maxVariants) {
        keep <- sort(sample(ncol(mat), maxVariants))
        if (returnGenotype || nrow(mat) > ncol(mat)) {
            mat <- mat[, keep]
        } else {
            mat <- mat[keep, keep]
        }
    }
    # Center and scale genotype matrices
    if (returnGenotype || nrow(mat) > ncol(mat)) {
        mat <- scale(mat)
        mat[is.na(mat)] <- 0
    }
    mat
}

# @noRd
.ldLoadIdInfo <- function(g, ldInfo, maxVariants) {
    ldPath <- ldInfo$LD_file[g]

    # Auto-detect format: genotype source or pre-computed block
    if (isGenotypeSource(ldPath)) {
        geno <- loadGenotypeRegion(ldPath)
        mat <- computeLd(geno)
    } else {
        # Pre-computed .cor.xz block
        snpFile <- if ("SNP_file" %in% colnames(ldInfo)) {
            ldInfo$SNP_file[g]
        } else {
            NULL # let processLdMatrix auto-detect .bim/.pvar/.pvar.zst
        }
        ld <- processLdMatrix(ldPath, snpFile)
        mat <- extractLdMatrix(ld)
    }

    if (!is.null(maxVariants) && ncol(mat) > maxVariants) {
        keep <- sort(sample(ncol(mat), maxVariants))
        mat <- mat[keep, keep]
    }
    mat
}

# Dispatch a single block load by the spec's source mode.
# @noRd
.ldLoadBlock <- function(spec, g) {
    switch(
        spec$mode,
        rList = .ldLoadRList(g, spec$rList, spec$maxVariants),
        xList = .ldLoadXList(g, spec$xList, spec$maxVariants),
        meta = .ldLoadRegionMeta(
            g,
            spec$ldMetaPath,
            spec$regions,
            spec$returnGenotype,
            spec$maxVariants
        ),
        info = .ldLoadIdInfo(g, spec$ldInfo, spec$maxVariants)
    )
}

#' Load one LD block from an ldLoader spec
#'
#' Given an \code{ldLoaderSpec} (from \code{\link{ldLoader}}) and a block index
#' \code{g}, load the corresponding LD correlation matrix (or the genotype
#' matrix, in region mode with \code{returnGenotype = TRUE}).
#'
#' @param spec An \code{ldLoaderSpec} object returned by \code{\link{ldLoader}}.
#' @param g Integer block index (1-based).
#' @return The LD correlation matrix or genotype matrix for block \code{g}.
#' @seealso \code{\link{ldLoader}}
#' @examples
#' spec <- ldLoader(rList = list(diag(10), diag(15)))
#' loadLdBlock(spec, 1)
#' @export
loadLdBlock <- function(spec, g) {
    if (!inherits(spec, "ldLoaderSpec")) {
        stop("`spec` must be an ldLoaderSpec (from ldLoader()).")
    }
    .ldLoadBlock(spec, g)
}


# =============================================================================
# LD correlation matrix from a dosage matrix
# -----------------------------------------------------------------------------
# Direct LD computation from an n x p dosage matrix. The internal backend
# uses Rfast::cora when available (else base cor()); optional snprelate /
# snpstats backends round-trip through a temp GDS or SnpMatrix. The
# population and GCTA methods match PLINK / GCTA conventions for missing
# data handling.
# =============================================================================

# --- computeLd method helpers -----------------------------------------------

# Non-sample methods only support the internal backend.
.computeLdRequireInternal <- function(backend) {
    if (backend != "internal") {
        stop(
            "backend '",
            backend,
            "' is only supported with method='sample'."
        )
    }
}

# Sample correlation (N-1 denominator) via the requested backend.
.computeLdSample <- function(X, backend) {
    if (backend == "snprelate") {
        return(.computeLdSnprelate(X))
    }
    if (backend == "snpstats") {
        return(.computeLdSnpstats(X))
    }
    # internal backend: Rfast::cora if available, else base cor(). Mean-impute
    # only when NAs exist (PLINK2 data typically has none).
    X_imp <- X
    if (anyNA(X_imp)) {
        colMeansX <- colMeans(X_imp, na.rm = TRUE)
        naPos <- which(is.na(X_imp), arr.ind = TRUE)
        X_imp[naPos] <- colMeansX[naPos[, 2]]
    }
    if (requireNamespace("Rfast", quietly = TRUE)) {
        # large=FALSE uses tcrossprod internally, ~40x faster than large=TRUE.
        Rfast::cora(X_imp, large = FALSE)
    } else {
        cor(X_imp)
    }
}

# Population variance (N denominator, GCTA-style; missing set to column mean 0).
.computeLdPopulation <- function(X, trimSamples) {
    if (trimSamples) {
        N_kept <- (nrow(X) %/% 4L) * 4L
        if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
    }
    N <- nrow(X)
    colMeansX <- colMeans(X, na.rm = TRUE)
    colVarsX <- colMeans(X^2, na.rm = TRUE) - colMeansX^2
    # Covariance divides by total N (GCTA convention); heterogeneous missingness
    # slightly deflates cross-column correlations.
    if (anyNA(X)) {
        naRates <- colMeans(is.na(X))
        if (max(naRates) - min(naRates) > 0.1) {
            warning(
                "Population LD method with heterogeneous missingness ",
                "(max NA rate ",
                round(max(naRates), 3),
                ", min ",
                round(min(naRates), 3),
                "): correlations may be biased. Consider using ",
                "method='sample' which handles missingness via mean imputation."
            )
        }
    }
    X_c <- sweep(X, 2, colMeansX)
    X_c[is.na(X_c)] <- 0
    covMat <- crossprod(X_c) / N
    sdVec <- sqrt(colVarsX)
    covMat / outer(sdVec, sdVec)
}

# GCTA per-pair missing-data covariance (matches DENTIST calcLDFromBfile_gcta):
# tracks per-pair non-missing counts and applies a correction term.
.gctaCovariance <- function(X, colMeansX, N, p) {
    notNa <- !is.na(X)
    X_zero <- X
    X_zero[is.na(X_zero)] <- 0
    pairCounts <- crossprod(notNa * 1.0)
    # E_i2[i,j] = pairSums[i,j] / N: mean of SNP i over samples where j is
    # observed; p x p, row i col j = sum of X_i where j non-missing, / N.
    pairSums <- crossprod(X_zero, notNa * 1.0)
    sum_XY <- crossprod(X_zero)
    E_i2 <- pairSums / N
    E_j2 <- t(E_i2)
    sum_XY /
        N +
        outer(colMeansX, colMeansX) * (pairCounts / N) -
        colMeansX * E_j2 -
        E_i2 * rep(colMeansX, each = p)
}

# GCTA LD: per-pair missing-data correction, then correlation.
.computeLdGcta <- function(X, trimSamples) {
    if (trimSamples) {
        N_kept <- (nrow(X) %/% 4L) * 4L
        if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
    }
    N <- nrow(X)
    p <- ncol(X)
    colMeansX <- colMeans(X, na.rm = TRUE)
    colVarsX <- colMeans(X^2, na.rm = TRUE) - colMeansX^2
    covMat <- .gctaCovariance(X, colMeansX, N, p)
    sdVec <- sqrt(colVarsX)
    sdOuter <- outer(sdVec, sdVec)
    R <- matrix(0.001, p, p)
    valid <- sdOuter > 0
    R[valid] <- covMat[valid] / sdOuter[valid]
    R
}

computeLd <- function(
    X,
    method = c("sample", "population", "gcta"),
    backend = c("internal", "snprelate", "snpstats"),
    trimSamples = FALSE,
    shrinkage = 0
) {
    if (is.null(X)) {
        stop("X must be provided.")
    }
    method <- match.arg(method)
    backend <- match.arg(backend)
    nms <- colnames(X)
    if (method == "sample") {
        R <- .computeLdSample(X, backend)
    } else if (method == "population") {
        .computeLdRequireInternal(backend)
        R <- .computeLdPopulation(X, trimSamples)
    } else {
        .computeLdRequireInternal(backend)
        R <- .computeLdGcta(X, trimSamples)
    }
    diag(R) <- 1.0
    R[is.na(R) | is.nan(R)] <- 0
    # Optional shrinkage toward identity (lassosum, Mak et al 2017).
    if (shrinkage > 0 && shrinkage <= 1) {
        R <- (1 - shrinkage) * R + shrinkage * diag(nrow(R))
    }
    colnames(R) <- rownames(R) <- nms
    R
}

#' Compute LD via SNPRelate (creates a temporary GDS file from the dosage
#' matrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix.
#' @noRd
.computeLdSnprelate <- function(X) {
    # nocov start
    if (!requireNamespace("SNPRelate", quietly = TRUE)) {
        stop("Package 'SNPRelate' is required for backend='snprelate'")
    }
    if (!requireNamespace("gdsfmt", quietly = TRUE)) {
        stop("Package 'gdsfmt' is required for backend='snprelate'")
    }
    # nocov end

    tmpGds <- tempfile(fileext = ".gds")
    on.exit(unlink(tmpGds), add = TRUE)

    # Round to integer dosage for GDS (0/1/2)
    X_int <- round(X)
    storage.mode(X_int) <- "integer"
    X_int[is.na(X_int)] <- 3L # GDS missing code

    snpIds <- colnames(X) %||% seq_len(ncol(X))
    sampleIds <- rownames(X) %||% seq_len(nrow(X))

    SNPRelate::snpgdsCreateGeno(
        tmpGds,
        genmat = X_int,
        sample.id = sampleIds,
        snp.id = snpIds,
        snp.chromosome = rep(1L, ncol(X)),
        snp.position = seq_len(ncol(X)),
        snpfirstdim = FALSE
    )

    gds <- SNPRelate::snpgdsOpen(tmpGds, readonly = TRUE)
    on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

    ldObj <- SNPRelate::snpgdsLDMat(
        gds,
        method = "corr",
        slide = -1,
        verbose = FALSE
    )
    ldObj$LD
}

#' Compute LD via snpStats (converts dosage matrix to SnpMatrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix (r, not r^2).
#' @noRd
.computeLdSnpstats <- function(X) {
    # nocov start
    if (!requireNamespace("snpStats", quietly = TRUE)) {
        stop("Package 'snpStats' is required for backend='snpstats'")
    }
    # nocov end

    # snpStats expects counts of the B allele as raw codes: 1=AA, 2=AB, 3=BB,
    # 0=NA pecotmr dosage is ALT count (0/1/2), so map: 0->1, 1->2, 2->3, NA->0
    X_raw <- round(X) + 1L
    X_raw[is.na(X) | X_raw < 1L] <- 0L
    X_raw[X_raw > 3L] <- 3L
    storage.mode(X_raw) <- "raw"
    sm <- new("SnpMatrix", X_raw)

    R <- as.matrix(snpStats::ld(sm, stats = "R", depth = ncol(X) - 1L))
    # snpStats::ld returns a sparse-like matrix; ensure full dense
    R[is.na(R)] <- 0
    diag(R) <- 1
    R
}
