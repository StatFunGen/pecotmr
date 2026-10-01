#' Deduplicate and sort genomic regions by chromosome and start position.
#' @importFrom dplyr distinct arrange
#' @noRd
orderDedupRegions <- function(df) {
    mutate(df, chrom = canonChrom(.data$chrom)) |>
        distinct(.data$chrom, .data$start, .keep_all = TRUE) |>
        arrange(chromOrder(.data$chrom), .data$start)
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
    chromData <- genomicData |> filter(.data$chrom == regionChrom)
    if (nrow(chromData) == 0) {
        msg <- glue("No data for chromosome {regionChrom}")
        abort(msg)
    }

    # Clamp query to available range
    regionStart <- max(regionStart, min(chromData$start))
    regionEnd <- min(regionEnd, max(chromData$end))

    startRow <- genomicData |>
        filter(
            .data$chrom == regionChrom,
            .data$start <= regionStart,
            .data$end > regionStart
        ) |>
        slice(1)
    endRow <- genomicData |>
        filter(
            .data$chrom == regionChrom,
            .data$start < regionEnd,
            .data$end >= regionEnd
        ) |>
        arrange(desc(.data$end)) |>
        slice(1)

    if (nrow(startRow) == 0 || nrow(endRow) == 0) {
        msg <- glue(
            "Region {regionChrom}:{regionStart}-{regionEnd} is not ",
            "covered by any rows in the LD metadata."
        )
        abort(msg)
    }
    list(startRow = startRow, endRow = endRow)
}

#' Validate that startRow..endRow fully covers [regionStart, regionEnd].
#' @noRd
validateSelectedRegion <- function(startRow, endRow, regionStart, regionEnd) {
    if (startRow$start > regionStart || endRow$end < regionEnd) {
        availStart <- startRow$start
        availEnd <- endRow$end
        msg <- glue(
            "Region {regionStart}-{regionEnd} is not fully covered by ",
            "the LD metadata (available: {availStart}-{availEnd})."
        )
        abort(msg)
    }
}

#' Extract values of a column for rows spanning the intersection range.
#' @noRd
extractFilePaths <- function(genomicData, intersectionRows, columnToExtract) {
    if (!is_in(columnToExtract, names(genomicData))) {
        msg <- glue("Column '{columnToExtract}' not found in genomic data.")
        abort(msg)
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
.findValidFilePath <- function(targetFilePath, referenceFilePath) {
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
    msg <- glue(
        "Both reference and target file paths do not work. Tried ",
        "paths: '{referenceFilePath}' and '{targetFullPath}'"
    )
    abort(msg)
}

# Vectorised .findValidFilePath over a vector of target paths.
# @noRd
.findValidFilePaths <- function(referenceFilePath, targetFilePaths) {
    map_chr(
        targetFilePaths,
        .findValidFilePath,
        referenceFilePath = referenceFilePath
    )
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
    parts <- str_split(genomicData$path, ",", simplify = TRUE)
    pathNames <- if (ncol(parts) == 2) {
        c("LD_file_path", "bim_file_path")
    } else {
        "LD_file_path"
    }
    filePath <- `names<-`(
        as_tibble(parts, .name_repair = "minimal"),
        pathNames
    )
    bind_cols(genomicData, filePath) |>
        select(-any_of("path"))
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
    bimPaths <- if (is_in("bim_file_path", names(genomicData))) {
        .findValidFilePaths(
            ldReferenceMetaFile,
            extractFilePaths(genomicData, intersectionRows, "bim_file_path")
        )
    } else {
        NULL
    }
    list(ldPaths = ldPaths, bimPaths = bimPaths)
}

# A metadata row recorded as start=0, end=0 covers the whole chromosome
# rather than an empty interval, so widen its end before any overlap test.
# @noRd
.ldMetaWidenWholeChrom <- function(meta) {
    mutate(
        meta,
        end = if_else(
            .data$start == 0 & .data$end == 0,
            Inf,
            as.numeric(.data$end)
        )
    )
}

getRegionalLdMeta <- function(
    ldReferenceMetaFile,
    region,
    completeCoverageRequired = FALSE
) {
    genomicData <- `names<-`(
        vroom(ldReferenceMetaFile),
        c("chrom", "start", "end", "path")
    ) |>
        .ldMetaWidenWholeChrom() |>
        orderDedupRegions() |>
        .regionalLdParsePaths()
    region <- `names<-`(
        parseRegion(region),
        c("chrom", "start", "end")
    ) |>
        orderDedupRegions()
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
    candidates <- str_c(ldFilePath, c(".bim", ".pvar", ".pvar.zst"))
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) {
        msg <- glue(
            "No variant file found for: {ldFilePath} ",
            "(tried .bim, .pvar, .pvar.zst)"
        )
        abort(msg)
    }
    found[1]
}

# Read + normalise the LD variant metadata (canonical chrom / variant id / GD).
.processLdVariants <- function(snpFilePath) {
    raw <- readVariantMetadata(snpFilePath)
    isPvar <- !is_in("gpos", names(raw))
    ldVariants <- mutate(
        raw,
        chrom = canonChrom(.data$chrom),
        variants = normalizeVariantId(.data$id)
    )
    if (!isPvar) {
        return(rename(ldVariants, GD = "gpos"))
    }
    # A .pvar carries no genetic distance, so GD and pos both take the
    # position parsed back out of the variant id.
    parsedPos <- map_int(ldVariants$variants, .ldVariantPos)
    rename(ldVariants, GD = "pos") |>
        mutate(GD = parsedPos, pos = parsedPos)
}

processLdMatrix <- function(ldFilePath, snpFilePath = NULL) {
    ldFileCon <- xzfile(ldFilePath)
    ldValues <- scan(ldFileCon, quiet = TRUE)
    close(ldFileCon)
    snpFilePath <- .processLdSnpFile(ldFilePath, snpFilePath)
    rawVariants <- .processLdVariants(snpFilePath)
    # Label and symmetrize the matrix.
    labelled <- `dimnames<-`(
        matrix(ldValues, ncol = sqrt(length(ldValues)), byrow = TRUE),
        list(rawVariants$variants, rawVariants$variants)
    )
    # Only one triangle is stored on disk; mirror whichever one is empty.
    lower <- lower.tri(labelled)
    upper <- upper.tri(labelled)
    ldMatrix <- if (all(labelled[lower] == 0)) {
        replace(labelled, lower, t(labelled)[lower])
    } else {
        replace(labelled, upper, t(labelled)[upper])
    }
    # Order variants by genomic position.
    posOrder <- order(map_int(rawVariants$variants, .ldVariantPos))
    ldVariants <- slice(rawVariants, posOrder)
    list(
        ldMatrix = ldMatrix[ldVariants$variants, ldVariants$variants],
        ldVariants = ldVariants
    )
}

#' Subset an LD matrix and variant info to a genomic region, optionally further
#' restricted to specific coordinates.
#' @importFrom dplyr mutate select inner_join
#' @noRd
extractLdForRegion <- function(ldMatrix, variants, region, extractCoordinates) {
    inRegion <- filter(
        variants,
        .data$chrom == region$chrom &
            .data$pos >= region$start &
            .data$pos <= region$end
    )
    extracted <- if (is.null(extractCoordinates)) {
        inRegion
    } else {
        .ldRestrictToCoordinates(inRegion, extractCoordinates)
    }
    mat <- ldMatrix[extracted$variants, extracted$variants, drop = FALSE]
    list(extractedLdMatrix = mat, extractedLdVariants = extracted)
}

# The region's variants restricted to an explicit (chrom, pos) list, keeping
# only the variant-info columns the callers read back.
# @noRd
.ldRestrictToCoordinates <- function(inRegion, extractCoordinates) {
    wanted <- extractCoordinates |>
        mutate(chrom = canonChrom(.data$chrom)) |>
        select("chrom", "pos")
    joined <- inRegion |>
        mutate(chrom = canonChrom(.data$chrom)) |>
        inner_join(wanted, by = c("chrom", "pos"))
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
        names(joined)
    )
    select(joined, all_of(keepCols))
}

# Concatenate per-block variant-id lists into one deduplicated vector, dropping
# a repeated boundary variant shared between adjacent blocks.
# @noRd
.ldMergeVariants <- function(variantList) {
    # Only the running tail decides whether a block's first id is a repeated
    # boundary variant, so the merge is a fold.
    reduce(variantList, .ldAppendBlockVariants, .init = character(0))
}

# @noRd
.ldAppendBlockVariants <- function(merged, v) {
    ids <- if (is.list(v) && !is.null(v$variants)) v$variants else v
    if (length(ids) == 0) {
        return(merged)
    }
    repeatsBoundary <- length(merged) > 0 && tail(merged, 1) == ids[1]
    c(merged, if (repeatsBoundary) ids[-1] else ids)
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

    # Deliberate preallocate-and-scatter: `combined` is variants x variants,
    # so building it by folding full-size copies would multiply the memory
    # this function needs by the number of blocks.
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
        msg <- glue(
            "returnGenotype=TRUE requires genotype files, not ",
            "pre-computed LD matrices."
        )
        abort(msg)
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
    full <- getCorrelation(result)
    corr <- if (is.null(full)) {
        NULL
    } else {
        full[-dupIdx, -dupIdx, drop = FALSE]
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
#' @param block Integer block index (or vector of them), for sources that
#'   carry no coordinates: an \code{ldInfo} table (the cTWAS \code{LD_map}
#'   shape -- an \code{LD_file} column and optionally \code{SNP_file}), a
#'   correlation or genotype matrix, or a list of either. Supply exactly one
#'   of \code{region} and \code{block}.
#' @param dropMonomorphic Logical. Drop variants with no variation, which
#'   carry no LD and make a correlation undefined. Needs allele frequencies,
#'   so it applies to genotype sources.
#' @param materializeGenotypes Logical. Read the dosages once and keep them on
#'   the returned object, so later access does no file I/O. Costs memory and
#'   saves repeated reads.
#' @param maxVariants Integer or \code{NULL}. Randomly thin any block larger
#'   than this, to bound memory.
#' @param seed Integer or \code{NULL}. Seeds the \code{maxVariants} draw,
#'   offset by the block index so blocks are independent yet reproducible, via
#'   a scoped \code{withr::local_seed} that leaves the session RNG alone.
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
#' # Several regions at once return one LdData per block.
#' length(loadLdMatrix(meta, region = rep("chr22:16000000-18000000", 2)))
#' @export
loadLdMatrix <- function(
    ldMetaFilePath,
    region = NULL,
    block = NULL,
    extractCoordinates = NULL,
    returnGenotype = FALSE,
    dropMonomorphic = FALSE,
    materializeGenotypes = FALSE,
    maxVariants = NULL,
    seed = NULL,
    nSample = NULL
) {
    .ldLoadValidateAddress(ldMetaFilePath, region, block)
    keys <- if (is.null(region)) block else region
    # Only an atomic vector addresses several blocks. `region` also accepts a
    # single data.frame / GRanges spec, whose length is its column count --
    # treating that as many keys would map over the columns.
    if (is.atomic(keys) && length(keys) > 1L) {
        return(map(
            keys,
            .loadLdOne,
            source = ldMetaFilePath,
            byRegion = is.null(block),
            extractCoordinates = extractCoordinates,
            returnGenotype = returnGenotype,
            dropMonomorphic = dropMonomorphic,
            materializeGenotypes = materializeGenotypes,
            maxVariants = maxVariants,
            seed = seed,
            nSample = nSample
        ))
    }
    .loadLdOne(
        keys,
        ldMetaFilePath,
        byRegion = is.null(block),
        extractCoordinates = extractCoordinates,
        returnGenotype = returnGenotype,
        dropMonomorphic = dropMonomorphic,
        materializeGenotypes = materializeGenotypes,
        maxVariants = maxVariants,
        seed = seed,
        nSample = nSample
    )
}

# A source is addressed either by genomic region (an LD meta file, which
# carries coordinates) or by block index (an `ldInfo` table or in-memory
# matrices, which do not). Exactly one addressing mode applies.
# @noRd
.ldLoadValidateAddress <- function(source, region, block) {
    if (is.null(region) && is.null(block)) {
        abort("loadLdMatrix: supply either `region` or `block`.")
    }
    if (!is.null(region) && !is.null(block)) {
        abort("loadLdMatrix: supply `region` or `block`, not both.")
    }
    if (!is.null(region) && !is.character(source)) {
        msg <- glue(
            "loadLdMatrix: `region` addresses an LD meta file, but ",
            "`ldMetaFilePath` is a {class(source)[[1L]]}. Coordinate-free ",
            "sources (an `ldInfo` table, a matrix, or a list of matrices) ",
            "are addressed with `block`."
        )
        abort(msg)
    }
    invisible(TRUE)
}

# One block, whatever the source. Everything that is not an LD meta file is
# coordinate-free and indexed by `block`.
# @noRd
.loadLdOne <- function(
    key,
    source,
    byRegion,
    extractCoordinates,
    returnGenotype,
    dropMonomorphic,
    materializeGenotypes,
    maxVariants,
    seed,
    nSample
) {
    loaded <- if (byRegion) {
        .loadLdFromMeta(
            source,
            key,
            extractCoordinates,
            returnGenotype,
            nSample
        )
    } else {
        .loadLdFromIndexed(source, key)
    }
    loaded |>
        .loadLdDedup() |>
        .ldApplyMonomorphic(dropMonomorphic) |>
        .ldApplySubsample(maxVariants, seed, key) |>
        .ldApplyMaterialize(materializeGenotypes)
}

# Coordinate-free sources, addressed by block index: an `ldInfo` table (the
# cTWAS LD_map shape -- LD_file plus optional SNP_file, no coordinates), a
# single in-memory matrix, or a list of them. Every source returns an LdData,
# so callers never have to branch on what they loaded from.
# @noRd
.loadLdFromIndexed <- function(source, block) {
    if (is.data.frame(source)) {
        return(.ldInfoBlock(source, block))
    }
    mat <- if (is.list(source)) source[[block]] else source
    if (!is.matrix(mat)) {
        msg <- glue(
            "loadLdMatrix: cannot address a {class(source)[[1L]]} by block. ",
            "Supply an LD meta file path, an `ldInfo` data.frame, a matrix, ",
            "or a list of matrices."
        )
        abort(msg)
    }
    .ldDataFromMatrix(mat, isGenotype = nrow(mat) > ncol(mat))
}

# One row of an `ldInfo` table. Genotype paths are read and correlated;
# pre-computed .cor.xz blocks are read directly.
# @noRd
.ldInfoBlock <- function(ldInfo, block) {
    if (!is_in("LD_file", colnames(ldInfo))) {
        abort("loadLdMatrix: an `ldInfo` table needs an `LD_file` column.")
    }
    ldPath <- as.character(ldInfo$LD_file)[block]
    if (isGenotypeSource(ldPath)) {
        geno <- loadGenotypeRegion(ldPath)
        return(.ldDataFromMatrix(geno, isGenotype = TRUE))
    }
    snpFile <- if (is_in("SNP_file", colnames(ldInfo))) {
        as.character(ldInfo$SNP_file)[block]
    } else {
        NULL # processLdMatrix auto-detects the .bim / .pvar companion
    }
    .ldDataFromProcessed(processLdMatrix(ldPath, snpFile))
}

# Wrap processLdMatrix()'s (matrix, variants) pair as an LdData so BOTH
# `ldInfo` sources return one type -- the contract .loadLdFromIndexed states
# and the post-load chain (dedup / monomorphic / subsample) relies on.
#
# Not .ldDataFromMatrix(): that one is for a bare matrix with no metadata and
# substitutes placeholder chrNA:1..n coordinates. A precomputed LD block comes
# with its .bim/.pvar, so the real chrom/pos/alleles are carried through.
# @noRd
.ldDataFromProcessed <- function(proc) {
    v <- proc$ldVariants
    n <- nrow(v)
    gr <- .refPanelToGranges(data.frame(
        chrom = as.character(v$chrom),
        pos = as.integer(v$pos),
        variant_id = as.character(v$variants),
        A1 = as.character(v$A1),
        A2 = as.character(v$A2),
        stringsAsFactors = FALSE
    ))
    LdData(
        correlation = proc$ldMatrix,
        variants = gr,
        blockMetadata = tibble(
            blockId = 1L,
            size = n,
            startIdx = 1L,
            endIdx = n
        ),
        nRef = 0L
    )
}

# Wrap a bare matrix as an LdData so every source returns one type. Variant
# identity comes from the dimnames when present.
# @noRd
.ldDataFromMatrix <- function(mat, isGenotype) {
    n <- if (isGenotype) ncol(mat) else nrow(mat)
    named <- if (isGenotype) colnames(mat) else rownames(mat)
    ids <- named %||% str_c("v", seq_len(n))
    gr <- S4Vectors::`mcols<-`(
        GRanges(
            seqnames = rep("chrNA", n),
            ranges = IRanges::IRanges(start = seq_len(n), width = 1L)
        ),
        value = S4Vectors::DataFrame(variant_id = ids)
    )
    LdData(
        correlation = if (isGenotype) NULL else mat,
        genotypeHandle = if (isGenotype) mat else NULL,
        variants = gr,
        blockMetadata = tibble(blockId = 1L, size = n),
        nRef = if (isGenotype) nrow(mat) else 0L
    )
}

# Drop variants with no variation. They carry no LD, and a zero-variance
# column makes a correlation undefined.
# @noRd
.ldApplyMonomorphic <- function(ld, dropMonomorphic) {
    if (!isTRUE(dropMonomorphic)) {
        return(ld)
    }
    refPanel <- getRefPanel(ld)
    if (is.null(refPanel) || !is_in("allele_freq", colnames(refPanel))) {
        return(ld)
    }
    p <- refPanel$allele_freq
    keep <- !is.na(p) & p > 0 & p < 1
    if (all(keep)) {
        return(ld)
    }
    .ldSubsetData(ld, which(keep))
}

# Randomly thin an oversized block, seeded per block so a run is reproducible
# without disturbing the session RNG.
# @noRd
.ldApplySubsample <- function(ld, maxVariants, seed, key) {
    if (is.null(maxVariants) || length(ld) <= maxVariants) {
        return(ld)
    }
    if (!is.null(seed)) {
        offset <- if (is.numeric(key)) as.integer(key) else 0L
        withr::local_seed(as.integer(seed) + offset)
    }
    .ldSubsetData(ld, sort(sample(length(ld), maxVariants)))
}

# Read the dosages once and keep them on the object, so later access does no
# file I/O.
# @noRd
.ldApplyMaterialize <- function(ld, materializeGenotypes) {
    if (!isTRUE(materializeGenotypes) || !hasGenotypes(ld)) {
        return(ld)
    }
    X <- getGenotypes(ld)
    if (!is.matrix(X)) {
        return(ld)
    }
    LdData(
        correlation = NULL,
        genotypeHandle = X,
        snpIdx = NULL,
        variants = getVariantInfo(ld),
        blockMetadata = getBlockMetadata(ld),
        nRef = getNRef(ld)
    )
}

# Narrow an LdData to a subset of its variants, keeping correlation and
# genotypes consistent with the ranges.
# @noRd
.ldSubsetData <- function(ld, idx) {
    full <- ld@correlation
    R <- if (!is.null(full) && is.matrix(full)) {
        full[idx, idx, drop = FALSE]
    } else {
        full
    }
    raw <- ld@genotypeHandle
    gh <- if (is.matrix(raw)) raw[, idx, drop = FALSE] else raw
    snpIdx <- if (is.matrix(raw) || is.null(ld@snpIdx)) {
        ld@snpIdx
    } else {
        ld@snpIdx[idx]
    }
    LdData(
        correlation = R,
        genotypeHandle = gh,
        snpIdx = if (is.matrix(gh)) NULL else snpIdx,
        variants = getVariantInfo(ld)[idx],
        blockMetadata = getBlockMetadata(ld),
        nRef = getNRef(ld)
    )
}

# @noRd
.loadLdFromMeta <- function(
    ldMetaFilePath,
    region,
    extractCoordinates,
    returnGenotype,
    nSample
) {
    source <- resolveLdSource(ldMetaFilePath)
    isGeno <- is_in(source$type, c("plink2", "plink1", "vcf", "gds"))
    # "auto": return X for genotype sources, R for pre-computed.
    if (identical(returnGenotype, "auto")) {
        returnGenotype <- isGeno
    }
    .loadLdDispatch(
        source,
        isGeno,
        region,
        extractCoordinates,
        returnGenotype,
        nSample
    )
}

# ---------- Internal: resolve LD source type ----------

#' @noRd
hasPlink2Files <- function(prefix) {
    file.exists(str_c(prefix, ".pgen")) &&
        (file.exists(str_c(prefix, ".pvar")) ||
            file.exists(str_c(prefix, ".pvar.zst"))) &&
        file.exists(str_c(prefix, ".psam"))
}

#' @noRd
hasPlink1Files <- function(prefix) {
    file.exists(str_c(prefix, ".bed")) &&
        file.exists(str_c(prefix, ".bim")) &&
        file.exists(str_c(prefix, ".fam"))
}

#' @noRd
isVcfPath <- function(path) {
    str_detect(path, "\\.(vcf|vcf\\.gz|bcf)$") && file.exists(path)
}

#' @noRd
isGdsPath <- function(path) {
    str_detect(path, "\\.gds$") && file.exists(path)
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
#' @importFrom checkmate checkFileExists
.resolveLdReadMeta <- function(path) {
    res <- checkFileExists(path, access = "r")
    if (!isTRUE(res)) {
        msg <- glue(
            "LD metadata file: {res}",
            "\n  Expected: a TSV file with columns chrom, start, end, path.",
            .trim = FALSE
        )
        abort(msg)
    }
    meta <- as.data.frame(vroom(path, show_col_types = FALSE, n_max = 1))
    if (ncol(meta) < 4) {
        msg <- glue(
            "LD metadata file must have at least 4 columns (chrom, ",
            "start, end, path): {path}"
        )
        abort(msg)
    }
    `colnames<-`(
        meta,
        replace(colnames(meta), seq_len(4), c("chrom", "start", "end", "path"))
    )
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
    rawPath <- str_remove(meta$path[1], ",.*$")
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
        msg <- glue(
            "Metadata has start=0, end=0 but path does not resolve to ",
            "genotype files: {resolved}",
            "\n  The 0:0 sentinel is only valid for whole-chromosome ",
            "genotype files.",
            .trim = FALSE
        )
        abort(msg)
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
    meta <- as.data.frame(vroom(metaPath, show_col_types = FALSE)) |>
        `colnames<-`(c("chrom", "start", "end", "path")) |>
        mutate(chrom = canonChrom(.data$chrom))
    queryChrom <- canonChrom(parsed$chrom)

    matching <- meta[meta$chrom == queryChrom, , drop = FALSE]
    if (nrow(matching) == 0) {
        msg <- glue(
            "No entry for chromosome {queryChrom} in metadata file: ",
            "{metaPath}"
        )
        abort(msg)
    }
    rawPath <- str_remove(matching$path[1], ",.*$")
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
    afreq <- readAfreq(genotypePath)
    alleleFreq <- if (is.null(afreq)) {
        colMeans(X, na.rm = TRUE) / 2
    } else {
        freqMatch <- match(variantInfo$id, afreq$id)
        nUnmatched <- sum(is.na(freqMatch))
        if (nUnmatched > 0) {
            nFreq <- length(freqMatch)
            msg <- glue(
                "{nUnmatched} out of {nFreq} variants have no allele ",
                "frequency in .afreq file."
            )
            warn(msg)
        }
        afreq$alt_freq[freqMatch]
    }
    mutate(
        parseVariantId(variantIds),
        variant_id = variantIds,
        allele_freq = alleleFreq,
        !!!.ldPanelSampleCols(alleleFreq, nSample)
    )
}

# The panel's per-variant variance / non-missing count, which are only
# derivable once the caller declares the sample size behind the frequencies.
# @noRd
.ldPanelSampleCols <- function(p, nSample) {
    if (is.null(nSample)) {
        return(list())
    }
    list(
        variance = 2 * p * (1 - p) * nSample / (nSample - 1),
        n_nomiss = nSample
    )
}

# Single-block metadata spanning the loaded region.
.loadLdGtBlockMeta <- function(variantInfo, variantIds) {
    positions <- variantInfo$pos
    tibble(
        blockId = 1L,
        chrom = as.character(variantInfo$chrom[1]),
        blockStart = min(positions),
        blockEnd = max(positions),
        size = length(variantIds),
        startIdx = 1L,
        endIdx = length(variantIds)
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
    handle <- .readGenotypeHandle(genotypePath)
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
    variantInfo <- result$variant_info
    variantIds <- normalizeVariantId(formatVariantId(
        variantInfo$chrom,
        variantInfo$pos,
        variantInfo$A2,
        variantInfo$A1
    ))
    X <- `colnames<-`(result$X, variantIds)
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
# subset of an `ldSketch` `GenotypeHandle`. Called directly by
# twasWeightsPipeline, fineMappingPipeline, jointEngine,
# causalInferencePipeline, colocboostPipeline and the coloc / summaryStatsQc
# paths, which differ only in their error message prefix and in whether
# variants absent from the panel raise an error or get silently dropped --
# both of which are arguments here rather than a per-pipeline wrapper.
#
# Arguments:
#   ldSketch    A GenotypeHandle.
#   variantIds  Character vector of SNP IDs to extract.
#   label       Error-message prefix, e.g. "twasWeightsPipeline".
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
# Normalise an LD-sketch input to the RangedSummarizedExperiment the slot
# holds. Accepts a bare GenotypeHandle, which is what every caller has passed
# until now, and wraps it -- so a sketch and QtlDataset's genotype experiment
# end up being the same object with different provenance rather than two
# spellings of one idea.
# @noRd
.asLdSketch <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    if (methods::is(x, "RangedSummarizedExperiment")) {
        return(x)
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(.genotypeExperiment(x))
    }
    abort(glue(
        "`ldSketch` must be a genotype panel -- a ",
        "RangedSummarizedExperiment or a GenotypeHandle (got ",
        "{class(x)[[1L]]})."
    ))
}

# The GenotypeHandle a sketch reads through. The handle is the DelayedArray
# seed of the dosage assay, so it stays reachable after the sketch is subset.
# Tolerates a bare handle, which is what objects serialised before the slot
# became an RSE still carry.
# @noRd
.ldSketchHandle <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(x)
    }
    .ghSeedHandle(
        DelayedArray::seed(SummarizedExperiment::assay(x, "dosage"))
    )
}

# Narrow a sketch to a variant subset. On an RSE this is ordinary subsetting,
# which brings the DelayedArray's index bookkeeping with it; a legacy bare
# handle falls back to the handle-level trim, which needs fileIdx to stay
# read-safe.
# @noRd
.ldSketchSubset <- function(x, keep) {
    if (is.null(x)) {
        return(NULL)
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(.subsetGenotypeHandle(x, keep))
    }
    x[keep, ]
}

# A sketch trimmed to zero variants, every other property preserved. Used when
# an entry set genuinely carries no variants: the object references no LD, so
# the correct retained panel is empty rather than the full genome-wide sketch,
# which on a real panel is ~135 MB of serialized dead weight per such entry.
# NULL-safe and idempotent.
# @noRd
.emptySketch <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(.emptyGenotypeHandle(x))
    }
    # A panel holds its handle TWICE -- as rowRanges and inside the dosage
    # assay's DelayedArray seed. `x[integer(0), ]` drops the rows but leaves
    # the seed carrying the whole panel's snpInfo, which is precisely the
    # weight this exists to shed. Rebuild from an emptied handle so both go.
    .genotypeExperiment(.emptyGenotypeHandle(.ldSketchHandle(x)))
}

# "format @ path" for show methods, which describe where a panel came from
# rather than what is in it.
# @noRd
.ldSketchLabel <- function(x) {
    handle <- .ldSketchHandle(x)
    glue("{getFormat(handle)} @ {getPath(handle)}")
}

# TRUE for anything that can serve as an LD panel. A bare GenotypeHandle
# still counts: that is what objects serialised before the slot became an RSE
# carry.
# @noRd
.ldIsPanel <- function(x) {
    methods::is(x, "RangedSummarizedExperiment") ||
        methods::is(x, "GenotypeHandle")
}

# The handle behind an already-open genotype source, or NULL when the value is
# not one. Callers that accept "a panel, a bare handle, or a spec to resolve"
# use this to separate the first two -- which are ready to use -- from the
# third. A panel is the public shape (readGenotypes()); a bare handle is what
# internal callers and objects serialised before the change still carry.
# @noRd
.openGenotypeHandle <- function(x) {
    if (!.ldIsPanel(x)) {
        return(NULL)
    }
    .ldSketchHandle(x)
}

# ---------- LD sketch: the panel's own view of itself ----------
#
# These read the panel through the SummarizedExperiment interface rather than
# reaching past it to the GenotypeHandle. The handle is the DelayedArray's
# seed -- the READ PATH -- and unwrapping it for metadata leaks that seed into
# callers that only ever wanted chrom/pos/alleles or a sample count, which the
# RSE already answers. Reaching for `.ldSketchHandle()` is now reserved for
# what genuinely lives on the handle: file format/path, chromosome shard
# routing, the `.afreq` sidecar, and building or pruning the seed itself.
#
# A bare GenotypeHandle -- what objects serialised before the slot became an
# RSE carry -- is projected exactly the way `.genotypeExperiment()` projects
# one, so both shapes answer identically.

# The panel's variant ranges: seqnames = chromosome (chr-prefixed), start =
# position, mcols = SNP / A1 / A2.
# @noRd
.ldSketchRanges <- function(x, label = "LD sketch") {
    # The chokepoint every panel accessor funnels through, so the guard lives
    # here rather than being repeated at each entry point. Without it a NULL
    # or non-panel sketch surfaces as "unable to find an inherited method for
    # 'getSnpInfo'" from whichever accessor happened to touch it first, which
    # says nothing about the LD reference being the problem. Callers that
    # validate with their own label (`.ldFromSketch`, the ctwas assembler) do
    # so first, so their message wins.
    .ldFromSketchValidate(x, label)
    if (methods::is(x, "RangedSummarizedExperiment")) {
        return(SummarizedExperiment::rowRanges(x))
    }
    .genotypeSnpRanges(x, normalizeVariantId(as.character(getSnpInfo(x)$SNP)))
}

# The panel's variant ids in the panel's OWN labelling. This is the mcols
# column, not `names()`: the names are normalized ids, while anything keyed on
# the panel FILE -- the `.afreq` sidecar, the ctwas snp_map -- has to see the
# panel's raw labels.
# @noRd
.ldSketchVariantIds <- function(x) {
    as.character(S4Vectors::mcols(.ldSketchRanges(x))$SNP)
}

# The panel's variant ids as MATCHING keys: the raw labels, with any id whose
# allele slot holds a tag rather than an allele (chr21:13988152:INS:T) re-
# rendered from the panel's own A1 / A2 columns, which are correct. Use this
# wherever a panel id is matched against summary statistics or another id set,
# and `.ldSketchVariantIds()` only where the file's literal label is wanted.
#
# The repair is what lets an untidily-named panel entry be MATCHED rather than
# discarded: harmonization re-keys a surviving variant to the reference-allele
# id, so the panel side has to speak the same form or every later id lookup
# misses the variant it just harmonized.
# @noRd
.ldSketchMatchIds <- function(x) {
    mc <- S4Vectors::mcols(.ldSketchRanges(x))
    ids <- as.character(mc$SNP)
    if (is.null(mc$A1) || is.null(mc$A2)) {
        return(ids)
    }
    .repairVariantIds(ids, as.character(mc$A2), as.character(mc$A1))
}

# Per-variant chromosome, canonical (no chr prefix).
# @noRd
.ldSketchChrom <- function(x) {
    canonChrom(as.character(GenomicRanges::seqnames(.ldSketchRanges(x))))
}

# @noRd
.ldSketchSampleIds <- function(x) {
    if (methods::is(x, "RangedSummarizedExperiment")) {
        return(as.character(colnames(x)))
    }
    as.character(getSampleIds(x))
}

# @noRd
.ldSketchNSamples <- function(x) {
    length(.ldSketchSampleIds(x))
}

# Dosage for a variant subset, samples x variants.
#
# Read through the assay, so the DelayedArray seed -- and with it the file
# ordering, the sharded routing and the empty-request guard in
# `extractBlockGenotypes()` -- stays the single reader. `extract_array()`
# delegates to that same function with `meanImpute = FALSE`, so imputation is
# layered on here rather than pushed down.
# @noRd
.ldSketchDosage <- function(x, snpIdx, meanImpute = TRUE) {
    # The other primitive that does not pass through `.ldSketchRanges()`.
    .ldFromSketchValidate(x, "LD sketch")
    if (!methods::is(x, "RangedSummarizedExperiment")) {
        return(.dosageMatrix(x, snpIdx, meanImpute = meanImpute))
    }
    dosage <- t(as.matrix(
        SummarizedExperiment::assay(x, "dosage")[snpIdx, , drop = FALSE]
    ))
    if (meanImpute) .qtlMeanImpute(dosage) else dosage
}

.ldFromSketchValidate <- function(ldSketch, label) {
    if (is.null(ldSketch)) {
        msg <- glue(
            "{label}: the SumStats/collection carries no ldSketch ",
            "(ldSketch = NULL); this step needs an LD reference."
        )
        abort(msg)
    }
    # A bare GenotypeHandle is still accepted: that is what objects
    # serialised before the slot became an RSE carry.
    if (!.ldIsPanel(ldSketch)) {
        msg <- glue(
            "{label}: ldSketch must be a genotype panel ",
            "(RangedSummarizedExperiment or GenotypeHandle)."
        )
        abort(msg)
    }
}

# Match requested ids to the panel; NULL if none match. Returns kept ids/order.
.ldFromSketchMatch <- function(ldSketch, variantIds, label, onMissing) {
    # Match by (chrom, pos, allele) tuple with an exact id-string fallback for
    # rsID panels; the caller's original ids and order are preserved.
    m <- matchVariants(
        variantIds,
        .ldSketchMatchIds(ldSketch),
        removeStrandAmbiguous = FALSE
    )
    nMissing <- length(variantIds) - length(m$idxA)
    if (nMissing > 0L && onMissing == "error") {
        msg <- glue(
            "{label}: {nMissing} variant id(s) not present in the LD ",
            "sketch panel."
        )
        abort(msg)
    }
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    o <- order(m$idxA) # restore the caller's requested order
    list(
        keptIds = variantIds[m$idxA[o]],
        idx = m$idxB[o],
        sign = m$sign[o]
    )
}

.ldFromSketch <- function(
    ldSketch,
    variantIds,
    label = ".ldFromSketch",
    onMissing = c("error", "drop")
) {
    .ldFromSketchValidate(ldSketch, label)
    onMissing <- arg_match(onMissing)
    matched <- .ldFromSketchMatch(ldSketch, variantIds, label, onMissing)
    if (is.null(matched)) {
        return(NULL)
    }
    ldMat <- computeLd(
        ldSketch,
        method = "sample",
        snpIdx = matched$idx
    )
    # The match above is allele-aware, so a panel entry whose alleles are
    # swapped relative to the caller's id still matches -- but the dosage it
    # returns counts the OTHER allele, which negates every correlation that
    # variant takes part in (r(2 - x, y) = -r(x, y)). Put the matrix back in
    # the caller's frame exactly as `.cbFlipPairToCanonical()` does:
    # LD_ij -> sign_i * sign_j * LD_ij. A no-op on a harmonized panel, where
    # every sign is +1.
    signed <- if (any(matched$sign < 0)) {
        ldMat * outer(matched$sign, matched$sign)
    } else {
        ldMat
    }
    named <- `dimnames<-`(signed, list(matched$keptIds, matched$keptIds))
    if (onMissing == "drop") {
        return(`attr<-`(named, "keptVariantIds", matched$keptIds))
    }
    named
}

# ---------- LD sketch: per-variant panel statistics and filtering ----------

# Per-variant minor allele frequency and missingness rate from panel dosage.
#
# The dosage must NOT be mean-imputed: imputation fills every hole, so the
# missingness rate would come back zero everywhere.
# @noRd
.panelVariantStats <- function(dosage) {
    nSamp <- nrow(dosage)
    nObs <- colSums(!is.na(dosage))
    af <- if_else(
        nObs > 0L,
        colSums(dosage, na.rm = TRUE) / (2 * nObs),
        NA_real_
    )
    list(
        af = af,
        maf = pmin(af, 1 - af),
        missRate = if (nSamp > 0L) {
            1 - nObs / nSamp
        } else {
            rep(0, ncol(dosage))
        }
    )
}

# A MAC cutoff is a MAF cutoff once expressed per panel sample, so the
# stricter of the two applies -- the same rule .qtlVariantFilters() uses on
# individual level data, so one number means the same thing on both paths.
# @noRd
.panelEffectiveMaf <- function(mafCutoff, macCutoff, nSamples) {
    macAsMaf <- if (nSamples > 0L) macCutoff / (2 * nSamples) else 0
    max(mafCutoff, macAsMaf)
}

# The .afreq stems to read for a panel: the handle's own stem, or one per
# chromosome the panel actually spans. A sharded handle lists every
# chromosome in its manifest, and reading the shards the summary statistics
# never touch is the cost per-chromosome sketch trimming exists to avoid.
# @noRd
.panelAfreqPrefixes <- function(handle) {
    chromPaths <- .genotypeChromPaths(handle)
    if (length(chromPaths) == 0L) {
        return(.genotypeReadPath(handle))
    }
    spanned <- unique(canonChrom(as.character(getSnpInfo(handle)$CHR)))
    paths <- unname(chromPaths[is_in(names(chromPaths), spanned)])
    map_chr(paths, .resolveGenotypeResourcePath)
}

# One shard's id/alt_freq pairs, or NULL when the sidecar is missing or does
# not carry them. A malformed sidecar must not take the whole filter down, so
# it degrades to "no frequencies here" and the caller falls back to dosage.
# @noRd
#' @importFrom rlang try_fetch
.panelAfreqTable <- function(prefix) {
    af <- try_fetch(readAfreq(prefix), error = function(cnd) NULL)
    if (is.null(af) || !all(is_in(c("id", "alt_freq"), colnames(af)))) {
        return(NULL)
    }
    select(af, all_of(c("id", "alt_freq")))
}

# Panel minor allele frequency for `variantIds` read from the PLINK2 .afreq
# sidecar, or NULL when the sidecar cannot answer for every id.
#
# Refusing a PARTIAL sidecar is deliberate: the dosage path treats a variant
# it cannot measure as a drop, so a fast path that quietly kept the ids
# .afreq happens to omit would disagree with it about the same variant.
# @noRd
.panelAfreqMaf <- function(handle, variantIds) {
    if (getFormat(handle) != "plink2") {
        return(NULL)
    }
    tbls <- compact(map(.panelAfreqPrefixes(handle), .panelAfreqTable))
    if (length(tbls) == 0L) {
        return(NULL)
    }
    afTbl <- list_rbind(tbls)
    idx <- match(variantIds, pull(afTbl, "id"))
    if (anyNA(idx)) {
        return(NULL)
    }
    altFreq <- as.numeric(pull(afTbl, "alt_freq"))[idx]
    pmin(altFreq, 1 - altFreq)
}

# TRUE for each matched panel variant that fails the cutoffs.
#
# Allele frequency alone answers MAF and MAC, so with no missingness cutoff
# set the .afreq sidecar decides the whole mask without materializing panel
# dosage -- which is what makes the filter affordable on a panel far larger
# than the analysed variant set. Missingness is not derivable from the
# sidecar (its observation count is an ALLELE count, whose ploidy the mask
# would have to guess), so an imissCutoff always takes the dosage path.
# @noRd
.panelDropMask <- function(ldSketch, matched, panelFilterArgs) {
    imissCutoff <- panelFilterArgs$imissCutoff
    if (imissCutoff >= 1) {
        # Keyed on the PANEL's own variant labels, not the caller's: .afreq
        # carries the .pvar ids, while `keptIds` is whatever id form the
        # caller asked in (normalized, most often).
        panelIds <- .ldSketchVariantIds(ldSketch)[matched$idx]
        # The handle stays for this one read: the sidecar is a file beside
        # the genotypes, which only the seed knows how to find.
        maf <- .panelAfreqMaf(.ldSketchHandle(ldSketch), panelIds)
        if (!is.null(maf)) {
            effMaf <- .panelEffectiveMaf(
                panelFilterArgs$mafCutoff,
                panelFilterArgs$macCutoff,
                .ldSketchNSamples(ldSketch)
            )
            return(is.na(maf) | maf < effMaf)
        }
    }
    dosage <- .ldSketchDosage(ldSketch, matched$idx, meanImpute = FALSE)
    stats <- .panelVariantStats(dosage)
    effMaf <- .panelEffectiveMaf(
        panelFilterArgs$mafCutoff,
        panelFilterArgs$macCutoff,
        nrow(dosage)
    )
    is.na(stats$maf) | stats$maf < effMaf | stats$missRate > imissCutoff
}

# The subset of `variantIds` whose LD-panel genotypes clear the cutoffs,
# returned in the caller's order.
#
# Variants absent from the panel are passed through untouched. Whether a
# missing variant is an error or is silently dropped belongs to `.ldFromSketch`
# (its `onMissing`), and deciding it a second time here would let the two
# disagree about the same variant.
# @noRd
.panelVariantFilter <- function(
    ldSketch,
    variantIds,
    panelFilterArgs = panelFilterConfig(),
    label = ".panelVariantFilter"
) {
    # .panelCutoffs answers NULL for a filter that would keep everything,
    # which is also the cheap exit here.
    cutoffs <- .panelCutoffs(panelFilterArgs)
    if (is.null(cutoffs)) {
        return(variantIds)
    }
    if (is.null(ldSketch) || length(variantIds) == 0L) {
        return(variantIds)
    }
    matched <- .ldFromSketchMatch(ldSketch, variantIds, label, "drop")
    if (is.null(matched)) {
        return(variantIds)
    }
    drop <- .panelDropMask(ldSketch, matched, cutoffs)
    variantIds[!is_in(variantIds, matched$keptIds[drop])]
}

# Row mask for the variants clearing the LD-panel cutoffs, reporting how many
# were dropped. A no-op when no cutoff is set.
# @noRd
.panelKeepMask <- function(variantIds, ldSketch, cutoffs, label) {
    if (is.null(ldSketch) || is.null(cutoffs)) {
        return(rep(TRUE, length(variantIds)))
    }
    kept <- .panelVariantFilter(ldSketch, variantIds, cutoffs, label = label)
    keep <- is_in(variantIds, kept)
    nDropped <- sum(!keep)
    if (nDropped > 0L) {
        msg <- glue(
            "{label}: dropped {nDropped} of {length(variantIds)} variant(s) ",
            "below the LD-panel MAF / MAC / missingness cutoffs."
        )
        inform(msg)
    }
    keep
}

# The panel-filter cutoffs a pipeline call carries, or NULL when none is set
# (so the filter short-circuits without touching the panel).
# @noRd
.panelCutoffs <- function(panelFilterArgs = panelFilterConfig()) {
    # NULL fields are how the pipelines spell "not set"; normalise them to the
    # no-op values so the short-circuit below is the only place that decides
    # whether a filter is worth running.
    maf <- panelFilterArgs$mafCutoff %||% 0
    mac <- panelFilterArgs$macCutoff %||% 0
    imiss <- panelFilterArgs$imissCutoff %||% 1
    if (maf <= 0 && mac <= 0 && imiss >= 1) {
        return(NULL)
    }
    list(mafCutoff = maf, macCutoff = mac, imissCutoff = imiss)
}
# ---------- LD sketch: cross-pipeline LD-panel compatibility check ----------

# Internal: assert that two LD sketches describe the same reference panel --
# the same samples, and the same allele coding (A1/A2 in the same orientation)
# on the variants they share. The SNP label is not compared, so a pure
# chr-prefix difference does not fail; an allele swap still does, since it
# means the two panels code their LD in opposite directions.
#
# The two panels' variant SETS need not agree, and after independent QC they
# normally do not: each object's sketch is trimmed to its own span and its own
# surviving variants, so one LD reference used for a QTL and a GWAS collection
# yields two differently trimmed sketches. Downstream LD lookups are all
# id-matched, so a partial overlap is reported (once per session) rather than
# refused; no overlap at all is an error. Shared by causalInferencePipeline,
# colocPipeline, qtlEnrichmentPipeline, ctwasPipeline, and colocboostPipeline.
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
        labelPart <- if (!is.null(label)) {
            glue("ldSketch on `{label}` is non-NULL ")
        } else {
            "qtl ldSketch is non-NULL "
        }
        msg <- glue(
            "{pipelineName}: {labelPart}but the GWAS ldSketch is NULL."
        )
        abort(msg)
    }
    FALSE
}

# Both must be genotype panels. Panel SIZE is deliberately NOT compared: two
# objects that share one LD sketch stop carrying identical variant sets as soon
# as they are QC'd separately, because each object's sketch is trimmed to its
# own position span at load (`.subsetSketchToRange`) and to its own surviving
# variants at the end of QC (`.subsetSketchToIds`). A QTL and a GWAS collection
# put through `summaryStatsQc` independently therefore diverge by construction,
# which is the normal case rather than a mistake. Every downstream LD lookup
# matches by variant id (`.ldFromSketch`), never by position in the panel, so
# divergent variant sets cost only the variants one side lacks. What would not
# be harmless -- two genuinely different reference panels -- is what
# `.ldSketchCheckContent` catches, by sample set and by allele coding.
.ldSketchCheckShape <- function(qtlLd, gwasLd, pipelineName, between) {
    if (!.ldIsPanel(qtlLd) || !.ldIsPanel(gwasLd)) {
        msg <- glue(
            "{pipelineName}: ldSketch slots{between} must both be genotype ",
            "panels for the cross-pipeline LD reference check."
        )
        abort(msg)
    }
}

# The (chrom, position, allele-pair) key the two panels are compared on. The
# allele pair is order-insensitive, so an A1/A2 swap keys to the SAME variant
# and is then reported as an allele-coding difference rather than looking like
# two unrelated variants that happen not to overlap. A panel carrying no allele
# columns falls back to chrom:position, which is all it can be keyed on.
# @noRd
.ldSketchVariantKeys <- function(x) {
    gr <- .ldSketchRanges(x)
    mc <- S4Vectors::mcols(gr)
    stem <- str_c(
        canonChrom(as.character(GenomicRanges::seqnames(gr))),
        ":",
        as.character(GenomicRanges::start(gr))
    )
    if (is.null(mc$A1) || is.null(mc$A2)) {
        return(stem)
    }
    # str_c propagates NA, which would make a variant with an unknown allele
    # match nothing on either side; an empty field keeps it comparable.
    rawA1 <- as.character(mc$A1)
    rawA2 <- as.character(mc$A2)
    a1 <- if_else(is.na(rawA1), "", rawA1)
    a2 <- if_else(is.na(rawA2), "", rawA2)
    str_c(stem, ":", if_else(a1 < a2, a1, a2), ":", if_else(a1 < a2, a2, a1))
}

# The panel's A1 (effect / counted) allele, which fixes the sign of every
# correlation the variant takes part in.
# @noRd
.ldSketchA1 <- function(x) {
    as.character(S4Vectors::mcols(.ldSketchRanges(x))$A1)
}

# A shared variant whose A1/A2 are swapped between the panels is the same
# variant coded in opposite directions, so the two panels' LD matrices disagree
# in sign wherever it appears. That is a different reference, not a trimming
# difference, and it is an error.
# @noRd
.ldSketchCheckAlleleCoding <- function(
    qtlLd,
    gwasLd,
    qIdx,
    gIdx,
    pipelineName,
    between
) {
    qA1 <- .ldSketchA1(qtlLd)[qIdx]
    gA1 <- .ldSketchA1(gwasLd)[gIdx]
    if (length(qA1) == 0L || length(gA1) == 0L) {
        return(invisible(NULL))
    }
    nSwapped <- sum(qA1 != gA1, na.rm = TRUE)
    if (nSwapped > 0L) {
        msg <- glue(
            "{pipelineName}: {nSwapped} variant(s) present in both ldSketch ",
            "panels carry swapped A1/A2 alleles{between}, so the two panels ",
            "code their LD in opposite directions; use the same ldSketch on ",
            "both."
        )
        abort(msg)
    }
    invisible(NULL)
}

# A partial overlap is the expected result of QC-ing the two sides separately,
# so it is reported once per session rather than on every call: ctwasPipeline
# and colocboostPipeline run this check once per region / per bundle.
# @noRd
.ldSketchReportOverlap <- function(nQ, nG, nShared, pipelineName, between) {
    if (nShared == nQ && nShared == nG) {
        return(invisible(NULL))
    }
    msg <- glue(
        "{pipelineName}: the two ldSketch panels share {nShared} variant(s) ",
        "of {nQ} (QTL side) and {nG} (GWAS side){between}. LD is looked up ",
        "per variant, so only the shared ones contribute. Differing variant ",
        "sets are expected when the two sides were QC'd separately."
    )
    warn(
        msg,
        .frequency = "once",
        .frequency_id = str_c("pecotmrLdSketchOverlap-", pipelineName)
    )
    invisible(NULL)
}

# Compare the two panels variant by variant. A partial overlap passes; what
# fails is no overlap at all (two unrelated panels) or a shared variant whose
# alleles are swapped.
# @noRd
.ldSketchCheckOverlap <- function(qtlLd, gwasLd, pipelineName, between) {
    qKey <- .ldSketchVariantKeys(qtlLd)
    gKey <- .ldSketchVariantKeys(gwasLd)
    # An emptied panel -- what a zero-variant object carries, since it
    # references no LD -- contradicts nothing, so there is nothing to compare.
    if (length(qKey) == 0L || length(gKey) == 0L) {
        return(invisible(NULL))
    }
    idx <- match(qKey, gKey)
    shared <- which(!is.na(idx))
    if (length(shared) == 0L) {
        msg <- glue(
            "{pipelineName}: the two ldSketch panels share no variant",
            "{between} ({length(qKey)} vs {length(gKey)} variants); they ",
            "describe different LD references."
        )
        abort(msg)
    }
    .ldSketchCheckAlleleCoding(
        qtlLd,
        gwasLd,
        shared,
        idx[shared],
        pipelineName,
        between
    )
    .ldSketchReportOverlap(
        length(qKey),
        length(gKey),
        length(shared),
        pipelineName,
        between
    )
}

# Panels must describe the same reference panel: the same samples, and the same
# allele coding on the variants they share. Their variant SETS need not agree
# -- see `.ldSketchCheckShape` for why.
.ldSketchCheckContent <- function(qtlLd, gwasLd, pipelineName, between) {
    qIds <- .ldSketchSampleIds(qtlLd)
    gIds <- .ldSketchSampleIds(gwasLd)
    if (!identical(qIds, gIds)) {
        msg <- glue(
            "{pipelineName}: ldSketch panels have different sample sets",
            "{between}; use the same ldSketch on both."
        )
        abort(msg)
    }
    .ldSketchCheckOverlap(qtlLd, gwasLd, pipelineName, between)
}

.requireMatchingLdSketches <- function(
    qtlLd,
    gwasLd,
    pipelineName,
    label = NULL,
    nullPolicy = c("qtl-required", "lenient")
) {
    nullPolicy <- arg_match(nullPolicy)
    if (.ldSketchNullGuard(qtlLd, gwasLd, pipelineName, label, nullPolicy)) {
        return(invisible(NULL))
    }
    between <- if (!is.null(label)) {
        glue(" between `{label}` and gwas inputs", .trim = FALSE)
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


# ---------- Internal: load LD from pre-computed blocks ----------

#' Load pre-computed LD from block-based metadata files.
#' @importFrom purrr map_chr map_int map_dbl
#' @noRd
# --- loadLdFromBlocks helpers -----------------------------------------------

# One block's record: the region-extracted matrix, its variants, and the
# chromosome it sits on (falling back to the requested region's when the
# block contributed nothing).
#
# A RECORD rather than a position in three parallel lists: the filter and the
# metadata builder below would otherwise have to keep four sequences in the
# same order by hand, and a slip there would silently attach one block's
# chromosome to another block's variants.
# @noRd
.loadLdOneBlock <- function(
    j,
    ldFilePaths,
    bimFilePaths,
    intersectedLdFiles,
    extractCoordinates
) {
    proc <- processLdMatrix(ldFilePaths[[j]], bimFilePaths[[j]])
    extracted <- extractLdForRegion(
        ldMatrix = proc$ldMatrix,
        variants = proc$ldVariants,
        region = intersectedLdFiles$region,
        extractCoordinates = extractCoordinates
    )
    blockVariants <- extracted$extractedLdVariants
    list(
        matrix = extracted$extractedLdMatrix,
        variants = blockVariants,
        chrom = if (nrow(blockVariants) > 0) {
            as.character(blockVariants$chrom[[1L]])
        } else {
            as.character(intersectedLdFiles$region$chrom)
        }
    )
}

# Load + region-extract each LD block, one record per block.
.loadLdBlocksLoop <- function(
    ldFilePaths,
    bimFilePaths,
    intersectedLdFiles,
    extractCoordinates
) {
    map(
        seq_along(ldFilePaths),
        .loadLdOneBlock,
        ldFilePaths = ldFilePaths,
        bimFilePaths = bimFilePaths,
        intersectedLdFiles = intersectedLdFiles,
        extractCoordinates = extractCoordinates
    )
}

# Drop blocks with no variants in the region (error if none remain).
.loadLdFilterEmpty <- function(blocks) {
    nonEmpty <- map_lgl(map(blocks, "variants"), .ldBlockHasVariants)
    if (!any(nonEmpty)) {
        abort("No variants found in any LD block for the specified region.")
    }
    if (any(!nonEmpty)) {
        nEmpty <- sum(!nonEmpty)
        msg <- glue(
            "Removing {nEmpty} empty LD block(s) with no variants in ",
            "the region."
        )
        inform(msg)
    }
    blocks[nonEmpty]
}

# Per-block metadata (id, chrom, span, size, index range in the merged matrix).
.loadLdBlockMetadata <- function(blocks, ldVariants) {
    variants <- map(blocks, "variants")
    blockVariants <- map(variants, "variants")
    blockPositions <- map(variants, "pos")
    tibble(
        blockId = seq_along(blocks),
        chrom = map_chr(blocks, "chrom"),
        blockStart = map_dbl(blockPositions, min),
        blockEnd = map_dbl(blockPositions, max),
        size = map_int(blockVariants, length),
        startIdx = map_int(
            blockVariants,
            .ldBlockStartIdx,
            ldVariants = ldVariants
        ),
        endIdx = map_int(
            blockVariants,
            .ldBlockEndIdx,
            ldVariants = ldVariants
        )
    )
}

# Build the reference panel: variant ids + merged per-variant annotations,
# deriving variance from nSample + allele_freq when it is otherwise absent.
# One annotation column aligned to the panel's variant order.
# @noRd
.ldAnnotationColumn <- function(col, mergedVariantList, ids) {
    mergedVariantList[[col]][match(ids, mergedVariantList$variants)]
}

.loadLdRefPanel <- function(ldMatrix, extractedLdVariantsList, nSample) {
    ids <- rownames(ldMatrix)
    mergedVariantList <- bind_rows(extractedLdVariantsList)
    annotations <- intersect(
        c("allele_freq", "variance", "n_nomiss"),
        colnames(mergedVariantList)
    )
    panel <- mutate(
        parseVariantId(ids),
        variant_id = ids,
        !!!set_names(
            map(
                annotations,
                .ldAnnotationColumn,
                mergedVariantList = mergedVariantList,
                ids = ids
            ),
            annotations
        )
    )
    needVar <- !is_in("variance", colnames(panel)) || all(is.na(panel$variance))
    if (
        is.null(nSample) ||
            !needVar ||
            !is_in("allele_freq", colnames(panel))
    ) {
        return(panel)
    }
    p <- panel$allele_freq
    mutate(
        panel,
        variance = 2 * p * (1 - p) * nSample / (nSample - 1),
        n_nomiss = nSample
    )
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
    kept <- .loadLdFilterEmpty(blocks)
    keptVariants <- map(kept, "variants")
    ldMatrix <- createLdMatrix(
        ldMatrices = map(kept, "matrix"),
        variants = keptVariants
    )
    ldVariants <- rownames(ldMatrix)
    blockMetadata <- .loadLdBlockMetadata(kept, ldVariants)
    refPanel <- .loadLdRefPanel(ldMatrix, keptVariants, nSample)
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
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' filterVariantsByLdReference(
#'   variantIds = c("chr22:16050000:A:G", "chr22:17000000:C:T"),
#'   ldReferenceMetaFile = meta)
#' @export
#' @importFrom checkmate assertCharacter assertFlag
filterVariantsByLdReference <- function(
    variantIds,
    ldReferenceMetaFile,
    keepIndel = TRUE
) {
    assertCharacter(variantIds, any.missing = FALSE)
    assertFlag(keepIndel)
    variantsDf <- parseVariantId(variantIds)

    # Derive region to scope the reference lookup
    regionDf <- variantsDf |>
        group_by(.data$chrom) |>
        summarise(start = min(.data$pos), end = max(.data$pos))

    # Use shared helper -- no genotype loading
    refInfo <- getRefVariantInfo(ldReferenceMetaFile, regionDf)
    refChrom <- canonChrom(refInfo$chrom)
    refKey <- str_c(refChrom, ":", refInfo$pos)

    variantKey <- str_c(variantsDf$chrom, ":", variantsDf$pos)
    inRef <- which(is_in(variantKey, refKey))
    keepIndices <- if (keepIndel) {
        inRef
    } else {
        intersect(inRef, which(isSnpAlleles(variantsDf$A1, variantsDf$A2)))
    }

    nDropped <- length(variantIds) - length(keepIndices)
    if (nDropped > 0) {
        nTotal <- length(variantIds)
        msg <- glue(
            "{nDropped} out of {nTotal} total variants dropped due to ",
            "absence on the reference LD panel."
        )
        inform(msg)
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
        abort("Empty or NULL LD matrix provided.")
    }
    wanted <- list(variantIds, variantIds)
    if (!identical(dimnames(combinedMatrix), wanted)) {
        return(`dimnames<-`(combinedMatrix, wanted))
    }
    combinedMatrix
}

# Drop blocks with invalid/out-of-range indices; renumber the survivors.
.partitionFilterBlocks <- function(blockMetadata, nVariants) {
    validBlocks <- map_lgl(
        seq_len(nrow(blockMetadata)),
        .ldBlockValid,
        blockMetadata = blockMetadata,
        nVariants = nVariants
    )
    if (!any(validBlocks)) {
        msg <- glue(
            "No valid LD blocks found. All block indices are out of ",
            "range or empty."
        )
        abort(msg)
    }
    if (any(!validBlocks)) {
        nInvalid <- sum(!validBlocks)
        msg <- glue(
            "Removing {nInvalid} LD block(s) with invalid or ",
            "out-of-range indices."
        )
        inform(msg)
        kept <- filter(blockMetadata, validBlocks)
        return(mutate(kept, blockId = seq_len(nrow(kept))))
    }
    blockMetadata
}

#' @importFrom checkmate assertClass
partitionLdMatrix <- function(
    ldData,
    mergeSmallBlocks = TRUE,
    minMergedBlockSize = 500,
    maxMergedBlockSize = 10000
) {
    assertClass(ldData, "LdData")
    variantIds <- getVariantIds(ldData)
    combinedMatrix <- .partitionValidateMatrix(
        getCorrelation(ldData),
        variantIds
    )
    rawBlocks <- getBlockMetadata(ldData)
    blocks <- .partitionFilterBlocks(
        if (is(rawBlocks, "GRanges")) as_tibble(rawBlocks) else rawBlocks,
        length(variantIds)
    )
    # Validate the block structure of the matrix (skip if only one block).
    if (nrow(blocks) > 1) {
        validateBlockStructure(combinedMatrix, blocks, variantIds)
    }
    merging <- mergeSmallBlocks &&
        any(blocks$size < minMergedBlockSize) &&
        nrow(blocks) > 1
    blockMetadata <- if (merging) {
        mergeBlocks(blocks, minMergedBlockSize, maxMergedBlockSize)
    } else {
        blocks
    }
    extractBlockMatrices(combinedMatrix, blockMetadata, variantIds)
}

# Every (i, j) block pair with i < j.
# @noRd
.ldUpperPairs <- function(nBlocks) {
    grid <- expand.grid(i = seq_len(nBlocks), j = seq_len(nBlocks))
    grid[grid$i < grid$j, , drop = FALSE]
}

# @noRd
.blockPairMessagesAt <- function(
    k,
    pairs,
    blockMetadata,
    matrix,
    variantIds,
    n
) {
    .blockPairMessages(
        pairs$i[[k]],
        pairs$j[[k]],
        blockMetadata,
        matrix,
        variantIds,
        n
    )
}

# Concatenate per-item message vectors, empty-safe.
# @noRd
.ldConcatChr <- function(pieces) {
    if (length(pieces) == 0L) {
        return(character(0))
    }
    as.character(list_c(pieces))
}

#' Validate that cross-block entries are zero (excluding boundary variants).
#' @noRd
validateBlockStructure <- function(matrix, blockMetadata, variantIds) {
    n <- length(variantIds)
    pairs <- .ldUpperPairs(nrow(blockMetadata))
    msgs <- .ldConcatChr(map(
        seq_len(nrow(pairs)),
        .blockPairMessagesAt,
        pairs = pairs,
        blockMetadata = blockMetadata,
        matrix = matrix,
        variantIds = variantIds,
        n = n
    ))
    if (length(msgs) > 0) {
        msgList <- str_flatten(msgs, collapse = "\n")
        msg <- glue(
            "Matrix lacks expected block structure:\n{msgList}",
            .trim = FALSE
        )
        abort(msg)
    }
}

# Overlap-consistency messages for one block pair (i, j): out-of-range indices,
# or non-zero cross-block correlation with boundary variants excluded. Returns a
# (possibly empty) character vector.
# @noRd
.blockPairMessages <- function(i, j, blockMetadata, matrix, variantIds, n) {
    si <- blockMetadata$startIdx[i]
    ei <- blockMetadata$endIdx[i]
    sj <- blockMetadata$startIdx[j]
    ej <- blockMetadata$endIdx[j]
    if (si > n || ei > n || sj > n || ej > n) {
        return(str_c(
            "Block indices out of range for blocks",
            i,
            "and",
            j,
            sep = " "
        ))
    }
    # Exclude boundary variants (potential overlaps)
    vi <- variantIds[si:(ei - 1)]
    vj <- variantIds[(sj + 1):ej]
    maxVal <- max(abs(matrix[vi, vj, drop = FALSE]))
    if (maxVal <= 1e-10) {
        return(character(0))
    }
    str_c(
        "Non-zero correlation between blocks",
        i,
        "and",
        j,
        "- max:",
        maxVal,
        sep = " "
    )
}

#' Can blocks `i` and `j` of `blockMetadata` merge (same chrom, combined size
#' within `maxSize`)? Indexes the columns directly rather than slicing rows.
#' @noRd
canMerge <- function(blockMetadata, i, j, maxSize) {
    blockMetadata$chrom[i] == blockMetadata$chrom[j] &&
        (blockMetadata$size[i] + blockMetadata$size[j]) <= maxSize
}

#' @noRd
mergeTwoBlocks <- function(blockMetadata, idx1, idx2) {
    if (idx1 > idx2) {
        tmp <- idx1
        idx1 <- idx2
        idx2 <- tmp
    }
    merged <- mutate(
        blockMetadata,
        endIdx = replace(.data$endIdx, idx1, blockMetadata$endIdx[idx2]),
        size = replace(
            .data$size,
            idx1,
            blockMetadata$size[idx1] + blockMetadata$size[idx2]
        )
    ) |>
        slice(-idx2)
    mutate(merged, blockId = seq_len(nrow(merged)))
}

#' Find blocks below minSize and identify the best neighbor to merge with.
#' @noRd
findMergeCandidates <- function(blockMetadata, minSize, maxSize) {
    found <- compact(map(
        seq_len(nrow(blockMetadata)),
        .ldMergeCandidateFor,
        blockMetadata = blockMetadata,
        minSize = minSize,
        maxSize = maxSize
    ))
    if (length(found) == 0L) {
        return(tibble(block_idx = integer(), merge_with = integer()))
    }
    bind_rows(found)
}

# Which neighbour block `i` should merge into (the smaller one when both
# qualify), or NULL when it is big enough already or neither neighbour can
# take it.
# @noRd
.ldMergeCandidateFor <- function(i, blockMetadata, minSize, maxSize) {
    if (blockMetadata$size[i] >= minSize) {
        return(NULL)
    }
    prevOk <- i > 1 && canMerge(blockMetadata, i, i - 1, maxSize)
    nextOk <- i < nrow(blockMetadata) &&
        canMerge(blockMetadata, i, i + 1, maxSize)
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
        return(NULL)
    }
    tibble(block_idx = i, merge_with = mergeWith)
}

#' Iteratively merge blocks below minSize with their smallest neighbor.
#' @noRd
mergeBlocks <- function(blockMetadata, minSize, maxSize) {
    if (nrow(blockMetadata) <= 1) {
        return(blockMetadata)
    }
    candidates <- findMergeCandidates(blockMetadata, minSize, maxSize)
    if (nrow(candidates) == 0) {
        return(blockMetadata)
    }
    # Each merge changes which blocks are still too small, so this is a fixed
    # point: recurse on the merged metadata rather than rebinding it.
    mergeBlocks(
        mergeTwoBlocks(
            blockMetadata,
            candidates$block_idx[1],
            candidates$merge_with[1]
        ),
        minSize,
        maxSize
    )
}

# Helper function to extract block matrices
# Extract one block's submatrix + mapping (NULL to skip empty/OOB blocks).
.extractOneBlock <- function(matrix, variantIds, startIdx, endIdx, i) {
    if (endIdx < startIdx) {
        return(NULL)
    }
    if (startIdx > length(variantIds) || endIdx > length(variantIds)) {
        msg <- glue(
            "Block {i} has indices outside the range of variantIds. ",
            "Skipping."
        )
        warn(msg)
        return(NULL)
    }
    blockVariants <- variantIds[startIdx:endIdx]
    list(
        matrix = matrix[blockVariants, blockVariants, drop = FALSE],
        mapping = tibble(
            variant_id = blockVariants,
            blockId = i
        )
    )
}

# @noRd
.ldExtractBlockAt <- function(i, matrix, variantIds, blockMetadata) {
    .extractOneBlock(
        matrix,
        variantIds,
        blockMetadata$startIdx[i],
        blockMetadata$endIdx[i],
        i
    )
}

# `x[[i]] <- v` in a loop never extends past the last assigned position, so
# trailing skipped blocks leave no entry at all.
# @noRd
.ldTrimTrailingNull <- function(xs) {
    filled <- which(!map_lgl(xs, is.null))
    if (length(filled) == 0L) {
        return(list())
    }
    xs[seq_len(max(filled))]
}

extractBlockMatrices <- function(matrix, blockMetadata, variantIds) {
    blocks <- map(
        seq_len(nrow(blockMetadata)),
        .ldExtractBlockAt,
        matrix = matrix,
        variantIds = variantIds,
        blockMetadata = blockMetadata
    )
    kept <- compact(blocks)
    mappings <- map(kept, "mapping")
    list(
        # A skipped block leaves a hole, as `ldMatrices[[i]] <- ...` did:
        # positions stay aligned with blockMetadata rows.
        ldMatrices = .ldTrimTrailingNull(map(blocks, "matrix")),
        variantIndices = if (length(mappings) == 0L) {
            tibble(variant_id = character(), blockId = integer())
        } else {
            bind_rows(mappings)
        },
        blockMetadata = blockMetadata
    )
}


# The PSD repair `method` asks for, or the matrix unchanged when it is
# already positive definite. Returns list(R, methodApplied).
# @noRd
.checkLdRepair <- function(R, eig, vals, method, isPd, shrinkage, p, rTol) {
    if (isPd) {
        return(list(R = R, methodApplied = "none"))
    }
    if (method == "shrink") {
        return(list(
            R = (1 - shrinkage) * R + shrinkage * diag(p),
            methodApplied = "shrink"
        ))
    }
    if (method != "eigenfix") {
        return(list(R = R, methodApplied = "none"))
    }
    # Negative eigenvalues raised to a small POSITIVE value, not zero: rTol
    # makes the result strictly positive definite, which the Cholesky-based
    # methods (PRS-CS, SDPR) require; exactly zero would be PSD but not PD.
    rebuilt <- eig$vectors %*% diag(pmax(vals, rTol)) %*% t(eig$vectors)
    # Restore exact symmetry and unit diagonal.
    list(
        R = `diag<-`((rebuilt + t(rebuilt)) / 2, 1),
        methodApplied = "eigenfix"
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
    method <- arg_match(method)
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

    fixed <- .checkLdRepair(R, eig, vals, method, isPd, shrinkage, p, rTol)
    Rout <- fixed$R
    methodApplied <- fixed$methodApplied

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
    # Keep the first member of each cluster and drop the rest -- which is
    # exactly the entries that repeat a cluster already seen.
    ind.delete <- which(duplicated(clusters))
    X.new <- X
    filter.id <- seq_len(p)
    if (length(ind.delete) > 0) {
        # drop = FALSE keeps the column names when a single column survives;
        # without it the result degrades to a vector and the names have to be
        # recomputed by index arithmetic afterwards.
        X.new <- as.matrix(X[, -ind.delete, drop = FALSE])
        filter.id <- filter.id[-ind.delete]
        if (verbose) {
            nDel <- length(ind.delete)
            msg <- glue(
                "ldPruneByCorrelation: pruned {nDel} of {p} columns at ",
                "|cor| > {corThres}"
            )
            inform(msg)
        }
    } else if (verbose) {
        msg <- glue(
            "ldPruneByCorrelation: no columns pruned at |cor| > {corThres}"
        )
        inform(msg)
    }
    list(X.new = X.new, filter.id = filter.id)
}

#' @title Options for the SNPRelate LD-Pruning Backend
#' @description Build a checked record of extra arguments for
#'   \code{SNPRelate::snpgdsLDpruning()}, the engine behind
#'   \code{ldPruneByCorrelation(backend = "snprelate")}.
#' @param ... Arguments for \code{SNPRelate::snpgdsLDpruning()} -- in
#'   practice the window controls \code{slide.max.bp} and
#'   \code{slide.max.n} (pruning is greedy within a window, so the window
#'   bounds which pairs are ever compared), \code{maf},
#'   \code{missing.rate}, \code{remove.monosnp} and \code{num.thread}.
#'   \code{gdsobj}, \code{method}, \code{ld.threshold} and \code{verbose}
#'   are supplied by pecotmr and refused. The temporary GDS pecotmr writes
#'   holds every sample and variant of \code{X} on one synthetic chromosome,
#'   so \code{sample.id}, \code{snp.id} and \code{autosome.only} are not
#'   meaningful selectors here.
#' @return A \code{MethodConfig} record for
#'   \code{ldPruneByCorrelation(methodArgs =)}.
#' @seealso \code{\link{ldPruneByCorrelation}}
#' @examples
#' ldPruningConfig(slide.max.bp = 1e6)
#' @export
ldPruningConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            gdsobj = "the temporary GDS pecotmr writes from `X`",
            method = "fixed at 'corr' by this backend",
            ld.threshold = "the caller's `corThres`",
            verbose = "the caller's `verbose`"
        ),
        "ldPruningConfig"
    )
    .newMethodConfig(
        "SNPRelate::snpgdsLDpruning",
        defaults = list(),
        extra = extra,
        label = "ldPruningConfig",
        engine = "ldPruning"
    )
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
#' @param methodArgs Extra arguments for
#'   \code{SNPRelate::snpgdsLDpruning()}, built with
#'   \code{\link{ldPruningConfig}} -- the window controls in particular.
#'   Only the \code{"snprelate"} backend has an engine to configure, so
#'   supplying options alongside \code{backend = "hclust"} is an error
#'   rather than a silent no-op.
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
    verbose = FALSE,
    methodArgs = ldPruningConfig()
) {
    backend <- arg_match(backend)
    .assertMethodConfig(methodArgs, "ldPruningConfig", "methodArgs")
    if (backend == "snprelate") {
        return(.ldPruneSnprelate(
            X,
            corThres = corThres,
            verbose = verbose,
            methodArgs = methodArgs
        ))
    }
    if (length(methodArgs) > 0L) {
        msg <- glue(
            "ldPruneByCorrelation: `methodArgs` configures ",
            "SNPRelate::snpgdsLDpruning(), which backend 'hclust' does not ",
            "call. Pass backend = 'snprelate', or drop `methodArgs`."
        )
        abort(msg)
    }
    .ldPruneHclust(X, corThres, verbose)
}

#' SNPRelate-based LD pruning helper
#' @noRd
.ldPruneSnprelateDeps <- function() {
    if (
        !requireNamespace("SNPRelate", quietly = TRUE) ||
            !requireNamespace("gdsfmt", quietly = TRUE)
    ) {
        msg <- glue(
            "Packages 'SNPRelate' and 'gdsfmt' are required for ",
            "backend='snprelate'."
        )
        abort(msg)
    }
}

# Write X (rounded to integer genotype codes) to a temporary GDS for SNPRelate.
.ldPruneSnprelateCreateGds <- function(tmpGds, X, snpNames, p) {
    genoInt <- `storage.mode<-`(round(X), "integer")
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

.ldPruneSnprelate <- function(X, corThres, verbose, methodArgs) {
    .ldPruneSnprelateDeps()
    p <- ncol(X)
    snpNames <- colnames(X) %||% str_c("snp", seq_len(p))
    tmpGds <- tempfile(fileext = ".gds")
    on.exit(unlink(tmpGds), add = TRUE)
    .ldPruneSnprelateCreateGds(tmpGds, X, snpNames, p)
    gds <- SNPRelate::snpgdsOpen(tmpGds, allow.duplicate = TRUE)
    on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)
    keepList <- exec(
        SNPRelate::snpgdsLDpruning,
        gds,
        method = "corr",
        ld.threshold = corThres,
        verbose = verbose,
        !!!as.list(methodArgs)
    )
    keepIds <- sort(unname(list_c(keepList)))
    X.new <- X[, keepIds, drop = FALSE]
    if (verbose) {
        nKept <- length(keepIds)
        msg <- glue(
            "ldPruneByCorrelation (snprelate): kept {nKept} of {p} ",
            "columns at |cor| > {corThres}"
        )
        inform(msg)
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
#'   the lowest-variance candidate), or \code{"responseCorrelation"} (remove
#'   the candidate whose correlation with \code{response} has the smallest
#'   magnitude).
#' @param response Numeric vector required when \code{strategy =
#'   "responseCorrelation"}; the outcome to correlate against.
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
.dropCollinearPickCor <- function(X, problematicCols, verbose, seed = NULL) {
    corMatrix <- `diag<-`(
        abs(cor(X[, problematicCols, drop = FALSE])),
        0
    )
    if (length(problematicCols) == 2) {
        if (!is.null(seed)) {
            withr::local_seed(seed)
        }
        colToRemove <- sample(problematicCols, 1)
        if (verbose) {
            inform(glue(
                "dropCollinearColumns: two candidates, randomly removing ",
                "{colToRemove}"
            ))
        }
        return(colToRemove)
    }
    colToRemove <- problematicCols[which.max(colSums(corMatrix))]
    if (verbose) {
        inform(glue(
            "dropCollinearColumns: highest sum |cor| -> removing ",
            "{colToRemove}"
        ))
    }
    colToRemove
}

# Choose which of >=2 collinear columns to drop, per the requested strategy.
.dropCollinearPick <- function(
    X,
    problematicCols,
    strategy,
    response,
    verbose,
    seed = NULL
) {
    if (strategy == "variance") {
        variances <- apply(X[, problematicCols, drop = FALSE], 2, var)
        colToRemove <- problematicCols[which.min(variances)]
        if (verbose) {
            inform(glue(
                "dropCollinearColumns: smallest variance -> removing ",
                "{colToRemove}"
            ))
        }
        return(colToRemove)
    }
    if (strategy == "correlation") {
        return(.dropCollinearPickCor(X, problematicCols, verbose, seed = seed))
    }
    if (is.null(response)) {
        abort(glue(
            "response must be supplied for strategy = ",
            "'responseCorrelation'"
        ))
    }
    corWithResponse <- apply(
        X[, problematicCols, drop = FALSE],
        2,
        cor,
        y = response
    )
    colToRemove <- problematicCols[which.min(abs(corWithResponse))]
    if (verbose) {
        inform(glue(
            "dropCollinearColumns: smallest |cor| with response -> ",
            "removing {colToRemove}"
        ))
    }
    colToRemove
}

dropCollinearColumns <- function(
    X,
    problematicCols,
    strategy = c("correlation", "variance", "responseCorrelation"),
    response = NULL,
    verbose = FALSE,
    seed = NULL
) {
    strategy <- arg_match(strategy)
    if (length(problematicCols) == 0) {
        return(X)
    }
    if (length(problematicCols) == 1) {
        colToRemove <- problematicCols[1]
        if (verbose) {
            msg <- glue(
                "dropCollinearColumns: removing single column {colToRemove}"
            )
            inform(msg)
        }
        return(X[, !is_in(colnames(X), colToRemove), drop = FALSE])
    }
    colToRemove <- .dropCollinearPick(
        X,
        problematicCols,
        strategy,
        response,
        verbose,
        seed = seed
    )
    X[, !is_in(colnames(X), colToRemove), drop = FALSE]
}

# Design matrix [1 | X | C] with the intercept + X columns named.
# @noRd
.ldBuildDesign <- function(X, C) {
    XD <- cbind(1, X, C)
    `colnames<-`(
        XD,
        replace(
            colnames(XD),
            seq_len(ncol(X) + 1L),
            c("Intercept", colnames(X))
        )
    )
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
    nms[is_in(nms, colnames(X))]
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
    Xtemp <- X[, !is_in(colnames(X), problematicColnames), drop = FALSE]
    tempDesign <- .ldBuildDesign(Xtemp, C)
    if (qr(tempDesign)$rank == ncol(tempDesign)) {
        if (verbose) {
            nCol <- length(problematicColnames)
            inform(glue(
                "enforceDesignFullRank: full rank after batch-removing ",
                "{nCol} column(s)"
            ))
        }
        return(FALSE)
    }
    if (verbose) {
        inform(glue(
            "enforceDesignFullRank: batch removal insufficient, ",
            "skipping to correlation-pruning fallback"
        ))
    }
    TRUE
}

# Iteratively drop collinear columns until the design is full rank.
.edfrIterate <- function(
    X,
    C,
    strategy,
    response,
    maxIterations,
    verbose,
    seed = NULL
) {
    # Deliberate iteration: each round drops collinear columns and re-tests
    # the rank, and `X` is a samples x variants genotype matrix -- recursing
    # would hold every intermediate copy alive on the stack.
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
            verbose = verbose,
            seed = seed
        )
        Xdesign <- .ldBuildDesign(X, C)
        matrixRank <- qr(Xdesign)$rank
        iteration <- iteration + 1L
        if (verbose) {
            nCol <- ncol(Xdesign)
            inform(glue(
                "enforceDesignFullRank: iter {iteration} rank ",
                "{matrixRank} / {nCol}"
            ))
        }
    }
    if (iteration == maxIterations) {
        warn(glue(
            "enforceDesignFullRank: maxIterations reached; design may ",
            "still be rank-deficient"
        ))
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
        inform("enforceDesignFullRank: applying ldPruneByCorrelation fallback")
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
            nCol <- ncol(Xdesign)
            msg <- glue(
                "enforceDesignFullRank: threshold {threshold} -> rank ",
                "{matrixRank} / {nCol}"
            )
            inform(msg)
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
#'   \code{strategy = "responseCorrelation"}.
#' @param maxIterations Integer. Hard cap on the iterative-prune loop. Default
#'   300.
#' @param corrThresholds Numeric vector of |cor| thresholds used for the
#'   \code{\link{ldPruneByCorrelation}} fallback, tried in order. Default
#'   \code{seq(0.75, 0.5, by = -0.05)}.
#' @param verbose Logical. If TRUE, print per-iteration progress. Default FALSE.
#' @param seed Integer or \code{NULL}. Seeds the tie-break draw on the rare
#'   occasion that exactly two equally-collinear columns are candidates for
#'   removal (\code{strategy = "correlation"}), via a scoped
#'   \code{withr::local_seed}, so the session RNG is left untouched. \code{NULL}
#'   (default) leaves the draw under the session RNG, so an outer
#'   \code{set.seed()} still governs it.
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
    strategy = c("correlation", "variance", "responseCorrelation"),
    response = NULL,
    maxIterations = 300L,
    corrThresholds = seq(0.75, 0.5, by = -0.05),
    verbose = FALSE,
    seed = NULL
) {
    strategy <- arg_match(strategy)
    originalColnames <- colnames(X)
    initialNcol <- ncol(X)
    Xdesign <- .ldBuildDesign(X, C)
    matrixRank <- qr(Xdesign)$rank
    if (verbose) {
        nCol <- ncol(Xdesign)
        msg <- glue(
            "enforceDesignFullRank: initial rank {matrixRank} / {nCol}"
        )
        inform(msg)
    }
    skipIterative <- .edfrCheckBatch(X, C, Xdesign, matrixRank, verbose)
    iterated <- if (skipIterative) {
        X
    } else {
        .edfrIterate(
            X,
            C,
            strategy,
            response,
            maxIterations,
            verbose,
            seed = seed
        )
    }
    reduced <- .edfrCorrelationFallback(iterated, C, corrThresholds, verbose)
    if (ncol(reduced) == 1L && initialNcol == 1L) {
        return(`colnames<-`(reduced, originalColnames))
    }
    reduced
}

# Require the bigsnpr/bigstatsr packages used for score-based LD clumping.
.ldClumpCheckDeps <- function() {
    if (!requireNamespace("bigsnpr", quietly = TRUE)) {
        abort("Package 'bigsnpr' is required.")
    }
    if (!requireNamespace("bigstatsr", quietly = TRUE)) {
        abort("Package 'bigstatsr' is required.")
    }
}

# Validate the clumping inputs (dimensions of score/chr/pos vs X).
#' @importFrom checkmate assertVector
.ldClumpValidate <- function(X, score, chr, pos) {
    # NOT assertMatrix: X may be a bigstatsr FBM, which is not a base matrix.
    if (ncol(X) < 1L) {
        abort("ldClumpByScore: X must have at least one column")
    }
    assertVector(score, len = ncol(X), null.ok = TRUE)
    assertVector(chr, len = ncol(X))
    assertVector(pos, len = ncol(X))
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
            inform("ldClumpByScore: single variant, skipping clumping")
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
        nKeep <- length(keep)
        nCol <- ncol(X)
        msg <- glue(
            "ldClumpByScore: {nKeep} / {nCol} variants retained at ",
            "r2 <= {r2}"
        )
        inform(msg)
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
#' @importFrom checkmate assertClass
#' @noRd
extractLdMatrix <- function(ld, wantGenotype = FALSE) {
    assertClass(ld, "LdData")
    if (wantGenotype && hasGenotypes(ld)) {
        return(getGenotypes(ld))
    }
    getCorrelation(ld)
}


# ---- Per-block helpers ------------------------------------------------------

# =============================================================================
# LD correlation matrix from a dosage matrix
# -----------------------------------------------------------------------------
# Direct LD computation from an n x p dosage matrix. The internal backend
# uses Rfast::cora when available (else base cor()); optional snprelate /
# snpstats backends round-trip through a temp GDS or SnpMatrix. The
# population and GCTA methods match PLINK / GCTA conventions for missing
# data handling.
# =============================================================================

# LD for a block of a genotype panel: read the dosages, then correlate them.
#
# The guards are the reason this is a function rather than two lines at the
# call site. A backend that cannot read the block returns NULL, and a
# single-variant block has no correlation structure; both mean "identity",
# and open-coding that in every caller is how one of them ends up omitting it.
# @noRd
.computeLdFromPanel <- function(
    panel,
    snpIdx,
    method,
    backend,
    trimSamples,
    shrinkage
) {
    idx <- if (is.null(snpIdx)) {
        seq_along(.ldSketchRanges(panel))
    } else {
        snpIdx
    }
    geno <- .ldSketchDosage(panel, idx)
    if (is.null(geno)) {
        return(diag(length(idx)))
    }
    if (ncol(geno) < 2L) {
        return(diag(length(idx)))
    }
    computeLd(
        geno,
        method = method,
        backend = backend,
        trimSamples = trimSamples,
        shrinkage = shrinkage
    )
}

# Compute LD without ever materialising the dosage matrix.
#
# The dosage matrix scales with samples x variants while the LD matrix scales
# with variants^2, so on a panel with many samples the input dwarfs the
# output: 500k samples x 10k variants is ~40 GB in, ~800 MB out. Reading the
# correlation straight off the GDS keeps the genotypes on disk, which is the
# difference between "slow" and "impossible" on a large block.
#
# It is NOT a drop-in for the in-memory path. SNPRelate applies its own
# missing-data policy, whereas the in-memory backends mean-impute before
# correlating. On complete data the two agree to ~2e-15; with missing calls
# they diverge, and the gap grows with the missing rate (max |difference|
# ~0.04 at 2% missing, ~0.11 at 10% on a synthetic panel). That is a choice
# of estimator, not a rounding difference, so it is opt-in rather than
# automatic.
# @noRd
.computeLdOnDisk <- function(X, method, backend, shrinkage, snpIdx = NULL) {
    if (backend != "snprelate") {
        abort(glue(
            "computeLd(onDisk = TRUE) is only available for ",
            "backend = 'snprelate' (got '{backend}'); the other backends ",
            "correlate a matrix that is already in memory."
        ))
    }
    if (method != "sample") {
        abort(glue(
            "computeLd(onDisk = TRUE) computes the sample correlation; ",
            "method = '{method}' has no on-disk implementation."
        ))
    }
    if (!.ldIsPanel(X)) {
        abort(glue(
            "computeLd(onDisk = TRUE) reads through a genotype panel, so ",
            "`X` must be one (got {class(X)[[1L]]}). A dosage matrix is ",
            "already in memory, which is what onDisk exists to avoid."
        ))
    }
    fmt <- getFormat(.ldSketchHandle(X))
    if (fmt != "gds") {
        abort(glue(
            "computeLd(onDisk = TRUE) needs a GDS-backed panel; this one is ",
            "'{fmt}'. Only GDS exposes an on-disk LD routine."
        ))
    }
    idx <- if (is.null(snpIdx)) {
        seq_along(.ldSketchRanges(X))
    } else {
        snpIdx
    }
    # The handle stays here on purpose: the GDS on-disk LD routine reads the
    # file directly, which is a seed-level operation with no assay equivalent.
    R <- .ldCleanCorrelation(.computeBlockLdGds(.ldSketchHandle(X), idx))
    .ldShrinkToIdentity(R, shrinkage)
}

# A raw estimator output made a usable correlation matrix: unit diagonal and
# no NA / NaN cells.
# @noRd
.ldCleanCorrelation <- function(R) {
    unitDiag <- `diag<-`(R, 1.0)
    replace(unitDiag, is.na(unitDiag) | is.nan(unitDiag), 0)
}

# Optional shrinkage toward the identity (lassosum, Mak et al 2017).
# @noRd
.ldShrinkToIdentity <- function(R, shrinkage) {
    if (shrinkage <= 0 || shrinkage > 1) {
        return(R)
    }
    (1 - shrinkage) * R + shrinkage * diag(nrow(R))
}

# --- computeLd method helpers -----------------------------------------------

# Non-sample methods only support the internal backend.
.computeLdRequireInternal <- function(backend) {
    if (backend != "internal") {
        msg <- glue(
            "backend '{backend}' is only supported with method='sample'."
        )
        abort(msg)
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
    X_imp <- if (anyNA(X)) .ldMeanImputeColumns(X) else X
    if (requireNamespace("Rfast", quietly = TRUE)) {
        # large=FALSE uses tcrossprod internally, ~40x faster than large=TRUE.
        Rfast::cora(X_imp, large = FALSE)
    } else {
        cor(X_imp)
    }
}

# Every missing cell replaced by its column mean.
# @noRd
.ldMeanImputeColumns <- function(X) {
    colMeansX <- colMeans(X, na.rm = TRUE)
    naPos <- which(is.na(X), arr.ind = TRUE)
    replace(X, naPos, colMeansX[naPos[, 2]])
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
            maxNa <- round(max(naRates), 3)
            minNa <- round(min(naRates), 3)
            msg <- glue(
                "Population LD method with heterogeneous missingness ",
                "(max NA rate {maxNa}, min {minNa}): correlations may be ",
                "biased. Consider using method='sample' which handles ",
                "missingness via mean imputation."
            )
            warn(msg)
        }
    }
    # Centering keeps the NA pattern (colMeansX is finite wherever a column
    # has data), so a missing cell contributes nothing to the crossprod.
    X_c <- replace(sweep(X, 2, colMeansX), is.na(X), 0)
    covMat <- crossprod(X_c) / N
    sdVec <- sqrt(colVarsX)
    covMat / outer(sdVec, sdVec)
}

# GCTA per-pair missing-data covariance (matches DENTIST calcLDFromBfile_gcta):
# tracks per-pair non-missing counts and applies a correction term.
.gctaCovariance <- function(X, colMeansX, N, p) {
    notNa <- !is.na(X)
    X_zero <- replace(X, !notNa, 0)
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
    # A zero-variance pair has no defined correlation; it keeps the 0.001
    # floor GCTA writes there.
    valid <- sdOuter > 0
    replace(matrix(0.001, p, p), valid, covMat[valid] / sdOuter[valid])
}

#' @title Compute an LD Correlation Matrix
#' @description Correlation between variants, from either a dosage matrix
#'   already in memory or a genotype panel read on demand.
#'
#'   Three backends compute the same sample correlation different ways.
#'   \code{"internal"} is a single BLAS crossprod (\code{Rfast::cora}, or
#'   \code{cor()} without it) and is both the fastest and the default;
#'   \code{"snprelate"} and \code{"snpstats"} exist to cross-check it and do
#'   strictly more work, since each converts the same matrix into its own
#'   representation first.
#' @param X A numeric dosage matrix (samples x variants), or a genotype panel
#'   -- a \code{RangedSummarizedExperiment} whose dosage assay reads through
#'   a genotype handle.
#' @param method Estimator: \code{"sample"} (default), \code{"population"}
#'   or \code{"gcta"}. The last two are \code{"internal"}-only.
#' @param backend \code{"internal"} (default), \code{"snprelate"} or
#'   \code{"snpstats"}.
#' @param trimSamples Logical; drop samples to a multiple of four, matching
#'   DENTIST's block handling. Ignored for \code{method = "sample"}.
#' @param shrinkage Numeric in (0, 1]; shrink towards the identity
#'   (lassosum, Mak et al. 2017). Zero (default) applies none.
#' @param onDisk Logical, default \code{FALSE}. Read the correlation
#'   straight off a GDS panel without materialising dosages. The dosage
#'   matrix scales with samples x variants while the result scales with
#'   variants^2, so on a wide panel the input dwarfs the output and keeping it
#'   on disk is what makes a large block feasible at all.
#'
#'   Not a drop-in substitute: SNPRelate applies its own missing-data policy
#'   where the in-memory backends mean-impute first. On complete data the two
#'   agree to ~2e-15, but with missing calls they diverge and the gap widens
#'   with the missing rate. Requires \code{backend = "snprelate"},
#'   \code{method = "sample"} and a GDS-backed panel.
#' @param snpIdx Optional variant selection when \code{X} is a panel; row
#'   indices into it. \code{NULL} (default) uses every variant. Not
#'   meaningful for a dosage matrix, which is already the block.
#' @return A numeric correlation matrix (variants x variants), with
#'   \code{dimnames} taken from the block where available.
#' @examples
#' data(qtlSumStatsExample)
#' panel <- getLdSketch(qtlSumStatsExample)
#' R <- computeLd(panel, snpIdx = 1:5)
#' dim(R)
#' # A dosage matrix already in memory is correlated directly.
#' X <- matrix(rbinom(200, 2, 0.3), nrow = 50, ncol = 4)
#' dim(computeLd(X))
#' @export
computeLd <- function(
    X,
    method = c("sample", "population", "gcta"),
    backend = c("internal", "snprelate", "snpstats"),
    trimSamples = FALSE,
    shrinkage = 0,
    onDisk = FALSE,
    snpIdx = NULL
) {
    if (is.null(X)) {
        abort("X must be provided.")
    }
    method <- arg_match(method)
    backend <- arg_match(backend)
    if (isTRUE(onDisk)) {
        return(.computeLdOnDisk(X, method, backend, shrinkage, snpIdx))
    }
    if (.ldIsPanel(X)) {
        return(.computeLdFromPanel(
            X,
            snpIdx,
            method,
            backend,
            trimSamples,
            shrinkage
        ))
    }
    if (!is.null(snpIdx)) {
        abort(glue(
            "`snpIdx` selects variants from a genotype panel; `X` is a ",
            "{class(X)[[1L]]}, which is already the block to correlate."
        ))
    }
    .computeLdMatrix(X, method, backend, trimSamples, shrinkage)
}

# Correlate a plain genotype block, once the panel and on-disk paths have been
# ruled out. Dimnames are carried across by hand: the estimators return a bare
# matrix, and callers key LD by variant name.
# @noRd
.computeLdMatrix <- function(X, method, backend, trimSamples, shrinkage) {
    nms <- colnames(X)
    raw <- if (method == "sample") {
        .computeLdSample(X, backend)
    } else if (method == "population") {
        .computeLdRequireInternal(backend)
        .computeLdPopulation(X, trimSamples)
    } else {
        .computeLdRequireInternal(backend)
        .computeLdGcta(X, trimSamples)
    }
    R <- .ldShrinkToIdentity(.ldCleanCorrelation(raw), shrinkage)
    `dimnames<-`(R, list(nms, nms))
}

#' Compute LD via SNPRelate (creates a temporary GDS file from the dosage
#' matrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix.
#' @noRd
.computeLdSnprelate <- function(X) {
    if (!requireNamespace("SNPRelate", quietly = TRUE)) {
        abort("Package 'SNPRelate' is required for backend='snprelate'")
    }
    if (!requireNamespace("gdsfmt", quietly = TRUE)) {
        abort("Package 'gdsfmt' is required for backend='snprelate'")
    }

    tmpGds <- tempfile(fileext = ".gds")
    on.exit(unlink(tmpGds), add = TRUE)

    # Round to integer dosage for GDS (0/1/2)
    # 3L is the GDS missing code.
    X_int <- replace(`storage.mode<-`(round(X), "integer"), is.na(X), 3L)

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
    if (!requireNamespace("snpStats", quietly = TRUE)) {
        abort("Package 'snpStats' is required for backend='snpstats'")
    }

    # snpStats expects counts of the B allele as raw codes: 1=AA, 2=AB, 3=BB,
    # 0=NA pecotmr dosage is ALT count (0/1/2), so map: 0->1, 1->2, 2->3, NA->0
    shifted <- round(X) + 1L
    coded <- pmin(replace(shifted, is.na(X) | shifted < 1L, 0L), 3L)
    sm <- new("SnpMatrix", `storage.mode<-`(coded, "raw"))

    raw <- as.matrix(snpStats::ld(sm, stats = "R", depth = ncol(X) - 1L))
    # snpStats::ld returns a sparse-like matrix; ensure full dense
    .ldCleanCorrelation(raw)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The integer BP position parsed from a "chrom:pos:a1:a2" variant id.
# @noRd
.ldVariantPos <- function(v) {
    as.integer(str_split(v, ":")[[1L]][2])
}

# TRUE when a per-block variant table has at least one row.
# @noRd
.ldBlockHasVariants <- function(v) {
    nrow(v) > 0
}

# The first index of a block's variants in the merged variant order.
# @noRd
.ldBlockStartIdx <- function(v, ldVariants) {
    min(match(v, ldVariants))
}

# The last index of a block's variants in the merged variant order.
# @noRd
.ldBlockEndIdx <- function(v, ldVariants) {
    max(match(v, ldVariants))
}

# TRUE when block `i`'s index range is present, finite, and within
# [1, nVariants].
# @noRd
.ldBlockValid <- function(i, blockMetadata, nVariants) {
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
}
