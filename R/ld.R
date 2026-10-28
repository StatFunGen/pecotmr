#' @include MethodParam.R
NULL

#' Deduplicate and sort genomic regions by chromosome and start position.
#' @importFrom dplyr distinct arrange
#' @noRd
orderDedupRegions <- function(df) {
    df |>
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

# Split the comma-joined path column into LD (+ optional bim) path columns.
# @noRd
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

#' Find LD blocks overlapping a query region from a metadata TSV file.
#'
#' @param ldReferenceMetaFile TSV with columns chrom, start, end, path. The path
#'   column may be comma-separated: "ld_file,bim_file".
#' @param region A \code{GRanges} window.
#' @param completeCoverageRequired If TRUE, error when the region extends beyond
#'   available LD blocks.
#' @return A list with: intersections (LD_file_paths, bimFilePaths), ldMetaData,
#'   and parsed region.
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom vroom vroom
#' @noRd
regionalLdMeta <- function(
    ldReferenceMetaFile,
    region,
    completeCoverageRequired = FALSE
) {
    genomicData <- `names<-`(
        vroom(ldReferenceMetaFile),
        c("chrom", "start", "end", "path")
    ) |>
        mutate(chrom = withChrPrefix(.data$chrom)) |>
        .ldMetaWidenWholeChrom() |>
        orderDedupRegions() |>
        .regionalLdParsePaths()
    region <- orderDedupRegions(.regionTable(region))
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

# ---------- Internal: resolve LD source type ----------

# Read + validate the first row of an LD metadata TSV (>=4 columns).
#' @importFrom checkmate checkFileExists
# @noRd
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

# Genotype source descriptor for the resolved path, or NULL when it is not a
# genotype file or PLINK prefix, its format told as readGenotypes() tells it.
# A PLINK source is described by its stem, which its readers take.
.resolveLdGenotypeType <- function(resolved, path) {
    format <- .genotypePathFormat(resolved)
    if (is.null(format)) {
        return(NULL)
    }
    dataPath <- if (str_starts(format, "plink")) {
        .genotypeStem(resolved)
    } else {
        resolved
    }
    list(type = format, dataPath = dataPath, metaPath = path)
}

#' Resolve an LD source metadata TSV to its genotype file type.
#'
#' The metadata TSV has columns: chrom, start, end, path, with one row per
#' chromosome (start=0, end=0) whose path is a per-chromosome genotype file or
#' prefix (PLINK2, PLINK1, VCF, GDS or BGEN). The actual region filter is
#' applied by the genotype loader. pecotmr does not read pre-computed LD from
#' disk, so a row pointing at anything else is an error.
#'
#' This function peeks at the first row to determine the data type. The actual
#' per-chromosome path is resolved later by
#' \code{resolveGenotypePathForRegion()} at load time.
#'
#' A path that is itself a genotype file or PLINK prefix is that source, with
#' no metadata.
#'
#' @param path Path to a metadata TSV file with columns chrom, start, end,
#'   path, or to a genotype file or PLINK prefix.
#' @return A list with:
#'   \item{type}{"plink2", "plink1", "vcf", "gds" or "bgen"}
#'   \item{dataPath}{Genotype path from first row (for type detection only;
#'   actual
#'     per-chromosome path is resolved at load time)}
#'   \item{metaPath}{The metadata TSV path, or NULL for a genotype source
#'     given directly}
#' @importFrom vroom vroom
#' @noRd
resolveLdSource <- function(path) {
    direct <- .resolveLdGenotypeType(path, NULL)
    if (!is.null(direct)) {
        return(direct)
    }
    meta <- .resolveLdReadMeta(path)
    # Strip the comma-separated bim path, then resolve relative to the meta dir.
    rawPath <- str_remove(meta$path[1], ",.*$")
    resolved <- file.path(dirname(path), rawPath)
    genoType <- .resolveLdGenotypeType(resolved, path)
    if (!is.null(genoType)) {
        return(genoType)
    }
    abort(c(
        glue(
            "The LD metadata's path does not resolve to genotype files: ",
            "{resolved}"
        ),
        i = str_c(
            "pecotmr reads LD only from genotypes (PLINK2, PLINK1, VCF, GDS, ",
            "BGEN); ",
            "wrap a pre-computed LD matrix in an LdMatrix instead."
        )
    ))
}

#' Resolve the correct genotype path for a given region from a metadata TSV.
#' Reads the TSV, finds the row matching the query region's chromosome, and
#' returns the resolved genotype file path or prefix.
#' @importFrom vroom vroom
#' @noRd
resolveGenotypePathForRegion <- function(metaPath, region) {
    queryChrom <- unique(.regionTable(region)$chrom)
    if (length(queryChrom) != 1L) {
        abort(glue(
            "`region` spans chromosomes {str_flatten(queryChrom, ', ')}; ",
            "a genotype file is resolved for one chromosome at a time."
        ))
    }
    meta <- as.data.frame(vroom(metaPath, show_col_types = FALSE)) |>
        `colnames<-`(c("chrom", "start", "end", "path")) |>
        mutate(chrom = withChrPrefix(.data$chrom))

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

# ---------- LD reference: what it is and how it is held ----------
#
# An LD reference is one of these, and every helper below answers for each:
#
#   a genotype panel  a RangedSummarizedExperiment over a GenotypeHandle,
#                     read lazily from disk -- the same shape as QtlDataset's
#                     genotype experiment. A bare GenotypeHandle, which objects
#                     serialised before the slot became an RSE still carry, is
#                     treated as one.
#   a GenotypeMixture several genotype panels whose LD is averaged, still
#                     lazy: each panel is read only for the variants asked
#                     about, and answers for its own samples.
#   an LdMatrix       a correlation matrix the caller built, held in memory.

# Normalise an LD-reference input to the shape the slot holds: a panel or an
# LdMatrix as given, a bare GenotypeHandle wrapped into a panel -- so a sketch
# and QtlDataset's genotype experiment end up being the same object with
# different provenance rather than two spellings of one idea.
# @noRd
.asLdRef <- function(x) {
    if (is.null(x) || .ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        return(x)
    }
    if (methods::is(x, "RangedSummarizedExperiment")) {
        return(x)
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(.genotypeExperiment(x))
    }
    abort(glue(
        "`ldReference` must be a genotype panel (see readGenotypes()), a ",
        "GenotypeMixture or an LdMatrix (got {class(x)[[1L]]})."
    ))
}

# TRUE for an LD reference held as a matrix rather than as genotypes.
# @noRd
.ldIsMatrixRef <- function(x) {
    methods::is(x, "LdMatrix")
}

# TRUE for a mixture of genotype panels.
# @noRd
.ldIsMixtureRef <- function(x) {
    methods::is(x, "GenotypeMixture")
}

# TRUE for an LD reference backed by genotypes: a panel or a mixture of them.
# @noRd
.ldIsGenotypeRef <- function(x) {
    .ldIsPanel(x) || .ldIsMixtureRef(x)
}

# TRUE for anything that can serve as an LD reference.
# @noRd
.ldIsReference <- function(x) {
    .ldIsGenotypeRef(x) || .ldIsMatrixRef(x)
}

# -----------------------------------------------------------------------------
# Per-row LD references
# -----------------------------------------------------------------------------
# A collection whose rows can sit on different LD references -- summary
# statistics, the fine-mapping results and TWAS weights fit from them, coloc
# results -- keeps each distinct reference ONCE, in a named list in its
# `ldReference` slot, and names each row's entry in an `ldReference` mcols
# column (NA for an LD-free row). Fifty rows on one sketch hold one copy of
# it, in memory and in a saved object alike, and the per-row key travels with
# the rows on subsetting the way every per-row value in mcols does.

# The per-tuple form of an `ldReference` argument: one reference for every
# tuple, or a list holding one per tuple. What ldReference() returns reads
# the same way -- an empty list is LD-free and a one-entry list is that
# reference for every tuple -- so a collection's references can be handed
# straight to the next one built from it.
# @noRd
.ldRefPerTuple <- function(ldReference, n) {
    asList <- is.list(ldReference) && !is.object(ldReference)
    if (asList && length(ldReference) <= 1L) {
        ldReference <- if (length(ldReference) == 0L) {
            NULL
        } else {
            ldReference[[1L]]
        }
        asList <- FALSE
    }
    if (!asList) {
        return(rep(list(.asLdRef(ldReference)), n))
    }
    if (length(ldReference) != n) {
        abort(glue(
            "`ldReference` must be one LD reference, or a list with one per ",
            "tuple (got {length(ldReference)} for {n} tuple(s))."
        ))
    }
    map(unname(ldReference), .asLdRef)
}

# Each distinct reference among `perRow` once, named ld1, ld2, ... in order
# of first use, and every row's key into them (NA for an LD-free row).
# @noRd
.ldRefDict <- function(perRow) {
    present <- which(!map_lgl(perRow, is.null))
    ids <- .ldRefIdentities(perRow[present])
    firstIdx <- which(!duplicated(ids))
    keys <- str_c("ld", seq_along(firstIdx))
    rowKeys <- replace(
        rep(NA_character_, length(perRow)),
        present,
        keys[match(ids, ids[firstIdx])]
    )
    refs <- if (length(firstIdx) == 0L) {
        list()
    } else {
        set_names(perRow[present[firstIdx]], keys)
    }
    list(refs = refs, keys = rowKeys)
}

# A content identity per reference, so equal references -- the same object
# handed to every row, or two equal copies -- collapse to one entry. Hashed
# once per distinct object rather than per row: a large LdMatrix shared by
# every row is read once, not once per row.
# @noRd
#' @importFrom rlang hash obj_address
.ldRefIdentities <- function(refs) {
    if (length(refs) == 0L) {
        return(character(0))
    }
    addresses <- map_chr(refs, obj_address)
    firstOfAddress <- !duplicated(addresses)
    hashes <- map_chr(refs[firstOfAddress], hash)
    hashes[match(addresses, addresses[firstOfAddress])]
}

# Each row's key into its collection's LD references (NA for an LD-free row).
# @noRd
.ldRefKeys <- function(x) {
    mcols(x, use.names = FALSE)$ldReference %||% rep(NA_character_, length(x))
}

# What a collection reports as its LD references: the entries some row names,
# so a subset does not report the ones its rows left behind. An object saved
# before the references became per-row holds one reference in the slot and no
# key column, and read as a list it would pass for LD-free; it is refused
# instead (one that held NULL was LD-free, and still reads that way). So is
# one saved before the slot was renamed, which holds it as `ldSketch`.
# @noRd
.ldRefsInUse <- function(x) {
    if (!methods::.hasSlot(x, "ldReference")) {
        .abortSavedLdRef(x)
    }
    refs <- x@ldReference
    if (is.null(refs)) {
        return(list())
    }
    if (!is.list(refs) || is.object(refs)) {
        .abortSavedLdRef(x)
    }
    .ldRefsUsed(refs, .ldRefKeys(x))
}

# @noRd
.abortSavedLdRef <- function(x) {
    abort(c(
        glue(
            "This {class(x)[[1L]]} was saved before LD references became ",
            "per-row, so its rows cannot be matched to its reference."
        ),
        i = "Bring it up to date with `x <- updateObject(x)`."
    ))
}

# The entries of `refs` some key in `keys` names.
# @noRd
.ldRefsUsed <- function(refs, keys) {
    used <- refs[is_in(names(refs), keys)]
    if (length(used) == 0L) list() else used
}

# Every row's LD reference, NULL for an LD-free row.
# @noRd
.rowLdRefs <- function(x) {
    refs <- ldReference(x)
    map(.ldRefKeys(x), .ldRefForKey, refs = refs)
}

# Result records tagged with the LD reference of the input row they came
# from, so the collection built from them keeps every row on its own
# reference.
# @noRd
.recordsWithLdRef <- function(records, ldReference) {
    map(records, list_assign, ldReference = ldReference)
}

# `x` with every row on `ldReference` -- one reference, or a list with one
# per row -- its references and its rows' keys rebuilt together. The slot is
# written first and the keys second, so the object is only validated once
# the two agree.
# @noRd
#' @importFrom methods slot<-
.withRowLdRefs <- function(x, ldReference) {
    dict <- .ldRefDict(.ldRefPerTuple(ldReference, length(x)))
    md <- mcols(x, use.names = FALSE)
    md[["ldReference"]] <- dict$keys
    S4Vectors::`mcols<-`(
        `slot<-`(x, "ldReference", check = FALSE, value = dict$refs),
        value = md
    )
}

# `object` as saved before its references became per-row -- one reference,
# or NULL, in the slot and no row key -- brought up to date: that reference
# becomes the slot's one entry and every row is keyed to it (every row is
# LD-free when it held none). An object saved before the slot was renamed
# holds that reference, a panel or NULL, as `ldSketch`; it is moved to
# `ldReference` first. Every reference is brought up to date too
# (.updateLdRefOrKeep()). A collection that keeps no LD references, such as
# a CtwasResult, is returned as it is.
# @noRd
.updateRowLdRefs <- function(object) {
    object <- .renameSavedSlots(object, c(ldSketch = "ldReference"))
    if (!methods::.hasSlot(object, "ldReference")) {
        return(object)
    }
    old <- methods::slot(object, "ldReference")
    if (!is.list(old) || is.object(old)) {
        return(.withRowLdRefs(object, .updateLdRefOrKeep(old)))
    }
    `slot<-`(
        object,
        "ldReference",
        check = FALSE,
        value = map(old, .updateLdRefOrKeep)
    )
}

# A collection's LD reference as saved, brought up to date (.updatePanel()).
# One whose files cannot be read again is kept as saved, with a warning,
# rather than failing the update: the rest of the collection is still worth
# having, and updating again once the files are back finishes the job.
# @noRd
#' @importFrom rlang try_fetch
.updateLdRefOrKeep <- function(ref) {
    try_fetch(
        .updatePanel(ref),
        pecotmrUnreadablePanel = function(cnd) {
            warn(
                "An LD reference is kept as saved: LD cannot be read from it.",
                parent = cnd
            )
            .asLdRef(ref)
        }
    )
}

# Row `i`'s LD reference, or NULL for an LD-free row.
# @noRd
.rowLdRef <- function(x, i) {
    .ldRefForKey(.ldRefKeys(x)[[i]], ldReference(x))
}

# @noRd
.ldRefForKey <- function(key, refs) {
    if (is.na(key)) NULL else refs[[key]]
}

# What is wrong with an LD-reference list and the keys pointing into it, or
# NULL. An entry's name is what a key points at, so an unnamed, blank or
# repeated one leaves a row pointing at nothing, or at two references.
# @noRd
.ldRefDictProblems <- function(refs, keys) {
    if (length(refs) > 0L) {
        nms <- names(refs)
        badNames <- is.null(nms) ||
            anyNA(nms) ||
            any(nms == "") ||
            anyDuplicated(nms) > 0L
        if (badNames) {
            return("ldReference entries must have unique, non-empty names")
        }
        bad <- nms[!map_lgl(refs, .ldIsReference)]
        if (length(bad) > 0L) {
            return(str_c(
                "ldReference entries that are not LD references: ",
                str_flatten(bad, ", ")
            ))
        }
    }
    unknown <- setdiff(keys[!is.na(keys)], names(refs))
    if (length(unknown) == 0L) {
        return(NULL)
    }
    str_c(
        "rows name LD reference(s) not in ldReference: ",
        str_flatten(unknown, ", ")
    )
}

# The GenotypeHandle a panel reads through, or NULL when it reads through
# none: an LdMatrix, a mixture, or a panel whose dosages come from another
# backend (a VCF, GDS or PLINK1 seed, or a matrix held in memory). The handle
# is the DelayedArray seed of the dosage assay, so it stays reachable after
# the panel is subset. Tolerates a bare handle, which is what objects
# serialised before the slot became an RSE still carry.
# @noRd
.ldRefHandle <- function(x) {
    if (is.null(x) || .ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        return(NULL)
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(x)
    }
    seed <- .panelSeed(x)
    if (!methods::is(seed, "GhSeed")) {
        return(NULL)
    }
    .ghSeedHandle(seed)
}

# The DelayedArray seed under a panel's dosage assay -- a GhSeed, or a
# GenotypeDosageSeed for a VCF, GDS or PLINK1 panel -- or NULL for an assay
# held in memory or combined from several seeds (a per-chromosome set of
# those).
# @noRd
#' @importFrom rlang try_fetch
.panelSeed <- function(panel) {
    dosage <- SummarizedExperiment::assay(
        panel,
        "dosage",
        withDimnames = FALSE
    )
    if (!methods::is(dosage, "DelayedArray")) {
        return(NULL)
    }
    try_fetch(DelayedArray::seed(dosage), error = function(cnd) NULL)
}

# ---------- LD reference: where its genotypes come from ----------
#
# A panel records when it is built where its genotypes come from --
# list(format, path, chromPaths) -- as metadata(panel)$genotypeSource. What
# needs the file reads it there: a sketch's .afreq sidecar, a label, ctwas's
# check that the payload exists, a per-chromosome manifest's shards. So a
# panel need not read through a GenotypeHandle: a VCF, GDS or PLINK1 panel
# reads through its own backend, and any RSE with a dosage assay can serve.

# The record of a handle's genotype source.
# @noRd
.genotypeSourceRecord <- function(handle) {
    list(
        format = genotypeFormat(handle),
        path = path(handle),
        chromPaths = .genotypeChromPaths(handle)
    )
}

# An LD reference's genotype source: the record its panel carries -- or, for
# a panel built before panels carried one and for a bare handle, its
# handle's -- or NULL for a panel held in memory, an LdMatrix or a mixture.
# @noRd
.ldRefSource <- function(x) {
    if (is.null(x) || .ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        return(NULL)
    }
    if (methods::is(x, "RangedSummarizedExperiment")) {
        recorded <- S4Vectors::metadata(x)$genotypeSource
        if (!is.null(recorded)) {
            return(recorded)
        }
    }
    handle <- .ldRefHandle(x)
    if (is.null(handle)) NULL else .genotypeSourceRecord(handle)
}

# @noRd
.ldRefFormat <- function(x) {
    .ldRefSource(x)$format
}

# The file an LD reference reads, resolved to a real path (a bundled
# example's portable reference resolves under the installed package), or
# NULL when it reads none.
# @noRd
.ldRefReadPath <- function(x) {
    p <- .ldRefSource(x)$path
    if (is.null(p)) NULL else .resolveGenotypeResourcePath(p)
}

# The per-chromosome files of a panel read from a manifest, named by
# chromosome; empty for a panel read from one file.
# @noRd
.ldRefChromPaths <- function(x) {
    .ldRefSource(x)$chromPaths %||% character(0)
}

# Narrow an LD reference to a variant subset. A panel is cut so that, where it
# can, its seed sheds the dropped variants too (.panelSubsetRows()); an
# LdMatrix subsets its matrix along with its variants; a legacy bare handle
# falls back to the handle-level trim, which needs fileIdx to stay read-safe.
# @noRd
.ldRefSubset <- function(x, keep) {
    if (is.null(x)) {
        return(NULL)
    }
    if (.ldIsMatrixRef(x)) {
        return(x[keep])
    }
    if (.ldIsMixtureRef(x)) {
        return(.gmRowsAt(x, keep))
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(.subsetGenotypeHandle(x, keep))
    }
    .panelSubsetRows(x, keep)
}

# A sketch trimmed to zero variants, every other property preserved. Used when
# an entry set genuinely carries no variants: the object references no LD, so
# the correct retained panel is empty rather than the full genome-wide sketch,
# which on a real panel is ~135 MB of serialized dead weight per such entry.
# NULL-safe and idempotent.
# @noRd
.emptyLdRef <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    if (.ldIsMatrixRef(x)) {
        return(x[integer(0)])
    }
    if (.ldIsMixtureRef(x)) {
        return(.gmRowsAt(x, integer(0)))
    }
    if (methods::is(x, "GenotypeHandle")) {
        return(.emptyGenotypeHandle(x))
    }
    .emptyPanel(x, keepSamples = FALSE)
}

# A panel cut to the variants `keep` marks (a mask or positions). A panel
# whose rows are exactly its seed's, cut to an ordered subset, has the seed
# cut too -- a GenotypeHandle trimmed, a VCF, GDS or PLINK1 seed cut to the
# records the kept rows read -- so the cut panel's seed does not carry the
# dropped variants' metadata (see .emptyPanel()); its rowRanges, samples,
# colData and recorded metadata are its own, cut alike, so no sidecar is read
# again. Any other cut -- of a panel that is itself a row subset, of one held
# in memory or whose handle predates fileIdx, or one that reorders or repeats
# rows -- is taken from its assay, which then reads only those rows.
# @noRd
.panelSubsetRows <- function(panel, keep) {
    keepIdx <- if (is.logical(keep)) which(keep) else as.integer(keep)
    dosage <- .panelCutDosage(panel, keepIdx)
    if (is.null(dosage)) {
        # `drop = FALSE`: a one-variant `keep` would otherwise collapse the
        # assay to a vector.
        return(panel[keep, , drop = FALSE])
    }
    SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = dosage[, colnames(panel), drop = FALSE]),
        rowRanges = SummarizedExperiment::rowRanges(panel)[keepIdx],
        colData = SummarizedExperiment::colData(panel),
        metadata = S4Vectors::metadata(panel)
    )
}

# The dosage assay of a panel cut to rows `keepIdx`, over a seed cut alike,
# or NULL when the seed cannot be cut so (see .panelSubsetRows()).
# @noRd
.panelCutDosage <- function(panel, keepIdx) {
    if (
        length(keepIdx) >= nrow(panel) ||
            is.unsorted(keepIdx, strictly = TRUE)
    ) {
        return(NULL)
    }
    seed <- .panelSeed(panel)
    if (
        methods::is(seed, "GenotypeDosageSeed") &&
            identical(rownames(panel), dimnames(seed)[[1L]])
    ) {
        return(DelayedArray::DelayedArray(.gdSeedRows(seed, keepIdx)))
    }
    handle <- .ldRefHandle(panel)
    if (is.null(handle) || !.panelRowsAreHandle(panel, handle)) {
        return(NULL)
    }
    genotypeDelayedArray(.subsetGenotypeHandle(handle, keepIdx))
}

# A panel's dimensions before any subsetting: its dosage seed's, or its own
# when it holds its dosages in memory.
# @noRd
.panelFullDims <- function(panel) {
    seed <- .panelSeed(panel)
    if (is.null(seed)) dim(panel) else dim(seed)
}

# TRUE when a panel's rows are its handle's own, in order, and the handle can
# be trimmed (it records each variant's file position).
# @noRd
.panelRowsAreHandle <- function(panel, handle) {
    si <- snpInfo(handle)
    is_in("fileIdx", names(si)) &&
        identical(rownames(panel), normalizeVariantId(as.character(si$SNP)))
}

# A panel with no variants, held in memory. `x[integer(0), ]` would drop the
# rows but keep the dosage seed -- a handle carrying every variant's
# metadata, or a backend over the whole file -- which is precisely the
# weight an empty panel exists to shed, so no seed is kept. Where it came
# from stays recorded (.ldRefSource()). `keepSamples` keeps its samples and
# their colData, as a mixture's sampleMap needs; otherwise they go too.
# @noRd
.emptyPanel <- function(panel, keepSamples) {
    cols <- if (keepSamples) seq_len(ncol(panel)) else integer(0)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            dosage = matrix(
                numeric(0),
                0L,
                length(cols),
                dimnames = list(NULL, colnames(panel)[cols])
            )
        ),
        rowRanges = SummarizedExperiment::rowRanges(panel)[integer(0)],
        colData = SummarizedExperiment::colData(panel)[cols, , drop = FALSE],
        metadata = list_assign(
            S4Vectors::metadata(panel),
            genotypeSource = .ldRefSource(panel)
        )
    )
}

# "format @ path" for show methods, which describe where a panel came from
# rather than what is in it. An LdMatrix came from no file, so it is
# described by its size instead.
# @noRd
.ldRefLabel <- function(x) {
    if (.ldIsMatrixRef(x)) {
        return(glue("LdMatrix of {length(x)} variants"))
    }
    if (.ldIsMixtureRef(x)) {
        return(glue(
            "GenotypeMixture of {length(.gmComponents(x))} panels: ",
            "{str_flatten(map_chr(.gmComponents(x), .ldRefLabel), '; ')}"
        ))
    }
    source <- .ldRefSource(x)
    if (is.null(source)) {
        return(glue("genotype panel of {ncol(x)} samples held in memory"))
    }
    glue("{source$format} @ {source$path}")
}

# The `show()` line for a collection's LD references: the one reference, or
# each by name when the rows use several; `none` when there are none.
# @noRd
.ldRefsLabel <- function(refs, none) {
    if (length(refs) == 0L) {
        return(none)
    }
    labels <- map_chr(refs, .ldRefLabel)
    if (length(refs) == 1L) {
        return(labels[[1L]])
    }
    str_c(
        length(refs),
        " references -- ",
        str_flatten(str_c(names(refs), ": ", labels), "; ")
    )
}

# TRUE for anything that can serve as an LD panel. A bare GenotypeHandle
# still counts: that is what objects serialised before the slot became an RSE
# carry.
# @noRd
.ldIsPanel <- function(x) {
    methods::is(x, "RangedSummarizedExperiment") ||
        methods::is(x, "GenotypeHandle")
}

# An already-open genotype source -- a panel, or a bare handle -- as given,
# or NULL when the value is not one. Callers that accept "a panel, a bare
# handle, or a spec to resolve" use this to separate the first two, which are
# ready to use, from the third. A panel is the public shape (readGenotypes());
# a bare handle is what internal callers and objects serialised before the
# change still carry.
# @noRd
.openGenotypeSource <- function(x) {
    if (.ldIsPanel(x)) x else NULL
}

# ---------- LD reference: its own view of itself ----------
#
# These read a panel through the SummarizedExperiment interface rather than
# reaching past it to the GenotypeHandle. The handle is the DelayedArray's
# seed -- the READ PATH -- and unwrapping it for metadata leaks that seed into
# callers that only ever wanted chrom/pos/alleles or a sample count, which the
# RSE already answers. Reaching for `.ldRefHandle()` is now reserved for
# what genuinely lives on the handle: file format/path, chromosome shard
# routing, the `.afreq` sidecar, and building or pruning the seed itself.
#
# A bare GenotypeHandle -- what objects serialised before the slot became an
# RSE carry -- is projected exactly the way `.genotypeExperiment()` projects
# one, so both shapes answer identically. An LdMatrix is itself the GRanges
# of its variants, with the same SNP / A1 / A2 columns.

# The reference's variant ranges: seqnames = chromosome (chr-prefixed), start =
# position, mcols = SNP / A1 / A2 (and AF on an LdMatrix that records it).
# @noRd
.ldRefRanges <- function(x, label = "LD reference") {
    # The chokepoint every panel accessor funnels through, so the guard lives
    # here rather than being repeated at each entry point. Without it a NULL
    # or non-panel sketch surfaces as "unable to find an inherited method for
    # 'snpInfo'" from whichever accessor happened to touch it first, which
    # says nothing about the LD reference being the problem. Callers that
    # validate with their own label (`.ldFromRef`, the ctwas assembler) do
    # so first, so their message wins.
    .ldFromRefValidate(x, label)
    if (.ldIsMatrixRef(x)) {
        return(as(x, "GRanges"))
    }
    if (.ldIsMixtureRef(x)) {
        return(.gmRanges(x))
    }
    if (methods::is(x, "RangedSummarizedExperiment")) {
        # A panel built outside pecotmr may spell its chromosomes otherwise.
        return(.ucscRanges(SummarizedExperiment::rowRanges(x)))
    }
    .genotypeSnpRanges(x, normalizeVariantId(as.character(snpInfo(x)$SNP)))
}

# How many variants an LD reference covers.
# @noRd
.ldRefNVariants <- function(x) {
    length(.ldRefRanges(x))
}

# The panel's variant ids in the panel's OWN labelling. This is the mcols
# column, not `names()`: the names are normalized ids, while anything keyed on
# the panel FILE -- the `.afreq` sidecar, the ctwas snp_map -- has to see the
# panel's raw labels.
# @noRd
.ldRefVariantIds <- function(x) {
    as.character(S4Vectors::mcols(.ldRefRanges(x))$SNP)
}

# The panel's variant ids as MATCHING keys: the raw labels, with any id whose
# allele slot holds a tag rather than an allele (chr21:13988152:INS:T) re-
# rendered from the panel's own A1 / A2 columns, which are correct. Use this
# wherever a panel id is matched against summary statistics or another id set,
# and `.ldRefVariantIds()` only where the file's literal label is wanted.
#
# The repair is what lets an untidily-named panel entry be MATCHED rather than
# discarded: harmonization re-keys a surviving variant to the reference-allele
# id, so the panel side has to speak the same form or every later id lookup
# misses the variant it just harmonized.
# @noRd
.ldRefMatchIds <- function(x) {
    mc <- S4Vectors::mcols(.ldRefRanges(x))
    ids <- as.character(mc$SNP)
    if (is.null(mc$A1) || is.null(mc$A2)) {
        return(ids)
    }
    .repairVariantIds(ids, as.character(mc$A2), as.character(mc$A1))
}

# Per-variant chromosome.
# @noRd
.ldRefChrom <- function(x) {
    as.character(GenomicRanges::seqnames(.ldRefRanges(x)))
}

# A panel's sample ids. An LdMatrix has no samples, so asking is a bug in the
# caller rather than something to answer with an empty vector.
# @noRd
.ldRefSampleIds <- function(x) {
    if (.ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        abort(
            "Only a single genotype panel has one sample set.",
            .internal = TRUE
        )
    }
    if (methods::is(x, "RangedSummarizedExperiment")) {
        return(as.character(colnames(x)))
    }
    as.character(sampleIds(x))
}

# A panel's column count -- its samples, or a sketch's pseudo-samples -- or
# the size an LdMatrix records (NA when it records none).
# @noRd
.ldRefNSamples <- function(x) {
    if (.ldIsMatrixRef(x)) {
        return(nSamples(x))
    }
    if (.ldIsMixtureRef(x)) {
        return(.gmPooledN(x, .ldRefNSamples))
    }
    length(.ldRefSampleIds(x))
}

# The number of real people behind an LD reference, for statistics that count
# them -- a MAC cutoff's allele counts, DENTIST's truncation rank. A genotype
# panel's columns are its samples. A stochastic sketch's columns are B
# pseudo-samples projected from n real ones (the projection cannot raise the
# rank of their LD above n), so its n is the one its projection record holds.
# An LdMatrix records its size, or NA.
# @noRd
.ldRefRealNSamples <- function(x) {
    if (.ldIsMatrixRef(x)) {
        return(nSamples(x))
    }
    if (.ldIsMixtureRef(x)) {
        return(sum(map_int(.gmComponents(x), .ldRefRealNSamples)))
    }
    .ldRefProjection(x)$nSamples %||% .ldRefNSamples(x)
}

# The number of independent samples a reference's LD is worth, for the
# estimators that correct for LD-estimation noise (ldsc's unbiased r^2, HDL's
# finite reference). A stochastic sketch's LD carries the noise of its n real
# samples and of its B pseudo-samples, so it is worth 1 / (1/n + 1/B). An
# LdMatrix records this directly.
# @noRd
.ldRefNoiseN <- function(x) {
    if (.ldIsMatrixRef(x)) {
        return(nSamples(x))
    }
    if (.ldIsMixtureRef(x)) {
        return(.gmPooledN(x, .ldRefNoiseN))
    }
    nCols <- .ldRefNSamples(x)
    n <- .ldRefProjection(x)$nSamples
    if (is.null(n)) {
        return(nCols)
    }
    as.integer(round(1 / (1 / n + 1 / nCols)))
}

# Each of the reference's variants' A1 frequency, or NULL when it records
# none. An LdMatrix and a stochastic sketch carry an `AF` column; a bare
# sketch handle (an object serialised before the slot became an RSE) answers
# from its sidecar. A panel of real genotypes measures frequencies from its
# dosages instead, so it answers NULL here.
# @noRd
.ldRefAf <- function(x) {
    if (.ldIsMixtureRef(x)) {
        return(.gmAf(x))
    }
    af <- S4Vectors::mcols(.ldRefRanges(x))$AF
    if (!is.null(af) || !methods::is(x, "GenotypeHandle")) {
        return(af)
    }
    record <- .ldRefProjection(x)
    if (is.null(record)) {
        return(NULL)
    }
    as.numeric(record$af$alt_freq)[match(.ldRefVariantIds(x), record$af$id)]
}

# The real sample size a step needs, or an error naming the step when the
# reference cannot say -- an LdMatrix built without `nSamples`.
# @noRd
.ldRefRequireNSamples <- function(x, step) {
    n <- .ldRefRealNSamples(x)
    if (is.na(n)) {
        abort(c(
            glue("{step} needs the LD reference's sample size."),
            i = "Record it with LdMatrix(nSamples = )."
        ))
    }
    n
}

# A stochastic sketch's projection record -- list(nSamples, af) -- or NULL for
# a panel of real genotypes and for an LdMatrix. A panel built by
# .genotypeExperiment() carries it from construction (see .withLdProjection);
# a bare handle, which only objects serialised before the slot became an RSE
# carry, is read from its sidecar on demand.
# @noRd
.ldRefProjection <- function(x) {
    if (is.null(x) || .ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        return(NULL)
    }
    if (methods::is(x, "RangedSummarizedExperiment")) {
        return(S4Vectors::metadata(x)$ldProjection)
    }
    .sketchProjectionRecord(.ldRefHandle(x))
}

# TRUE for a stochastic sketch: a panel whose columns are pseudo-samples
# projected from real ones, so its dosages are not genotypes and say nothing
# about allele frequency or missingness.
# @noRd
.ldRefIsProjection <- function(x) {
    !is.null(.ldRefProjection(x))
}

# Read a stochastic sketch's sidecar: its real sample count -- the largest
# OBS_CT, halved (OBS_CT counts alleles, two per diploid sample; variants with
# missing calls count fewer) -- and each variant's real A1 frequency. NULL
# when the handle is not a sketch, i.e. carries no PLINK2-style .afreq beside
# it (in any genotype format) with U_MIN / U_MAX and OBS_CT.
# @noRd
.sketchProjectionRecord <- function(handle) {
    if (is.null(handle)) {
        return(NULL)
    }
    .afreqProjectionRecord(.genotypeSourceRecord(handle), snpInfo(handle)$CHR)
}

# .sketchProjectionRecord() for a genotype source and the chromosomes its
# panel spans, for a panel that reads through no handle.
# @noRd
.afreqProjectionRecord <- function(source, chroms) {
    tbls <- compact(map(.afreqPrefixes(source, chroms), .afreqIfStochastic))
    if (length(tbls) == 0L) {
        return(NULL)
    }
    af <- list_rbind(tbls)
    list(
        nSamples = as.integer(max(as.numeric(af$obs_ct), na.rm = TRUE) %/% 2),
        af = select(af, all_of(c("id", "alt_freq")))
    )
}

# Record a stochastic sketch's real sample count and allele frequencies on the
# panel built over it: the count in metadata, the frequencies as an `AF`
# column of rowRanges (the frequency of A1, which PLINK2 writes as ALT). Done
# once, at construction, so nothing downstream re-reads the sidecar -- or,
# worse, measures frequencies off the projected dosages.
# @noRd
.withLdProjection <- function(panel, record) {
    if (is.null(record)) {
        return(panel)
    }
    rr <- SummarizedExperiment::rowRanges(panel)
    mc <- S4Vectors::mcols(rr)
    mc$AF <- as.numeric(record$af$alt_freq)[
        match(as.character(mc$SNP), record$af$id)
    ]
    withAf <- SummarizedExperiment::`rowRanges<-`(
        panel,
        value = S4Vectors::`mcols<-`(rr, value = mc)
    )
    S4Vectors::`metadata<-`(
        withAf,
        value = c(
            S4Vectors::metadata(withAf),
            list(ldProjection = list(nSamples = record$nSamples))
        )
    )
}

# .withLdProjection() for a panel that reads through no handle -- a VCF or
# GDS panel, or a per-chromosome set of them -- its sidecar found from the
# genotype source the panel records.
# @noRd
.withSourceProjection <- function(panel) {
    .withLdProjection(
        panel,
        .afreqProjectionRecord(.ldRefSource(panel), .ldRefChrom(panel))
    )
}

# The .afreq at `prefix` when it is a stochastic sketch's (it carries U_MIN /
# U_MAX) and records OBS_CT, else NULL. A malformed sidecar reads as absent,
# as it does for the panel filter.
# @noRd
#' @importFrom rlang try_fetch
.afreqIfStochastic <- function(prefix) {
    af <- try_fetch(readAfreq(prefix), error = function(cnd) NULL)
    if (is.null(af) || !all(is_in(c("u_min", "obs_ct"), colnames(af)))) {
        return(NULL)
    }
    af
}

# Dosage for a variant subset, samples x variants.
#
# Read through the assay, so the DelayedArray seed -- and with it the file
# ordering, the sharded routing and the empty-request guard in
# `extractBlockGenotypes()` -- stays the single reader. `extract_array()`
# delegates to that same function with `meanImpute = FALSE`, so imputation is
# layered on here rather than pushed down.
# @noRd
.ldRefDosage <- function(x, snpIdx, meanImpute = TRUE) {
    # The other primitive that does not pass through `.ldRefRanges()`.
    .ldFromRefValidate(x, "LD reference")
    if (.ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        abort(
            "Only a single genotype panel has dosages to read.",
            .internal = TRUE
        )
    }
    if (!methods::is(x, "RangedSummarizedExperiment")) {
        return(.dosageMatrix(x, snpIdx, meanImpute = meanImpute))
    }
    dosage <- t(as.matrix(
        SummarizedExperiment::assay(x, "dosage")[snpIdx, , drop = FALSE]
    ))
    if (meanImpute) .qtlMeanImpute(dosage) else dosage
}

# Require a non-NULL LD reference: a genotype panel (a
# RangedSummarizedExperiment, or a bare GenotypeHandle) or an LdMatrix.
# @noRd
.ldFromRefValidate <- function(ldReference, label) {
    if (is.null(ldReference)) {
        msg <- glue(
            "{label}: the SumStats/collection carries no ldReference ",
            "(ldReference = NULL); this step needs an LD reference."
        )
        abort(msg)
    }
    # A bare GenotypeHandle is still accepted: that is what objects
    # serialised before the slot became an RSE carry.
    if (!.ldIsReference(ldReference)) {
        msg <- glue(
            "{label}: ldReference must be a genotype panel ",
            "(RangedSummarizedExperiment or GenotypeHandle), a ",
            "GenotypeMixture or an LdMatrix."
        )
        abort(msg)
    }
}

# Match requested ids to the panel; NULL if none match. Returns kept ids/order.
.ldFromRefMatch <- function(ldReference, variantIds, label, onMissing) {
    # Match by (chrom, pos, allele) tuple with an exact id-string fallback for
    # rsID panels; the caller's original ids and order are preserved.
    m <- matchVariants(
        variantIds,
        .ldRefMatchIds(ldReference),
        removeStrandAmbiguous = FALSE
    )
    nMissing <- length(variantIds) - length(m$idxA)
    if (nMissing > 0L && onMissing == "error") {
        msg <- glue(
            "{label}: {nMissing} variant id(s) not present in the LD ",
            "reference."
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

# The LD matrix over a specified variant subset of an LD reference: computed
# from a genotype panel's dosages, or read out of an LdMatrix. Called directly
# by
# twasWeightsPipeline, fineMappingPipeline, jointEngine,
# causalInferencePipeline, colocboostPipeline and the coloc / summaryStatsQc
# paths, which differ only in their error message prefix and in whether
# variants absent from the panel raise an error or get silently dropped --
# both of which are arguments here rather than a per-pipeline wrapper.
#
# Arguments:
#   ldReference A genotype panel (RangedSummarizedExperiment, or a bare
#               GenotypeHandle) or an LdMatrix.
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
# @noRd
.ldFromRef <- function(
    ldReference,
    variantIds,
    label = ".ldFromRef",
    onMissing = c("error", "drop")
) {
    .ldFromRefValidate(ldReference, label)
    onMissing <- arg_match(onMissing)
    matched <- .ldFromRefMatch(ldReference, variantIds, label, onMissing)
    if (is.null(matched)) {
        return(NULL)
    }
    ldMat <- .ldRefLdAt(ldReference, matched$idx)
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

# The factor of a genotype reference's LD over `variantIds` (see
# .ldRefFactor()), matched and allele-aligned exactly as .ldFromRef() does:
# columns follow the caller's ids and are named by them, and a column whose
# panel alleles are swapped is negated, so crossprod(F) / (nrow(F) - 1) is
# .ldFromRef()'s matrix. An LdMatrix has no factor, so callers send it down
# the matrix path instead.
# @noRd
.ldFactorFromRef <- function(
    ldReference,
    variantIds,
    label = ".ldFactorFromRef",
    onMissing = c("error", "drop")
) {
    .ldFromRefValidate(ldReference, label)
    onMissing <- arg_match(onMissing)
    matched <- .ldFromRefMatch(ldReference, variantIds, label, onMissing)
    if (is.null(matched)) {
        return(NULL)
    }
    factor <- .ldRefFactor(ldReference, matched$idx)
    signed <- if (any(matched$sign < 0)) {
        sweep(factor, 2L, matched$sign, `*`)
    } else {
        factor
    }
    named <- `dimnames<-`(signed, list(NULL, matched$keptIds))
    if (onMissing == "drop") {
        return(`attr<-`(named, "keptVariantIds", matched$keptIds))
    }
    named
}

# The LD among the reference's variants at `idx`, in the reference's own
# allele coding: the sample correlation of a panel's dosages, or the
# corresponding block of an LdMatrix.
# @noRd
.ldRefLdAt <- function(ldReference, idx) {
    if (.ldIsMatrixRef(ldReference)) {
        return(unname(ldMatrix(ldReference)[idx, idx, drop = FALSE]))
    }
    if (.ldIsMixtureRef(ldReference)) {
        return(.gmLdAt(ldReference, idx))
    }
    computeLd(ldReference, method = "sample", snpIdx = idx)
}

# A factor F of a genotype reference's LD over `idx`: a matrix with
# zero-mean columns whose crossprod(F) / (nrow(F) - 1) is exactly the LD
# .ldRefLdAt() returns, so every computation on that LD -- quadratic forms,
# products, columns, eigenvectors -- can run on F without forming the p x p
# matrix. Its rows behave as samples, so an engine that correlates them
# itself recovers the same LD. A panel's factor is its standardized dosages;
# a mixture's stacks its panels' (.gmFactor()).
# @noRd
.ldRefFactor <- function(ldReference, idx) {
    if (.ldIsMixtureRef(ldReference)) {
        return(.gmFactor(ldReference, idx))
    }
    .ldFactorFromDosage(.ldRefDosage(ldReference, idx, meanImpute = TRUE))
}

# A panel's mean-imputed, standardized dosages over `idx`, with a monomorphic
# column left at zero: the blocks a mixture's factor is stacked from.
# @noRd
.ldRefStandardized <- function(ldReference, idx) {
    .ldStandardize(.ldRefDosage(ldReference, idx, meanImpute = TRUE))
}

# The factor of a dosage matrix's (samples x variants) LD, under
# .ldRefFactor()'s contract. A monomorphic variant is independent of every
# other with unit variance in computeLd()'s LD, while its standardized
# column is all zero; two padding rows of +/- sqrt((m - 1) / 2) restore its
# variance without correlating it with anything or moving any column's mean,
# and the other rows are rescaled for the m rows the factor then has.
# @noRd
.ldFactorFromDosage <- function(X) {
    standardized <- .ldStandardize(X)
    dead <- which(colSums(standardized^2) == 0)
    if (length(dead) == 0L) {
        return(standardized)
    }
    n <- nrow(standardized)
    m <- n + 2L * length(dead)
    rbind(
        standardized * sqrt((m - 1) / (n - 1)),
        .ldPaddingRows(dead, ncol(standardized), sqrt((m - 1) / 2))
    )
}

# Mean-imputed columns, centred and scaled to unit variance (n - 1
# denominator), as scale() would make them, except that a constant column
# stays zero rather than turning into NaN.
# @noRd
.ldStandardize <- function(X) {
    imputed <- if (anyNA(X)) .ldMeanImputeColumns(X) else X
    centred <- sweep(imputed, 2L, colMeans(imputed))
    sds <- sqrt(colSums(centred^2) / (nrow(centred) - 1))
    sweep(centred, 2L, replace(sds, sds == 0, 1), `/`)
}

# Two rows per column in `cols`, holding +value and -value in that column
# and zero elsewhere. They add 2 value^2 to the column's sum of squares and
# nothing to its mean or to any crossproduct.
# @noRd
.ldPaddingRows <- function(cols, p, value) {
    rows <- matrix(0, 2L * length(cols), p)
    i <- seq_along(cols)
    rows[cbind(2L * i - 1L, cols)] <- value
    rows[cbind(2L * i, cols)] <- -value
    rows
}

# The eigendecomposition of the LD R = crossprod(F) / (m - 1), from the thin
# SVD of its factor F (m x p) rather than from R: the eigenvalues d^2 / (m - 1)
# and the right singular vectors, for the min(m, p) directions the SVD
# returns, largest first. R's other p - min(m, p) eigenvalues are exactly
# zero; `nNull` counts them, for an estimator that counts every direction.
# @noRd
.ldFactorEigen <- function(ldFactor) {
    s <- svd(ldFactor, nu = 0L)
    list(
        values = s$d^2 / (nrow(ldFactor) - 1),
        vectors = s$v,
        nNull = ncol(ldFactor) - length(s$d)
    )
}

# ---------- LD rows a stripe at a time ----------
#
# An LD source is one block's LD in any of three forms: the matrix itself,
# list(R = ); a factor of it, list(F = ) with R = crossprod(F) / (m - 1); or
# an eigendecomposition, list(vectors = , values = ) with R = V diag(values)
# V'. The helpers below read R's rows from any of them, so a computation over
# R's rows never needs the matrix whole.

# How many of R's rows a stripe holds.
# @noRd
.ldStripeSize <- function() 256L

# The number of variants an LD source covers.
# @noRd
.ldSourceSize <- function(ld) {
    if (!is.null(ld$R)) {
        return(nrow(ld$R))
    }
    if (!is.null(ld$F)) {
        return(ncol(ld$F))
    }
    nrow(ld$vectors)
}

# Rows `rows` of an LD source's R, without forming the rest.
# @noRd
.ldSourceRows <- function(ld, rows) {
    if (!is.null(ld$R)) {
        return(ld$R[rows, , drop = FALSE])
    }
    if (!is.null(ld$F)) {
        return(
            crossprod(ld$F[, rows, drop = FALSE], ld$F) / (nrow(ld$F) - 1)
        )
    }
    scaled <- sweep(ld$vectors[rows, , drop = FALSE], 2L, ld$values, "*")
    tcrossprod(scaled, ld$vectors)
}

# (R o R) A -- R squared element-wise, times A -- for an LD source: column a
# of the result is the LD score sum_k r^2_jk a_k. The element-wise square
# does not inherit a factor's low rank (R has rank below m, R o R is in
# general full rank), so no product through the factor serves every column
# of A at once: one Gram matrix F diag(a) F' per column would cost m^2 p for
# each. Instead R's rows are formed a stripe at a time, squared, and
# multiplied by all of A, which costs p^2 (m + ncol(A)) and holds one stripe.
# A stored matrix is used as it is.
# @noRd
.ldSquaredProduct <- function(ld, A) {
    if (!is.null(ld$R)) {
        return((ld$R^2) %*% A)
    }
    p <- .ldSourceSize(ld)
    starts <- seq.int(1L, p, by = .ldStripeSize())
    stripes <- map(starts, .ldSquaredProductStripe, ld = ld, A = A, p = p)
    exec(rbind, !!!stripes)
}

# @noRd
.ldSquaredProductStripe <- function(start, ld, A, p) {
    rows <- seq.int(start, min(start + .ldStripeSize() - 1L, p))
    .ldSourceRows(ld, rows)^2 %*% A
}

# `ldReference` cut to the variants `variantIds` names (matched as
# .ldFromRef() matches them), for a caller that reads the same block's LD
# several times: a panel's dosage is read from its files once, here, and held
# in memory, rather than read again by every reader. An LdMatrix is cut to
# those variants, and a GenotypeMixture is cut with its panels left to read
# as they are needed.
# @noRd
.ldRefInMemory <- function(ldReference, variantIds, label) {
    .ldFromRefValidate(ldReference, label)
    matched <- .ldFromRefMatch(ldReference, variantIds, label, "error")
    rows <- sort(unique(matched$idx))
    if (!methods::is(ldReference, "RangedSummarizedExperiment")) {
        if (methods::is(ldReference, "GenotypeHandle")) {
            return(ldReference)
        }
        return(ldReference[rows])
    }
    block <- ldReference[rows, ]
    SummarizedExperiment::`assay<-`(
        block,
        "dosage",
        withDimnames = FALSE,
        value = as.matrix(SummarizedExperiment::assay(
            block,
            "dosage",
            withDimnames = FALSE
        ))
    )
}

# b' R b for each column b of B, from an LD matrix or a factor of it:
# ||F b||^2 / (m - 1), never forming R.
# @noRd
.ldQuadForms <- function(ld, B) {
    if (!is.null(ld$R)) {
        return(colSums(B * (ld$R %*% B)))
    }
    colSums((ld$F %*% B)^2) / (nrow(ld$F) - 1)
}

# ---------- LD reference: per-variant statistics and filtering ----------

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

# The .afreq stems to read for a genotype source: its file's, or one per
# chromosome of `chroms` (those its panel spans) for a per-chromosome
# manifest. A manifest lists every chromosome, and reading the shards the
# summary statistics never touch is the cost per-chromosome sketch trimming
# exists to avoid. None for a panel that reads no file.
# @noRd
.afreqPrefixes <- function(source, chroms) {
    if (is.null(source)) {
        return(character(0))
    }
    chromPaths <- source$chromPaths %||% character(0)
    if (length(chromPaths) == 0L) {
        return(.genotypeStem(.resolveGenotypeResourcePath(source$path)))
    }
    spanned <- unique(as.character(chroms))
    paths <- unname(chromPaths[is_in(names(chromPaths), spanned)])
    .genotypeStem(map_chr(paths, .resolveGenotypeResourcePath))
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

# Panel minor allele frequency for `variantIds` read from the PLINK2-style
# .afreq sidecar beside the panel's file, in any genotype format, or NULL
# when the sidecar cannot answer for every id.
#
# Refusing a PARTIAL sidecar is deliberate: the dosage path treats a variant
# it cannot measure as a drop, so a fast path that quietly kept the ids
# .afreq happens to omit would disagree with it about the same variant.
# @noRd
.panelAfreqMaf <- function(source, chroms, variantIds) {
    tbls <- compact(map(.afreqPrefixes(source, chroms), .panelAfreqTable))
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

# TRUE for each matched variant of the LD reference that fails the cutoffs.
# The statistics come from `.ldRefVariantStats()`, which knows where each
# kind of reference keeps them truthfully; missingness is only asked for when
# a missingness cutoff is set.
# @noRd
.panelDropMask <- function(ldReference, matched, panelFilterParam) {
    imissCutoff <- panelFilterParam$imissCutoff
    stats <- .ldRefVariantStats(
        ldReference,
        matched$idx,
        needMissing = imissCutoff < 1,
        step = "The LD-panel filter"
    )
    effMaf <- .panelEffectiveMaf(
        panelFilterParam$mafCutoff,
        panelFilterParam$macCutoff,
        .ldRefMacSamples(ldReference, panelFilterParam$macCutoff)
    )
    fails <- is.na(stats$maf) | stats$maf < effMaf
    if (is.null(stats$missRate)) {
        return(fails)
    }
    fails | stats$missRate > imissCutoff
}

# Per-variant MAF -- and missingness, when `needMissing` -- for the
# reference's variants at `idx`, from wherever the reference records them
# truthfully:
#
#   genotype panel  the .afreq sidecar when it answers for every variant and
#                   no missingness is wanted (no dosage need be read, which is
#                   what keeps the filter affordable on a panel far larger
#                   than the analysed variant set); the dosages otherwise
#   sketch          the real frequencies its projection record holds. Its
#                   dosages are projections taken after mean-filling, so they
#                   carry neither frequencies nor missingness
#   LdMatrix        its AF column; it has no genotypes to be missing
#
# `step` names the caller in the error when the reference cannot answer.
# @noRd
.ldRefVariantStats <- function(ldReference, idx, needMissing, step) {
    if (.ldIsMixtureRef(ldReference)) {
        return(.gmVariantStats(ldReference, idx, needMissing, step))
    }
    if (.ldIsMatrixRef(ldReference) || .ldRefIsProjection(ldReference)) {
        return(.ldRefRecordedStats(ldReference, idx, needMissing, step))
    }
    if (!needMissing) {
        # Keyed on the PANEL's own variant labels: .afreq carries the .pvar
        # ids. The handle stays for this one read, since the sidecar is a
        # file beside the genotypes, which only the seed knows how to find.
        panelIds <- .ldRefVariantIds(ldReference)[idx]
        maf <- .panelAfreqMaf(
            .ldRefSource(ldReference),
            .ldRefChrom(ldReference),
            panelIds
        )
        if (!is.null(maf)) {
            return(list(maf = maf, missRate = NULL))
        }
    }
    .panelVariantStats(.ldRefDosage(ldReference, idx, meanImpute = FALSE))
}

# Frequencies from a reference that records them rather than measuring them:
# an LdMatrix or a stochastic sketch. Neither can say how often a genotype is
# missing, so a missingness cutoff is refused rather than silently skipped.
# @noRd
.ldRefRecordedStats <- function(ldReference, idx, needMissing, step) {
    what <- if (.ldIsMatrixRef(ldReference)) {
        "an LdMatrix, which has no genotypes"
    } else {
        "a stochastic sketch, whose dosages are projections"
    }
    if (needMissing) {
        abort(c(
            glue(
                "{step} needs genotype missingness, and the LD reference is ",
                "{what}."
            ),
            i = "Leave its missingness cutoff unset."
        ))
    }
    af <- .ldRefAf(ldReference)
    if (is.null(af)) {
        abort(c(
            glue(
                "{step} needs allele frequencies, and the LD reference ",
                "records none."
            ),
            i = "Give the LdMatrix's variants an `AF` column."
        ))
    }
    af <- as.numeric(af)[idx]
    list(maf = pmin(af, 1 - af), missRate = NULL)
}

# Per-variant genotype variance at `idx`, for scaling non-standardized weights
# to the correlation scale. A panel of real genotypes measures it on its
# (mean-imputed) dosages. A sketch or an LdMatrix has no genotypes to measure,
# so it gives 2p(1 - p) from its recorded A1 frequency -- the dosage variance
# under Hardy-Weinberg -- or NA when it records none.
# @noRd
.ldRefVariance <- function(ldReference, idx) {
    if (.ldIsMixtureRef(ldReference)) {
        return(.gmVariance(ldReference, idx))
    }
    if (.ldIsMatrixRef(ldReference) || .ldRefIsProjection(ldReference)) {
        af <- .ldRefAf(ldReference)
        if (is.null(af)) {
            return(rep(NA_real_, length(idx)))
        }
        af <- as.numeric(af)[idx]
        return(2 * af * (1 - af))
    }
    geno <- .ldRefDosage(ldReference, idx, meanImpute = TRUE)
    apply(geno, 2, stats::var, na.rm = TRUE)
}

# The real sample count a MAC cutoff divides by, or 0 when no MAC cutoff is
# set (which .panelEffectiveMaf() reads as "no MAC term"). A mixture of
# panels has no single count, so a MAC cutoff on one is refused.
# @noRd
.ldRefMacSamples <- function(ldReference, macCutoff) {
    if (macCutoff <= 0) {
        return(0L)
    }
    # A mixture's statistics are its panels' worst (see .gmVariantStats()), so
    # its counts are judged against its smallest panel.
    if (.ldIsMixtureRef(ldReference)) {
        return(min(map_int(.gmComponents(ldReference), .ldRefRealNSamples)))
    }
    if (.ldIsMatrixRef(ldReference) && !is.null(mixtureWeights(ldReference))) {
        abort(c(
            str_c(
                "A MAC cutoff counts one panel's alleles, and the LdMatrix ",
                "is a mixture."
            ),
            i = "Use `mafCutoff` instead."
        ))
    }
    .ldRefRequireNSamples(ldReference, "A MAC cutoff")
}

# The subset of `variantIds` whose LD-panel genotypes clear the cutoffs,
# returned in the caller's order.
#
# Variants absent from the panel are passed through untouched. Whether a
# missing variant is an error or is silently dropped belongs to `.ldFromRef`
# (its `onMissing`), and deciding it a second time here would let the two
# disagree about the same variant.
# @noRd
.panelVariantFilter <- function(
    ldReference,
    variantIds,
    panelFilterParam = PanelFilterParam(),
    label = ".panelVariantFilter"
) {
    # .panelCutoffs answers NULL for a filter that would keep everything,
    # which is also the cheap exit here.
    cutoffs <- .panelCutoffs(panelFilterParam)
    if (is.null(cutoffs)) {
        return(variantIds)
    }
    if (is.null(ldReference) || length(variantIds) == 0L) {
        return(variantIds)
    }
    matched <- .ldFromRefMatch(ldReference, variantIds, label, "drop")
    if (is.null(matched)) {
        return(variantIds)
    }
    drop <- .panelDropMask(ldReference, matched, cutoffs)
    variantIds[!is_in(variantIds, matched$keptIds[drop])]
}

# Row mask for the variants clearing the LD-panel cutoffs, reporting how many
# were dropped. A no-op when no cutoff is set.
# @noRd
.panelKeepMask <- function(variantIds, ldReference, cutoffs, label) {
    if (is.null(ldReference) || is.null(cutoffs)) {
        return(rep(TRUE, length(variantIds)))
    }
    kept <- .panelVariantFilter(ldReference, variantIds, cutoffs, label = label)
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
.panelCutoffs <- function(panelFilterParam = PanelFilterParam()) {
    # NULL fields are how the pipelines spell "not set"; normalise them to the
    # no-op values so the short-circuit below is the only place that decides
    # whether a filter is worth running.
    maf <- panelFilterParam$mafCutoff %||% 0
    mac <- panelFilterParam$macCutoff %||% 0
    imiss <- panelFilterParam$imissCutoff %||% 1
    if (maf <= 0 && mac <= 0 && imiss >= 1) {
        return(NULL)
    }
    list(mafCutoff = maf, macCutoff = mac, imissCutoff = imiss)
}
# ---------- LD reference: cross-pipeline compatibility check ----------

# Internal: assert that two LD references describe the same reference panel --
# the same samples (for genotype panels) or the same LD (for LdMatrix
# objects), and the same allele coding (A1/A2 in the same orientation) on the
# variants they share. The SNP label is not compared, so a pure
# chr-prefix difference does not fail; an allele swap still does, since it
# means the two panels code their LD in opposite directions.
#
# The two panels' variant SETS need not agree, and after independent QC they
# normally do not: each object's sketch is trimmed to its own span and its own
# surviving variants, so one LD reference used for a QTL and a GWAS collection
# yields two differently trimmed sketches. Downstream LD lookups are all
# id-matched, so a partial overlap is reported (once per session) rather than
# refused; no overlap at all is an error. Shared by causalInferencePipeline,
# colocPipeline, qtlEnrichmentPipeline and ctwasPipeline, each over the rows
# it combines (.requireMatchingRowLdRefs()).
#
# NULL handling:
#   nullPolicy = "qtl-required" (default): a NULL qtlLd skips the check; a
#     non-NULL qtlLd with a NULL gwasLd is an error. Used by cip / coloc
#     (enloc included) / ctwas.
#   nullPolicy = "lenient": a NULL on either side skips the check. Used by
#     qtlEnrichmentPipeline.
#
# `label` is the human-readable name of the QTL-side input (e.g.
# "twasWeights" or "fineMappingResult"); it is woven into the error
# messages when provided. `pipelineName` prefixes every error so the
# failure source remains discoverable.
# --- .requireMatchingLdRefs helpers -------------------------------------

# Null handling: TRUE if the caller should return early (qtl NULL, or gwas NULL
# under a lenient policy); errors when a non-NULL qtl faces a NULL gwas sketch.
.ldRefNullGuard <- function(qtlLd, gwasLd, pipelineName, label, nullPolicy) {
    if (is.null(qtlLd)) {
        return(TRUE)
    }
    if (is.null(gwasLd)) {
        if (nullPolicy == "lenient") {
            return(TRUE)
        }
        labelPart <- if (!is.null(label)) {
            glue("ldReference on `{label}` is non-NULL ")
        } else {
            "qtl ldReference is non-NULL "
        }
        msg <- glue(
            "{pipelineName}: {labelPart}but the GWAS ldReference is NULL."
        )
        abort(msg)
    }
    FALSE
}

# Both must be genotype panels. Panel SIZE is deliberately NOT compared: two
# objects that share one LD sketch stop carrying identical variant sets as soon
# as they are QC'd separately, because each object's sketch is trimmed to its
# own position span at load (`.subsetLdRefToRange`) and to its own surviving
# variants at the end of QC (`.subsetLdRefToIds`). A QTL and a GWAS collection
# put through `summaryStatsQc` independently therefore diverge by construction,
# which is the normal case rather than a mistake. Every downstream LD lookup
# matches by variant id (`.ldFromRef`), never by position in the panel, so
# divergent variant sets cost only the variants one side lacks. What would not
# be harmless -- two genuinely different reference panels -- is what
# `.ldRefCheckContent` catches, by sample set and by allele coding.
.ldRefCheckShape <- function(qtlLd, gwasLd, pipelineName, between) {
    if (!.ldIsReference(qtlLd) || !.ldIsReference(gwasLd)) {
        msg <- glue(
            "{pipelineName}: ldReference slots{between} must both be genotype ",
            "panels or LdMatrix objects for the cross-pipeline LD reference ",
            "check."
        )
        abort(msg)
    }
    kinds <- c(.ldRefKind(qtlLd), .ldRefKind(gwasLd))
    if (kinds[[1L]] != kinds[[2L]]) {
        msg <- glue(
            "{pipelineName}: one side's ldReference is {kinds[[1L]]} and ",
            "the other's {kinds[[2L]]}{between}; use the same LD reference ",
            "on both."
        )
        abort(msg)
    }
}

# What kind of LD reference `x` is, in words, for messages that compare two.
# @noRd
.ldRefKind <- function(x) {
    if (.ldIsMatrixRef(x)) {
        return("an LdMatrix")
    }
    if (.ldIsMixtureRef(x)) {
        return("a GenotypeMixture")
    }
    "a genotype panel"
}

# The (chrom, position, allele-pair) key the two panels are compared on. The
# allele pair is order-insensitive, so an A1/A2 swap keys to the SAME variant
# and is then reported as an allele-coding difference rather than looking like
# two unrelated variants that happen not to overlap. A panel carrying no allele
# columns falls back to chrom:position, which is all it can be keyed on.
# @noRd
.ldRefVariantKeys <- function(x) {
    gr <- .ldRefRanges(x)
    mc <- S4Vectors::mcols(gr)
    stem <- str_c(
        as.character(GenomicRanges::seqnames(gr)),
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
.ldRefA1 <- function(x) {
    as.character(S4Vectors::mcols(.ldRefRanges(x))$A1)
}

# A shared variant whose A1/A2 are swapped between the panels is the same
# variant coded in opposite directions, so the two panels' LD matrices disagree
# in sign wherever it appears. That is a different reference, not a trimming
# difference, and it is an error.
# @noRd
.ldRefCheckAlleleCoding <- function(
    qtlLd,
    gwasLd,
    qIdx,
    gIdx,
    pipelineName,
    between
) {
    qA1 <- .ldRefA1(qtlLd)[qIdx]
    gA1 <- .ldRefA1(gwasLd)[gIdx]
    if (length(qA1) == 0L || length(gA1) == 0L) {
        return(invisible(NULL))
    }
    nSwapped <- sum(qA1 != gA1, na.rm = TRUE)
    if (nSwapped > 0L) {
        msg <- glue(
            "{pipelineName}: {nSwapped} variant(s) present in both ",
            "ldReference panels carry swapped A1/A2 alleles{between}, so the ",
            "two panels code their LD in opposite directions; use the same ",
            "ldReference on both."
        )
        abort(msg)
    }
    invisible(NULL)
}

# A partial overlap is the expected result of QC-ing the two sides separately,
# so it is reported once per session rather than on every call: ctwasPipeline
# and colocboostPipeline run this check once per region / per bundle.
# @noRd
.ldRefReportOverlap <- function(nQ, nG, nShared, pipelineName, between) {
    if (nShared == nQ && nShared == nG) {
        return(invisible(NULL))
    }
    msg <- glue(
        "{pipelineName}: the two ldReference panels share {nShared} ",
        "variant(s) of {nQ} (QTL side) and {nG} (GWAS side){between}. LD is ",
        "looked up per variant, so only the shared ones contribute. ",
        "Differing variant sets are expected when the two sides were QC'd ",
        "separately."
    )
    warn(
        msg,
        .frequency = "once",
        .frequency_id = str_c("pecotmrLdReferenceOverlap-", pipelineName)
    )
    invisible(NULL)
}

# Compare the two panels variant by variant. A partial overlap passes; what
# fails is no overlap at all (two unrelated panels) or a shared variant whose
# alleles are swapped.
# @noRd
.ldRefCheckOverlap <- function(qtlLd, gwasLd, pipelineName, between) {
    qKey <- .ldRefVariantKeys(qtlLd)
    gKey <- .ldRefVariantKeys(gwasLd)
    # An emptied panel -- what a zero-variant object carries, since it
    # references no LD -- contradicts nothing, so there is nothing to compare.
    if (length(qKey) == 0L || length(gKey) == 0L) {
        return(invisible(NULL))
    }
    idx <- match(qKey, gKey)
    shared <- which(!is.na(idx))
    if (length(shared) == 0L) {
        msg <- glue(
            "{pipelineName}: the two ldReference panels share no variant",
            "{between} ({length(qKey)} vs {length(gKey)} variants); they ",
            "describe different LD references."
        )
        abort(msg)
    }
    .ldRefCheckAlleleCoding(
        qtlLd,
        gwasLd,
        shared,
        idx[shared],
        pipelineName,
        between
    )
    .ldRefReportOverlap(
        length(qKey),
        length(gKey),
        length(shared),
        pipelineName,
        between
    )
    invisible(list(qIdx = shared, gIdx = idx[shared]))
}

# Two LdMatrix objects describe the same LD reference when they agree on the
# LD among the variants they share -- what two QC'd subsets of one matrix do,
# and two unrelated matrices do not. Checked on an evenly spaced subset of the
# shared variants, as validity checks symmetry, so the cost does not grow with
# the overlap.
# @noRd
.ldRefCheckSameLd <- function(qtlLd, gwasLd, shared, pipelineName, between) {
    if (is.null(shared)) {
        return(invisible(NULL))
    }
    if (!.ldmAgreesOn(qtlLd, gwasLd, shared$qIdx, shared$gIdx)) {
        msg <- glue(
            "{pipelineName}: the two LdMatrix references disagree on the LD ",
            "among the variants they share{between}; they describe different ",
            "LD references."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Two mixtures describe the same LD reference when they average the same
# panels (by their samples) with the same weights.
# @noRd
.ldRefCheckSameMixture <- function(qtlLd, gwasLd, pipelineName, between) {
    if (!.gmSameMixture(qtlLd, gwasLd)) {
        msg <- glue(
            "{pipelineName}: the two GenotypeMixture references{between} ",
            "average different panels or weight them differently; use the ",
            "same LD reference on both."
        )
        abort(msg)
    }
    invisible(NULL)
}

# TRUE when two LdMatrix objects hold the same LD among the variants at
# `xIdx` in one and `yIdx` in the other (matched pairs, same allele coding).
# Compared on an evenly spaced subset, so the cost does not grow with the
# overlap.
# @noRd
.ldmAgreesOn <- function(x, y, xIdx, yIdx, tol = 1e-6) {
    n <- length(xIdx)
    if (n == 0L) {
        return(TRUE)
    }
    probe <- unique(round(seq(1, n, length.out = min(n, 200L))))
    xR <- ldMatrix(x)[xIdx[probe], xIdx[probe], drop = FALSE]
    yR <- ldMatrix(y)[yIdx[probe], yIdx[probe], drop = FALSE]
    max(abs(xR - yR)) <= tol
}

# Panels must describe the same reference panel: the same samples, and the same
# allele coding on the variants they share. Two LdMatrix objects must agree on
# the LD they share instead, having no samples to compare. Their variant SETS
# need not agree -- see `.ldRefCheckShape` for why.
.ldRefCheckContent <- function(qtlLd, gwasLd, pipelineName, between) {
    if (.ldIsMatrixRef(qtlLd)) {
        shared <- .ldRefCheckOverlap(qtlLd, gwasLd, pipelineName, between)
        return(
            .ldRefCheckSameLd(qtlLd, gwasLd, shared, pipelineName, between)
        )
    }
    if (.ldIsMixtureRef(qtlLd)) {
        .ldRefCheckSameMixture(qtlLd, gwasLd, pipelineName, between)
        return(.ldRefCheckOverlap(qtlLd, gwasLd, pipelineName, between))
    }
    qIds <- .ldRefSampleIds(qtlLd)
    gIds <- .ldRefSampleIds(gwasLd)
    if (!identical(qIds, gIds)) {
        msg <- glue(
            "{pipelineName}: ldReference panels have different sample sets",
            "{between}; use the same ldReference on both."
        )
        abort(msg)
    }
    .ldRefCheckOverlap(qtlLd, gwasLd, pipelineName, between)
}

.requireMatchingLdRefs <- function(
    qtlLd,
    gwasLd,
    pipelineName,
    label = NULL,
    nullPolicy = c("qtl-required", "lenient")
) {
    nullPolicy <- arg_match(nullPolicy)
    if (.ldRefNullGuard(qtlLd, gwasLd, pipelineName, label, nullPolicy)) {
        return(invisible(NULL))
    }
    between <- if (!is.null(label)) {
        glue(" between `{label}` and gwas inputs", .trim = FALSE)
    } else {
        ""
    }
    .ldRefCheckShape(qtlLd, gwasLd, pipelineName, between)
    .ldRefCheckContent(qtlLd, gwasLd, pipelineName, between)
    invisible(NULL)
}

# .requireMatchingLdRefs() for the rows an analysis actually combines: QTL row
# qtlRows[k] with GWAS row gwasRows[k]. Rows on different references are fine
# as long as no combined pair puts two of them together. Each distinct pair of
# references is checked once, and an error names the first rows that combine
# them.
# @noRd
#' @importFrom purrr pwalk
.requireMatchingRowLdRefs <- function(
    qtl,
    gwas,
    qtlRows,
    gwasRows,
    pipelineName,
    label = NULL,
    nullPolicy = c("qtl-required", "lenient")
) {
    nullPolicy <- arg_match(nullPolicy)
    pairs <- distinct(
        tibble(
            qtlKey = .ldRefKeys(qtl)[qtlRows],
            gwasKey = .ldRefKeys(gwas)[gwasRows],
            qtlRow = qtlRows,
            gwasRow = gwasRows
        ),
        .data$qtlKey,
        .data$gwasKey,
        .keep_all = TRUE
    )
    pwalk(
        pairs,
        .requireMatchingPairLdRefs,
        qtlRefs = ldReference(qtl),
        gwasRefs = ldReference(gwas),
        pipelineName = pipelineName,
        label = label %||% "qtl input",
        nullPolicy = nullPolicy
    )
    invisible(NULL)
}

# @noRd
.requireMatchingPairLdRefs <- function(
    qtlKey,
    gwasKey,
    qtlRow,
    gwasRow,
    qtlRefs,
    gwasRefs,
    pipelineName,
    label,
    nullPolicy
) {
    .requireMatchingLdRefs(
        .ldRefForKey(qtlKey, qtlRefs),
        .ldRefForKey(gwasKey, gwasRefs),
        pipelineName = pipelineName,
        label = glue("{label} row {qtlRow} (with gwas row {gwasRow})"),
        nullPolicy = nullPolicy
    )
}

# .requireMatchingRowLdRefs() over the QTL / GWAS row pairs whose regions
# overlap: the pairs an analysis that brings the two sides together region by
# region combines.
# @noRd
.requireOverlapLdRefs <- function(
    qtl,
    gwas,
    pipelineName,
    label = NULL,
    nullPolicy = c("qtl-required", "lenient")
) {
    pairs <- .overlappingRowPairs(qtl, gwas)
    .requireMatchingRowLdRefs(
        qtl,
        gwas,
        pairs$qtlRows,
        pairs$gwasRows,
        pipelineName = pipelineName,
        label = label,
        nullPolicy = nullPolicy
    )
}

# The (QTL row, GWAS row) pairs whose element spans overlap. Both sides are
# narrowed to the chromosomes they share first: GenomicRanges warns about
# sequence levels the two have no part of in common, which would be noise for
# collections on different chromosomes rather than a sign of anything wrong.
# @noRd
.overlappingRowPairs <- function(qtl, gwas) {
    q <- .rowSpans(qtl)
    g <- .rowSpans(gwas)
    shared <- intersect(
        as.character(GenomicRanges::seqnames(q$spans)),
        as.character(GenomicRanges::seqnames(g$spans))
    )
    if (length(shared) == 0L) {
        return(list(qtlRows = integer(0), gwasRows = integer(0)))
    }
    qOn <- is_in(as.character(GenomicRanges::seqnames(q$spans)), shared)
    gOn <- is_in(as.character(GenomicRanges::seqnames(g$spans)), shared)
    hits <- GenomicRanges::findOverlaps(
        .onSeqlevels(q$spans[qOn], shared),
        .onSeqlevels(g$spans[gOn], shared),
        ignore.strand = TRUE
    )
    list(
        qtlRows = q$rows[qOn][S4Vectors::queryHits(hits)],
        gwasRows = g$rows[gOn][S4Vectors::subjectHits(hits)]
    )
}

# Which of `query`'s ranges overlap any of `subject`'s, compared on the
# chromosomes the two share -- GenomicRanges warns about sequence levels two
# objects have none of in common, which is noise for ranges that simply sit
# on different chromosomes.
# @noRd
.overlapsOnShared <- function(query, subject) {
    query <- .locusRanges(query)
    subject <- .locusRanges(subject)
    qChr <- as.character(GenomicRanges::seqnames(query))
    sChr <- as.character(GenomicRanges::seqnames(subject))
    shared <- intersect(qChr, sChr)
    hit <- rep(FALSE, length(query))
    if (length(shared) == 0L) {
        return(hit)
    }
    qOn <- is_in(qChr, shared)
    replace(
        hit,
        qOn,
        IRanges::overlapsAny(
            .onSeqlevels(query[qOn], shared),
            .onSeqlevels(subject[is_in(sChr, shared)], shared)
        )
    )
}

# Each non-empty row's span, with the row it belongs to.
# @noRd
.rowSpans <- function(x) {
    rr <- range(x)
    list(
        spans = .locusRanges(unlist(rr, use.names = FALSE)),
        rows = rep(seq_along(rr), lengths(rr))
    )
}

# Where `gr`'s ranges sit, for comparing collections' ranges: strand and
# metadata dropped, the genome build kept, which is what stops ranges of two
# builds from meeting.
# @noRd
#' @importFrom BiocGenerics unstrand
.locusRanges <- function(gr) {
    GenomicRanges::granges(unstrand(gr))
}

# The one LD reference rows `rows` of `x` share, for an analysis that reads
# their statistics together against one LD matrix. Rows on different
# references are refused rather than read from one of them; `what` names the
# rows in the error.
# @noRd
.rowsSharedLdRef <- function(x, rows, what) {
    keys <- unique(.ldRefKeys(x)[rows])
    if (length(keys) > 1L) {
        abort(c(
            glue("{what} use {length(keys)} different LD references."),
            i = str_c(
                "Their statistics are read together against one LD matrix; ",
                "give them one LD reference."
            )
        ))
    }
    .rowLdRef(x, rows[[1L]])
}

# `gr` with exactly the sequence levels `levels`.
# @noRd
.onSeqlevels <- function(gr, levels) {
    GenomeInfoDb::keepSeqlevels(gr, levels, pruning.mode = "coarse")
}

# ---------- Variant lookups against an LD reference on disk ----------

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

    # The span of the variants on each chromosome scopes the reference lookup.
    spans <- variantsDf |>
        group_by(.data$chrom) |>
        summarise(start = min(.data$pos), end = max(.data$pos))

    # Use shared helper -- no genotype loading
    refInfo <- refVariantInfo(
        ldReferenceMetaFile,
        GenomicRanges::GRanges(
            spans$chrom,
            IRanges::IRanges(spans$start, spans$end)
        )
    )
    refKey <- str_c(refInfo$chrom, ":", refInfo$pos)

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

# ---------- LD matrix checks and repair ----------

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

# Single-linkage LD pruning: keep one representative per |cor| cluster.
# Single-linkage clusters cut at 1 - corThres are the connected components
# of the graph joining columns with |cor| >= corThres, so they are found from
# the columns' correlations a stripe of rows at a time
# (ldThresholdComponents()), never the full matrix.
.ldPruneHclust <- function(X, corThres, verbose) {
    p <- ncol(X)
    clusters <- ldThresholdComponents(
        .ldStandardize(X),
        corThres,
        .ldStripeSize()
    )
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
#' @return A \code{MethodOptions} record for
#'   \code{ldPruneByCorrelation(methodArgs =)}.
#' @seealso \code{\link{ldPruneByCorrelation}}
#' @examples
#' LdPruningOptions(slide.max.bp = 1e6)
#' @export
LdPruningOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            gdsobj = "the temporary GDS pecotmr writes from `X`",
            method = "fixed at 'corr' by this backend",
            ld.threshold = "the caller's `corThres`",
            verbose = "the caller's `verbose`"
        ),
        "LdPruningOptions"
    )
    .newMethodOptions(
        "SNPRelate::snpgdsLDpruning",
        defaults = list(),
        extra = extra,
        label = "LdPruningOptions",
        engine = "ldPruning"
    )
}

#' Prune columns by pairwise correlation (LD-style prune)
#'
#' Performs LD pruning using one of two backends. The default \code{"hclust"}
#' backend groups the columns as single-linkage hierarchical clustering on the
#' distance (1 - |cor|), cut at 1 - \code{corThres}, would, and keeps the
#' first column of each group. Those groups are the connected components
#' of the graph joining columns with |cor| >= \code{corThres}, found from the
#' correlations a stripe of rows at a time, so the full correlation matrix is
#' never formed. The \code{"snprelate"} backend delegates
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
#'     \item{\code{"hclust"}}{Single-linkage groups, found as the connected
#'       components of the |cor| >= \code{corThres} graph.}
#'     \item{\code{"snprelate"}}{Requires \pkg{SNPRelate} and \pkg{gdsfmt}.
#'       Creates a temporary GDS file and runs
#'       \code{SNPRelate::snpgdsLDpruning(method = "corr")}. GDS stores
#'       genotypes as the integer codes 0, 1 and 2, so fractional dosages are
#'       refused rather than rounded.}
#'   }
#' @param verbose Logical. If TRUE, print progress messages. Default FALSE.
#' @param methodArgs Extra arguments for
#'   \code{SNPRelate::snpgdsLDpruning()}, built with
#'   \code{\link{LdPruningOptions}} -- the window controls in particular.
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
#' @importFrom stats cor
#' @export
ldPruneByCorrelation <- function(
    X,
    corThres = 0.8,
    backend = c("hclust", "snprelate"),
    verbose = FALSE,
    methodArgs = LdPruningOptions()
) {
    backend <- arg_match(backend)
    .assertMethodOptions(methodArgs, "LdPruningOptions", "methodArgs")
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

# Write X, integer genotype codes, to a temporary GDS for SNPRelate.
.ldPruneSnprelateCreateGds <- function(tmpGds, X, snpNames, p) {
    genoInt <- `storage.mode<-`(X, "integer")
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
    .ldRequireIntegerGenotypes(X, "snprelate", "hclust")
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

# Correlation strategy: drop the most-connected column (random tie-break at 2).
# @noRd
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
    # Engine-call exception: bigstatsr::FBM.code256 builds a file-backed
    # matrix -- a data structure the clumping step needs, not an engine.
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
#'   chromosome for each variant, in any style (\code{"chr1"} or \code{1}).
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
    # Engine-call exception: bigsnpr::snp_clumping is an LD algorithm, and
    # ld.R is the LD module that owns it -- neither a pipeline nor a wrapper,
    # and the call's only consumer. A one-function ldWrapper.R would separate
    # the call from the module whose job it is.
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
        infos.chr = as.integer(chromOrder(chr)),
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
# The guard is the reason this is a function rather than two lines at the
# call site: a single-variant block has no correlation structure, so its LD
# is the identity, and open-coding that in every caller is how one of them
# ends up omitting it.
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
        seq_along(.ldRefRanges(panel))
    } else {
        snpIdx
    }
    geno <- .ldRefDosage(panel, idx)
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
    fmt <- .ldRefFormat(X) %||% "held in memory"
    if (fmt != "gds") {
        abort(glue(
            "computeLd(onDisk = TRUE) needs a GDS-backed panel; this one is ",
            "'{fmt}'. Only GDS exposes an on-disk LD routine."
        ))
    }
    if (.withGdsFile(.ldRefReadPath(X), .gdsFileFormat) != "SNP_ARRAY") {
        abort(glue(
            "computeLd(onDisk = TRUE) reads SNPRelate's SNP GDS files; ",
            "{.ldRefReadPath(X)} is a SeqArray GDS file. Use ",
            "computeLd(onDisk = FALSE), which reads its dosages."
        ))
    }
    idx <- if (is.null(snpIdx)) {
        seq_along(.ldRefRanges(X))
    } else {
        snpIdx
    }
    # The GDS on-disk LD routine reads the file directly, which has no assay
    # equivalent: it is handed the file and the variants' own ids.
    R <- .ldCleanCorrelation(
        .computeBlockLdGds(.ldRefReadPath(X), .ldRefVariantIds(X)[idx])
    )
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
    # Engine-call exception: Rfast::cora stands in for base cor() (see the
    # backend note below) -- a numeric primitive, not an analysis engine.
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
#'   representation first. Both store genotypes as the integer codes 0, 1 and
#'   2, so they refuse fractional dosages (imputed, or a stochastic sketch's)
#'   rather than round them; \code{"internal"} correlates dosages as they
#'   are.
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
#' panel <- ldReference(qtlSumStatsExample)[[1L]]
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

# SNPRelate (GDS) and snpStats (SnpMatrix) hold a genotype as a 2-bit code,
# so a dosage outside 0, 1 and 2 -- an imputed or projected one -- would be
# rounded, and its LD changed, without a word. They take genotype codes only;
# dosages go to the backend that correlates them as they are.
# @noRd
.ldRequireIntegerGenotypes <- function(X, backend, instead) {
    bad <- !is.na(X) & !is_in(X, c(0, 1, 2))
    if (!any(bad)) {
        return(invisible(TRUE))
    }
    abort(c(
        glue(
            "backend = \"{backend}\" stores genotypes as the codes 0, 1 and ",
            "2, and {sum(bad)} value(s) here are not one of them."
        ),
        i = glue(
            "Rounding them would change their LD. Use backend = ",
            "\"{instead}\", which takes dosages as they are."
        )
    ))
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

    .ldRequireIntegerGenotypes(X, "snprelate", "internal")
    tmpGds <- tempfile(fileext = ".gds")
    on.exit(unlink(tmpGds), add = TRUE)

    # Integer genotype codes 0/1/2, with 3L, the GDS missing code, for NA.
    X_int <- replace(`storage.mode<-`(X, "integer"), is.na(X), 3L)

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
    .ldRequireIntegerGenotypes(X, "snpstats", "internal")

    # snpStats expects counts of the B allele as raw codes: 1=AA, 2=AB, 3=BB,
    # 0=NA pecotmr dosage is ALT count (0/1/2), so map: 0->1, 1->2, 2->3, NA->0
    coded <- replace(X + 1L, is.na(X), 0L)
    sm <- new("SnpMatrix", `storage.mode<-`(coded, "raw"))

    raw <- as.matrix(snpStats::ld(sm, stats = "R", depth = ncol(X) - 1L))
    # snpStats::ld returns a sparse-like matrix; ensure full dense
    .ldCleanCorrelation(raw)
}

#' @rdname PanelFilterParam
#' @aliases PanelFilterParam-class
#' @exportClass PanelFilterParam
setClass(
    "PanelFilterParam",
    contains = "MethodParam",
    slots = c(
        mafCutoff = "numeric",
        macCutoff = "numeric",
        imissCutoff = "numeric"
    )
)

#' @title LD Reference Panel Filtering Options
#' @description Which of an LD panel's variants are kept before the panel is
#'   used, as one checked bundle. Narrower than
#'   \code{\link{GenotypeFilterParam}}: a panel is not a study's genotypes, so
#'   there is no sample restriction and no variance cutoff, and the
#'   missingness default is \code{1} (no filter) rather than \code{0}.
#' @param mafCutoff Minor-allele-frequency floor. Default \code{0}.
#' @param macCutoff Minor-allele-count floor; the stricter of this and
#'   \code{mafCutoff} applies, using the panel's real sample count -- for a
#'   stochastic LD sketch, the one its \code{.afreq} records, not its
#'   pseudo-samples. Default \code{0}.
#' @param imissCutoff Per-variant missingness ceiling. Default \code{1}, which
#'   filters nothing and lets the allele-frequency sidecar be read instead of
#'   materializing dosage.
#' @return A \code{PanelFilterParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{GenotypeFilterParam}},
#'   \code{\link{SumstatsFilterParam}}
#' @examples
#' PanelFilterParam(mafCutoff = 0.001)
#' @export
PanelFilterParam <- function(
    mafCutoff = 0,
    macCutoff = 0,
    imissCutoff = 1
) {
    new(
        "PanelFilterParam",
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        imissCutoff = imissCutoff
    )
}
