# =============================================================================
# GenotypeHandle S4 class
# -----------------------------------------------------------------------------
# Lazy file handle for genotype data (PLINK2, BGEN). Opens the file for
# metadata only (sample IDs and SNP info); dosage extraction is deferred until
# extractBlockGenotypes() is called. A VCF, GDS or PLINK1 file is read as a
# panel directly, through a GenotypeDosageSeed (genotypeIo.R), not through a
# handle.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Genotype File Handle
#' @description S4 container holding a path + format + metadata for lazy
#'   genotype access. Supports PLINK2 (.pgen/.pvar/.psam) and BGEN (.bgen
#'   with its .bgi index); a VCF, GDS or PLINK1 file is read as a panel by
#'   \code{readGenotypes()} without one.
#' @slot path Character, file path. For a one-file-per-chromosome handle this is
#'   the chrom-meta file path (a display/provenance value); the per-chromosome
#'   payload files live in \code{chromPaths}.
#' @slot format Character, \code{"plink2"} or \code{"bgen"}.
#' @slot snpInfo data.frame, SNP metadata read from the index/sidecar, its
#'   \code{CHR} in UCSC style whatever the file writes. For a sharded handle
#'   this is the union across chromosomes (row-bound in the order the shards
#'   were supplied), so \code{snpIdx} stays a single global index space.
#' @slot nSamples Integer, number of samples.
#' @slot sampleIds Character vector of sample identifiers.
#' @slot pgenPtr Opaque pointer for PLINK2 reader state (NULL otherwise).
#' @slot chromPaths Named character vector mapping each chromosome, in UCSC
#'   style (e.g. \code{"chr21"}, \code{"chrX"}), to the per-chromosome
#'   payload path/prefix. Empty
#'   (\code{character(0)}) for a single-file handle; non-empty marks a
#'   one-file-per-chromosome (sharded) handle whose extraction is routed by
#'   chromosome.
#' @details Internal: the handle is the seed layer behind a genotype
#'   panel. Public code obtains a panel from \code{readGenotypes()},
#'   which is a \code{RangedSummarizedExperiment}.
#' @keywords internal
setClass(
    "GenotypeHandle",
    representation(
        path = "character",
        format = "character",
        snpInfo = "data.frame",
        nSamples = "integer",
        sampleIds = "character",
        pgenPtr = "ANY",
        chromPaths = "character"
    ),
    prototype = prototype(
        chromPaths = character(0)
    ),
    validity = function(object) {
        validFormats <- c("plink2", "bgen")
        nm <- names(object@chromPaths)
        badChromPaths <- length(object@chromPaths) > 0L &&
            (is.null(nm) || any(str_length(nm) == 0L) || anyDuplicated(nm))
        errors <- c(
            if (length(object@path) != 1L) {
                "'path' must be a single character string"
            },
            if (!is_in(object@format, validFormats)) {
                str_c(
                    "'format' must be one of: ",
                    str_flatten(validFormats, ", ")
                )
            },
            if (badChromPaths) {
                str_c(
                    "'chromPaths' must be a uniquely-named ",
                    "character vector (names = chromosomes)"
                )
            }
        )
        if (length(errors) == 0) TRUE else errors
    }
)

#' @rdname show-methods
#' @export
setMethod("show", "GenotypeHandle", function(object) {
    cat(glue("GenotypeHandle [{object@format}]\n", .trim = FALSE))
    chromPaths <- .genotypeChromPaths(object)
    if (length(chromPaths) > 0L) {
        cat(glue("  Chrom-meta: {object@path}\n", .trim = FALSE))
        cat(glue(
            "  {length(chromPaths)} per-chromosome files: ",
            "{str_flatten(names(chromPaths), ', ')}\n",
            .trim = FALSE
        ))
    } else {
        cat(glue("  Path: {object@path}\n", .trim = FALSE))
    }
    cat(glue(
        "  {object@nSamples} samples, {nrow(object@snpInfo)} SNPs\n",
        .trim = FALSE
    ))
})

#' @title Create a GenotypeHandle Object
#' @description Construct a \code{GenotypeHandle} from one of several input
#'   forms. Exactly one of the following must be specified:
#'   \describe{
#'     \item{\code{path}}{A single file path with a recognized extension:
#'       \code{.bgen} or \code{.pgen}. Format is auto-detected from the
#'       extension. A BGEN file needs its \code{.bgi} index beside it, and
#'       its sample ids in the file or in a \code{.sample} file beside it. A
#'       VCF, GDS or PLINK1 file is refused: it is read as a panel by
#'       \code{readGenotypes()}.}
#'     \item{\code{plink2Prefix}}{A path prefix; the constructor appends
#'       \code{.pgen}, \code{.pvar} (or \code{.pvar.zst}), and \code{.psam}
#'       to locate the triplet.}
#'     \item{\code{pgen} + \code{pvar} + \code{psam}}{Explicit PLINK2
#'       triplet. The three files must share a stem; if they don't, use
#'       \code{plink2Prefix} or arrange symlinks at a common stem.}
#'     \item{\code{genoMeta}}{One genotype file per chromosome. Either a path
#'       to a whitespace/TSV meta file whose first column is the chromosome
#'       (\code{#chr}) and second column the per-chromosome payload
#'       (a \code{.pgen}/\code{.bgen} file or a PLINK2 prefix;
#'       relative paths resolve against the meta file's directory), or a
#'       named character vector mapping chromosome to payload. All shards
#'       must share one format and identical sample IDs in the same order.
#'       Extraction is routed to the correct file by the requested region's
#'       chromosome.}
#'   }
#'   The constructor opens the file for metadata only (sample IDs and SNP
#'   info); dosage extraction is deferred until \code{extractBlockGenotypes()}
#'   is called.
#'
#' @param path Single file path (.bgen/.pgen), or \code{NULL}.
#' @param plink2Prefix Path prefix for a PLINK2 triplet, or \code{NULL}.
#' @param pgen,pvar,psam Explicit paths to the PLINK2 triplet, or all
#'   \code{NULL}.
#' @param ldMeta Path to an LD-meta TSV file (columns \code{chrom},
#'   \code{start}, \code{end}, \code{path}; the \code{path} column may be
#'   comma-separated as \code{ld_file,bim_file}). Requires \code{block}. The
#'   constructor resolves the row covering \code{block}, then delegates to the
#'   appropriate file-based handler. When the resolved row points at PLINK2
#'   or BGEN files, the corresponding reader is used; a row pointing
#'   at a pre-computed LD matrix (\code{.cor.xz}) is refused, as pecotmr reads
#'   no pre-computed LD from disk -- wrap such a matrix in an
#'   \code{\link{LdMatrix}} instead.
#' @param block The \code{GRanges} window to look up in \code{ldMeta}; the
#'   whole file of the row covering it is opened.
#' @param genoMeta One-file-per-chromosome specification: a path to a
#'   \code{#chr,path} meta file or a named character vector (names =
#'   chromosomes, values = payload paths/prefixes). Optionally pass
#'   \code{format} to force a single backend for every shard.
#' @param chroms Optional character vector of chromosomes. With \code{genoMeta}
#'   (a sharded, one-file-per-chromosome panel) only the shards whose chromosome
#'   is listed are read, skipping the rest -- an I/O optimisation when the panel
#'   is genome-wide but only a few chromosomes are needed. Chromosome labels may
#'   be in any style (\code{"chr1"}/\code{"1"} match, \code{23}/\code{X}
#'   etc.). If none of the requested chromosomes are present in the panel, every
#'   shard is read (so the caller's own absence check can report the mismatch).
#'   Only meaningful with \code{genoMeta}; supplying it with any other source is
#'   an error (a single-file panel has no per-chromosome shards to skip).
#' @param format Character, \code{"plink2"} or \code{"bgen"}, forcing the
#'   backend of \code{path} or of every \code{genoMeta} shard. If NULL,
#'   inferred from the file extension.
#' @return A \code{GenotypeHandle} object.
#' @keywords internal
#' @rdname GenotypeHandle-class
GenotypeHandle <- function(
    path = NULL,
    plink2Prefix = NULL,
    pgen = NULL,
    pvar = NULL,
    psam = NULL,
    ldMeta = NULL,
    block = NULL,
    genoMeta = NULL,
    chroms = NULL,
    format = NULL
) {
    flags <- .ghValidateArgs(
        bed = NULL,
        bim = NULL,
        fam = NULL,
        pgen = pgen,
        pvar = pvar,
        psam = psam,
        ldMeta = ldMeta,
        block = block
    )
    sources <- .ghResolveSources(
        path = path,
        plink1Prefix = NULL,
        plink2Prefix = plink2Prefix,
        ldMeta = ldMeta,
        genoMeta = genoMeta,
        chroms = chroms,
        flags = flags
    )
    switch(
        names(sources)[which(sources)],
        path = .readGenotypeHandle(path, format = format),
        plink2Prefix = .makePlink2Handle(plink2Prefix),
        plink2Triplet = .genotypeHandleFromPlink2Triplet(pgen, pvar, psam),
        ldMeta = .genotypeHandleFromLdMeta(ldMeta, block),
        genoMeta = .genotypeHandleFromChromMeta(genoMeta, chroms, format)
    )
}

# Validate the bed/bim/fam + pgen/pvar/psam triplet completeness and the
# ldMeta<->block coupling. Returns list(bedComplete, pgenComplete).
# @noRd
.ghValidateArgs <- function(bed, bim, fam, pgen, pvar, psam, ldMeta, block) {
    bedComplete <- .ghTrioComplete(bed, bim, fam, "bed/bim/fam")
    pgenComplete <- .ghTrioComplete(pgen, pvar, psam, "pgen/pvar/psam")
    if (!is.null(ldMeta) && is.null(block)) {
        abort("`ldMeta` requires `block`, the GRanges window to look up.")
    }
    if (is.null(ldMeta) && !is.null(block)) {
        abort("`block` is only meaningful when `ldMeta` is supplied.")
    }
    if (!is.null(block)) {
        .assertRegionGRanges(block, "`block`")
    }
    list(bedComplete = bedComplete, pgenComplete = pgenComplete)
}

# TRUE when a file triplet is fully supplied; error when partially supplied.
# @noRd
.ghTrioComplete <- function(a, b, c, label) {
    given <- !is.null(a) || !is.null(b) || !is.null(c)
    complete <- !is.null(a) && !is.null(b) && !is.null(c)
    if (given && !complete) {
        msg <- glue(
            "If specifying the {label} triplet, all three must be provided."
        )
        abort(msg)
    }
    complete
}

# Build the exactly-one-source indicator vector; error unless exactly one input
# source was supplied, and gate `chroms` to the genoMeta path. Because exactly
# one entry is TRUE, GenotypeHandle()'s switch on its name is total.
# @noRd
.ghResolveSources <- function(
    path,
    plink1Prefix,
    plink2Prefix,
    ldMeta,
    genoMeta,
    chroms,
    flags
) {
    sources <- c(
        path = !is.null(path),
        plink1Prefix = !is.null(plink1Prefix),
        plink2Prefix = !is.null(plink2Prefix),
        plink1Triplet = flags$bedComplete,
        plink2Triplet = flags$pgenComplete,
        ldMeta = !is.null(ldMeta),
        genoMeta = !is.null(genoMeta)
    )
    if (sum(sources) != 1L) {
        nSrc <- sum(sources)
        msg <- glue(
            "Exactly one of `path`, `plink1Prefix`, `plink2Prefix`, the ",
            "bed/bim/fam triplet, the pgen/pvar/psam triplet, `ldMeta`, or ",
            "`genoMeta` must be specified (got {nSrc})."
        )
        abort(msg)
    }
    if (!is.null(chroms) && !sources[["genoMeta"]]) {
        msg <- glue(
            "`chroms` restricts which per-chromosome shards are read and is ",
            "only supported with `genoMeta` (a single-file panel has no ",
            "shards to skip)."
        )
        abort(msg)
    }
    sources
}

.genotypeHandleFromLdMeta <- function(ldMeta, region) {
    .ghLdPathToHandle(.ghResolveLdPath(ldMeta, region))
}

# Resolve the single genotype-payload path for `region` from an LD meta,
# erroring on no-coverage, multi-row spans, or a pre-computed .cor(.xz) matrix.
# @noRd
.ghResolveLdPath <- function(ldMeta, region) {
    ldPaths <- regionalLdMeta(ldMeta, region)$intersections$LD_file_paths
    if (length(ldPaths) == 0L) {
        regionStr <- str_flatten(as.character(region), ", ")
        msg <- glue(
            "ldMeta: no LD-meta row covers region {regionStr} ",
            "in {ldMeta}."
        )
        abort(msg)
    }
    if (length(ldPaths) > 1L) {
        regionStr <- str_flatten(as.character(region), ", ")
        msg <- glue(
            "ldMeta: region {regionStr} spans multiple LD-meta rows; ",
            "only a region within one row is resolved. Restrict the region ",
            "to a single LD block."
        )
        abort(msg)
    }
    ldPath <- ldPaths[[1L]]
    if (str_detect(ldPath, regex("\\.cor(\\.xz)?$", ignore_case = TRUE))) {
        regionStr <- str_flatten(as.character(region), ", ")
        msg <- glue(
            "ldMeta: the LD-meta row for region {regionStr} points at a ",
            "pre-computed correlation matrix ({ldPath}). pecotmr does not ",
            "read pre-computed LD from disk; wrap the matrix in an LdMatrix ",
            "instead. An LD-meta row names genotypes (VCF/GDS/BGEN/PLINK)."
        )
        abort(msg)
    }
    ldPath
}

# The handle of the genotype file or PLINK prefix an LD-meta row names, for a
# format a GenotypeHandle reads.
# @noRd
.ghLdPathToHandle <- function(ldPath) {
    .readGenotypeHandle(
        ldPath,
        format = .requireGenotypeFormat(ldPath, "ldMeta")
    )
}

# The panel of the genotype file or PLINK prefix an LD-meta row names, in any
# format.
# @noRd
.ghLdPathToPanel <- function(ldPath, vcfArgs = list()) {
    .readGenotypePanel(
        ldPath,
        format = .requireGenotypeFormat(ldPath, "ldMeta"),
        vcfArgs = vcfArgs
    )
}

# The stem the three files of a PLINK1 triplet share, which is how its panel
# is read.
# @noRd
#' @importFrom checkmate assertString
.plink1TripletStem <- function(bed, bim, fam) {
    assertString(bed)
    assertString(bim)
    assertString(fam)
    stems <- c(
        bed = .genotypeStem(bed),
        bim = .genotypeStem(bim),
        fam = .genotypeStem(fam)
    )
    if (n_distinct(stems) != 1L) {
        stemList <- str_c("  ", names(stems), ": ", stems, collapse = "\n")
        msg <- glue(
            "`bed`, `bim`, and `fam` must share a common path stem. ",
            "Got:\n{stemList}\n",
            "If your files are at different paths, either rename them to ",
            "share a stem or arrange symlinks at a common prefix and pass ",
            "`plink1Prefix` instead.",
            .trim = FALSE
        )
        abort(msg)
    }
    unname(stems[1L])
}

.genotypeHandleFromPlink2Triplet <- function(pgen, pvar, psam) {
    for (f in list(pgen = pgen, pvar = pvar, psam = psam)) {
        if (!is.character(f) || length(f) != 1L) {
            abort("Each of `pgen`, `pvar`, `psam` must be a single file path.")
        }
    }
    stems <- c(
        pgen = .genotypeStem(pgen),
        pvar = .genotypeStem(pvar),
        psam = .genotypeStem(psam)
    )
    if (n_distinct(stems) != 1L) {
        stemList <- str_c("  ", names(stems), ": ", stems, collapse = "\n")
        msg <- glue(
            "`pgen`, `pvar`, and `psam` must share a common path stem. ",
            "Got:\n{stemList}\n",
            "If your files are at different paths, either rename them to ",
            "share a stem or arrange symlinks at a common prefix and pass ",
            "`plink2Prefix` instead.",
            .trim = FALSE
        )
        abort(msg)
    }
    .makePlink2Handle(unname(stems[1L]))
}

# ---------------------------------------------------------------------------
# One-file-per-chromosome (sharded) handle support
# ---------------------------------------------------------------------------

# Chromosome labels are spelled in UCSC style (withChrPrefix(), variantId.R)
# as a meta or a shard is read, so the @chromPaths routing keys and the
# @snpInfo$CHR lookups compare exactly (and X/Y/M survive).

# Per-chromosome shard map, tolerant of GenotypeHandle objects deserialized
# from before the `chromPaths` slot existed (e.g. an RDS saved by an older
# pecotmr). Such objects have no `chromPaths` slot, so a direct `@` access
# errors; treat them as single-file handles.
#' @keywords internal
#' @importFrom purrr possibly
.genotypeChromPaths <- function(handle) {
    # A handle with no chromosome paths is an ordinary outcome, not a fault.
    possibly(chromPaths, otherwise = character(0))(handle)
}

# Case-insensitive match of the first of `aliases` present in `cols`; falls
# back to the `default` column position when none of the aliases are present.
.metaMatchCol <- function(cols, aliases, default) {
    hit <- which(is_in(str_to_lower(cols), str_to_lower(aliases)))
    if (length(hit) > 0L) hit[[1L]] else default
}

# Enforce one payload per chromosome key: LD sketches / sharded genotypes are
# genome-wide or per-chromosome, never sub-chromosomal, so a chromosome mapping
# to two distinct payloads is an error. Returns the data.frame unchanged.
.metaCheckUniqueChrom <- function(df) {
    byChr <- split(df$path, df$chrom)
    multi <- names(byChr)[map_lgl(byChr, .ghMultiplePayloads)]
    if (length(multi) > 0L) {
        multiStr <- str_flatten(multi, ", ")
        msg <- glue(
            "genoMeta: chromosome(s) {multiStr} map to ",
            "multiple genotype payloads; each chromosome must map to ",
            "exactly one payload (no sub-chromosomal blocks)."
        )
        abort(msg)
    }
    df
}

# Parse the genoMeta input into a data.frame(chrom, path). Accepts either a
# path to a meta file (whitespace- or tab-delimited, with header) or a named
# character vector (names = chromosomes). In a meta file the chromosome and
# path columns are matched BY NAME -- chromosome from #chr/#chrom/chr/chrom and
# path from path/payload/prefix/genotype -- so extra columns (e.g. legacy
# start/end) and non-standard column order are tolerated; the first two columns
# are used only when the header carries none of those aliases. Relative payload
# paths are resolved against the meta file's own directory.
#' @keywords internal
.parseChromMeta <- function(genoMeta) {
    if (.isChromMetaFile(genoMeta)) {
        return(.parseChromMetaFile(genoMeta))
    }
    if (
        is.character(genoMeta) &&
            length(genoMeta) >= 1L &&
            !is.null(names(genoMeta))
    ) {
        return(.metaCheckUniqueChrom(tibble(
            chrom = withChrPrefix(names(genoMeta)),
            path = unname(genoMeta)
        )))
    }
    msg <- glue(
        "genoMeta: expected a path to a `#chr,path` meta ",
        "file or a named character vector (names = chromosomes, ",
        "values = paths)."
    )
    abort(msg)
}

# TRUE when `genoMeta` is a single existing-file path (vs a named vector).
# @noRd
.isChromMetaFile <- function(genoMeta) {
    is.character(genoMeta) &&
        length(genoMeta) == 1L &&
        is.null(names(genoMeta)) &&
        file.exists(genoMeta)
}

# Parse a `#chr,path` meta file into a (chrom, path) data.frame, resolving
# relative payload paths against the meta file's directory.
# @noRd
.parseChromMetaFile <- function(genoMeta) {
    meta <- read_table(
        genoMeta,
        col_names = TRUE,
        col_types = cols(.default = col_character()),
        comment = ""
    )
    if (ncol(meta) < 2L) {
        msg <- glue(
            "genoMeta: meta file '{genoMeta}' must have ",
            "at least 2 columns (chromosome, path)."
        )
        abort(msg)
    }
    chromCol <- .metaMatchCol(
        names(meta),
        c("#chr", "#chrom", "chr", "chrom"),
        default = 1L
    )
    pathCol <- .metaMatchCol(
        names(meta),
        c("path", "payload", "prefix", "genotype"),
        default = 2L
    )
    base <- dirname(normalizePath(genoMeta))
    .metaCheckUniqueChrom(tibble(
        chrom = withChrPrefix(meta[[chromCol]]),
        path = .metaResolvePaths(as.character(meta[[pathCol]]), base)
    ))
}

# Resolve relative payload paths against `base` (absolute/existing kept as-is).
# @noRd
.metaResolvePaths <- function(pth, base) {
    map_chr(pth, .ghResolvePath, base = base)
}

# The format of one shard payload: `format` when it is forced, else told by
# the file extension, falling back to probing a PLINK prefix's sidecars.
# Nothing is read, so a set's shards can be checked to share one format
# before any of them is.
# @noRd
.shardFormat <- function(p, format = NULL) {
    format %||% .requireGenotypeFormat(p, "genoMeta")
}

# Build a single-file GenotypeHandle for one shard payload, of a format a
# handle reads. `format` (optional) forces a backend; otherwise it is detected
# as .shardFormat() detects it.
#' @keywords internal
.resolveGenotypeShard <- function(p, format = NULL) {
    .readGenotypeHandle(p, format = .shardFormat(p, format))
}

# The panel of a per-chromosome meta. A VCF, GDS or PLINK1 set is its
# shards' panels stacked (.stackDosagePanels()); a set in a format a handle
# reads, one sharded handle's panel. `chroms` (optional) restricts the read
# to the shards for those chromosomes, as for .genotypeHandleFromChromMeta().
# @noRd
.chromMetaPanel <- function(genoMeta, chroms = NULL, format = NULL) {
    parsed <- .parseChromMeta(genoMeta)
    if (nrow(parsed) == 0L) {
        abort(
            "genoMeta: no chromosomes found in the meta input."
        )
    }
    selected <- .chromMetaSelect(parsed, chroms)
    formats <- map_chr(selected$path, .shardFormat, format = format)
    shared <- .chromMetaCheckFormats(formats)
    if (!is_in(shared, .DOSAGE_PANEL_FORMATS)) {
        return(.genotypeExperiment(
            .genotypeHandleFromChromMeta(genoMeta, chroms, format)
        ))
    }
    panels <- map2(selected$path, formats, .readDosagePanel)
    .withSourceProjection(.stackDosagePanels(
        panels,
        shared,
        .chromMetaPath(genoMeta),
        selected$path
    ))
}

# Assemble a sharded handle from a per-chromosome meta. Reads each shard's
# metadata via the existing single-file readers, validates a single shared
# format and identical sample IDs (same order, required for cross-shard
# cbind), and row-binds the per-shard snpInfo into one global index space.
# `chroms` (optional) restricts the read to the shards for those chromosomes:
# the other per-chromosome files are never opened, which is the I/O win when a
# genome-wide panel backs summary statistics on only a few chromosomes.
#' @keywords internal
.genotypeHandleFromChromMeta <- function(
    genoMeta,
    chroms = NULL,
    format = NULL
) {
    parsed <- .parseChromMeta(genoMeta)
    if (nrow(parsed) == 0L) {
        abort(
            "genoMeta: no chromosomes found in the meta input."
        )
    }
    selected <- .chromMetaSelect(parsed, chroms)
    shards <- map(selected$path, .resolveGenotypeShard, format = format)
    sharedFormat <- .chromMetaCheckFormats(map_chr(shards, .ghFormat))
    .chromMetaCheckSamples(shards, selected)
    unifiedSnpInfo <- bind_rows(map(shards, .ghSnpInfo))
    new(
        "GenotypeHandle",
        path = .chromMetaPath(genoMeta),
        format = sharedFormat,
        snpInfo = unifiedSnpInfo,
        nSamples = nSamples(shards[[1L]]),
        sampleIds = sampleIds(shards[[1L]]),
        pgenPtr = NULL,
        chromPaths = .chromMetaPaths(shards)
    )
}

# Restrict a parsed meta to the requested chromosomes' shards so the rest are
# never read. If nothing matches (requested chromosomes absent from the panel)
# fall back to every shard, leaving the caller's own containment check to
# produce the usual diagnostic instead of a confusing empty-handle error here.
# @noRd
.chromMetaSelect <- function(parsed, chroms) {
    if (is.null(chroms)) {
        return(parsed)
    }
    keep <- is_in(parsed$chrom, withChrPrefix(chroms))
    if (any(keep)) filter(parsed, keep) else parsed
}

# Require a single shared format across shards' `formats`; returns it.
# @noRd
.chromMetaCheckFormats <- function(formats) {
    if (n_distinct(formats) != 1L) {
        formatStr <- str_flatten(unique(formats), ", ")
        msg <- glue(
            "genoMeta: all per-chromosome files must share ",
            "one format; got: {formatStr}."
        )
        abort(msg)
    }
    formats[[1L]]
}

# Require identical sample IDs (same order) across shards for cross-shard cbind.
# @noRd
.chromMetaCheckSamples <- function(shards, parsed) {
    sample0 <- shards[[1L]]@sampleIds
    for (i in seq_along(shards)[-1L]) {
        if (!identical(shards[[i]]@sampleIds, sample0)) {
            mismatchPath <- parsed$path[[i]]
            msg <- glue(
                "genoMeta: all per-chromosome files must ",
                "have identical sample IDs in the same order (mismatch ",
                "at '{mismatchPath}')."
            )
            abort(msg)
        }
    }
    invisible(NULL)
}

# One shard's chromosome -> path mapping.
# @noRd
.chromShardPaths <- function(shard) {
    chroms <- unique(snpInfo(shard)$CHR)
    set_names(rep(path(shard), length(chroms)), chroms)
}

# A chromosome may only come from one file; two shards claiming it means the
# panel is not split cleanly per chromosome.
# @noRd
.checkChromNotSplit <- function(ch, count) {
    if (count <= 1L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "genoMeta: chromosome '{ch}' appears ",
        "in more than one per-chromosome file."
    )
    abort(msg)
}

# Map each chromosome to its shard path, erroring if a chromosome spans files.
# @noRd
.chromMetaPaths <- function(shards) {
    if (length(shards) == 0L) {
        return(character(0))
    }
    chromPaths <- list_c(map(shards, .chromShardPaths))
    counts <- table(names(chromPaths))
    walk2(names(counts), as.integer(counts), .checkChromNotSplit)
    chromPaths
}

# The handle's recorded path: the normalized meta-file path, or a placeholder
# for a named-vector meta.
# @noRd
.chromMetaPath <- function(genoMeta) {
    if (.isChromMetaFile(genoMeta)) normalizePath(genoMeta) else "<chrom-meta>"
}

#' @rdname snpInfo
#' @keywords internal
setMethod("snpInfo", "GenotypeHandle", function(x) x@snpInfo)

#' @rdname genotypeFormat
#' @keywords internal
setMethod("genotypeFormat", "GenotypeHandle", function(x) x@format)

# Adopts BiocGenerics' `path` generic rather than defining a pecotmr one:
# a file-backed handle is exactly what that generic is for. Deliberately not
# exported -- GenotypeHandle is internal, so only package code dispatches
# here. The generic's first formal is `object`, not `x`.
#' @importFrom BiocGenerics path
#' @noRd
setMethod("path", "GenotypeHandle", function(object, ...) object@path)

#' @rdname chromPaths
#' @keywords internal
setMethod("chromPaths", "GenotypeHandle", function(x) x@chromPaths)

# A stem-based genotype format is named by a companion file rather than by
# the stem itself: plink1 writes <stem>.bed, plink2 writes <stem>.pgen.
# Returns the stem when one of those is bundled, NULL otherwise -- a
# single-file format (a .gds) is the stem, and the caller handles that.
#
# plink2 was missing here, so a .pgen panel could not be referenced
# portably at all and had to bake an absolute build-machine path.
# @noRd
.resolveBundledStem <- function(stem) {
    exts <- c(".bed", ".pgen")
    paths <- map_chr(exts, .bundledExtdataPath, stem = stem)
    hit <- which(str_length(paths) > 0L)
    if (length(hit) == 0L) {
        return(NULL)
    }
    first <- hit[[1L]]
    str_remove(paths[[first]], str_c("\\", exts[[first]], "$"))
}

# @noRd
.bundledExtdataPath <- function(ext, stem) {
    system.file("extdata", str_c(stem, ext), package = "pecotmr")
}

# Resolve a portable bundled-resource reference of the form
# "pecotmr://extdata/<stem>" to a concrete filesystem path via
# `system.file()`; any ordinary path is returned unchanged. This lets the
# packaged example objects (which cannot bake an install-specific absolute
# path into their serialized GenotypeHandle) carry a portable reference that
# resolves at extraction time, on whatever machine the package is installed.
.resolveGenotypeResourcePath <- function(path) {
    if (length(path) != 1L || is.na(path)) {
        return(path)
    }
    m <- str_match(path, "^pecotmr://extdata/(.+)$")[1L, ]
    if (is.na(m[[1L]])) {
        return(path)
    }
    stem <- m[[2L]]
    byStem <- .resolveBundledStem(stem)
    if (!is.null(byStem)) {
        return(byStem)
    }
    resolved <- system.file("extdata", stem, package = "pecotmr")
    if (str_length(resolved) == 0L) {
        msg <- glue(
            "cannot resolve bundled genotype resource '{stem}' under ",
            "inst/extdata/ (package 'pecotmr')."
        )
        abort(msg)
    }
    resolved
}

# The concrete filesystem path/stem to open for a handle's data. Resolves a
# bundled-resource reference; ordinary paths pass through. Used by the block
# extractors and LD-panel keying, so `path()` keeps returning the stored
# (portable) value for display/provenance while file access resolves it.
.genotypeReadPath <- function(handle) {
    .resolveGenotypeResourcePath(path(handle))
}

#' @rdname sampleIds
#' @keywords internal
setMethod("sampleIds", "GenotypeHandle", function(x) x@sampleIds)

#' @rdname pgenPtr
#' @keywords internal
setMethod("pgenPtr", "GenotypeHandle", function(x) x@pgenPtr)

# Shares the public `nSamples` page: the generic's @details already records
# that the handle is the seed layer behind a panel, and a @keywords internal
# here would mark that whole merged page internal.
#' @rdname nSamples
setMethod("nSamples", "GenotypeHandle", function(x) x@nSamples)

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# TRUE when one chromosome's payload set has more than one distinct path.
# @noRd
.ghMultiplePayloads <- function(p) {
    n_distinct(p) > 1L
}

# Resolve a payload path against `base` (absolute/existing paths pass through).
# @noRd
.ghResolvePath <- function(p, base) {
    if (str_detect(p, "^(/|[A-Za-z]:)") || file.exists(p)) {
        p
    } else {
        file.path(base, p)
    }
}

# The snpInfo slot of one genotype shard.
# @noRd
.ghSnpInfo <- function(h) {
    snpInfo(h)
}

# The format slot of one genotype shard.
# @noRd
.ghFormat <- function(h) {
    genotypeFormat(h)
}


# =============================================================================
# GhSeed -- a DelayedArray backend over a GenotypeHandle
#
# `extractBlockGenotypes(handle, snpIdx)` already reads ONLY the requested
# variants, so it is already the lazy interface a DelayedArray seed needs; this
# is a thin adapter, not a new reader. Wrapping it means a genotype matrix can
# sit in a SummarizedExperiment assay -- and so in a MultiAssayExperiment --
# without materialising, which is what lets the whole panel be described
# without being read.
#
# Orientation: pecotmr works in samples x variants, but Bioconductor assays are
# features x samples, and `assay(se, "dosage")` is already variants x samples.
# The seed therefore reports variants as rows. Getting this backwards returns a
# plausible-looking matrix and a silently transposed LD matrix, so the tests
# pin it against genotypes().
# =============================================================================

#' DelayedArray seed over a GenotypeHandle
#'
#' Lets genotype dosages back a \code{DelayedArray} that reads through the
#' handle on demand, so an assay can describe a panel far larger than memory.
#'
#' @slot handle The \code{\link{GenotypeHandle}} to read through.
#' @slot dm Integer length-2 dimensions (variants, samples).
#' @slot dn List of dimnames (variant ids, sample ids).
#' @name GhSeed-class
#' @keywords internal
setClass(
    "GhSeed",
    representation(handle = "GenotypeHandle", dm = "integer", dn = "list")
)

#' @rdname GhSeed-class
setMethod("dim", "GhSeed", function(x) x@dm)

#' @rdname GhSeed-class
setMethod("dimnames", "GhSeed", function(x) {
    # The sample axis is derived from the handle, not stored a second time.
    # handle@sampleIds is the single copy; storing dn[[2]] as well would
    # duplicate ~10k anonymous projection-dimension names in every serialised
    # sketch (see genotypeDelayedArray()). sampleIds() is character(0) on an
    # emptied handle, so this stays consistent with dm = c(nVar, 0).
    list(x@dn[[1L]], as.character(sampleIds(x@handle)))
})

# The handle a seed reads through. Internal accessor: GhSeed is not exported,
# so this adds no public surface, and it keeps the `@` next to the class.
# @noRd
.ghSeedHandle <- function(x) x@handle

#' @rdname GhSeed-class
#' @importFrom S4Arrays extract_array
setMethod("extract_array", "GhSeed", function(x, index) {
    # A NULL index means "everything along this dimension", not "nothing".
    # Reading it as nothing yields an empty block that looks exactly like a
    # legitimately empty region.
    dm <- dim(x)
    i <- if (is.null(index[[1L]])) seq_len(dm[[1L]]) else index[[1L]]
    j <- if (is.null(index[[2L]])) seq_len(dm[[2L]]) else index[[2L]]
    if (length(i) == 0L || length(j) == 0L) {
        return(matrix(numeric(0), length(i), length(j)))
    }
    se <- extractBlockGenotypes(
        .ghSeedHandle(x),
        snpIdx = i,
        meanImpute = FALSE
    )
    as.matrix(SummarizedExperiment::assay(se, "dosage"))[, j, drop = FALSE]
})

#' Genotype dosages as a lazily-read DelayedArray
#'
#' Wraps a \code{\link{GenotypeHandle}} so its dosages can back a
#' \code{SummarizedExperiment} assay without being read. Nothing is read when
#' the array is built, nor when it is placed in a
#' \code{SummarizedExperiment} or \code{MultiAssayExperiment}: only the blocks
#' an operation actually touches are fetched.
#'
#' The result is \strong{variants x samples}, the Bioconductor assay
#' orientation, which is the transpose of what \code{\link{genotypes}}
#' returns.
#'
#' @param handle A \code{\link{GenotypeHandle}}.
#' @return A \code{DelayedMatrix} of dosages.
#' @keywords internal
genotypeDelayedArray <- function(handle) {
    if (!methods::is(handle, "GenotypeHandle")) {
        msg <- glue(
            "`handle` must be a GenotypeHandle (got {class(handle)[[1L]]})."
        )
        abort(msg)
    }
    si <- snpInfo(handle)
    seed <- methods::new(
        "GhSeed",
        handle = handle,
        dm = c(nrow(si), as.integer(nSamples(handle))),
        # Only the variant axis is stored; the sample axis is derived from the
        # handle in the dimnames() method above, so the ~10k sample ids are not
        # duplicated between dn[[2]] and handle@sampleIds.
        dn = list(
            normalizeVariantId(as.character(si$SNP)),
            NULL
        )
    )
    DelayedArray::DelayedArray(seed)
}
