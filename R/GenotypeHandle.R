# =============================================================================
# GenotypeHandle S4 class
# -----------------------------------------------------------------------------
# Lazy file handle for genotype data (PLINK1, PLINK2, VCF/BCF, GDS). Opens
# the file for metadata only (sample IDs and SNP info); dosage extraction
# is deferred until extractBlockGenotypes() is called.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Genotype File Handle
#' @description S4 container holding a path + format + metadata for lazy
#'   genotype access. Supports PLINK1 (.bed/.bim/.fam), PLINK2
#'   (.pgen/.pvar/.psam), VCF/BCF, and GDS.
#' @slot path Character, file path. For a one-file-per-chromosome handle this is
#'   the chrom-meta file path (a display/provenance value); the per-chromosome
#'   payload files live in \code{chromPaths}.
#' @slot format Character, one of \code{"plink1"}, \code{"plink2"},
#'   \code{"vcf"}, \code{"gds"}.
#' @slot snpInfo data.frame, SNP metadata read from the index/sidecar. For a
#'   sharded handle this is the union across chromosomes (row-bound in the order
#'   the shards were supplied), so \code{snpIdx} stays a single global index
#'   space.
#' @slot nSamples Integer, number of samples.
#' @slot sampleIds Character vector of sample identifiers.
#' @slot pgenPtr Opaque pointer for PLINK2 reader state (NULL otherwise).
#' @slot chromPaths Named character vector mapping canonical chromosome (e.g.
#'   \code{"21"}, \code{"X"}) to the per-chromosome payload path/prefix. Empty
#'   (\code{character(0)}) for a single-file handle; non-empty marks a
#'   one-file-per-chromosome (sharded) handle whose extraction is routed by
#'   chromosome.
#' @export
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
        errors <- character()
        if (length(object@path) != 1L) {
            errors <- c(errors, "'path' must be a single character string")
        }
        valid_formats <- c("gds", "vcf", "plink1", "plink2")
        if (!is_in(object@format, valid_formats)) {
            errors <- c(
                errors,
                str_c(
                    "'format' must be one of: ",
                    str_flatten(valid_formats, ", ")
                )
            )
        }
        if (length(object@chromPaths) > 0L) {
            nm <- names(object@chromPaths)
            if (
                is.null(nm) ||
                    any(str_length(nm) == 0L) ||
                    anyDuplicated(nm)
            ) {
                errors <- c(
                    errors,
                    str_c(
                        "'chromPaths' must be a uniquely-named ",
                        "character vector (names = chromosomes)"
                    )
                )
            }
        }
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
#'       \code{.vcf}, \code{.vcf.gz}, \code{.vcf.bgz}, \code{.bcf}, or
#'       \code{.gds}. Format is auto-detected from the extension.}
#'     \item{\code{plink1Prefix}}{A path prefix; the constructor appends
#'       \code{.bed}, \code{.bim}, and \code{.fam} to locate the triplet.}
#'     \item{\code{plink2Prefix}}{A path prefix; the constructor appends
#'       \code{.pgen}, \code{.pvar} (or \code{.pvar.zst}), and \code{.psam}
#'       to locate the triplet.}
#'     \item{\code{bed} + \code{bim} + \code{fam}}{Explicit PLINK1 triplet.
#'       The three files must share a stem; if they don't, use
#'       \code{plink1Prefix} or arrange symlinks at a common stem.}
#'     \item{\code{pgen} + \code{pvar} + \code{psam}}{Explicit PLINK2
#'       triplet. Same shared-stem requirement.}
#'     \item{\code{genoMeta}}{One genotype file per chromosome. Either a path
#'       to a whitespace/TSV meta file whose first column is the chromosome
#'       (\code{#chr}) and second column the per-chromosome payload
#'       (a \code{.bed}/\code{.pgen}/\code{.vcf[.gz]}/\code{.bcf}/\code{.gds}
#'       file or a PLINK prefix; relative paths resolve against the meta
#'       file's directory), or a named character vector mapping chromosome to
#'       payload. All shards must share one format and identical sample IDs in
#'       the same order. Extraction is routed to the correct file by the
#'       requested region's chromosome.}
#'   }
#'   The constructor opens the file for metadata only (sample IDs and SNP
#'   info); dosage extraction is deferred until \code{extractBlockGenotypes()}
#'   is called.
#'
#' @param path Single file path (.vcf/.vcf.gz/.vcf.bgz/.bcf/.gds), or
#'   \code{NULL}.
#' @param plink1Prefix Path prefix for a PLINK1 triplet, or \code{NULL}.
#' @param plink2Prefix Path prefix for a PLINK2 triplet, or \code{NULL}.
#' @param bed,bim,fam Explicit paths to the PLINK1 triplet, or all \code{NULL}.
#' @param pgen,pvar,psam Explicit paths to the PLINK2 triplet, or all
#'   \code{NULL}.
#' @param ldMeta Path to an LD-meta TSV file (columns \code{chrom},
#'   \code{start}, \code{end}, \code{path}; the \code{path} column may be
#'   comma-separated as \code{ld_file,bim_file}). Requires \code{region}. The
#'   constructor resolves the row covering \code{region}, then delegates to the
#'   appropriate file-based handler. When the resolved row points at PLINK1 /
#'   PLINK2 / VCF / GDS files, the corresponding reader is used; \code{.cor.xz}
#'   (pre-computed LD-matrix) rows are not supported here -- use
#'   \code{\link{loadLdMatrix}} for that case.
#' @param region Region specification for \code{ldMeta} lookup:
#'   \code{"chr:start-end"} string or a one-row data.frame with \code{chrom},
#'   \code{start}, \code{end}.
#' @param genoMeta One-file-per-chromosome specification: a path to a
#'   \code{#chr,path} meta file or a named character vector (names =
#'   chromosomes, values = payload paths/prefixes). Optionally pass
#'   \code{format} via \code{...} to force a single backend for every shard.
#' @param chroms Optional character vector of chromosomes. With \code{genoMeta}
#'   (a sharded, one-file-per-chromosome panel) only the shards whose chromosome
#'   is listed are read, skipping the rest -- an I/O optimisation when the panel
#'   is genome-wide but only a few chromosomes are needed. Chromosome labels are
#'   compared canonically (\code{"chr1"}/\code{"1"} match, \code{23}/\code{X}
#'   etc.). If none of the requested chromosomes are present in the panel, every
#'   shard is read (so the caller's own absence check can report the mismatch).
#'   Only meaningful with \code{genoMeta}; supplying it with any other source is
#'   an error (a single-file panel has no per-chromosome shards to skip).
#' @param ... Additional arguments forwarded to the format-specific reader.
#' @return A \code{GenotypeHandle} object.
#' @examples
#' ex <- function(f) system.file("extdata", f, package = "pecotmr")
#' GenotypeHandle(bed = ex("toy_ref.bed"), bim = ex("toy_ref.bim"),
#'   fam = ex("toy_ref.fam"))
#' @export
GenotypeHandle <- function(
    path = NULL,
    plink1Prefix = NULL,
    plink2Prefix = NULL,
    bed = NULL,
    bim = NULL,
    fam = NULL,
    pgen = NULL,
    pvar = NULL,
    psam = NULL,
    ldMeta = NULL,
    region = NULL,
    genoMeta = NULL,
    chroms = NULL,
    ...
) {
    p <- as.list(environment())
    flags <- .ghValidateArgs(p)
    sources <- .ghResolveSources(p, flags)
    if (sources[["path"]]) {
        return(readGenotypes(path, ...))
    }
    if (sources[["plink1Prefix"]]) {
        return(.makePlink1Handle(plink1Prefix, ...))
    }
    if (sources[["plink2Prefix"]]) {
        return(.makePlink2Handle(plink2Prefix, ...))
    }
    if (sources[["plink1Triplet"]]) {
        return(.genotypeHandleFromPlink1Triplet(bed, bim, fam, ...))
    }
    if (sources[["plink2Triplet"]]) {
        return(.genotypeHandleFromPlink2Triplet(pgen, pvar, psam, ...))
    }
    if (sources[["ldMeta"]]) {
        return(.genotypeHandleFromLdMeta(ldMeta, region, ...))
    }
    # nSources == 1 is enforced above, so the only remaining source is genoMeta.
    .genotypeHandleFromChromMeta(genoMeta, chroms = chroms, ...)
}

# Validate the bed/bim/fam + pgen/pvar/psam triplet completeness and the
# ldMeta<->region coupling. Returns list(bedComplete, pgenComplete).
# @noRd
.ghValidateArgs <- function(p) {
    bedComplete <- .ghTrioComplete(p$bed, p$bim, p$fam, "bed/bim/fam")
    pgenComplete <- .ghTrioComplete(p$pgen, p$pvar, p$psam, "pgen/pvar/psam")
    if (!is.null(p$ldMeta) && is.null(p$region)) {
        msg <- glue(
            "`ldMeta` requires a `region` (a 'chr:start-end' string or a ",
            "one-row data.frame with chrom/start/end)."
        )
        abort(msg)
    }
    if (is.null(p$ldMeta) && !is.null(p$region)) {
        abort("`region` is only meaningful when `ldMeta` is supplied.")
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
# source was supplied, and gate `chroms` to the genoMeta path.
# @noRd
.ghResolveSources <- function(p, flags) {
    sources <- c(
        path = !is.null(p$path),
        plink1Prefix = !is.null(p$plink1Prefix),
        plink2Prefix = !is.null(p$plink2Prefix),
        plink1Triplet = flags$bedComplete,
        plink2Triplet = flags$pgenComplete,
        ldMeta = !is.null(p$ldMeta),
        genoMeta = !is.null(p$genoMeta)
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
    if (!is.null(p$chroms) && !sources[["genoMeta"]]) {
        msg <- glue(
            "`chroms` restricts which per-chromosome shards are read and is ",
            "only supported with `genoMeta` (a single-file panel has no ",
            "shards to skip)."
        )
        abort(msg)
    }
    sources
}

.genotypeHandleFromLdMeta <- function(ldMeta, region, ...) {
    ldPath <- .ghResolveLdPath(ldMeta, region)
    .ghLdPathToHandle(ldPath, ...)
}

# Resolve the single genotype-payload path for `region` from an LD meta,
# erroring on no-coverage, multi-row spans, or a pre-computed .cor(.xz) matrix.
# @noRd
.ghResolveLdPath <- function(ldMeta, region) {
    ldPaths <- getRegionalLdMeta(ldMeta, region)$intersections$LD_file_paths
    if (length(ldPaths) == 0L) {
        regionStr <- deparse(region)
        msg <- glue(
            "GenotypeHandle: no LD-meta row covers region {regionStr} ",
            "in {ldMeta}."
        )
        abort(msg)
    }
    if (length(ldPaths) > 1L) {
        regionStr <- deparse(region)
        msg <- glue(
            "GenotypeHandle: region {regionStr} spans multiple LD-meta ",
            "rows; the GenotypeHandle constructor only resolves single-row ",
            "regions. Use loadLdMatrix() for multi-row regions, or restrict ",
            "the region to a single LD block."
        )
        abort(msg)
    }
    ldPath <- ldPaths[[1L]]
    if (str_detect(ldPath, regex("\\.cor(\\.xz)?$", ignore_case = TRUE))) {
        regionStr <- deparse(region)
        msg <- glue(
            "GenotypeHandle: the LD-meta row for region {regionStr} points ",
            "at a pre-computed correlation matrix ({ldPath}). Use ",
            "loadLdMatrix() / loadLdSketch() for .cor.xz inputs; ",
            "GenotypeHandle accepts only genotype payloads (VCF/GDS/PLINK)."
        )
        abort(msg)
    }
    ldPath
}

# Dispatch a resolved genotype path to the format-specific reader by extension.
# @noRd
.ghLdPathToHandle <- function(ldPath, ...) {
    lower <- str_to_lower(ldPath)
    if (str_detect(lower, "\\.vcf(\\.b?gz)?$") || str_ends(lower, "\\.bcf")) {
        return(readGenotypes(ldPath, format = "vcf", ...))
    }
    if (str_ends(lower, "\\.gds")) {
        return(readGenotypes(ldPath, format = "gds", ...))
    }
    if (str_ends(lower, "\\.bed")) {
        return(.makePlink1Handle(
            str_remove(ldPath, regex("\\.bed$", ignore_case = TRUE)),
            ...
        ))
    }
    if (str_ends(lower, "\\.pgen")) {
        return(.makePlink2Handle(
            str_remove(ldPath, regex("\\.pgen$", ignore_case = TRUE)),
            ...
        ))
    }
    msg <- glue(
        "GenotypeHandle: unsupported LD-meta file extension on ",
        "'{ldPath}'. Expected one of ",
        ".vcf/.vcf.gz/.vcf.bgz/.bcf/.gds/.bed/.pgen."
    )
    abort(msg)
}

.genotypeHandleFromPlink1Triplet <- function(bed, bim, fam, ...) {
    for (f in list(bed = bed, bim = bim, fam = fam)) {
        if (!is.character(f) || length(f) != 1L) {
            abort("Each of `bed`, `bim`, `fam` must be a single file path.")
        }
    }
    stems <- c(
        bed = file_path_sans_ext(bed),
        bim = file_path_sans_ext(bim),
        fam = file_path_sans_ext(fam)
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
    .makePlink1Handle(unname(stems[1L]), ...)
}

.genotypeHandleFromPlink2Triplet <- function(pgen, pvar, psam, ...) {
    for (f in list(pgen = pgen, pvar = pvar, psam = psam)) {
        if (!is.character(f) || length(f) != 1L) {
            abort("Each of `pgen`, `pvar`, `psam` must be a single file path.")
        }
    }
    pvarStem <- str_remove(pvar, regex("\\.zst$", ignore_case = TRUE))
    stems <- c(
        pgen = file_path_sans_ext(pgen),
        pvar = file_path_sans_ext(pvarStem),
        psam = file_path_sans_ext(psam)
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
    .makePlink2Handle(unname(stems[1L]), ...)
}

# ---------------------------------------------------------------------------
# One-file-per-chromosome (sharded) handle support
# ---------------------------------------------------------------------------

# Chromosome labels are canonicalized via canonChrom() (variantId.R) -- the
# single chromosome normalizer used for @chromPaths routing keys and for
# @snpInfo$CHR lookups, so both sides stay consistent (and X/Y/MT survive).

# Per-chromosome shard map, tolerant of GenotypeHandle objects deserialized
# from before the `chromPaths` slot existed (e.g. an RDS saved by an older
# pecotmr). Such objects have no `chromPaths` slot, so a direct `@` access
# errors; treat them as single-file handles.
#' @keywords internal
.genotypeChromPaths <- function(handle) {
    tryCatch(getChromPaths(handle), error = function(e) character(0))
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
            "GenotypeHandle(genoMeta): chromosome(s) {multiStr} map to ",
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
            chrom = names(genoMeta),
            path = unname(genoMeta)
        )))
    }
    msg <- glue(
        "GenotypeHandle(genoMeta): expected a path to a `#chr,path` meta ",
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
            "GenotypeHandle(genoMeta): meta file '{genoMeta}' must have ",
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
        chrom = as.character(meta[[chromCol]]),
        path = .metaResolvePaths(as.character(meta[[pathCol]]), base)
    ))
}

# Resolve relative payload paths against `base` (absolute/existing kept as-is).
# @noRd
.metaResolvePaths <- function(pth, base) {
    map_chr(pth, .ghResolvePath, base = base)
}

# Build a single-file GenotypeHandle for one shard payload, dispatching to the
# right reader. `format` (optional) forces a backend; otherwise it is detected
# from the file extension, falling back to PLINK prefix probing.
#' @keywords internal
.resolveGenotypeShard <- function(p, format = NULL) {
    lower <- str_to_lower(p)
    if (!is.null(format)) {
        if (format == "plink1") {
            return(.makePlink1Handle(p))
        }
        if (format == "plink2") {
            return(.makePlink2Handle(p))
        }
        return(readGenotypes(p, format = format))
    }
    if (str_detect(lower, "\\.vcf(\\.b?gz)?$") || str_ends(lower, "\\.bcf")) {
        return(readGenotypes(p, format = "vcf"))
    }
    if (str_ends(lower, "\\.gds")) {
        return(readGenotypes(p, format = "gds"))
    }
    if (str_ends(lower, "\\.bed")) {
        return(.makePlink1Handle(
            str_remove(p, regex("\\.bed$", ignore_case = TRUE))
        ))
    }
    if (str_ends(lower, "\\.pgen")) {
        return(.makePlink2Handle(
            str_remove(p, regex("\\.pgen$", ignore_case = TRUE))
        ))
    }
    # No recognized extension: treat as a PLINK prefix, probe for the sidecar.
    if (file.exists(str_c(p, ".bed"))) {
        return(.makePlink1Handle(p))
    }
    if (file.exists(str_c(p, ".pgen"))) {
        return(.makePlink2Handle(p))
    }
    msg <- glue(
        "GenotypeHandle(genoMeta): cannot determine genotype format for ",
        "'{p}'. Use a recognized extension ",
        "(.bed/.pgen/.vcf[.gz]/.bcf/.gds), a PLINK prefix, or pass ",
        "`format=`."
    )
    abort(msg)
}

# Assemble a sharded handle from a per-chromosome meta. Reads each shard's
# metadata via the existing single-file readers, validates a single shared
# format and identical sample IDs (same order, required for cross-shard
# cbind), and row-binds the per-shard snpInfo into one global index space.
# `chroms` (optional) restricts the read to the shards for those chromosomes:
# the other per-chromosome files are never opened, which is the I/O win when a
# genome-wide panel backs summary statistics on only a few chromosomes.
#' @keywords internal
.genotypeHandleFromChromMeta <- function(genoMeta, chroms = NULL, ...) {
    format <- list(...)$format
    parsed <- .parseChromMeta(genoMeta)
    if (nrow(parsed) == 0L) {
        abort(
            "GenotypeHandle(genoMeta): no chromosomes found in the meta input."
        )
    }
    parsed <- .chromMetaSelect(parsed, chroms)
    shards <- map(parsed$path, .resolveGenotypeShard, format = format)
    sharedFormat <- .chromMetaCheckFormats(shards)
    .chromMetaCheckSamples(shards, parsed)
    unifiedSnpInfo <- bind_rows(map(shards, .ghSnpInfo))
    new(
        "GenotypeHandle",
        path = .chromMetaPath(genoMeta),
        format = sharedFormat,
        snpInfo = unifiedSnpInfo,
        nSamples = getNSamples(shards[[1L]]),
        sampleIds = getSampleIds(shards[[1L]]),
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
    keep <- is_in(
        canonChrom(as.character(parsed$chrom)),
        canonChrom(as.character(chroms))
    )
    if (any(keep)) filter(parsed, keep) else parsed
}

# Require a single shared format across shards; returns it.
# @noRd
.chromMetaCheckFormats <- function(shards) {
    formats <- map_chr(shards, .ghFormat)
    if (n_distinct(formats) != 1L) {
        formatStr <- str_flatten(unique(formats), ", ")
        msg <- glue(
            "GenotypeHandle(genoMeta): all per-chromosome files must share ",
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
                "GenotypeHandle(genoMeta): all per-chromosome files must ",
                "have identical sample IDs in the same order (mismatch ",
                "at '{mismatchPath}')."
            )
            abort(msg)
        }
    }
    invisible(NULL)
}

# Map each chromosome to its shard path, erroring if a chromosome spans files.
# (Sequential uniqueness accumulation -- kept as a loop.)
# @noRd
.chromMetaPaths <- function(shards) {
    chromPaths <- character(0)
    for (i in seq_along(shards)) {
        for (ch in unique(canonChrom(getSnpInfo(shards[[i]])$CHR))) {
            if (is_in(ch, names(chromPaths))) {
                msg <- glue(
                    "GenotypeHandle(genoMeta): chromosome '{ch}' appears ",
                    "in more than one per-chromosome file."
                )
                abort(msg)
            }
            chromPaths[[ch]] <- getPath(shards[[i]])
        }
    }
    chromPaths
}

# The handle's recorded path: the normalized meta-file path, or a placeholder
# for a named-vector meta.
# @noRd
.chromMetaPath <- function(genoMeta) {
    if (.isChromMetaFile(genoMeta)) normalizePath(genoMeta) else "<chrom-meta>"
}

#' @rdname getSnpInfo
#' @export
setMethod("getSnpInfo", "GenotypeHandle", function(x) x@snpInfo)

#' @rdname getFormat
#' @export
setMethod("getFormat", "GenotypeHandle", function(x) x@format)

#' @rdname getPath
#' @export
setMethod("getPath", "GenotypeHandle", function(x) x@path)

#' @rdname getChromPaths
#' @export
setMethod("getChromPaths", "GenotypeHandle", function(x, ...) x@chromPaths)

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
    bed <- system.file("extdata", str_c(stem, ".bed"), package = "pecotmr")
    if (str_length(bed) > 0L) {
        return(str_remove(bed, "\\.bed$"))
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
# extractors and LD-panel keying, so `getPath()` keeps returning the stored
# (portable) value for display/provenance while file access resolves it.
.genotypeReadPath <- function(handle) {
    .resolveGenotypeResourcePath(getPath(handle))
}

#' @rdname getSampleIds
#' @export
setMethod("getSampleIds", "GenotypeHandle", function(x) x@sampleIds)

#' @rdname getPgenPtr
#' @export
setMethod("getPgenPtr", "GenotypeHandle", function(x) x@pgenPtr)

#' @rdname getNSamples
#' @export
setMethod("getNSamples", "GenotypeHandle", function(x) x@nSamples)

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
    getSnpInfo(h)
}

# The format slot of one genotype shard.
# @noRd
.ghFormat <- function(h) {
    getFormat(h)
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
# pin it against getGenotypes().
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
setMethod("dimnames", "GhSeed", function(x) x@dn)

#' @rdname GhSeed-class
#' @importFrom S4Arrays extract_array
setMethod("extract_array", "GhSeed", function(x, index) {
    # A NULL index means "everything along this dimension", not "nothing".
    # Reading it as nothing yields an empty block that looks exactly like a
    # legitimately empty region.
    i <- if (is.null(index[[1L]])) seq_len(x@dm[[1L]]) else index[[1L]]
    j <- if (is.null(index[[2L]])) seq_len(x@dm[[2L]]) else index[[2L]]
    if (length(i) == 0L || length(j) == 0L) {
        return(matrix(numeric(0), length(i), length(j)))
    }
    se <- extractBlockGenotypes(x@handle, snpIdx = i, meanImpute = FALSE)
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
#' orientation, which is the transpose of what \code{\link{getGenotypes}}
#' returns.
#'
#' @param handle A \code{\link{GenotypeHandle}}.
#' @return A \code{DelayedMatrix} of dosages.
#' @examples
#' data(qtlDatasetExample)
#' da <- genotypeDelayedArray(getGenotypeHandle(qtlDatasetExample))
#' dim(da)
#' @export
genotypeDelayedArray <- function(handle) {
    if (!methods::is(handle, "GenotypeHandle")) {
        msg <- glue(
            "`handle` must be a GenotypeHandle (got {class(handle)[[1L]]})."
        )
        abort(msg)
    }
    si <- getSnpInfo(handle)
    seed <- methods::new(
        "GhSeed",
        handle = handle,
        dm = c(nrow(si), as.integer(getNSamples(handle))),
        dn = list(
            normalizeVariantId(as.character(si$SNP)),
            as.character(getSampleIds(handle))
        )
    )
    DelayedArray::DelayedArray(seed)
}
