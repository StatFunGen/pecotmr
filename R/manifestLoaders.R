# =============================================================================
# Manifest loaders
# -----------------------------------------------------------------------------
# Reader functions that build the individual-level / summary-statistic S4
# containers from an on-disk (or in-memory) manifest:
#
#   loadQtlDatasetFromManifest         -> QtlDataset
#   loadQtlSumStatsFromManifest        -> QtlSumStats
#   loadGwasSumStatsFromManifest       -> GwasSumStats
#   loadMultiStudyQtlDatasetFromManifest -> MultiStudyQtlDataset
#
# These are pure readers: they never run summaryStatsQc() (the loaded
# sumstats objects carry qcInfo = list()) and never filter genotypes (the
# QtlDataset QC knobs are stored as lazy slots, applied at extraction time).
#
# A manifest is either a data.frame or a path to a delimited file. A ".csv"
# extension is read with readr::read_csv; anything else is read with
# readr::read_tsv. Both snake_case and camelCase column names are accepted at
# the parsing boundary and canonicalised to camelCase; snake_case is never
# propagated internally.
# =============================================================================

#' @include GenotypeHandle.R QtlDataset.R gwasSumStats.R qtlSumStats.R
#' @include MultiStudyQtlDataset.R variantId.R sumstatsQc.R
NULL

# =============================================================================
# Manifest ingestion + column canonicalisation
# =============================================================================

# Read a manifest into a base data.frame. A data.frame is passed through; a
# single path is read by extension (.csv -> read_csv, else read_tsv).
.readManifest <- function(manifest) {
    if (is.data.frame(manifest)) {
        return(as.data.frame(
            manifest,
            stringsAsFactors = FALSE,
            check.names = FALSE
        ))
    }
    if (!is.character(manifest) || length(manifest) != 1L) {
        stop("`manifest` must be a data.frame or a single file path.")
    }
    if (!file.exists(manifest)) {
        stop("manifest file not found: ", manifest)
    }
    parsed <- if (grepl("\\.csv$", manifest, ignore.case = TRUE)) {
        readr::read_csv(manifest, show_col_types = FALSE, progress = FALSE)
    } else {
        readr::read_tsv(manifest, show_col_types = FALSE, progress = FALSE)
    }
    as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)
}

# Directory a manifest's relative paths resolve against: the file's own
# directory when the manifest was a path, otherwise NULL (paths used as-is).
.manifestBase <- function(manifest) {
    if (
        is.character(manifest) &&
            length(manifest) == 1L &&
            file.exists(manifest)
    ) {
        dirname(normalizePath(manifest))
    } else {
        NULL
    }
}

# Resolve a possibly-relative path against a manifest base directory.
.resolveRel <- function(p, base) {
    if (is.null(p) || length(p) != 1L || is.na(p) || !nzchar(p)) {
        return(p)
    }
    if (grepl("^(/|[A-Za-z]:|~)", p) || is.null(base) || file.exists(p)) {
        return(p)
    }
    file.path(base, p)
}

# Rename accepted aliases to the canonical camelCase name and validate that all
# `required` canonical columns are present. `aliases` is a named list keyed by
# canonical name; each value is the character vector of accepted source names
# (including the canonical name itself). The first alias present wins.
.canonManifestCols <- function(df, aliases, required, label) {
    for (canon in names(aliases)) {
        if (canon %in% names(df)) {
            next
        }
        hit <- intersect(aliases[[canon]], names(df))
        if (length(hit) >= 1L) {
            names(df)[match(hit[[1L]], names(df))] <- canon
        }
    }
    missingCols <- setdiff(required, names(df))
    if (length(missingCols) > 0L) {
        stop(
            label,
            " manifest is missing required column(s): ",
            paste(missingCols, collapse = ", ")
        )
    }
    df
}

# Resolve a single collection-level scalar (e.g. genome) that may be supplied
# both as a constant manifest column and as a function argument. The column (if
# present) must be constant; a supplied argument must agree with it.
.reconcileScalar <- function(colValues, arg, name, required = TRUE) {
    colVal <- NULL
    if (!is.null(colValues)) {
        v <- unique(as.character(colValues[
            !is.na(colValues) &
                nzchar(as.character(colValues))
        ]))
        if (length(v) > 1L) {
            stop(
                "`",
                name,
                "` column must be constant across the manifest; got: ",
                paste(v, collapse = ", ")
            )
        }
        if (length(v) == 1L) colVal <- v
    }
    if (!is.null(arg) && !is.null(colVal)) {
        if (!identical(as.character(arg), colVal)) {
            stop(
                "`",
                name,
                "` argument ('",
                arg,
                "') disagrees with the manifest ",
                "column value ('",
                colVal,
                "')."
            )
        }
        return(arg)
    }
    if (!is.null(arg)) {
        return(arg)
    }
    if (!is.null(colVal)) {
        return(colVal)
    }
    if (required) {
        stop(
            "`",
            name,
            "` must be provided as an argument or a manifest column."
        )
    }
    NULL
}

# =============================================================================
# Genotype-format detection + LD-sketch resolution
# =============================================================================

# Build a single-file / prefix GenotypeHandle, detecting the backend from the
# path. Recognised: .vcf[.gz|.bgz]/.bcf/.gds single files, .bed/.pgen single
# files, and bare PLINK1/PLINK2 prefixes (probed by sidecar existence).
.detectGenotypeFormat <- function(path) {
    lower <- tolower(path)
    if (
        grepl("\\.vcf(\\.b?gz)?$", lower) ||
            endsWith(lower, ".bcf") ||
            endsWith(lower, ".gds")
    ) {
        return(GenotypeHandle(path = path))
    }
    if (endsWith(lower, ".bed")) {
        return(GenotypeHandle(
            plink1Prefix = sub("\\.bed$", "", path, ignore.case = TRUE)
        ))
    }
    if (endsWith(lower, ".pgen")) {
        return(GenotypeHandle(
            plink2Prefix = sub("\\.pgen$", "", path, ignore.case = TRUE)
        ))
    }
    if (file.exists(paste0(path, ".bed"))) {
        return(GenotypeHandle(plink1Prefix = path))
    }
    if (file.exists(paste0(path, ".pgen"))) {
        return(GenotypeHandle(plink2Prefix = path))
    }
    stop(
        "Could not determine genotype format for '",
        path,
        "' (expected .bed/.pgen/.vcf[.gz]/.bcf/.gds or a PLINK prefix)."
    )
}

# Turn an ldSketch specification into a GenotypeHandle. Accepts a prebuilt
# handle, a named chrom -> path vector (genoMeta), or a single string that is
# either a genotype file/prefix or a chrom-sharded genoMeta meta file.
# `chroms` (optional) restricts a chrom-sharded panel to those chromosomes so
# the other per-chromosome shards are never read; it is ignored for a single
# genotype file (which has no shards to skip).
.resolveLdSketch <- function(spec, chroms = NULL) {
    if (is.null(spec)) {
        return(NULL)
    }
    if (methods::is(spec, "GenotypeHandle")) {
        return(spec)
    }
    if (is.character(spec) && !is.null(names(spec))) {
        return(GenotypeHandle(genoMeta = spec, chroms = chroms))
    }
    if (is.character(spec) && length(spec) == 1L) {
        lower <- tolower(spec)
        looksLikeGenotype <- grepl("\\.vcf(\\.b?gz)?$", lower) ||
            endsWith(lower, ".bcf") ||
            endsWith(lower, ".gds") ||
            endsWith(lower, ".bed") ||
            endsWith(lower, ".pgen") ||
            file.exists(paste0(spec, ".bed")) ||
            file.exists(paste0(spec, ".pgen"))
        if (looksLikeGenotype) {
            return(.detectGenotypeFormat(spec))
        }
        # Otherwise treat it as a chrom-sharded genoMeta meta file.
        return(GenotypeHandle(genoMeta = spec, chroms = chroms))
    }
    stop(
        "`ldSketch` must be a GenotypeHandle, a genotype path/prefix, or a ",
        "genoMeta spec (named chrom->path vector or meta-file path)."
    )
}

# Resolve the LD sketch SPEC from the (argument, ldSketchPath column) pair,
# without reading any genotype metadata yet. Returns a prebuilt GenotypeHandle
# (when the caller passed one), or the reconciled spec (a path/prefix or a
# genoMeta path/named-vector) for later materialisation. Deferring the read to
# `.materializeLdSketch` lets the loader first learn which chromosomes the
# summary statistics cover and skip the panel's other per-chromosome shards.
.resolveLdSketchInput <- function(df, ldSketch, base) {
    if (methods::is(ldSketch, "GenotypeHandle")) {
        return(ldSketch)
    }
    colSpec <- NULL
    if (!is.null(df$ldSketchPath)) {
        v <- unique(as.character(df$ldSketchPath[
            !is.na(df$ldSketchPath) &
                nzchar(as.character(df$ldSketchPath))
        ]))
        if (length(v) > 1L) {
            stop("`ldSketchPath` column must be constant across the manifest.")
        }
        if (length(v) == 1L) colSpec <- .resolveRel(v, base)
    }
    spec <- ldSketch
    if (
        !is.null(spec) &&
            !is.null(colSpec) &&
            !identical(as.character(spec), as.character(colSpec))
    ) {
        stop(
            "`ldSketch` argument disagrees with the manifest ",
            "`ldSketchPath` column."
        )
    }
    spec <- spec %||% colSpec
    if (is.null(spec)) {
        stop(
            "`ldSketch` must be provided as an argument or an ",
            "`ldSketchPath` column."
        )
    }
    spec
}

# Materialise a resolved ldSketch spec into a GenotypeHandle, reading only the
# shards for `chroms` when the spec is a chrom-sharded panel. A spec that is
# already a GenotypeHandle (the caller passed a prebuilt handle) is returned
# unchanged: its snpInfo is already in memory, so there is no read to restrict.
.materializeLdSketch <- function(spec, chroms) {
    if (is.null(spec)) {
        return(NULL)
    }
    if (methods::is(spec, "GenotypeHandle")) {
        return(spec)
    }
    .resolveLdSketch(spec, chroms = chroms)
}

# Canonical chromosomes present across a list of sumstats entry GRanges. NULL /
# empty entries contribute nothing; NA seqnames are dropped. Always returns a
# character vector (character(0) when nothing is present, never NULL).
.entriesChroms <- function(entries) {
    ch <- as.character(unlist(
        lapply(entries, function(gr) {
            if (is.null(gr) || length(gr) == 0L) {
                return(character(0))
            }
            canonChrom(as.character(GenomicRanges::seqnames(gr)))
        }),
        use.names = FALSE
    ))
    unique(ch[!is.na(ch)])
}

# =============================================================================
# LD-sketch containment check (validation only; no filtering)
# =============================================================================

# Canonical chromosomes present in an LD sketch: for a sharded handle these are
# the chromPaths names; for a single-file handle the unique snpInfo chromosomes.
.ldSketchChroms <- function(ldSketch) {
    if (length(ldSketch@chromPaths) > 0L) {
        canonChrom(names(ldSketch@chromPaths))
    } else {
        unique(canonChrom(as.character(ldSketch@snpInfo$CHR)))
    }
}

# Two-tier containment check for one entry GRanges against the LD sketch:
#   * chromosome tier (hard error): any entry chromosome absent from the sketch.
#   * variant tier: match entry variants against the sketch's variants on the
#     shared chromosomes. Zero overlap -> error; overlap fraction below
#     `minLdOverlapWarn` -> warning. Never filters.
.checkLdContainment <- function(ldSketch, entryGr, minLdOverlapWarn, label) {
    if (length(entryGr) == 0L) {
        return(invisible(NULL))
    }
    reqChrom <- .ldContainmentCheckChroms(entryGr, ldSketch, label)
    snpInfo <- ldSketch@snpInfo
    if (is.null(snpInfo) || nrow(snpInfo) == 0L) {
        warning(
            label,
            ": LD sketch carries no variant metadata; ",
            "skipping the variant-overlap check."
        )
        return(invisible(NULL))
    }
    .ldContainmentCheckOverlap(
        entryGr,
        snpInfo,
        reqChrom,
        minLdOverlapWarn,
        label
    )
}

# The entry's chromosomes must all be present in the LD sketch; returns them.
# @noRd
.ldContainmentCheckChroms <- function(entryGr, ldSketch, label) {
    reqChrom <- unique(canonChrom(
        as.character(GenomicRanges::seqnames(entryGr))
    ))
    availChrom <- .ldSketchChroms(ldSketch)
    missingChrom <- setdiff(reqChrom, availChrom)
    if (length(missingChrom) > 0L) {
        stop(
            label,
            ": chromosome(s) ",
            paste(missingChrom, collapse = ", "),
            " are absent from the LD sketch (available: ",
            paste(availChrom, collapse = ", "),
            ")."
        )
    }
    reqChrom
}

# Variant-overlap check: error on zero overlap, warn below minLdOverlapWarn.
# @noRd
.ldContainmentCheckOverlap <- function(
    entryGr,
    snpInfo,
    reqChrom,
    minLdOverlapWarn,
    label
) {
    keep <- canonChrom(as.character(snpInfo$CHR)) %in% reqChrom
    snpInfo <- snpInfo[keep, , drop = FALSE]
    entryVid <- formatVariantId(
        chrom = canonChrom(as.character(GenomicRanges::seqnames(entryGr))),
        pos = GenomicRanges::start(entryGr),
        A2 = as.character(entryGr$A2),
        A1 = as.character(entryGr$A1)
    )
    sketchVid <- formatVariantId(
        chrom = canonChrom(as.character(snpInfo$CHR)),
        pos = as.integer(snpInfo$BP),
        A2 = as.character(snpInfo$A2),
        A1 = as.character(snpInfo$A1)
    )
    nOverlap <- length(matchVariants(entryVid, sketchVid)$idxA)
    if (nOverlap == 0L) {
        stop(
            label,
            ": none of the ",
            length(entryGr),
            " summary-statistic ",
            "variants are present in the LD sketch."
        )
    }
    frac <- nOverlap / length(entryGr)
    if (frac < minLdOverlapWarn) {
        warning(sprintf(
            paste0(
                "%s: only %.1f%% of summary-statistic variants (%d/%d) are ",
                "present in the LD sketch."
            ),
            label,
            100 * frac,
            nOverlap,
            length(entryGr)
        ))
    }
    invisible(NULL)
}

# =============================================================================
# Summary-statistic file readers
# =============================================================================

# Coerce a region argument (chr:start-end string, GRanges, or one-row
# data.frame with chrom/start/end) into a GRanges.
.asGRegion <- function(region) {
    if (methods::is(region, "GRanges")) {
        return(region)
    }
    if (is.data.frame(region)) {
        return(GenomicRanges::GRanges(
            seqnames = as.character(region$chrom),
            ranges = IRanges::IRanges(
                start = as.integer(region$start),
                end = as.integer(region$end)
            )
        ))
    }
    if (is.character(region) && length(region) == 1L) {
        m <- regmatches(
            region,
            regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", region)
        )[[1L]]
        if (length(m) != 4L) {
            stop(
                "`region` must be a 'chr:start-end' string, a GRanges, or a ",
                "one-row data.frame with chrom/start/end."
            )
        }
        return(GenomicRanges::GRanges(
            seqnames = m[[2L]],
            ranges = IRanges::IRanges(
                start = as.integer(gsub(",", "", m[[3L]])),
                end = as.integer(gsub(",", "", m[[4L]]))
            )
        ))
    }
    stop(
        "`region` must be a 'chr:start-end' string, a GRanges, or a ",
        "one-row data.frame with chrom/start/end."
    )
}

# Remap a region's seqnames to the seqnames actually present in a tabix file
# (tolerant of chr-prefix differences). Returns NULL when no seqname matches.
.matchRegionSeqnames <- function(gr, available) {
    want <- as.character(GenomicRanges::seqnames(gr))
    mapped <- vapply(
        want,
        function(w) {
            hit <- available[canonChrom(available) == canonChrom(w)]
            if (length(hit) >= 1L) hit[[1L]] else NA_character_
        },
        character(1)
    )
    if (anyNA(mapped)) {
        return(NULL)
    }
    GenomicRanges::GRanges(
        seqnames = mapped,
        ranges = IRanges::IRanges(
            start = GenomicRanges::start(gr),
            end = GenomicRanges::end(gr)
        )
    )
}

# Read the region of a bgzipped + tabix-indexed delimited file into a
# data.frame, using the file's last header (#...) line for column names.
.readTabixRegion <- function(path, region) {
    tf <- Rsamtools::TabixFile(path)
    hdr <- Rsamtools::headerTabix(tf)$header
    colLine <- if (length(hdr) > 0L) {
        sub("^#", "", hdr[[length(hdr)]])
    } else {
        NULL
    }
    gr <- .matchRegionSeqnames(.asGRegion(region), Rsamtools::seqnamesTabix(tf))
    emptyDf <- if (!is.null(colLine)) {
        as.data.frame(
            readr::read_tsv(
                I(colLine),
                show_col_types = FALSE,
                progress = FALSE
            ),
            check.names = FALSE
        )[0, , drop = FALSE]
    } else {
        data.frame()
    }
    if (is.null(gr)) {
        return(emptyDf)
    }
    lines <- unlist(Rsamtools::scanTabix(tf, param = gr), use.names = FALSE)
    if (length(lines) == 0L) {
        return(emptyDf)
    }
    txt <- paste(c(colLine, lines), collapse = "\n")
    # Read every column as character; .resolveSumstatCols() coerces each field
    # to its target type. This prevents readr from type-guessing an allele/ID
    # column
    # that is uniformly T/F within the region as logical (the allele "T" ->
    # TRUE),
    # which would corrupt the variant id and break LD-containment matching.
    as.data.frame(
        readr::read_tsv(
            I(txt),
            show_col_types = FALSE,
            progress = FALSE,
            col_types = readr::cols(.default = readr::col_character())
        ),
        check.names = FALSE
    )
}

# Default source-column aliases used when no explicit columnMapping is given
# (mirrors gwasSumStats_construct.R).
.sumstatColumnAliases <- list(
    chrom = c("chrom", "#chrom", "chr", "#chr"),
    pos = c("pos", "position", "BP", "bp"),
    variant_id = c("variant_id", "SNP", "snp", "rsid"),
    A1 = c("A1", "a1"),
    A2 = c("A2", "a2"),
    Z = c("z", "Z"),
    N = c("n_sample", "N", "n"),
    N_CASE = c("N_CASE", "n_case", "ncase"),
    N_CONTROL = c("N_CONTROL", "n_control", "ncontrol"),
    BETA = c("beta", "BETA"),
    SE = c("se", "SE"),
    P = c("p", "P", "pvalue", "pval"),
    MAF = c("maf", "MAF", "effect_allele_frequency"),
    INFO = c("info", "INFO")
)

# Accepted columnMapping *key* spellings per canonical field. A mapping file
# uses the xqtl-protocol standard-key names (e.g. `z`, `n_sample`), which
# differ from the canonical mcol names (`Z`, `N`); this bridges the two.
.sumstatMappingKeys <- list(
    chrom = c("chrom", "chr"),
    pos = c("pos", "position"),
    variant_id = c("variant_id", "snp", "SNP"),
    A1 = c("A1", "a1"),
    A2 = c("A2", "a2"),
    Z = c("z", "Z"),
    N = c("n_sample", "N", "n"),
    N_CASE = c("n_case", "N_CASE", "ncase"),
    N_CONTROL = c("n_control", "N_CONTROL", "ncontrol"),
    BETA = c("beta", "BETA"),
    SE = c("se", "SE"),
    P = c("p", "P", "pvalue"),
    MAF = c("maf", "MAF", "effect_allele_frequency"),
    INFO = c("info", "INFO")
)

# Read a columnMapping specification into a named character vector mapping a
# standard key (chrom/pos/variant_id/...) to the source column name. Accepts a
# named list/vector, or a path to a YAML file of `standardName: sourceName`
# entries (the xqtl-protocol column-mapping format).
.readColumnMapping <- function(columnMapping) {
    if (is.null(columnMapping)) {
        return(NULL)
    }
    if (
        is.list(columnMapping) ||
            (is.character(columnMapping) &&
                !is.null(names(columnMapping)))
    ) {
        return(vapply(columnMapping, as.character, character(1)))
    }
    if (is.character(columnMapping) && length(columnMapping) == 1L) {
        if (!file.exists(columnMapping)) {
            stop("columnMapping file not found: ", columnMapping)
        }
        mapping <- yaml::read_yaml(columnMapping)
        if (
            !is.list(mapping) ||
                is.null(names(mapping)) ||
                any(!nzchar(names(mapping)))
        ) {
            stop(
                "columnMapping file '",
                columnMapping,
                "' must be a YAML mapping of standardName: sourceName."
            )
        }
        vapply(mapping, as.character, character(1))
    } else {
        stop(
            "`columnMapping` must be a named list/vector or a path to a YAML ",
            "mapping file."
        )
    }
}

# Resolve one canonical sumstats field to its source column: honor an explicit
# columnMapping (under any accepted key spelling), else fall back to known
# column
# aliases; NA when unresolved.
# @noRd
.resolveSumstatKey <- function(key, df, mapping, label) {
    # Look the field up in the mapping under any accepted standard-key spelling
    # (`z`/`Z`, `n_sample`/`N`, ...). `intersect` avoids the "subscript out of
    # bounds" a named character vector throws for an absent `[[` key.
    if (!is.null(mapping)) {
        cand <- intersect(.sumstatMappingKeys[[key]], names(mapping))
        if (length(cand) >= 1L) {
            src <- mapping[[cand[[1L]]]]
            if (!(src %in% names(df))) {
                stop(
                    label,
                    ": columnMapping['",
                    cand[[1L]],
                    "'] = '",
                    src,
                    "' is not a column in the sumstats file."
                )
            }
            return(src)
        }
    }
    intersect(.sumstatColumnAliases[[key]], names(df))[1L]
}

# Standardise the columns of a raw sumstats data.frame into the canonical
# schema expected by .dfToEntryGranges: chrom, pos, SNP, A1, A2, Z, N (+
# optional N_CASE, N_CONTROL, BETA, SE, P, MAF, INFO). N is required UNLESS
# per-variant N_CASE + N_CONTROL are both present (from which summaryStatsQc
# derives the effective sample size), OR `allowNoN = TRUE` because the caller
# has a study-level scalar (nCase/nControl or nSample) that summaryStatsQc will
# fill N from. Z is required UNLESS BETA + SE are both present, in which case
# the Wald z (Z = BETA / SE) is derived; a supplied Z always takes precedence.
# `mapping` (named std -> source) overrides the default aliases per key.
.resolveSumstatCols <- function(df, columnMapping, label, allowNoN = FALSE) {
    # Strip a leading '#' from the first column name so a '#CHR'-style header
    # (common in GWAS TSVs) resolves the same on plain-text and tabix reads.
    if (ncol(df) > 0L) {
        names(df)[1L] <- sub("^#", "", names(df)[1L])
    }
    mapping <- .readColumnMapping(columnMapping)
    resolved <- .resolveSumstatRequired(df, mapping, label)
    z <- .resolveSumstatZ(df, mapping, label)
    n <- .resolveSumstatN(df, mapping, label, allowNoN)
    .buildSumstatOut(df, resolved, z, n, mapping, label)
}

# Resolve + validate the always-required columns (chrom/pos/variant_id/A1/A2).
# @noRd
.resolveSumstatRequired <- function(df, mapping, label) {
    required <- c("chrom", "pos", "variant_id", "A1", "A2")
    resolved <- set_names(
        map(required, .resolveSumstatKey, df, mapping, label),
        required
    )
    missingKeys <- required[map_lgl(
        resolved,
        function(x) is.na(x) || is.null(x)
    )]
    if (length(missingKeys) > 0L) {
        stop(
            label,
            ": sumstats file is missing required field(s): ",
            paste(missingKeys, collapse = ", "),
            " (supply a columnMapping if the source names differ)."
        )
    }
    resolved
}

# Resolve the z source: a Z column, OR BETA + SE (Wald z = beta/se). At least
# one required; Z takes precedence. Returns list(hasZ, zSrc, betaSrc, seSrc).
# @noRd
.resolveSumstatZ <- function(df, mapping, label) {
    zSrc <- .resolveSumstatKey("Z", df, mapping, label)
    betaSrc <- .resolveSumstatKey("BETA", df, mapping, label)
    seSrc <- .resolveSumstatKey("SE", df, mapping, label)
    hasZ <- !is.na(zSrc) && !is.null(zSrc)
    hasBetaSe <- !is.na(betaSrc) &&
        !is.null(betaSrc) &&
        !is.na(seSrc) &&
        !is.null(seSrc)
    if (!hasZ && !hasBetaSe) {
        stop(
            label,
            ": sumstats file needs a z field (z/Z), or beta + se to derive ",
            "the Wald z-score beta/se (supply a columnMapping if the ",
            "source names ",
            "differ)."
        )
    }
    list(hasZ = hasZ, zSrc = zSrc, betaSrc = betaSrc, seSrc = seSrc)
}

# Resolve the sample-size source: an N column, OR N_CASE + N_CONTROL counts.
# Returns list(hasN, hasCounts, nSrc, ncaseSrc, ncontrolSrc).
# @noRd
.resolveSumstatN <- function(df, mapping, label, allowNoN) {
    nSrc <- .resolveSumstatKey("N", df, mapping, label)
    ncaseSrc <- .resolveSumstatKey("N_CASE", df, mapping, label)
    ncontrolSrc <- .resolveSumstatKey("N_CONTROL", df, mapping, label)
    hasN <- !is.na(nSrc) && !is.null(nSrc)
    hasCounts <- !is.na(ncaseSrc) &&
        !is.null(ncaseSrc) &&
        !is.na(ncontrolSrc) &&
        !is.null(ncontrolSrc)
    if (!hasN && !hasCounts && !allowNoN) {
        stop(
            label,
            ": sumstats file needs an N field (n_sample/N/n), or ",
            "N_CASE + N_CONTROL columns for the effective sample size, or a ",
            "study-level scalar (nCase/nControl or nSample) in the manifest ",
            "(supply a columnMapping if the source names differ)."
        )
    }
    list(
        hasN = hasN,
        hasCounts = hasCounts,
        nSrc = nSrc,
        ncaseSrc = ncaseSrc,
        ncontrolSrc = ncontrolSrc
    )
}

# Assemble the standardized sumstat data.frame (chrom/pos/SNP/A1/A2/Z, optional
# N / N_CASE / N_CONTROL, and any present BETA/SE/P/MAF/INFO).
# @noRd
.buildSumstatOut <- function(df, resolved, z, n, mapping, label) {
    out <- data.frame(
        chrom = as.character(df[[resolved$chrom]]),
        pos = as.integer(df[[resolved$pos]]),
        SNP = as.character(df[[resolved$variant_id]]),
        A1 = as.character(df[[resolved$A1]]),
        A2 = as.character(df[[resolved$A2]]),
        stringsAsFactors = FALSE
    )
    out$Z <- if (z$hasZ) {
        as.numeric(df[[z$zSrc]])
    } else {
        as.numeric(df[[z$betaSrc]]) / as.numeric(df[[z$seSrc]])
    }
    if (n$hasN) {
        out$N <- as.numeric(df[[n$nSrc]])
    }
    if (n$hasCounts) {
        out$N_CASE <- as.numeric(df[[n$ncaseSrc]])
        out$N_CONTROL <- as.numeric(df[[n$ncontrolSrc]])
    }
    for (key in c("BETA", "SE", "P", "MAF", "INFO")) {
        src <- .resolveSumstatKey(key, df, mapping, label)
        if (!is.na(src) && !is.null(src)) {
            out[[key]] <- as.numeric(df[[src]])
        }
    }
    out
}

# Default GWAS-VCF FORMAT tag mapping (canonical stat -> FORMAT field).
.gwasVcfFormatDefaults <- list(
    BETA = "ES",
    SE = "SE",
    P = "LP",
    N = "SS",
    MAF = "EAF"
)

# Pick the FORMAT sample column of a GWAS-VCF to read. Defaults to the single
# sample when there is exactly one; otherwise `sampleSelect` must name it.
.pickVcfSample <- function(samples, sampleSelect, label) {
    if (!is.null(sampleSelect)) {
        if (!(sampleSelect %in% samples)) {
            stop(
                label,
                ": sampleSelect '",
                sampleSelect,
                "' is not a sample in the VCF (have: ",
                paste(samples, collapse = ", "),
                ")."
            )
        }
        return(sampleSelect)
    }
    if (length(samples) == 1L) {
        return(samples[[1L]])
    }
    stop(
        label,
        ": the VCF has ",
        length(samples),
        " samples (",
        paste(samples, collapse = ", "),
        "); pass `sampleSelect` to choose the study column."
    )
}

# Read column `tag` from the GWAS-VCF geno matrix for one sample as numeric, or
# NULL when the tag is absent/unmapped.
# @noRd
.vcfGetField <- function(tag, geno, sample) {
    if (!is.null(tag) && tag %in% names(geno)) {
        as.numeric(geno[[tag]][, sample])
    } else {
        NULL
    }
}

# Convert a VCF (GWAS-VCF format) into the canonical sumstats data.frame. The
# effect allele (A1) is ALT; the other allele (A2) is REF. Stats come from the
# per-study FORMAT fields ES/SE/LP/SS/EAF, with Z = ES / SE.
.vcfToSumstatDf <- function(vcf, sampleSelect, formatMapping, label) {
    fmap <- modifyList(
        .gwasVcfFormatDefaults,
        if (is.null(formatMapping)) list() else as.list(formatMapping)
    )
    rr <- SummarizedExperiment::rowRanges(vcf)
    altList <- VariantAnnotation::alt(vcf)
    alt <- vapply(
        as.list(altList),
        function(a) as.character(a)[[1L]],
        character(1)
    )
    ref <- as.character(VariantAnnotation::ref(vcf))
    geno <- VariantAnnotation::geno(vcf)
    sample <- .pickVcfSample(colnames(vcf), sampleSelect, label)
    es <- .vcfGetField(fmap$BETA, geno, sample)
    se <- .vcfGetField(fmap$SE, geno, sample)
    lp <- .vcfGetField(fmap$P, geno, sample)
    ss <- .vcfGetField(fmap$N, geno, sample)
    eaf <- .vcfGetField(fmap$MAF, geno, sample)
    if (is.null(es) || is.null(se)) {
        stop(
            label,
            ": GWAS-VCF must carry ES and SE FORMAT fields to derive Z ",
            "(configure via `formatMapping`)."
        )
    }
    if (is.null(ss)) {
        stop(label, ": GWAS-VCF must carry an SS (sample size) FORMAT field.")
    }
    out <- data.frame(
        chrom = as.character(GenomicRanges::seqnames(rr)),
        pos = as.integer(GenomicRanges::start(rr)),
        SNP = names(rr),
        A1 = alt,
        A2 = ref,
        Z = es / se,
        N = ss,
        BETA = es,
        SE = se,
        stringsAsFactors = FALSE
    )
    if (!is.null(lp)) {
        out$P <- 10^(-lp)
    }
    if (!is.null(eaf)) {
        out$MAF <- pmin(eaf, 1 - eaf)
    }
    out
}

# Read one GWAS-VCF sumstats file (region-restricted only when bgzipped +
# tabix-indexed) into the canonical sumstats data.frame.
.readSumStatsVcf <- function(path, region, sampleSelect, formatMapping, label) {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        stop(
            label,
            ": reading VCF sumstats requires the 'VariantAnnotation' ",
            "package; please install it."
        )
    }
    hasTbi <- file.exists(paste0(path, ".tbi"))
    if (!is.null(region) && hasTbi) {
        tf <- Rsamtools::TabixFile(path)
        gr <- .matchRegionSeqnames(
            .asGRegion(region),
            Rsamtools::seqnamesTabix(tf)
        )
        vcf <- if (is.null(gr)) {
            VariantAnnotation::readVcf(path)[0, ]
        } else {
            VariantAnnotation::readVcf(tf, param = gr)
        }
    } else {
        if (!is.null(region) && !hasTbi) {
            warning(
                label,
                ": VCF '",
                path,
                "' is not bgzipped + tabix-indexed; ",
                "region ignored, reading the whole file."
            )
        }
        vcf <- VariantAnnotation::readVcf(path)
    }
    .vcfToSumstatDf(vcf, sampleSelect, formatMapping, label)
}

# Read one delimited-text sumstats file. Region reads only when a <path>.tbi
# sidecar exists; otherwise the whole file is read and any region is ignored.
.readSumStatsText <- function(
    path,
    region,
    columnMapping,
    label,
    allowNoN = FALSE
) {
    hasTbi <- file.exists(paste0(path, ".tbi"))
    raw <- if (!is.null(region) && hasTbi) {
        .readTabixRegion(path, region)
    } else {
        if (!is.null(region) && !hasTbi) {
            warning(
                label,
                ": '",
                path,
                "' is not tabix-indexed; ",
                "region ignored, reading the whole file."
            )
        }
        # Read every column as character; .resolveSumstatCols() coerces each
        # field to its target type. This prevents readr from type-guessing an
        # allele/ID column that is uniformly T/F as logical (the allele "T" ->
        # TRUE).
        as.data.frame(
            readr::read_tsv(
                path,
                show_col_types = FALSE,
                progress = FALSE,
                col_types = readr::cols(.default = readr::col_character())
            ),
            check.names = FALSE
        )
    }
    .resolveSumstatCols(raw, columnMapping, label, allowNoN = allowNoN)
}

# Dispatch a sumstats file to the text or GWAS-VCF reader. BCF is rejected.
.readSumStatsFile <- function(
    path,
    region,
    columnMapping,
    sampleSelect,
    formatMapping,
    label,
    allowNoN = FALSE
) {
    lower <- tolower(path)
    if (endsWith(lower, ".bcf")) {
        stop(
            label,
            ": BCF sumstats are not supported; convert '",
            path,
            "' to a bgzipped, tabix-indexed VCF."
        )
    }
    if (grepl("\\.vcf(\\.b?gz)?$", lower)) {
        # GWAS-VCF always carries an SS (sample-size) FORMAT field, so the study
        # scalar never has to stand in for a per-variant N here.
        .readSumStatsVcf(path, region, sampleSelect, formatMapping, label)
    } else {
        .readSumStatsText(
            path,
            region,
            columnMapping,
            label,
            allowNoN = allowNoN
        )
    }
}

# Read one sumstats source into an entry GRanges (canonical mcols).
.loadSumStatsEntry <- function(
    path,
    region,
    columnMapping,
    sampleSelect,
    formatMapping,
    label,
    allowNoN = FALSE
) {
    df <- .readSumStatsFile(
        path,
        region,
        columnMapping,
        sampleSelect,
        formatMapping,
        label,
        allowNoN = allowNoN
    )
    .dfToEntryGranges(df)
}

# =============================================================================
# Phenotype SummarizedExperiment builder (individual-level QTL)
# =============================================================================

# Read a covariate TSV into a samples-as-rows numeric matrix. The first column
# holds row identifiers; `transpose = TRUE` handles QTLtools-format inputs
# (covariates as rows, samples as columns).
.readCovariateMatrix <- function(path, transpose = FALSE) {
    raw <- as.data.frame(
        readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
        check.names = FALSE
    )
    rn <- as.character(raw[[1L]])
    m <- as.matrix(raw[, -1L, drop = FALSE])
    rownames(m) <- rn
    storage.mode(m) <- "double"
    if (transpose) {
        m <- t(m)
    }
    m
}

# Build one per-context SummarizedExperiment from a bgzipped BED phenotype
# file (QTLtools layout: #chr/start/end/ID + sample columns) and an optional
# per-context covariate TSV (attached as colData).
.buildContextSe <- function(bedPath, covPath, transposeCov) {
    bed <- .readPhenotypeBed(bedPath)
    cd <- .buildContextColData(covPath, transposeCov, bed$expr, bedPath)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = cd$expr),
        rowRanges = bed$rr,
        colData = cd$cd
    )
}

# Read a phenotype BED into an expression matrix (genes x samples) + rowRanges.
# @noRd
.readPhenotypeBed <- function(bedPath) {
    bed <- as.data.frame(
        readr::read_tsv(
            bedPath,
            show_col_types = FALSE,
            progress = FALSE,
            comment = ""
        ),
        check.names = FALSE
    )
    chrCol <- intersect(c("#chr", "#chrom", "chrom", "chr"), names(bed))[1L]
    startCol <- intersect(c("start", "Start"), names(bed))[1L]
    endCol <- intersect(c("end", "End"), names(bed))[1L]
    geneCol <- intersect(c("gene_id", "ID", "phenotype_id"), names(bed))[1L]
    if (any(is.na(c(chrCol, startCol, endCol, geneCol)))) {
        stop(
            "phenotype BED is missing one of chrom/start/end/gene_id columns: ",
            bedPath
        )
    }
    meta <- bed[, c(chrCol, startCol, endCol, geneCol)]
    names(meta) <- c("chrom", "start", "end", "gene_id")
    sampleCols <- setdiff(
        names(bed),
        c(
            "#chr",
            "#chrom",
            "chrom",
            "chr",
            "start",
            "Start",
            "end",
            "End",
            "gene_id",
            "ID",
            "phenotype_id",
            "strand",
            "Strand"
        )
    )
    expr <- as.matrix(bed[, sampleCols, drop = FALSE])
    storage.mode(expr) <- "double"
    rownames(expr) <- meta$gene_id
    rr <- GenomicRanges::GRanges(
        seqnames = meta$chrom,
        ranges = IRanges::IRanges(start = meta$start + 1L, end = meta$end)
    )
    names(rr) <- meta$gene_id
    list(expr = expr, rr = rr)
}

# Per-context colData from an optional covariate file (intersecting samples with
# the expression matrix, which is subset to the shared set). Returns list(expr,
# cd).
# @noRd
.buildContextColData <- function(covPath, transposeCov, expr, bedPath) {
    if (is.null(covPath) || !nzchar(covPath)) {
        return(list(
            expr = expr,
            cd = S4Vectors::DataFrame(row.names = colnames(expr))
        ))
    }
    pcov <- .readCovariateMatrix(covPath, transposeCov)
    common <- intersect(rownames(pcov), colnames(expr))
    if (length(common) == 0L) {
        stop(
            "No shared samples between phenotype and covariate file: ",
            bedPath
        )
    }
    list(
        expr = expr[, common, drop = FALSE],
        cd = S4Vectors::DataFrame(
            pcov[common, , drop = FALSE],
            row.names = common
        )
    )
}

# Aliases for the QtlDataset (phenotype) manifest.
.qtlDatasetManifestAliases <- list(
    context = c("context", "cond", "condition"),
    phenotypePath = c(
        "phenotypePath",
        "phenotype_path",
        "path",
        "phenotypeFile",
        "phenotype_file"
    ),
    covariatePath = c("covariatePath", "covariate_path", "cov_path", "covPath"),
    study = c("study", "study_id", "studyId"),
    genotypePath = c(
        "genotypePath",
        "genotype_path",
        "genotypePrefix",
        "genotype_prefix",
        "genotypeFile"
    ),
    genotypeCovariatePath = c(
        "genotypeCovariatePath",
        "genotype_covariate_path",
        "genotype_covariates",
        "genotypeCovariates"
    )
)

# Build one QtlDataset from an already-canonicalised set of manifest rows for a
# single study. `genotypesOverride` / `genoCovOverride` short-circuit the
# genotypePath / genotypeCovariatePath columns when supplied as arguments.
.buildQtlDatasetFromRows <- function(
    rows,
    study,
    base,
    transposeCov,
    qc,
    genotypesOverride = NULL,
    genoCovOverride = NULL
) {
    contexts <- unique(as.character(rows$context))
    phenotypes <- set_names(
        map(contexts, function(cx) {
            .buildOneContextSe(cx, rows, study, base, transposeCov)
        }),
        contexts
    )
    genoHandle <- .resolveGenoHandle(rows, study, base, genotypesOverride)
    genoCov <- .resolveGenoCov(
        rows,
        study,
        base,
        transposeCov,
        genoCovOverride
    )
    do.call(
        QtlDataset,
        c(
            list(
                study = study,
                genotypes = genoHandle,
                phenotypes = phenotypes,
                genotypeCovariates = genoCov
            ),
            qc
        )
    )
}

# The SummarizedExperiment for one context: exactly one phenotypePath (+ at most
# one covariatePath) across its manifest rows.
# @noRd
.buildOneContextSe <- function(cx, rows, study, base, transposeCov) {
    sub <- rows[as.character(rows$context) == cx, , drop = FALSE]
    pths <- unique(sub$phenotypePath[!is.na(sub$phenotypePath)])
    if (length(pths) != 1L) {
        stop(
            "Context '",
            cx,
            "' (study '",
            study,
            "') must reference exactly one phenotypePath; got: ",
            paste(pths, collapse = ", ")
        )
    }
    covPath <- NULL
    if ("covariatePath" %in% names(sub)) {
        covs <- unique(sub$covariatePath[
            !is.na(sub$covariatePath) & nzchar(as.character(sub$covariatePath))
        ])
        if (length(covs) > 1L) {
            stop(
                "Context '",
                cx,
                "' (study '",
                study,
                "') references multiple covariatePath values: ",
                paste(covs, collapse = ", ")
            )
        }
        if (length(covs) == 1L) {
            covPath <- .resolveRel(covs[[1L]], base)
        }
    }
    .buildContextSe(.resolveRel(pths[[1L]], base), covPath, transposeCov)
}

# The genotype handle: an override handle/path, else the study's single
# genotypePath (auto-detecting the format).
# @noRd
.resolveGenoHandle <- function(rows, study, base, genotypesOverride) {
    if (methods::is(genotypesOverride, "GenotypeHandle")) {
        return(genotypesOverride)
    }
    spec <- if (is.character(genotypesOverride)) {
        genotypesOverride
    } else {
        v <- unique(as.character(rows$genotypePath[!is.na(rows$genotypePath)]))
        if (length(v) != 1L) {
            stop(
                "study '",
                study,
                "' must reference exactly one genotypePath; got: ",
                paste(v, collapse = ", ")
            )
        }
        v[[1L]]
    }
    .detectGenotypeFormat(.resolveRel(spec, base))
}

# The genotype covariate matrix: an override matrix/path, else the study's
# single genotypeCovariatePath, else an empty matrix.
# @noRd
.resolveGenoCov <- function(rows, study, base, transposeCov, genoCovOverride) {
    if (is.matrix(genoCovOverride)) {
        return(genoCovOverride)
    }
    spec <- if (is.character(genoCovOverride)) {
        genoCovOverride
    } else if ("genotypeCovariatePath" %in% names(rows)) {
        v <- unique(as.character(rows$genotypeCovariatePath[
            !is.na(rows$genotypeCovariatePath) &
                nzchar(as.character(rows$genotypeCovariatePath))
        ]))
        if (length(v) > 1L) {
            stop(
                "study '",
                study,
                "' references multiple genotypeCovariatePath values."
            )
        }
        if (length(v) == 1L) v[[1L]] else NULL
    } else {
        NULL
    }
    if (is.null(spec)) {
        return(matrix(numeric(0), nrow = 0, ncol = 0))
    }
    .readCovariateMatrix(.resolveRel(spec, base), transposeCov)
}

# =============================================================================
# Public loaders
# =============================================================================

#' @title Load a QtlDataset from a manifest
#' @description Build a single-study \code{\link{QtlDataset}} from a manifest
#'   describing one row per QTL context (or per trait x context). The manifest
#'   is a data.frame or a path to a delimited file (\code{.csv} -> CSV, else
#'   TSV). Column names may be snake_case or camelCase.
#' @param manifest A data.frame or a path to a manifest file. Recognised columns
#'   (canonical camelCase; snake_case aliases accepted): \code{context}
#'   (required), \code{phenotypePath} (required, bgzipped BED),
#'   \code{covariatePath} (optional, per-context covariates), and the
#'   single-valued \code{study} / \code{genotypePath} /
#'   \code{genotypeCovariatePath}.
#' @param study Study identifier; reconciled with a \code{study} column.
#' @param genotypes A \code{\link{GenotypeHandle}} or a genotype path/prefix;
#'   reconciled with a \code{genotypePath} column.
#' @param genotypeCovariates A numeric matrix (samples x covariates) or a path
#'   to a covariate TSV; reconciled with a \code{genotypeCovariatePath} column.
#' @param scaleResiduals,mafCutoff,macCutoff,xvarCutoff Pass-through
#'   \code{\link{QtlDataset}} arguments (stored as lazy QC slots).
#' @param imissCutoff,keepSamples,keepVariants,keepIndel Pass-through
#'   \code{\link{QtlDataset}} arguments (stored as lazy QC slots).
#' @param transposeCovariates Transpose covariate TSVs (QTLtools layout) before
#'   treating them as samples-as-rows.
#' @return A \code{QtlDataset} object.
#' @examples
#' d <- system.file("extdata", "qtl_mini", package = "pecotmr")
#' manifest <- data.frame(context = "brain",
#'   phenotypePath = file.path(d, "example_geneexpr.bed.gz"),
#'   study = "s1", genotypePath = file.path(d, "example.chr22"))
#' loadQtlDatasetFromManifest(manifest = manifest, study = "s1")
#' @export
loadQtlDatasetFromManifest <- function(
    manifest,
    study = NULL,
    genotypes = NULL,
    genotypeCovariates = NULL,
    scaleResiduals = TRUE,
    mafCutoff = 0,
    macCutoff = 0,
    xvarCutoff = 0,
    imissCutoff = 0,
    keepSamples = character(0),
    keepVariants = character(0),
    keepIndel = TRUE,
    transposeCovariates = FALSE
) {
    base <- .manifestBase(manifest)
    df <- .canonManifestCols(
        .readManifest(manifest),
        .qtlDatasetManifestAliases,
        required = c("context", "phenotypePath"),
        label = "QtlDataset"
    )
    study <- .reconcileScalar(df$study, study, "study")
    qc <- list(
        scaleResiduals = scaleResiduals,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        xvarCutoff = xvarCutoff,
        imissCutoff = imissCutoff,
        keepSamples = keepSamples,
        keepVariants = keepVariants,
        keepIndel = keepIndel
    )
    .buildQtlDatasetFromRows(
        df,
        study,
        base,
        transposeCovariates,
        qc,
        genotypesOverride = genotypes,
        genoCovOverride = genotypeCovariates
    )
}

# Aliases shared by the sumstats manifests.
.gwasSumStatsManifestAliases <- list(
    study = c("study", "study_id", "studyId"),
    sumStatsPath = c(
        "sumStatsPath",
        "sumstats_path",
        "path",
        "file_path",
        "gwas_tsv",
        "filePath"
    ),
    columnMapping = c(
        "columnMapping",
        "column_mapping",
        "column_mapping_file",
        "columnMappingFile"
    ),
    nCase = c("nCase", "n_case"),
    nControl = c("nControl", "n_control"),
    nSample = c("nSample", "n_sample"),
    varY = c("varY", "var_y"),
    genome = c("genome"),
    ldSketchPath = c("ldSketchPath", "ld_sketch_path", "ld_sketch")
)

.qtlSumStatsManifestAliases <- list(
    study = c("study", "study_id", "studyId"),
    context = c("context", "cond", "condition"),
    trait = c("trait", "gene_id", "geneId", "molecular_trait_id"),
    sumStatsPath = c(
        "sumStatsPath",
        "sumstats_path",
        "path",
        "file_path",
        "filePath"
    ),
    columnMapping = c(
        "columnMapping",
        "column_mapping",
        "column_mapping_file",
        "columnMappingFile"
    ),
    nSample = c("nSample", "n_sample"),
    varY = c("varY", "var_y"),
    genome = c("genome"),
    ldSketchPath = c("ldSketchPath", "ld_sketch_path", "ld_sketch")
)

# TRUE if study `i` has usable case/control counts OR a total-N scalar, so a
# missing per-variant N is acceptable for that study.
# @noRd
.manifestHasStudyScalar <- function(i, nCaseCol, nControlCol, nSampleCol) {
    ccOk <- !is.null(nCaseCol) &&
        !is.null(nControlCol) &&
        is.finite(nCaseCol[[i]]) &&
        is.finite(nControlCol[[i]])
    nOk <- !is.null(nSampleCol) && is.finite(nSampleCol[[i]])
    isTRUE(ccOk) || isTRUE(nOk)
}

#' @title Load a GwasSumStats collection from a manifest
#' @description Build a \code{\link{GwasSumStats}} from a manifest with one row
#'   per study. No QC is run (the result carries \code{qcInfo = list()}). Each
#'   sumstats file needs a \code{z} column, or \code{beta}+\code{se} from which
#'   the Wald z (\code{z = beta/se}) is derived when \code{z} is absent (a
#'   supplied \code{z} takes precedence).
#' @param manifest A data.frame or path. Columns (snake_case aliases accepted):
#'   \code{study} (required, unique), \code{sumStatsPath} (required),
#'   \code{columnMapping} (optional), \code{nCase} / \code{nControl} (optional),
#'   \code{nSample} (optional study-level total N), \code{varY} (optional), and
#'   the single-valued \code{genome} / \code{ldSketchPath}. When a row supplies
#'   a study-level scalar (\code{nCase} + \code{nControl}, or \code{nSample}),
#'   its sumstats file need not carry a per-variant \code{N} column;
#'   \code{\link{summaryStatsQc}} fills \code{N} from the scalar.
#' @param genome Genome build; reconciled with a \code{genome} column.
#' @param ldSketch A \code{\link{GenotypeHandle}} or a spec
#'   (path/prefix/genoMeta); reconciled with an \code{ldSketchPath} column.
#' @param region Optional \code{chr:start-end} string, GRanges, or one-row
#'   data.frame restricting the variants read (honoured only for tabix-indexed
#'   text and bgzipped+tabixed VCF; ignored with a warning otherwise).
#' @param minLdOverlapWarn Warn when the fraction of sumstats variants present
#'   in the LD sketch falls below this (default 0.5).
#' @param columnMapping Optional default column mapping (a named list/vector, or
#'   a path to a YAML file of \code{standardName: sourceName} entries) applied
#'   when a row has no \code{columnMapping}.
#' @param sampleSelect Optional GWAS-VCF FORMAT sample (study) column to read.
#' @param formatMapping Optional GWAS-VCF FORMAT tag mapping (canonical stat ->
#'   FORMAT field), overriding ES/SE/LP/SS/EAF defaults.
#' @return A \code{GwasSumStats} object.
#' @examples
#' gwasTsv <- system.file("extdata", "manifests",
#'   "protocol_example.twas.gwas_sumstats.chr22.tsv.gz", package = "pecotmr")
#' ldStem <- file.path(system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr"), "protocol_example.LD.chr22")
#' manifest <- data.frame(study = "g1", sumStatsPath = gwasTsv)
#' loadGwasSumStatsFromManifest(manifest = manifest, genome = "hg38",
#'   ldSketch = ldStem, region = "chr22:10000000-19000000")
#' @export
loadGwasSumStatsFromManifest <- function(
    manifest,
    genome = NULL,
    ldSketch = NULL,
    region = NULL,
    minLdOverlapWarn = 0.5,
    columnMapping = NULL,
    sampleSelect = NULL,
    formatMapping = NULL
) {
    base <- .manifestBase(manifest)
    df <- .canonManifestCols(
        .readManifest(manifest),
        .gwasSumStatsManifestAliases,
        required = c("study", "sumStatsPath"),
        label = "GwasSumStats"
    )
    if (anyDuplicated(as.character(df$study))) {
        stop("GwasSumStats manifest `study` values must be unique.")
    }
    genome <- .reconcileScalar(df$genome, genome, "genome")
    ldSketchSpec <- .resolveLdSketchInput(df, ldSketch, base)
    ns <- .gwasStudyScalars(df)
    entries <- map(seq_len(nrow(df)), function(i) {
        .loadGwasEntry(
            i,
            df,
            base,
            region,
            columnMapping,
            sampleSelect,
            formatMapping,
            ns
        )
    })
    # Materialise the LD sketch reading only the chromosomes the sumstats cover,
    # then run the deferred per-study containment checks and trim to range.
    ldSketch <- .materializeLdSketch(ldSketchSpec, .entriesChroms(entries))
    .checkGwasLdContainment(ldSketch, entries, df, minLdOverlapWarn)
    ldSketch <- .subsetSketchToRange(ldSketch, entries)
    do.call(GwasSumStats, .gwasSumStatsArgs(df, entries, genome, ldSketch, ns))
}

# Study-level sample-size scalars from the manifest (n_case + n_control, or
# n_sample). When present a row need not carry a per-variant N.
# @noRd
.gwasStudyScalars <- function(df) {
    list(
        nCase = if ("nCase" %in% names(df)) as.numeric(df$nCase) else NULL,
        nControl = if ("nControl" %in% names(df)) {
            as.numeric(df$nControl)
        } else {
            NULL
        },
        nSample = if ("nSample" %in% names(df)) as.numeric(df$nSample) else NULL
    )
}

# Load one GwasSumStats manifest row's entry (resolving a per-row columnMapping
# override), allowing a missing per-variant N when a study scalar is present.
# @noRd
.loadGwasEntry <- function(
    i,
    df,
    base,
    region,
    columnMapping,
    sampleSelect,
    formatMapping,
    ns
) {
    label <- paste0("GwasSumStats[study=", df$study[[i]], "]")
    mapping <- if (
        "columnMapping" %in%
            names(df) &&
            !is.na(df$columnMapping[[i]]) &&
            nzchar(as.character(df$columnMapping[[i]]))
    ) {
        .resolveRel(as.character(df$columnMapping[[i]]), base)
    } else {
        columnMapping
    }
    .loadSumStatsEntry(
        .resolveRel(as.character(df$sumStatsPath[[i]]), base),
        region,
        mapping,
        sampleSelect,
        formatMapping,
        label,
        allowNoN = .manifestHasStudyScalar(
            i,
            ns$nCase,
            ns$nControl,
            ns$nSample
        )
    )
}

# Per-study LD-containment checks (deferred until the sketch is materialised).
# @noRd
.checkGwasLdContainment <- function(ldSketch, entries, df, minLdOverlapWarn) {
    for (i in seq_len(nrow(df))) {
        .checkLdContainment(
            ldSketch,
            entries[[i]],
            minLdOverlapWarn,
            paste0("GwasSumStats[study=", df$study[[i]], "]")
        )
    }
}

# GwasSumStats() constructor args, including any present study scalars + varY.
# @noRd
.gwasSumStatsArgs <- function(df, entries, genome, ldSketch, ns) {
    args <- list(
        study = as.character(df$study),
        entry = entries,
        genome = genome,
        ldSketch = ldSketch
    )
    if (!is.null(ns$nCase)) {
        args$nCase <- ns$nCase
    }
    if (!is.null(ns$nControl)) {
        args$nControl <- ns$nControl
    }
    if (!is.null(ns$nSample)) {
        args$nSample <- ns$nSample
    }
    if ("varY" %in% names(df)) {
        args$varY <- as.numeric(df$varY)
    }
    args
}

# Build the entry list for a QtlSumStats manifest. `allowNoN` is an optional
# per-row logical vector: when TRUE for a row, its sumstats file need not carry
# a per-variant N (a study-level nSample scalar fills it later in
# summaryStatsQc). NULL (default) requires a per-variant N on every row. The LD
# containment check is deferred to the caller so the sketch can first be read
# restricted to the chromosomes these entries cover.
.loadQtlSumStatsEntries <- function(
    df,
    base,
    region,
    columnMapping,
    sampleSelect,
    formatMapping,
    allowNoN = NULL
) {
    lapply(seq_len(nrow(df)), function(i) {
        label <- paste0(
            "QtlSumStats[",
            df$study[[i]],
            "/",
            df$context[[i]],
            "/",
            df$trait[[i]],
            "]"
        )
        mapping <- if (
            "columnMapping" %in%
                names(df) &&
                !is.na(df$columnMapping[[i]]) &&
                nzchar(as.character(df$columnMapping[[i]]))
        ) {
            .resolveRel(as.character(df$columnMapping[[i]]), base)
        } else {
            columnMapping
        }
        .loadSumStatsEntry(
            .resolveRel(as.character(df$sumStatsPath[[i]]), base),
            region,
            mapping,
            sampleSelect,
            formatMapping,
            label,
            allowNoN = !is.null(allowNoN) && isTRUE(allowNoN[[i]])
        )
    })
}

#' @title Load a QtlSumStats collection from a manifest
#' @description Build a \code{\link{QtlSumStats}} from a manifest with one row
#'   per \code{(study, context, trait)} tuple. No QC is run. Each sumstats file
#'   needs a \code{z} column, or \code{beta}+\code{se} from which the Wald z
#'   (\code{z = beta/se}) is derived when \code{z} is absent (a supplied
#'   \code{z} takes precedence).
#' @param manifest A data.frame or path. Columns (snake_case aliases accepted):
#'   \code{study}, \code{context}, \code{trait} (required), \code{sumStatsPath}
#'   (required), \code{columnMapping} (optional), \code{nSample} (optional
#'   tuple-level total N), \code{varY} (optional), and the single-valued
#'   \code{genome} / \code{ldSketchPath}. When a row supplies \code{nSample},
#'   its sumstats file need not carry a per-variant \code{N} column;
#'   \code{\link{summaryStatsQc}} fills \code{N} from the scalar. (Unlike the
#'   GWAS loader, there are no \code{nCase}/\code{nControl} columns: molecular
#'   QTL traits are quantitative.)
#' @param genome Genome build; reconciled with a \code{genome} column.
#' @param ldSketch A \code{\link{GenotypeHandle}} or spec; reconciled with an
#'   \code{ldSketchPath} column.
#' @param region,minLdOverlapWarn,columnMapping,sampleSelect,formatMapping As
#'   for \code{\link{loadGwasSumStatsFromManifest}}.
#' @return A \code{QtlSumStats} object.
#' @examples
#' tsv <- system.file("extdata", "manifests",
#'   "protocol_example.twas.gwas_sumstats.chr22.tsv.gz", package = "pecotmr")
#' ldStem <- file.path(system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr"), "protocol_example.LD.chr22")
#' manifest <- data.frame(study = "s1", context = "brain",
#'   trait = "ENSG1", sumStatsPath = tsv)
#' loadQtlSumStatsFromManifest(manifest = manifest, genome = "hg38",
#'   ldSketch = ldStem, region = "chr22:10000000-19000000")
#' @export
loadQtlSumStatsFromManifest <- function(
    manifest,
    genome = NULL,
    ldSketch = NULL,
    region = NULL,
    minLdOverlapWarn = 0.5,
    columnMapping = NULL,
    sampleSelect = NULL,
    formatMapping = NULL
) {
    base <- .manifestBase(manifest)
    df <- .canonManifestCols(
        .readManifest(manifest),
        .qtlSumStatsManifestAliases,
        required = c("study", "context", "trait", "sumStatsPath"),
        label = "QtlSumStats"
    )
    genome <- .reconcileScalar(df$genome, genome, "genome")
    ldSketchSpec <- .resolveLdSketchInput(df, ldSketch, base)

    # Tuple-level total-N scalar (from the manifest). When a row carries a
    # usable nSample the sumstats file need not supply a per-variant N:
    # summaryStatsQc fills N from the scalar. Gate the entry reader's N check
    # per row on it.
    nSampleCol <- if ("nSample" %in% names(df)) as.numeric(df$nSample) else NULL
    allowNoN <- if (!is.null(nSampleCol)) is.finite(nSampleCol) else NULL

    entries <- .loadQtlSumStatsEntries(
        df,
        base,
        region,
        columnMapping,
        sampleSelect,
        formatMapping,
        allowNoN = allowNoN
    )

    # Materialise the LD sketch reading only the chromosomes the summary stats
    # cover (skips other shards of a chrom-sharded panel), check containment,
    # then trim its snpInfo to the summary stats' per-chromosome position span.
    ldSketch <- .materializeLdSketch(ldSketchSpec, .entriesChroms(entries))
    .qtlSumStatsCheckContainment(ldSketch, entries, df, minLdOverlapWarn)
    ldSketch <- .subsetSketchToRange(ldSketch, entries)

    do.call(
        QtlSumStats,
        .qtlSumStatsArgs(df, entries, genome, ldSketch, nSampleCol)
    )
}

# Per-row LD-containment check (no-op when the sketch is NULL).
# @noRd
.qtlSumStatsCheckContainment <- function(
    ldSketch,
    entries,
    df,
    minLdOverlapWarn
) {
    if (is.null(ldSketch)) {
        return(invisible(NULL))
    }
    for (i in seq_len(nrow(df))) {
        .checkLdContainment(
            ldSketch,
            entries[[i]],
            minLdOverlapWarn,
            paste0(
                "QtlSumStats[",
                df$study[[i]],
                "/",
                df$context[[i]],
                "/",
                df$trait[[i]],
                "]"
            )
        )
    }
    invisible(NULL)
}

# Assemble the QtlSumStats() argument list from a canonicalised manifest.
# @noRd
.qtlSumStatsArgs <- function(df, entries, genome, ldSketch, nSampleCol) {
    args <- list(
        study = as.character(df$study),
        context = as.character(df$context),
        trait = as.character(df$trait),
        entry = entries,
        genome = genome,
        ldSketch = ldSketch
    )
    if (!is.null(nSampleCol)) {
        args$nSample <- nSampleCol
    }
    if ("varY" %in% names(df)) {
        args$varY <- as.numeric(df$varY)
    }
    args
}

#' @title Load a MultiStudyQtlDataset from manifests
#' @description Build a \code{\link{MultiStudyQtlDataset}} from a QtlDataset
#'   manifest (individual-level studies) and, optionally, a QtlSumStats manifest
#'   (summary-only studies).
#' @param qtlDatasetsManifest A data.frame or path. The QtlDataset schema with
#'   \code{study} and \code{genotypePath} \strong{mandatory} (grouping key =
#'   \code{study}); one \code{QtlDataset} is built per study group.
#' @param sumStatsManifest Optional QtlSumStats manifest (see
#'   \code{\link{loadQtlSumStatsFromManifest}}).
#' @param genome,ldSketch,region,minLdOverlapWarn Passed to the QtlSumStats
#'   loader for \code{sumStatsManifest}.
#' @param columnMapping,sampleSelect,formatMapping Passed to the QtlSumStats
#'   loader for \code{sumStatsManifest}.
#' @param transposeCovariates Transpose covariate TSVs (QTLtools layout).
#' @param scaleResiduals,mafCutoff,macCutoff,xvarCutoff Pass-through
#'   \code{\link{QtlDataset}} arguments applied to every study.
#' @param imissCutoff,keepSamples,keepVariants,keepIndel Pass-through
#'   \code{\link{QtlDataset}} arguments applied to every study.
#' @return A \code{MultiStudyQtlDataset} object.
#' @examples
#' d <- system.file("extdata", "qtl_mini", package = "pecotmr")
#' manifest <- data.frame(study = c("s1", "s2"), context = "brain",
#'   phenotypePath = file.path(d, "example_geneexpr.bed.gz"),
#'   genotypePath = file.path(d, "example.chr22"))
#' loadMultiStudyQtlDatasetFromManifest(qtlDatasetsManifest = manifest)
#' @export
loadMultiStudyQtlDatasetFromManifest <- function(
    qtlDatasetsManifest,
    sumStatsManifest = NULL,
    genome = NULL,
    ldSketch = NULL,
    region = NULL,
    minLdOverlapWarn = 0.5,
    columnMapping = NULL,
    sampleSelect = NULL,
    formatMapping = NULL,
    transposeCovariates = FALSE,
    scaleResiduals = TRUE,
    mafCutoff = 0,
    macCutoff = 0,
    xvarCutoff = 0,
    imissCutoff = 0,
    keepSamples = character(0),
    keepVariants = character(0),
    keepIndel = TRUE
) {
    qc <- .msqQcArgs(
        scaleResiduals,
        mafCutoff,
        macCutoff,
        xvarCutoff,
        imissCutoff,
        keepSamples,
        keepVariants,
        keepIndel
    )
    qtlDatasets <- .msqBuildDatasets(
        qtlDatasetsManifest,
        transposeCovariates,
        qc
    )
    sumStats <- .msqLoadSumStats(
        sumStatsManifest,
        genome,
        ldSketch,
        region,
        minLdOverlapWarn,
        columnMapping,
        sampleSelect,
        formatMapping
    )
    MultiStudyQtlDataset(qtlDatasets = qtlDatasets, sumStats = sumStats)
}

# Bundle the per-study QC / sample-filter arguments into a list.
# @noRd
.msqQcArgs <- function(
    scaleResiduals,
    mafCutoff,
    macCutoff,
    xvarCutoff,
    imissCutoff,
    keepSamples,
    keepVariants,
    keepIndel
) {
    list(
        scaleResiduals = scaleResiduals,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        xvarCutoff = xvarCutoff,
        imissCutoff = imissCutoff,
        keepSamples = keepSamples,
        keepVariants = keepVariants,
        keepIndel = keepIndel
    )
}

# Read the qtlDatasets manifest and build one QtlDataset per study.
# @noRd
.msqBuildDatasets <- function(qtlDatasetsManifest, transposeCovariates, qc) {
    base <- .manifestBase(qtlDatasetsManifest)
    df <- .canonManifestCols(
        .readManifest(qtlDatasetsManifest),
        .qtlDatasetManifestAliases,
        required = c("study", "genotypePath", "context", "phenotypePath"),
        label = "MultiStudyQtlDataset (qtlDatasetsManifest)"
    )
    studies <- unique(as.character(df$study))
    setNames(
        lapply(studies, function(s) {
            .buildQtlDatasetFromRows(
                df[as.character(df$study) == s, , drop = FALSE],
                s,
                base,
                transposeCovariates,
                qc
            )
        }),
        studies
    )
}

# Load the optional companion sumStats manifest (NULL when none supplied).
# @noRd
.msqLoadSumStats <- function(
    sumStatsManifest,
    genome,
    ldSketch,
    region,
    minLdOverlapWarn,
    columnMapping,
    sampleSelect,
    formatMapping
) {
    if (is.null(sumStatsManifest)) {
        return(NULL)
    }
    loadQtlSumStatsFromManifest(
        sumStatsManifest,
        genome = genome,
        ldSketch = ldSketch,
        region = region,
        minLdOverlapWarn = minLdOverlapWarn,
        columnMapping = columnMapping,
        sampleSelect = sampleSelect,
        formatMapping = formatMapping
    )
}
