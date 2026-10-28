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
# sumstats objects carry qcOptions = list()) and never filter genotypes (that is
# the pipelines' GenotypeFilterParam, applied per extracted block).
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

# Read a manifest into a tibble. A data.frame/tibble is passed through; a
# single path is read by extension (.csv -> read_csv, else read_tsv).
#' @importFrom checkmate assertString
#' @importFrom checkmate assertFileExists
.readManifest <- function(manifest) {
    if (is.data.frame(manifest)) {
        return(as_tibble(manifest, .name_repair = "minimal"))
    }
    assertString(manifest, .var.name = "manifest (data.frame or file path)")
    assertFileExists(manifest, access = "r", .var.name = "manifest file")
    if (str_detect(manifest, regex("\\.csv$", ignore_case = TRUE))) {
        readr::read_csv(manifest, show_col_types = FALSE, progress = FALSE)
    } else {
        readr::read_tsv(manifest, show_col_types = FALSE, progress = FALSE)
    }
}

# The windows a BED file holds, as a GRanges. BED's columns are positional --
# chromosome, start and end, then a name for the window -- and its start is
# 0-based, so it is shifted by one here. Blank, comment, `track` and
# `browser` lines are skipped, and so is a header line naming the columns,
# which LD-block files often have.
# @noRd
.readBedWindows <- function(path, what) {
    lines <- discard(
        readr::read_lines(path, skip_empty_rows = TRUE, progress = FALSE),
        str_detect,
        "^(#|track|browser)"
    )
    record <- str_detect(lines, "^[^\t]+\t[0-9]+\t[0-9]+(\t|$)")
    if (length(lines) > 0L && !record[[1L]]) {
        lines <- lines[-1L]
        record <- record[-1L]
    }
    if (length(lines) == 0L || !all(record)) {
        abort(glue(
            "{what}: '{path}' is not a BED file with windows in it; each ",
            "line must give a chromosome, a start and an end, tab-separated."
        ))
    }
    fields <- str_split(lines, "\t", simplify = TRUE)
    exec(
        GenomicRanges::GRanges,
        seqnames = fields[, 1L],
        ranges = IRanges::IRanges(
            start = as.integer(fields[, 2L]) + 1L,
            end = as.integer(fields[, 3L])
        ),
        !!!compact(list(blockId = if (ncol(fields) > 3L) fields[, 4L]))
    )
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
    if (is.null(p) || length(p) != 1L || is.na(p) || str_length(p) == 0L) {
        return(p)
    }
    if (str_detect(p, "^(/|[A-Za-z]:|~)") || is.null(base) || file.exists(p)) {
        return(p)
    }
    file.path(base, p)
}

#' @importFrom checkmate assertNames
# The existing column an alias set points at, or NULL when the canonical name
# is already present or nothing matches.
# @noRd
.manifestAliasHit <- function(canon, aliases, df) {
    if (is_in(canon, names(df))) {
        return(NULL)
    }
    hit <- intersect(aliases[[canon]], names(df))
    if (length(hit) == 0L) {
        return(NULL)
    }
    hit[[1L]]
}

# Rename accepted aliases to the canonical camelCase name and validate that all
# `required` canonical columns are present. `aliases` is a named list keyed by
# canonical name; each value is the character vector of accepted source names
# (including the canonical name itself). The first alias present wins.
# @noRd
.canonManifestCols <- function(df, aliases, required, label) {
    # Each alias resolves against the ORIGINAL names, so renaming one cannot
    # change what another matches.
    renames <- compact(set_names(
        map(names(aliases), .manifestAliasHit, aliases = aliases, df = df),
        names(aliases)
    ))
    df <- rename(df, !!!renames)
    assertNames(
        names(df),
        must.include = required,
        what = "colnames",
        .var.name = str_c(label, " manifest")
    )
    df
}

# Resolve a single collection-level scalar (e.g. genome) that may be supplied
# both as a constant manifest column and as a function argument. The column (if
# present) must be constant; a supplied argument must agree with it.
.reconcileScalar <- function(colValues, arg, name, required = TRUE) {
    colVal <- if (is.null(colValues)) {
        NULL
    } else {
        v <- unique(as.character(colValues[
            !is.na(colValues) &
                str_length(as.character(colValues)) > 0L
        ]))
        if (length(v) > 1L) {
            msg <- glue(
                "`{name}` column must be constant across the manifest; ",
                "got: {str_flatten(v, ', ')}"
            )
            abort(msg)
        }
        if (length(v) == 1L) v else NULL
    }
    if (!is.null(arg) && !is.null(colVal)) {
        if (!identical(as.character(arg), colVal)) {
            msg <- glue(
                "`{name}` argument ('{arg}') disagrees with the manifest ",
                "column value ('{colVal}')."
            )
            abort(msg)
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
        msg <- glue(
            "`{name}` must be provided as an argument or a manifest column."
        )
        abort(msg)
    }
    NULL
}

# The build a summary-statistics load names: the one given, checked against
# any its GWAS-VCF files record on ##reference, or else theirs. Files that
# disagree, or that disagree with the one given, are refused.
# @noRd
.manifestBuild <- function(genome, recorded, caller) {
    if (length(recorded) > 1L) {
        abort(glue(
            "{caller}: the GWAS-VCF files record different genome builds ",
            "({str_flatten(recorded, ', ')})."
        ))
    }
    if (is.null(genome)) {
        if (length(recorded) == 1L) {
            return(recorded)
        }
        abort("`genome` must be provided as an argument or a manifest column.")
    }
    if (length(recorded) == 1L && !identical(as.character(genome), recorded)) {
        abort(glue(
            "{caller}: `genome` is '{genome}', but the GWAS-VCF files ",
            "record '{recorded}'."
        ))
    }
    genome
}

# The build the manifest's rows are on: the `genome` argument or column,
# reconciled with the builds its GWAS-VCF files record.
# @noRd
.manifestGenome <- function(df, genome, base, caller) {
    .manifestBuild(
        .reconcileScalar(df[["genome"]], genome, "genome", required = FALSE),
        .manifestVcfBuilds(df, base),
        caller
    )
}

# The genome builds the manifest's GWAS-VCF files record, read from their
# headers alone.
# @noRd
.manifestVcfBuilds <- function(df, base) {
    paths <- map_chr(as.character(df$sumStatsPath), .resolveRel, base = base)
    vcfs <- paths[str_detect(str_to_lower(paths), "\\.vcf(\\.b?gz)?$")]
    if (length(vcfs) == 0L) {
        return(character(0))
    }
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        return(character(0))
    }
    unique(discard(map_chr(vcfs, .vcfHeaderBuild), is.na))
}

# @noRd
.vcfHeaderBuild <- function(path) {
    .vcfMetaBuild(VariantAnnotation::scanVcfHeader(path))
}

# =============================================================================
# Genotype-format detection + LD-sketch resolution
# =============================================================================

# The genotype panel of a single file or PLINK prefix, its format told as
# readGenotypes() tells it (.genotypePathFormat()).
.readGenotypePath <- function(path) {
    .readGenotypePanel(
        path,
        format = .requireGenotypeFormat(path, "manifest")
    )
}

# Turn an ldReference specification into a genotype panel, or pass an LdMatrix
# or a mixture through. Accepts a panel (or a bare handle), a named chrom ->
# path vector (genoMeta), or a single string that is either a genotype
# file/prefix or a chrom-sharded genoMeta meta file. `chroms` (optional)
# restricts a chrom-sharded panel to those chromosomes so the other
# per-chromosome shards are never read; it is ignored for a single genotype
# file (which has no shards to skip).
.resolveLdRef <- function(spec, chroms = NULL) {
    if (is.null(spec)) {
        return(NULL)
    }
    open <- .openLdRef(spec)
    if (!is.null(open)) {
        return(open)
    }
    if (is.character(spec) && !is.null(names(spec))) {
        return(.chromMetaPanel(spec, chroms))
    }
    if (is.character(spec) && length(spec) == 1L) {
        if (!is.null(.genotypePathFormat(spec))) {
            return(.readGenotypePath(spec))
        }
        # Otherwise treat it as a chrom-sharded genoMeta meta file.
        return(.chromMetaPanel(spec, chroms))
    }
    msg <- glue(
        "`ldReference` must be a genotype panel, a genotype path/prefix, a ",
        "genoMeta spec (named chrom->path vector or meta-file path), or an ",
        "LdMatrix."
    )
    abort(msg)
}

# An LD reference that is already open -- a genotype panel or a bare handle,
# an LdMatrix or a mixture, as it stands -- or NULL for a spec still to be
# resolved.
# @noRd
.openLdRef <- function(x) {
    if (.ldIsMatrixRef(x) || .ldIsMixtureRef(x)) {
        return(x)
    }
    .openGenotypeSource(x)
}

# The LD sketch SPEC of every manifest row, without reading any genotype
# metadata yet: the `ldReference` argument for all of them, or each row's own
# `ldReferencePath` (resolved against the manifest's directory). Rows may name
# different references; where a row's path and the argument are both given
# they must agree. An already-open panel passed as the argument is returned as
# it is, for every row. Deferring the reads to .manifestRowLdRefs() lets the
# loader first learn which chromosomes each reference's rows cover and skip
# the panel's other per-chromosome shards.
# @noRd
.resolveLdRefSpecs <- function(df, ldReference, base) {
    n <- nrow(df)
    open <- .openLdRef(ldReference)
    if (!is.null(open)) {
        return(rep(list(open), n))
    }
    paths <- .manifestLdRefPaths(df, base)
    given <- !is.na(paths)
    if (!is.null(ldReference)) {
        agrees <- map_lgl(paths, identical, as.character(ldReference))
        if (any(given & !agrees)) {
            msg <- glue(
                "`ldReference` argument disagrees with the manifest ",
                "`ldReferencePath` column."
            )
            abort(msg)
        }
        return(rep(list(ldReference), n))
    }
    if (!any(given)) {
        msg <- glue(
            "`ldReference` must be provided as an argument or an ",
            "`ldReferencePath` column."
        )
        abort(msg)
    }
    if (!all(given)) {
        abort(glue(
            "every manifest row needs an LD reference: row(s) ",
            "{str_flatten(which(!given), ', ')} have no `ldReferencePath`, ",
            "and no `ldReference` argument was given."
        ))
    }
    as.list(paths)
}

# Each manifest row's own `ldReferencePath`, resolved against the manifest's
# directory; NA where the column is absent or the cell blank. is_in guard,
# not `!is.null(df$ldReferencePath)`: on a tibble `$` for an absent column
# NULLs but WARNS ("Unknown or uninitialised column").
# @noRd
.manifestLdRefPaths <- function(df, base) {
    if (!is_in("ldReferencePath", names(df))) {
        return(rep(NA_character_, nrow(df)))
    }
    raw <- as.character(df$ldReferencePath)
    keep <- which(!is.na(raw) & str_length(raw) > 0L)
    replace(
        rep(NA_character_, length(raw)),
        keep,
        map_chr(raw[keep], .resolveRel, base = base)
    )
}

# The LD reference of every manifest row. Each distinct spec is read once,
# only for the chromosomes its rows cover, checked for containment against
# its rows with `check` (.checkGwasLdContainment() or
# .qtlSumStatsCheckContainment()), and trimmed to their span; rows on one
# spec share that one reference, so the collection stores it once.
# @noRd
.manifestRowLdRefs <- function(specs, entries, df, minLdOverlapWarn, check) {
    keys <- map_chr(specs, .ldSpecKey)
    groups <- split(seq_along(specs), factor(keys, levels = unique(keys)))
    refs <- map(
        groups,
        .manifestGroupLdRef,
        specs = specs,
        entries = entries,
        df = df,
        minLdOverlapWarn = minLdOverlapWarn,
        check = check
    )
    unname(refs[match(keys, names(groups))])
}

# One spec's reference for the manifest rows naming it. The containment
# checks sit between the read and the trim because containment is only
# answerable once the sketch is read, and trimming would hide a miss.
# @noRd
.manifestGroupLdRef <- function(
    rows,
    specs,
    entries,
    df,
    minLdOverlapWarn,
    check
) {
    rowEntries <- entries[rows]
    materialized <- .materializeLdRef(
        specs[[rows[[1L]]]],
        .entriesChroms(rowEntries)
    )
    check(materialized, rowEntries, slice(df, rows), minLdOverlapWarn)
    .subsetLdRefToRange(materialized, rowEntries)
}

# A key telling LD reference specs apart: a path, or per-chromosome paths, by
# its text; an already-open reference by the object itself.
# @noRd
.ldSpecKey <- function(spec) {
    if (!is.character(spec)) {
        return(obj_address(spec))
    }
    str_flatten(str_c(names(spec) %||% "", "=", spec), "|")
}

# Materialise a resolved ldReference spec into a genotype panel, reading only
# the shards for `chroms` when the spec is a chrom-sharded panel. A spec that is
# already open (the caller passed a panel) is returned unchanged: its variants
# are already in memory, so there is no read to restrict.
.materializeLdRef <- function(spec, chroms) {
    if (is.null(spec)) {
        return(NULL)
    }
    open <- .openLdRef(spec)
    if (!is.null(open)) {
        return(open)
    }
    .resolveLdRef(spec, chroms = chroms)
}

# The chromosomes present across a list of sumstats entry GRanges. NULL /
# empty entries contribute nothing; NA seqnames are dropped. Always returns a
# character vector (character(0) when nothing is present, never NULL).
.entriesChroms <- function(entries) {
    # as.character() is load-bearing: an all-empty list concatenates to NULL,
    # and this helper promises character(0).
    ch <- as.character(unname(list_c(map(entries, .mlEntryChroms))))
    unique(ch[!is.na(ch)])
}

# =============================================================================
# LD-sketch containment check (validation only; no filtering)
# =============================================================================

# The chromosomes present in an LD reference: for a panel read from a
# per-chromosome manifest these are the manifest's chromosomes; otherwise the
# chromosomes of its variants.
.ldRefChroms <- function(ldReference) {
    chromPaths <- .ldRefChromPaths(ldReference)
    if (length(chromPaths) > 0L) {
        return(names(chromPaths))
    }
    unique(.ldRefChrom(ldReference))
}

# Two-tier containment check for one entry GRanges against the LD sketch:
#   * chromosome tier (hard error): any entry chromosome absent from the sketch.
#   * variant tier: match entry variants against the sketch's variants on the
#     shared chromosomes. Zero overlap -> error; overlap fraction below
#     `minLdOverlapWarn` -> warning. Never filters.
.checkLdContainment <- function(ldReference, entryGr, minLdOverlapWarn, label) {
    if (length(entryGr) == 0L) {
        return(invisible(NULL))
    }
    reqChrom <- .ldContainmentCheckChroms(entryGr, ldReference, label)
    sketchGr <- .ldRefRanges(ldReference)
    if (length(sketchGr) == 0L) {
        msg <- glue(
            "{label}: the LD reference carries no variant metadata; ",
            "skipping the variant-overlap check."
        )
        warn(msg)
        return(invisible(NULL))
    }
    .ldContainmentCheckOverlap(
        entryGr,
        sketchGr,
        reqChrom,
        minLdOverlapWarn,
        label
    )
}

# The entry's chromosomes must all be present in the LD sketch; returns them.
# @noRd
.ldContainmentCheckChroms <- function(entryGr, ldReference, label) {
    reqChrom <- unique(as.character(GenomicRanges::seqnames(entryGr)))
    availChrom <- .ldRefChroms(ldReference)
    missingChrom <- setdiff(reqChrom, availChrom)
    if (length(missingChrom) > 0L) {
        msg <- glue(
            "{label}: chromosome(s) {str_flatten(missingChrom, ', ')} ",
            "are absent from the LD reference ",
            "(available: {str_flatten(availChrom, ', ')})."
        )
        abort(msg)
    }
    reqChrom
}

# Variant-overlap check: error on zero overlap, warn below minLdOverlapWarn.
# @noRd
.ldContainmentCheckOverlap <- function(
    entryGr,
    sketchGr,
    reqChrom,
    minLdOverlapWarn,
    label
) {
    keep <- is_in(as.character(GenomicRanges::seqnames(sketchGr)), reqChrom)
    # Same matching relation `.ldFromRef()` will use when the LD is
    # actually loaded, palindromes included -- otherwise this counts an A/T
    # variant as absent, under-reports the overlap and warns about a panel
    # that will in fact match.
    nOverlap <- length(
        matchVariants(
            .grVariantIds(entryGr),
            .grVariantIds(sketchGr[keep]),
            removeStrandAmbiguous = FALSE
        )$idxA
    )
    if (nOverlap == 0L) {
        msg <- glue(
            "{label}: none of the {length(entryGr)} summary-statistic ",
            "variants are present in the LD reference."
        )
        abort(msg)
    }
    frac <- nOverlap / length(entryGr)
    if (frac < minLdOverlapWarn) {
        pct <- sprintf("%.1f", 100 * frac)
        msg <- glue(
            "{label}: only {pct}% of summary-statistic variants ",
            "({nOverlap}/{length(entryGr)}) are present in the LD reference."
        )
        warn(msg)
    }
    invisible(NULL)
}

# =============================================================================
# Summary-statistic file readers
# =============================================================================

# Remap a region's seqnames to the seqnames actually present in a tabix file
# (tolerant of chr-prefix differences), merging windows that overlap so no line
# is read twice. A window on a chromosome the file does not hold reads nothing,
# so it is left out; NULL when no window matches. A genome-wide set of windows
# against a one-chromosome file reads that chromosome's windows.
.matchRegionSeqnames <- function(gr, available) {
    want <- as.character(GenomicRanges::seqnames(gr))
    mapped <- map_chr(want, .mlMatchSeqname, available = available)
    ok <- !is.na(mapped)
    if (!any(ok)) {
        return(NULL)
    }
    GenomicRanges::reduce(GenomicRanges::GRanges(
        seqnames = mapped[ok],
        ranges = IRanges::IRanges(
            start = GenomicRanges::start(gr)[ok],
            end = GenomicRanges::end(gr)[ok]
        )
    ))
}

# Read the region of a bgzipped + tabix-indexed delimited file into a
# data.frame, using the file's last header (#...) line for column names.
.readTabixRegion <- function(path, region) {
    tf <- Rsamtools::TabixFile(path)
    hdr <- Rsamtools::headerTabix(tf)$header
    colLine <- if (length(hdr) > 0L) {
        str_remove(hdr[[length(hdr)]], "^#")
    } else {
        NULL
    }
    gr <- .matchRegionSeqnames(region, Rsamtools::seqnamesTabix(tf))
    emptyDf <- if (!is.null(colLine)) {
        # A header-only read yields a 0-row tibble with the right columns.
        readr::read_tsv(I(colLine), show_col_types = FALSE, progress = FALSE)
    } else {
        tibble()
    }
    if (is.null(gr)) {
        return(emptyDf)
    }
    lines <- unname(list_c(Rsamtools::scanTabix(tf, param = gr)))
    if (length(lines) == 0L) {
        return(emptyDf)
    }
    txt <- str_flatten(c(colLine, lines), "\n")
    # Read every column as character; .resolveSumstatCols() coerces each field
    # to its target type. This prevents readr from type-guessing an allele/ID
    # column
    # that is uniformly T/F within the region as logical (the allele "T" ->
    # TRUE),
    # which would corrupt the variant id and break LD-containment matching.
    readr::read_tsv(
        I(txt),
        show_col_types = FALSE,
        progress = FALSE,
        col_types = readr::cols(.default = readr::col_character())
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
    # bhat/sebhat (and the betahat/sebetahat spelling) are what a cis-QTL scan
    # writes -- the pipeline's own TensorQTL nominal output uses bhat/sebhat.
    BETA = c("beta", "BETA", "bhat", "betahat"),
    SE = c("se", "SE", "sebhat", "sebetahat"),
    P = c("p", "P", "pvalue", "pval"),
    # AF: the DIRECTIONAL effect-allele frequency, exported as top_loci$af.
    # Declared only via an explicit `af`/`AF` mapping key, or a source column
    # literally named `af`. A column named effect_allele_frequency is NOT
    # auto-trusted as directional -- it must be mapped `af: <col>`.
    AF = c("af"),
    # MAF: a DIRECTIONLESS minor-allele frequency, QC only, never exported as
    # af. effect_allele_frequency stays here so an unmapped such column is
    # still used for filtering without being read as directional.
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
    # `af:` (directional effect-allele frequency) vs `maf:` (directionless,
    # QC only). A study may declare either, both, or neither.
    AF = c("af", "AF"),
    MAF = c("maf", "MAF", "effect_allele_frequency"),
    INFO = c("info", "INFO")
)

# Read a columnMapping specification into a named character vector mapping a
# standard key (chrom/pos/variant_id/...) to the source column name. Accepts a
# named list/vector, or a path to a YAML file of `standardName: sourceName`
# entries (the xqtl-protocol column-mapping format).
#' @importFrom checkmate assertFileExists
.readColumnMapping <- function(columnMapping) {
    if (is.null(columnMapping)) {
        return(NULL)
    }
    if (
        is.list(columnMapping) ||
            (is.character(columnMapping) &&
                !is.null(names(columnMapping)))
    ) {
        return(map_chr(columnMapping, as.character))
    }
    if (is.character(columnMapping) && length(columnMapping) == 1L) {
        assertFileExists(
            columnMapping,
            access = "r",
            .var.name = "columnMapping file"
        )
        mapping <- yaml::read_yaml(columnMapping)
        if (
            !is.list(mapping) ||
                is.null(names(mapping)) ||
                any(str_length(names(mapping)) == 0L)
        ) {
            msg <- glue(
                "columnMapping file '{columnMapping}' must be a YAML ",
                "mapping of standardName: sourceName."
            )
            abort(msg)
        }
        map_chr(mapping, as.character)
    } else {
        msg <- glue(
            "`columnMapping` must be a named list/vector or a path to a YAML ",
            "mapping file."
        )
        abort(msg)
    }
}

# Resolve one canonical sumstats field to its source column: honor an explicit
# columnMapping (under any accepted key spelling), else fall back to known
# column
# aliases; NA when unresolved.
# @noRd
.resolveSumstatKey <- function(
    key,
    df,
    mapping,
    label,
    claimed = character()
) {
    # Look the field up in the mapping under any accepted standard-key spelling
    # (`z`/`Z`, `n_sample`/`N`, ...). `intersect` avoids the "subscript out of
    # bounds" a named character vector throws for an absent `[[` key.
    if (!is.null(mapping)) {
        cand <- intersect(.sumstatMappingKeys[[key]], names(mapping))
        if (length(cand) >= 1L) {
            src <- mapping[[cand[[1L]]]]
            if (!is_in(src, names(df))) {
                msg <- glue(
                    "{label}: columnMapping['{cand[[1L]]}'] = '{src}' is ",
                    "not a column in the sumstats file."
                )
                abort(msg)
            }
            return(src)
        }
    }
    # Auto-detect never steals a source column another key explicitly claimed,
    # so `af: effect_allele_frequency` is not ALSO auto-read as MAF.
    intersect(
        .sumstatColumnAliases[[key]],
        setdiff(names(df), claimed)
    )[1L]
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
    unhashed <- if (ncol(df) == 0L) {
        df
    } else {
        `names<-`(
            df,
            replace(names(df), 1L, str_remove(names(df)[1L], "^#"))
        )
    }
    mapping <- .readColumnMapping(columnMapping)
    resolved <- .resolveSumstatRequired(unhashed, mapping, label)
    z <- .resolveSumstatZ(unhashed, mapping, label)
    n <- .resolveSumstatN(unhashed, mapping, label, allowNoN)
    .buildSumstatOut(unhashed, resolved, z, n, mapping, label)
}

# Resolve + validate the always-required columns (chrom/pos/variant_id/A1/A2).
# @noRd
.resolveSumstatRequired <- function(df, mapping, label) {
    required <- c("chrom", "pos", "variant_id", "A1", "A2")
    resolved <- set_names(
        map(required, .resolveSumstatKey, df, mapping, label),
        required
    )
    missingKeys <- required[map_lgl(resolved, .mlIsNaOrNull)]
    if (length(missingKeys) > 0L) {
        msg <- glue(
            "{label}: sumstats file is missing required field(s): ",
            "{str_flatten(missingKeys, ', ')} ",
            "(supply a columnMapping if the source names differ)."
        )
        abort(msg)
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
        msg <- glue(
            "{label}: sumstats file needs a z field (z/Z), or beta + se to ",
            "derive the Wald z-score beta/se (supply a columnMapping if the ",
            "source names differ)."
        )
        abort(msg)
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
        msg <- glue(
            "{label}: sumstats file needs an N field (n_sample/N/n), or ",
            "N_CASE + N_CONTROL columns for the effective sample size, or a ",
            "study-level scalar (nCase/nControl or nSample) in the manifest ",
            "(supply a columnMapping if the source names differ)."
        )
        abort(msg)
    }
    list(
        hasN = hasN,
        hasCounts = hasCounts,
        nSrc = nSrc,
        ncaseSrc = ncaseSrc,
        ncontrolSrc = ncontrolSrc
    )
}

# The source columns named by explicit mappings, so an auto-detected key never
# re-reads a column another key already claimed.
# @noRd
.sumstatClaimedCols <- function(mapping) {
    if (is.null(mapping)) {
        return(character())
    }
    unname(as.character(mapping))
}

# The identity columns every summary-statistic table carries.
# @noRd
.sumstatIdentityCols <- function(df, resolved) {
    tibble(
        chrom = as.character(df[[resolved$chrom]]),
        pos = as.integer(df[[resolved$pos]]),
        SNP = as.character(df[[resolved$variant_id]]),
        A1 = as.character(df[[resolved$A1]]),
        A2 = as.character(df[[resolved$A2]])
    )
}

# N, or the case / control pair, whichever the manifest supplied.
# @noRd
.sumstatSampleSizeCols <- function(df, n) {
    c(
        if (n$hasN) list(N = as.numeric(df[[n$nSrc]])) else list(),
        if (n$hasCounts) {
            list(
                N_CASE = as.numeric(df[[n$ncaseSrc]]),
                N_CONTROL = as.numeric(df[[n$ncontrolSrc]])
            )
        } else {
            list()
        }
    )
}

# The optional statistic columns that `df` actually carries, keyed by the
# canonical name each will be written under.
# @noRd
.sumstatOptionalSources <- function(df, mapping, label, claimed) {
    optionalKeys <- c("BETA", "SE", "P", "AF", "MAF", "INFO")
    sources <- set_names(
        map(
            optionalKeys,
            .resolveSumstatKeyFor,
            df = df,
            mapping = mapping,
            label = label,
            claimed = claimed
        ),
        optionalKeys
    )
    keep(sources, .sumstatSourceFound)
}

# Assemble the standardized sumstat data.frame (chrom/pos/SNP/A1/A2/Z, optional
# N / N_CASE / N_CONTROL, and any present BETA/SE/P/MAF/INFO).
# @noRd
.buildSumstatOut <- function(df, resolved, z, n, mapping, label) {
    claimed <- .sumstatClaimedCols(mapping)
    present <- .sumstatOptionalSources(df, mapping, label, claimed)
    result <- mutate(
        .sumstatIdentityCols(df, resolved),
        Z = if (z$hasZ) {
            as.numeric(df[[z$zSrc]])
        } else {
            as.numeric(df[[z$betaSrc]]) / as.numeric(df[[z$seSrc]])
        },
        # Whether the SOURCE supplied z, which is what GWAS-VCF's `EZ` means.
        # Recorded here because this is the only place that still knows: past
        # it, a supplied z and a Wald z computed from beta/se are both just
        # `Z`.
        z_supplied = z$hasZ,
        !!!.sumstatSampleSizeCols(df, n),
        !!!map(present, .sumstatNumericColumn, df = df)
    )
    .warnAfProvenance(
        result,
        .resolveSumstatKey("AF", df, mapping, label, claimed),
        label
    )
    result
}

# @noRd
.resolveSumstatKeyFor <- function(key, df, mapping, label, claimed) {
    .resolveSumstatKey(key, df, mapping, label, claimed)
}

# @noRd
.sumstatSourceFound <- function(src) {
    !is.null(src) && !is.na(src)
}

# @noRd
.sumstatNumericColumn <- function(src, df) {
    as.numeric(df[[src]])
}

# An exported af of NA must never be silent, and the two ways to get one are
# different problems: the study never declared a directional frequency, or it
# declared one whose values are unusable.
# @noRd
.warnAfProvenance <- function(out, afSrc, label) {
    if (is.na(afSrc) || is.null(afSrc)) {
        warn(glue(
            "{label}: no effect-allele frequency declared (map ",
            "`af: <col>` to export a directional af); top_loci$af will be ",
            "NA. A directionless `maf`/`FRQ` is used for QC only, never as af."
        ))
        return(invisible(NULL))
    }
    # An empty result has no values to be unusable: all() of none is TRUE.
    if (length(out$AF) > 0L && all(is.na(out$AF))) {
        warn(glue(
            "{label}: effect-allele frequency `af` was declared but its ",
            "values are all missing/unusable; top_loci$af will be NA."
        ))
    }
    invisible(NULL)
}

# Default GWAS-VCF FORMAT tag mapping (canonical stat -> FORMAT field).
# Candidate FORMAT keys per field, in priority order: the GWAS-VCF reserved
# key first, then spellings found in the wild.
#
# `N` used to be `"SS"` alone -- which is NOT a GWAS-VCF key; the spec's
# sample-size field is `NS`. The writer emitted `SS` too, so the two
# round-tripped with each other while neither conformed, and fixing the writer
# would have broken this reader on pecotmr's own output. `SS` stays as a
# fallback because files in the wild use it.
# @noRd
.gwasVcfFormatDefaults <- list(
    BETA = "ES",
    SE = "SE",
    # `EZ` is the z the SOURCE supplied -- the spec defines it as the z
    # "provided if it was used to derive the ES and SE fields". Reading it was
    # missing entirely, so a z-derived study round-tripped as though pecotmr
    # had computed its z from ES/SE.
    Z = "EZ",
    P = "LP",
    N = c("NS", "SS"),
    # `AF` is the spec's ALTERNATIVE-allele frequency: DIRECTIONAL, and the
    # one the posterior view publishes as `af`. It used to be read as a
    # fallback for the MAF slot and then folded with pmin(af, 1 - af), so an
    # AF of 0.8 came back as a MAF of 0.2 with no `AF` column at all --
    # pecotmr's own VCF output did not round-trip, and the direction the
    # standardized-to-raw weight conversion needs was silently gone. The
    # mirror of the released writer bug where `AF` held MAF.
    AF = "AF",
    # `EAF` stays with MAF deliberately, matching the delimited-table loader:
    # a column calling itself an effect-allele frequency is not auto-trusted
    # as directional, because the effect allele it refers to is not stated.
    MAF = c("MAF", "EAF")
)

# Pick the FORMAT sample column of a GWAS-VCF to read. Defaults to the single
# sample when there is exactly one; otherwise `sampleSelect` must name it.
.pickVcfSample <- function(samples, sampleSelect, label) {
    if (!is.null(sampleSelect)) {
        if (!is_in(sampleSelect, samples)) {
            msg <- glue(
                "{label}: sampleSelect '{sampleSelect}' is not a sample in ",
                "the VCF (have: {str_flatten(samples, ', ')})."
            )
            abort(msg)
        }
        return(sampleSelect)
    }
    if (length(samples) == 1L) {
        return(samples[[1L]])
    }
    msg <- glue(
        "{label}: the VCF has {length(samples)} samples ",
        "({str_flatten(samples, ', ')}); pass `sampleSelect` to ",
        "choose the study column."
    )
    abort(msg)
}

# Read column `tag` from the GWAS-VCF geno matrix for one sample as numeric, or
# NULL when the tag is absent/unmapped.
#
# `tag` may name SEVERAL candidate keys in priority order (the spec's key then
# the spellings found in the wild), so the first one the file actually carries
# wins. A single string still works -- it is a length-1 candidate list. The
# result is unnamed: geno matrix columns carry variant-id element names that
# the canonical sumstat tibble must not keep (data.frame stripped them).
# @noRd
.vcfGetField <- function(tag, geno, sample) {
    hit <- intersect(as.character(tag %||% character(0)), names(geno))
    if (length(hit) == 0L) {
        return(NULL)
    }
    unname(as.numeric(geno[[hit[[1L]]]][, sample]))
}

# Convert a VCF (GWAS-VCF format) into the canonical sumstats data.frame. The
# effect allele (A1) is ALT; the other allele (A2) is REF. Stats come from the
# per-study FORMAT fields ES/SE/LP/SS/EAF, with Z = ES / SE.
# @noRd
.vcfToSumstatDf <- function(vcf, sampleSelect, formatMapping, label) {
    fmap <- list_modify(
        .gwasVcfFormatDefaults,
        !!!compact(
            if (is.null(formatMapping)) list() else as.list(formatMapping)
        )
    )
    rr <- SummarizedExperiment::rowRanges(vcf)
    altList <- VariantAnnotation::alt(vcf)
    sample <- .pickVcfSample(colnames(vcf), sampleSelect, label)
    fld <- map(
        fmap[c("BETA", "SE", "Z", "P", "N", "AF", "MAF")],
        .vcfGetField,
        geno = VariantAnnotation::geno(vcf),
        sample = sample
    )
    .vcfRequireSumstatFields(fld, label)
    z <- .vcfSumstatZ(fld$BETA, fld$SE, fld$Z, length(rr))
    out <- tibble(
        chrom = as.character(GenomicRanges::seqnames(rr)),
        pos = unname(as.integer(GenomicRanges::start(rr))),
        SNP = names(rr),
        A1 = .mlFirstAlts(altList),
        A2 = unname(as.character(VariantAnnotation::ref(vcf))),
        Z = z$z,
        z_supplied = z$supplied,
        N = fld$N
    )
    mutate(out, !!!.vcfOptionalSumstatCols(fld))
}

# A GWAS-VCF must yield a z (supplied, or derivable from ES and SE) and a
# sample size. `fld` is the per-field read of `.vcfToSumstatDf()`.
# @noRd
.vcfRequireSumstatFields <- function(fld, label) {
    if ((is.null(fld$BETA) || is.null(fld$SE)) && is.null(fld$Z)) {
        msg <- glue(
            "{label}: GWAS-VCF must carry either a supplied z (EZ) or ES and ",
            "SE to derive one from (configure via `formatMapping`)."
        )
        abort(msg)
    }
    if (is.null(fld$N)) {
        msg <- glue(
            "{label}: GWAS-VCF must carry a sample-size FORMAT field ",
            "(NS, or SS in older files; configure via `formatMapping`)."
        )
        abort(msg)
    }
}

# The optional sumstat columns a GWAS-VCF may carry, NULL ones dropped. MAF is
# the file's own directionless MAF when it has one; else folded from AF, which
# is QC-usable and keeps `AF` intact alongside it rather than consuming it.
# @noRd
.vcfOptionalSumstatCols <- function(fld) {
    eaf <- fld$MAF
    af <- fld$AF
    compact(list(
        BETA = fld$BETA,
        SE = fld$SE,
        P = if (!is.null(fld$P)) 10^(-fld$P),
        AF = af,
        MAF = if (!is.null(eaf)) {
            pmin(eaf, 1 - eaf)
        } else if (!is.null(af)) {
            pmin(af, 1 - af)
        }
    ))
}

# The z and its provenance, per variant. A supplied `EZ` wins where it has a
# value; elsewhere the Wald z is derived from ES / SE. Per VARIANT rather than
# per file because a field declared in the header can be absent from
# individual records, and publishing a derived z as a supplied one would tell
# a reader that ES and SE came from it.
# @noRd
.vcfSumstatZ <- function(es, se, ez, n) {
    derived <- if (is.null(es) || is.null(se)) {
        rep(NA_real_, n)
    } else {
        es / se
    }
    if (is.null(ez)) {
        return(list(z = derived, supplied = rep(FALSE, n)))
    }
    list(
        z = if_else(is.na(ez), derived, ez),
        supplied = !is.na(ez)
    )
}

# Read one GWAS-VCF sumstats file (region-restricted only when bgzipped +
# tabix-indexed) into the canonical sumstats data.frame.
.readSumStatsVcf <- function(path, region, sampleSelect, formatMapping, label) {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        msg <- glue(
            "{label}: reading VCF sumstats requires the ",
            "'VariantAnnotation' package."
        )
        abort(msg)
    }
    hasTbi <- file.exists(str_c(path, ".tbi"))
    if (!is.null(region) && hasTbi) {
        tf <- Rsamtools::TabixFile(path)
        gr <- .matchRegionSeqnames(
            region,
            Rsamtools::seqnamesTabix(tf)
        )
        vcf <- if (is.null(gr)) {
            VariantAnnotation::readVcf(path)[0, ]
        } else {
            VariantAnnotation::readVcf(tf, param = gr)
        }
    } else {
        if (!is.null(region) && !hasTbi) {
            msg <- glue(
                "{label}: VCF '{path}' is not bgzipped + tabix-indexed, so ",
                "the whole file is read; `which` still selects the variants ",
                "kept."
            )
            warn(msg)
        }
        vcf <- VariantAnnotation::readVcf(path)
    }
    .vcfToSumstatDf(vcf, sampleSelect, formatMapping, label)
}

# ---------------------------------------------------------------------------
# Raw-table prefilters. Applied to a sumstats table as read, before
# `.resolveSumstatCols()` maps it onto the canonical schema.
# ---------------------------------------------------------------------------

# Keep only the rows belonging to one trait. A cis-QTL scan writes every gene's
# window into a single table (tensorqtl does), so without this a manifest row
# pulls in every other gene too; restricting by `which` only helps when the
# file is tabix-indexed, and cis windows of neighbouring genes overlap anyway.
# @noRd
.filterSumstatsByTrait <- function(df, traitColumn, trait, label) {
    if (is.null(traitColumn) || is.null(trait) || is.na(trait)) {
        return(df)
    }
    if (!is_in(traitColumn, names(df))) {
        cols <- str_flatten(names(df), ", ")
        msg <- glue(
            "{label}: traitColumn '{traitColumn}' is not a column of the ",
            "sumstats file (has: {cols})."
        )
        abort(msg)
    }
    keep <- !is.na(df[[traitColumn]]) &
        as.character(df[[traitColumn]]) == as.character(trait)
    if (!any(keep)) {
        msg <- glue(
            "{label}: no rows with {traitColumn} == '{trait}' in the sumstats ",
            "file."
        )
        abort(msg)
    }
    df[keep, , drop = FALSE]
}

# Populate A1/A2 from the variant id when the file carries no allele columns
# (a cis-QTL scan reports effect + se against the genotype's counted allele and
# leaves the alleles inside the id).
#
# The allele ORDER within an id is a property of whoever wrote it and cannot be
# recovered from the string: a .pvar writes REF:ALT, which is pecotmr's
# canonical A2:A1, while a PLINK .bim writes A1:A2. So the caller declares it
# rather than the reader guessing -- a wrong guess is silent, not loud, because
# summaryStatsQc() "corrects" the apparent mismatch by sign- and strand-flipping
# every variant against the LD panel.
# @noRd
.deriveAllelesFromVariantId <- function(df, order, mapping, label) {
    if (is.null(order) || identical(order, "none")) {
        return(df)
    }
    alleleCols <- map(
        c("A1", "A2"),
        .resolveSumstatKey,
        df = df,
        mapping = mapping,
        label = label
    )
    if (!all(map_lgl(alleleCols, .mlIsNaOrNull))) {
        return(df) # a real allele column always wins
    }
    idCol <- .resolveSumstatKey("variant_id", df, mapping, label)
    if (.mlIsNaOrNull(idCol)) {
        return(df) # the missing-field error is raised later
    }
    parsed <- suppressWarnings(parseVariantId(as.character(df[[idCol]])))
    bad <- sum(is.na(parsed$A1) | is.na(parsed$A2))
    if (bad > 0L) {
        n <- nrow(df)
        msg <- glue(
            "{label}: variantIdAlleles = '{order}' was requested but {bad} of ",
            "{n} variant id(s) carry no allele pair."
        )
        abort(msg)
    }
    # parseVariantId() reads an id as chr:pos:A2:A1 (pecotmr's canonical
    # order), so "A2A1" is a straight take and "A1A2" is the swap.
    if (identical(order, "A2A1")) {
        return(mutate(df, A1 = parsed$A1, A2 = parsed$A2))
    }
    mutate(df, A1 = parsed$A2, A2 = parsed$A1)
}

# @noRd
.applySumstatPrefilters <- function(df, prefilters, mapping, label) {
    if (is.null(prefilters)) {
        return(df)
    }
    df <- .filterSumstatsByTrait(
        df,
        prefilters$traitColumn,
        prefilters$trait,
        label
    )
    .deriveAllelesFromVariantId(
        df,
        prefilters$variantIdAlleles,
        mapping,
        label
    )
}

# Read one delimited-text sumstats file. Region reads only when a <path>.tbi
# sidecar exists; otherwise the whole file is read, with a warning, and the
# region selects the variants kept further on.
# @noRd
.readSumStatsText <- function(
    path,
    region,
    columnMapping,
    label,
    allowNoN = FALSE,
    prefilters = NULL
) {
    hasTbi <- file.exists(str_c(path, ".tbi"))
    raw <- if (!is.null(region) && hasTbi) {
        .readTabixRegion(path, region)
    } else {
        if (!is.null(region) && !hasTbi) {
            msg <- glue(
                "{label}: '{path}' is not tabix-indexed, so the whole file ",
                "is read; `which` still selects the variants kept."
            )
            warn(msg)
        }
        # Read every column as character; .resolveSumstatCols() coerces each
        # field to its target type. This prevents readr from type-guessing an
        # allele/ID column that is uniformly T/F as logical (the allele "T" ->
        # TRUE).
        readr::read_tsv(
            path,
            show_col_types = FALSE,
            progress = FALSE,
            col_types = readr::cols(.default = readr::col_character())
        )
    }
    filtered <- .applySumstatPrefilters(
        raw,
        prefilters,
        .readColumnMapping(columnMapping),
        label
    )
    .resolveSumstatCols(filtered, columnMapping, label, allowNoN = allowNoN)
}

# Dispatch a sumstats file to the text or GWAS-VCF reader. BCF is rejected.
.readSumStatsFile <- function(
    path,
    region,
    columnMapping,
    sampleSelect,
    formatMapping,
    label,
    allowNoN = FALSE,
    prefilters = NULL
) {
    lower <- str_to_lower(path)
    if (str_ends(lower, "\\.bcf")) {
        msg <- glue(
            "{label}: BCF sumstats are not supported; convert '{path}' to a ",
            "bgzipped, tabix-indexed VCF."
        )
        abort(msg)
    }
    if (str_detect(lower, "\\.vcf(\\.b?gz)?$")) {
        # GWAS-VCF always carries an SS (sample-size) FORMAT field, so the study
        # scalar never has to stand in for a per-variant N here.
        .readSumStatsVcf(path, region, sampleSelect, formatMapping, label)
    } else {
        .readSumStatsText(
            path,
            region,
            columnMapping,
            label,
            allowNoN = allowNoN,
            prefilters = prefilters
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
    allowNoN = FALSE,
    prefilters = NULL
) {
    df <- .readSumStatsFile(
        path,
        region,
        columnMapping,
        sampleSelect,
        formatMapping,
        label,
        allowNoN = allowNoN,
        prefilters = prefilters
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
    raw <- readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "minimal"
    )
    m <- `rownames<-`(
        `storage.mode<-`(as.matrix(select(raw, -1L)), "double"),
        as.character(pull(raw, 1L))
    )
    if (transpose) {
        return(t(m))
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

# Column names (and their case/alias variants) that are metadata, not
# expression samples, in a phenotype BED.
# @noRd
.bedReservedCols <- function() {
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
}

# Read a phenotype BED into an expression matrix (genes x samples) + rowRanges,
# their chromosomes in UCSC style whatever the file writes (1, chr1).
# @noRd
.readPhenotypeBed <- function(bedPath) {
    bed <- readr::read_tsv(
        bedPath,
        show_col_types = FALSE,
        progress = FALSE,
        comment = "",
        name_repair = "minimal"
    )
    chrCol <- intersect(c("#chr", "#chrom", "chrom", "chr"), names(bed))[1L]
    startCol <- intersect(c("start", "Start"), names(bed))[1L]
    endCol <- intersect(c("end", "End"), names(bed))[1L]
    geneCol <- intersect(c("gene_id", "ID", "phenotype_id"), names(bed))[1L]
    if (any(is.na(c(chrCol, startCol, endCol, geneCol)))) {
        msg <- glue(
            "phenotype BED is missing one of chrom/start/end/gene_id ",
            "columns: {bedPath}"
        )
        abort(msg)
    }
    meta <- `names<-`(
        select(bed, all_of(c(chrCol, startCol, endCol, geneCol))),
        c("chrom", "start", "end", "gene_id")
    )
    sampleCols <- setdiff(names(bed), .bedReservedCols())
    expr <- `rownames<-`(
        `storage.mode<-`(as.matrix(select(bed, all_of(sampleCols))), "double"),
        meta$gene_id
    )
    rr <- `names<-`(
        GenomicRanges::GRanges(
            seqnames = withChrPrefix(meta$chrom),
            ranges = IRanges::IRanges(start = meta$start + 1L, end = meta$end)
        ),
        meta$gene_id
    )
    list(expr = expr, rr = rr)
}

# Per-context colData from an optional covariate file (intersecting samples with
# the expression matrix, which is subset to the shared set). Returns list(expr,
# cd).
# @noRd
.buildContextColData <- function(covPath, transposeCov, expr, bedPath) {
    if (is.null(covPath) || str_length(covPath) == 0L) {
        return(list(
            expr = expr,
            cd = S4Vectors::DataFrame(row.names = colnames(expr))
        ))
    }
    pcov <- .readCovariateMatrix(covPath, transposeCov)
    common <- intersect(rownames(pcov), colnames(expr))
    if (length(common) == 0L) {
        msg <- glue(
            "No shared samples between phenotype and covariate file: {bedPath}"
        )
        abort(msg)
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
    studyName,
    base,
    transposeCov,
    genotypesOverride = NULL,
    genoCovOverride = NULL
) {
    contexts <- unique(as.character(rows$context))
    phenotypes <- set_names(
        map(
            contexts,
            .buildOneContextSe,
            rows = rows,
            studyName = studyName,
            base = base,
            transposeCov = transposeCov
        ),
        contexts
    )
    genotypes <- .resolveGenotypes(rows, studyName, base, genotypesOverride)
    genoCov <- .resolveGenoCov(
        rows,
        studyName,
        base,
        transposeCov,
        genoCovOverride
    )
    QtlDataset(
        studyName = studyName,
        genotypes = genotypes,
        phenotypes = phenotypes,
        genotypeCovariates = genoCov
    )
}

# The SummarizedExperiment for one context: exactly one phenotypePath (+ at most
# one covariatePath) across its manifest rows.
# @noRd
.buildOneContextSe <- function(cx, rows, studyName, base, transposeCov) {
    sub <- rows[as.character(rows$context) == cx, , drop = FALSE]
    pths <- unique(sub$phenotypePath[!is.na(sub$phenotypePath)])
    if (length(pths) != 1L) {
        msg <- glue(
            "Context '{cx}' (study '{studyName}') must reference exactly one ",
            "phenotypePath; got: {str_flatten(pths, ', ')}"
        )
        abort(msg)
    }
    covPath <- .oneContextCovariatePath(sub, cx, studyName, base)
    .buildContextSe(.resolveRel(pths[[1L]], base), covPath, transposeCov)
}

# The single covariate path a context declares, or NULL when the manifest has
# no covariatePath column or the context leaves it blank.
# @noRd
.oneContextCovariatePath <- function(sub, cx, studyName, base) {
    if (!is_in("covariatePath", names(sub))) {
        return(NULL)
    }
    covs <- unique(sub$covariatePath[
        !is.na(sub$covariatePath) &
            str_length(as.character(sub$covariatePath)) > 0L
    ])
    if (length(covs) > 1L) {
        abort(glue(
            "Context '{cx}' (study '{studyName}') references multiple ",
            "covariatePath values: {str_flatten(covs, ', ')}"
        ))
    }
    if (length(covs) == 0L) {
        return(NULL)
    }
    .resolveRel(covs[[1L]], base)
}

# The genotype panel: an override panel (or bare handle) as given, an
# override path, else the study's single genotypePath (auto-detecting the
# format).
# @noRd
.resolveGenotypes <- function(rows, studyName, base, genotypesOverride) {
    open <- .openGenotypeSource(genotypesOverride)
    if (!is.null(open)) {
        return(open)
    }
    spec <- if (is.character(genotypesOverride)) {
        genotypesOverride
    } else {
        v <- unique(as.character(rows$genotypePath[!is.na(rows$genotypePath)]))
        if (length(v) != 1L) {
            msg <- glue(
                "study '{studyName}' must reference exactly one genotypePath; ",
                "got: {str_flatten(v, ', ')}"
            )
            abort(msg)
        }
        v[[1L]]
    }
    .readGenotypePath(.resolveRel(spec, base))
}

# The genotype covariate matrix: an override matrix/path, else the study's
# single genotypeCovariatePath, else an empty matrix.
# @noRd
.resolveGenoCov <- function(
    rows,
    studyName,
    base,
    transposeCov,
    genoCovOverride
) {
    if (is.matrix(genoCovOverride)) {
        return(genoCovOverride)
    }
    spec <- if (is.character(genoCovOverride)) {
        genoCovOverride
    } else if (is_in("genotypeCovariatePath", names(rows))) {
        v <- unique(as.character(rows$genotypeCovariatePath[
            !is.na(rows$genotypeCovariatePath) &
                str_length(as.character(rows$genotypeCovariatePath)) > 0L
        ]))
        if (length(v) > 1L) {
            msg <- glue(
                "study '{studyName}' references multiple ",
                "genotypeCovariatePath ",
                "values."
            )
            abort(msg)
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
#' @param studyName Study identifier; reconciled with a \code{study} column.
#' @param genotypes A genotype panel (see \code{\link{readGenotypes}}) or
#'   a genotype path/prefix; reconciled with a \code{genotypePath} column.
#' @param genotypeCovariates A numeric matrix (samples x covariates) or a path
#'   to a covariate TSV; reconciled with a \code{genotypeCovariatePath} column.
#' @param transposeCovariates Transpose covariate TSVs (QTLtools layout) before
#'   treating them as samples-as-rows.
#' @return A \code{QtlDataset} object.
#' @examples
#' d <- system.file("extdata", "qtl_mini", package = "pecotmr")
#' manifest <- data.frame(context = "brain",
#'   phenotypePath = file.path(d, "example_geneexpr.bed.gz"),
#'   study = "s1", genotypePath = file.path(d, "example.chr22"))
#' loadQtlDatasetFromManifest(manifest = manifest, studyName = "s1")
#' @importFrom stringr str_ends
#' @importFrom checkmate assertString assertFlag
#' @export
loadQtlDatasetFromManifest <- function(
    manifest,
    studyName = NULL,
    genotypes = NULL,
    genotypeCovariates = NULL,
    transposeCovariates = FALSE
) {
    assertString(studyName, null.ok = TRUE)
    assertFlag(transposeCovariates)
    base <- .manifestBase(manifest)
    df <- .canonManifestCols(
        .readManifest(manifest),
        .qtlDatasetManifestAliases,
        required = c("context", "phenotypePath"),
        label = "QtlDataset"
    )
    # df[["study"]] (not df$study): study is optional here (it may be passed as
    # the `study` arg instead), and a tibble `$` on an absent column warns.
    studyName <- .reconcileScalar(df[["study"]], studyName, "study")
    .buildQtlDatasetFromRows(
        df,
        studyName,
        base,
        transposeCovariates,
        genotypesOverride = genotypes,
        genoCovOverride = genotypeCovariates
    )
}

# Aliases shared by the sumstats manifests: the column names a manifest may
# spell its BED-style window coordinates with -- BED's own (`chrom`,
# `chromStart`, `chromEnd`), the short `start` / `end`, and a leading `#` on
# the header's first column.
# @noRd
.BED_ALIASES <- list(
    chrom = c("chrom", "#chrom", "chr", "#chr"),
    start = c("start", "chromStart"),
    end = c("end", "chromEnd")
)

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
    ldReferencePath = c("ldReferencePath", "ld_reference_path", "ld_reference")
) |>
    c(.BED_ALIASES)

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
    ldReferencePath = c("ldReferencePath", "ld_reference_path", "ld_reference"),
    traitStart = c("traitStart", "trait_start"),
    traitEnd = c("traitEnd", "trait_end")
) |>
    c(.BED_ALIASES)

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
#'   per study. No QC is run (the result carries empty \code{qcOptions()}).
#'   Each sumstats file needs a \code{z} column, or \code{beta}+\code{se} from
#'   which the Wald z (\code{z = beta/se}) is derived when \code{z} is absent
#'   (a supplied \code{z} takes precedence).
#' @param manifest A data.frame or path. Columns (snake_case aliases accepted):
#'   \code{study} (required, unique), \code{sumStatsPath} (required),
#'   \code{columnMapping} (optional), \code{nCase} / \code{nControl} (optional),
#'   \code{nSample} (optional study-level total N), \code{varY} (optional), the
#'   single-valued \code{genome}, and \code{ldReferencePath}, which may differ
#'   between studies (each distinct reference is read once, for the
#'   chromosomes its studies cover, and stored once). When a row
#'   supplies a study-level scalar (\code{nCase} + \code{nControl}, or
#'   \code{nSample}), its sumstats file need not carry a per-variant \code{N}
#'   column;
#'   \code{\link{summaryStatsQc}} fills \code{N} from the scalar. Optional
#'   BED-style \code{chrom}/\code{start}/\code{end} columns, with a 0-based
#'   \code{start} as in BED, give each study its own window, read and kept as
#'   with \code{which}; a manifest that has them cannot also be given
#'   \code{which}.
#' @param genome Genome build; reconciled with a \code{genome} column.
#' @param ldReference A genotype panel (see \code{\link{readGenotypes}}) or a
#'   spec (path/prefix/genoMeta), used for every study; where a row also names
#'   an \code{ldReferencePath} the two must agree. Omit it to give each study
#'   the reference its \code{ldReferencePath} names.
#' @param which Optional analysis windows shared by every study, typically
#'   LD blocks: a \code{GRanges}, or the path of a BED file, read as BED
#'   with a 0-based start. Only the variants inside them are read, when the
#'   file is tabix-indexed text or a bgzipped + tabixed VCF (otherwise the
#'   whole file is read, with a warning), and each study is split into one
#'   element per window it has variants in, each window recorded as that
#'   element's \code{region}, as \code{\link{GwasSumStats}} does with its
#'   \code{region}; a named \code{GRanges} names the blocks. Without
#'   windows a genome-wide file splits
#'   into one element per chromosome; split by LD block, it has the
#'   granularity \code{\link{assembleCtwasInputs}} needs.
#' @param minLdOverlapWarn Warn when the fraction of sumstats variants present
#'   in the LD reference falls below this (default 0.5).
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
#'   ldReference = ldStem,
#'   which = GenomicRanges::GRanges("chr22:10000000-19000000"))
#' @export
loadGwasSumStatsFromManifest <- function(
    manifest,
    genome = NULL,
    ldReference = NULL,
    which = NULL,
    minLdOverlapWarn = 0.5,
    columnMapping = NULL,
    sampleSelect = NULL,
    formatMapping = NULL
) {
    base <- .manifestBase(manifest)
    df <- .gwasSumStatsManifest(manifest)
    genome <- .manifestGenome(df, genome, base, "loadGwasSumStatsFromManifest")
    ldReferenceSpecs <- .resolveLdRefSpecs(df, ldReference, base)
    ns <- .gwasStudyScalars(df)
    windows <- .manifestWindows(df, which, "GwasSumStats")
    entries <- map(
        seq_len(nrow(df)),
        .loadGwasEntry,
        df = df,
        base = base,
        windows = windows,
        columnMapping = columnMapping,
        sampleSelect = sampleSelect,
        formatMapping = formatMapping,
        ns = ns
    )
    ldReference <- .manifestRowLdRefs(
        ldReferenceSpecs,
        entries,
        df,
        minLdOverlapWarn,
        check = .checkGwasLdContainment
    )
    # Without windows a genome-wide file splits by chromosome, which is too
    # coarse for cTWAS; with them -- typically LD blocks -- each study becomes
    # one element per window. The constructor does the splitting either way.
    gwasArgs <- list_assign(
        .gwasSumStatsArgs(df, entries, genome, ldReference, ns),
        region = .gwasManifestRegion(windows)
    )
    exec(GwasSumStats, !!!gwasArgs)
}

# The windows a manifest load reads and keeps: the `which` argument, shared
# by every row, or the manifest's BED-style chrom/start/end columns, one
# window per row. Both is an error, being two answers to one question.
# @noRd
.manifestWindows <- function(df, which, label) {
    perRow <- .manifestRowWindows(df, label)
    if (!is.null(perRow) && !is.null(which)) {
        abort(glue(
            "{label} manifest gives each row a window in its ",
            "chrom/start/end columns, so `which` cannot give them too."
        ))
    }
    list(
        perRow = perRow,
        shared = if (!is.null(which)) .asRegionRanges(which, "`which`")
    )
}

# The window a manifest row reads: its own, else the shared ones.
# @noRd
.manifestWindowFor <- function(i, windows) {
    if (!is.null(windows$perRow)) windows$perRow[i] else windows$shared
}

# The GwasSumStats `region` a load builds: the shared windows, or a list with
# each study's own.
# @noRd
.gwasManifestRegion <- function(windows) {
    if (is.null(windows$perRow)) {
        return(windows$shared)
    }
    unname(as.list(GenomicRanges::split(
        windows$perRow,
        seq_along(windows$perRow)
    )))
}

# The window each manifest row gives in its chrom/start/end columns, or NULL
# when it has no such columns.
# @noRd
.manifestRowWindows <- function(df, label) {
    if (!any(is_in(c("start", "end"), names(df)))) {
        return(NULL)
    }
    .manifestRanges(df, "start", "end", label, "window")
}

# The trait position each manifest row gives in its traitStart / traitEnd
# columns, on the row's `chrom`, or NULL when it has no such columns. A trait
# sits on the chromosome its window does, so there is no second chromosome
# column.
# @noRd
.manifestTraitPos <- function(df, label) {
    if (!any(is_in(c("traitStart", "traitEnd"), names(df)))) {
        return(NULL)
    }
    .manifestRanges(df, "traitStart", "traitEnd", label, "trait position")
}

# Ranges read from a manifest's BED-style columns: `chrom` plus a start and
# an end column, the start 0-based as in BED. Every row must fill them in.
# The chromosomes come back in UCSC style, as every range read from a file
# does.
# @noRd
.manifestRanges <- function(df, startCol, endCol, label, what) {
    need <- c("chrom", startCol, endCol)
    absent <- setdiff(need, names(df))
    if (length(absent) > 0L) {
        abort(glue(
            "{label} manifest: a {what} needs the columns ",
            "{str_flatten(need, ', ')}; missing {str_flatten(absent, ', ')}."
        ))
    }
    empty <- which(
        is.na(df$chrom) | is.na(df[[startCol]]) | is.na(df[[endCol]])
    )
    if (length(empty) > 0L) {
        abort(glue(
            "{label} manifest: row(s) ",
            "{str_flatten(utils::head(empty, 5L), ', ')} ",
            "leave their {what} empty; give every row one, or none."
        ))
    }
    GenomicRanges::GRanges(
        withChrPrefix(df$chrom),
        IRanges::IRanges(
            start = as.integer(df[[startCol]]) + 1L,
            end = as.integer(df[[endCol]])
        )
    )
}

# A GwasSumStats manifest's rows, its columns under their canonical names, the
# required ones present and one row per study.
# @noRd
.gwasSumStatsManifest <- function(manifest) {
    df <- .canonManifestCols(
        .readManifest(manifest),
        .gwasSumStatsManifestAliases,
        required = c("study", "sumStatsPath"),
        label = "GwasSumStats"
    )
    if (anyDuplicated(as.character(df$study))) {
        abort("GwasSumStats manifest `study` values must be unique.")
    }
    df
}

# Study-level sample-size scalars from the manifest (n_case + n_control, or
# n_sample). When present a row need not carry a per-variant N.
# @noRd
.gwasStudyScalars <- function(df) {
    list(
        nCase = if (is_in("nCase", names(df))) as.numeric(df$nCase) else NULL,
        nControl = if (is_in("nControl", names(df))) {
            as.numeric(df$nControl)
        } else {
            NULL
        },
        nSample = if (is_in("nSample", names(df))) {
            as.numeric(df$nSample)
        } else {
            NULL
        }
    )
}

# Load one GwasSumStats manifest row's entry (resolving a per-row columnMapping
# override), allowing a missing per-variant N when a study scalar is present.
# @noRd
.loadGwasEntry <- function(
    i,
    df,
    base,
    windows,
    columnMapping,
    sampleSelect,
    formatMapping,
    ns
) {
    label <- str_c("GwasSumStats[study=", df$study[[i]], "]")
    mapping <- if (
        is_in("columnMapping", names(df)) &&
            !is.na(df$columnMapping[[i]]) &&
            str_length(as.character(df$columnMapping[[i]])) > 0L
    ) {
        .resolveRel(as.character(df$columnMapping[[i]]), base)
    } else {
        columnMapping
    }
    .loadSumStatsEntry(
        .resolveRel(as.character(df$sumStatsPath[[i]]), base),
        .manifestWindowFor(i, windows),
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
.checkGwasLdContainment <- function(
    ldReference,
    entries,
    df,
    minLdOverlapWarn
) {
    for (i in seq_len(nrow(df))) {
        .checkLdContainment(
            ldReference,
            entries[[i]],
            minLdOverlapWarn,
            str_c("GwasSumStats[study=", df$study[[i]], "]")
        )
    }
}

# GwasSumStats() constructor args, including any present study scalars + varY.
# @noRd
.gwasSumStatsArgs <- function(df, entries, genome, ldReference, ns) {
    c(
        list(
            study = as.character(df$study),
            entry = entries,
            genome = genome,
            ldReference = ldReference
        ),
        compact(list(
            nCase = ns$nCase,
            nControl = ns$nControl,
            nSample = ns$nSample,
            varY = if (is_in("varY", names(df))) as.numeric(df$varY)
        ))
    )
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
    windows,
    columnMapping,
    sampleSelect,
    formatMapping,
    allowNoN = NULL,
    traitColumn = NULL,
    variantIdAlleles = "none"
) {
    map(
        seq_len(nrow(df)),
        .loadQtlSumStatsEntry,
        df = df,
        base = base,
        windows = windows,
        columnMapping = columnMapping,
        sampleSelect = sampleSelect,
        formatMapping = formatMapping,
        allowNoN = allowNoN,
        traitColumn = traitColumn,
        variantIdAlleles = variantIdAlleles
    )
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
#'   tuple-level total N), \code{varY} (optional), the single-valued
#'   \code{genome}, and \code{ldReferencePath}, which may differ between rows
#'   (each distinct reference is read once, for the chromosomes its rows
#'   cover, and stored once). When a row supplies \code{nSample},
#'   its sumstats file need not carry a per-variant \code{N} column;
#'   \code{\link{summaryStatsQc}} fills \code{N} from the scalar. (Unlike the
#'   GWAS loader, there are no \code{nCase}/\code{nControl} columns: molecular
#'   QTL traits are quantitative.) Optional BED-style
#'   \code{chrom}/\code{start}/\code{end} columns, with a 0-based
#'   \code{start} as in BED, give each tuple its own window, read and kept as
#'   with \code{which}; a manifest that has them cannot also be given
#'   \code{which}. Optional \code{traitStart}/\code{traitEnd} columns, on
#'   the row's \code{chrom} and 0-based the same way, give the trait's
#'   position (\code{traitPos}).
#' @param genome Genome build; reconciled with a \code{genome} column.
#' @param ldReference A genotype panel (see \code{\link{readGenotypes}}) or
#'   spec, used for every row; where a row also names an
#'   \code{ldReferencePath} the two must agree. Omit it to give each row the
#'   reference its \code{ldReferencePath} names.
#' @param which Optional analysis window shared by every tuple: a
#'   \code{GRanges}, or the path of a BED file, read as BED with a 0-based
#'   start. One window only, since a tuple has one: tuples with windows
#'   of their own give them in the manifest's columns. Only the variants
#'   inside it are read, when the file is tabix-indexed (otherwise the whole
#'   file is read, with a warning), and it is each tuple's \code{region} (see
#'   \code{\link{QtlSumStats}}).
#' @param minLdOverlapWarn,columnMapping,sampleSelect,formatMapping As for
#'   \code{\link{loadGwasSumStatsFromManifest}}.
#' @param traitColumn Optional column of the sumstats file naming the trait of
#'   each row (e.g. \code{"molecular_trait_id"} in a cis-QTL scan, which writes
#'   every gene into one table). When set, each manifest row keeps only the
#'   rows whose \code{traitColumn} equals that row's \code{trait}. Without it
#'   a multi-trait file loads whole, since \code{which} restriction needs a
#'   tabix index and neighbouring cis windows overlap regardless.
#' @param variantIdAlleles Where to find the alleles when the file has no
#'   \code{A1}/\code{A2} columns -- a cis-QTL scan reports effect + se against
#'   the counted allele and leaves the alleles inside the variant id. The order
#'   within an id cannot be recovered from the string (a \code{.pvar} writes
#'   \code{REF:ALT}, i.e. pecotmr's canonical \code{A2:A1}; a PLINK
#'   \code{.bim} writes \code{A1:A2}), so declare it: \code{"none"} (default,
#'   require real columns), \code{"A2A1"}, or \code{"A1A2"}. An explicit
#'   allele column always wins. Declaring the wrong order fails silently --
#'   \code{\link{summaryStatsQc}} will sign- and strand-flip every variant to
#'   "correct" the apparent mismatch against the panel.
#' @return A \code{QtlSumStats} object.
#' @examples
#' tsv <- system.file("extdata", "manifests",
#'   "protocol_example.twas.gwas_sumstats.chr22.tsv.gz", package = "pecotmr")
#' ldStem <- file.path(system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr"), "protocol_example.LD.chr22")
#' manifest <- data.frame(study = "s1", context = "brain",
#'   trait = "ENSG1", sumStatsPath = tsv)
#' loadQtlSumStatsFromManifest(manifest = manifest, genome = "hg38",
#'   ldReference = ldStem,
#'   which = GenomicRanges::GRanges("chr22:10000000-19000000"))
#' @export
loadQtlSumStatsFromManifest <- function(
    manifest,
    genome = NULL,
    ldReference = NULL,
    which = NULL,
    minLdOverlapWarn = 0.5,
    columnMapping = NULL,
    sampleSelect = NULL,
    formatMapping = NULL,
    traitColumn = NULL,
    variantIdAlleles = c("none", "A2A1", "A1A2")
) {
    variantIdAlleles <- arg_match(variantIdAlleles)
    base <- .manifestBase(manifest)
    df <- .qtlSumStatsManifest(manifest)
    genome <- .manifestGenome(df, genome, base, "loadQtlSumStatsFromManifest")
    ldReferenceSpecs <- .resolveLdRefSpecs(df, ldReference, base)
    nSample <- .qtlManifestNSample(df)
    windows <- .qtlManifestWindows(df, which)
    entries <- .loadQtlSumStatsEntries(
        df,
        base,
        windows,
        columnMapping,
        sampleSelect,
        formatMapping,
        allowNoN = nSample$allowNoN,
        traitColumn = traitColumn,
        variantIdAlleles = variantIdAlleles
    )
    ldReference <- .manifestRowLdRefs(
        ldReferenceSpecs,
        entries,
        df,
        minLdOverlapWarn,
        check = .qtlSumStatsCheckContainment
    )
    exec(
        QtlSumStats,
        !!!.qtlSumStatsArgs(
            df,
            entries,
            genome,
            ldReference,
            nSample$nSampleCol
        ),
        region = windows$perRow %||% windows$shared,
        traitPos = .manifestTraitPos(df, "QtlSumStats")
    )
}

# A QTL load's windows (.manifestWindows()). A tuple has one window, so the
# `which` argument, which every tuple shares, must be one; tuples with
# windows of their own give them in the manifest's chrom/start/end columns.
# @noRd
.qtlManifestWindows <- function(df, which) {
    windows <- .manifestWindows(df, which, "QtlSumStats")
    if (length(windows$shared) > 1L) {
        abort(glue(
            "loadQtlSumStatsFromManifest(): `which` must be one window, ",
            "which every tuple shares (got {length(windows$shared)}). A ",
            "tuple with a window of its own gives it in the manifest's ",
            "chrom/start/end columns."
        ))
    }
    windows
}

# A QtlSumStats manifest's rows, its columns under their canonical names and
# the required ones present.
# @noRd
.qtlSumStatsManifest <- function(manifest) {
    .canonManifestCols(
        .readManifest(manifest),
        .qtlSumStatsManifestAliases,
        required = c("study", "context", "trait", "sumStatsPath"),
        label = "QtlSumStats"
    )
}

# The tuple-level total-N scalar a manifest may carry, plus the per-row
# "this tuple needs no per-variant N" flags derived from it: where a row has a
# usable nSample, summaryStatsQc fills N from it, so the entry reader's N check
# is relaxed for that row.
# @noRd
.qtlManifestNSample <- function(df) {
    nSampleCol <- if (is_in("nSample", names(df))) {
        as.numeric(df$nSample)
    } else {
        NULL
    }
    list(
        nSampleCol = nSampleCol,
        allowNoN = if (!is.null(nSampleCol)) is.finite(nSampleCol) else NULL
    )
}

# Per-row LD-containment check (no-op when the sketch is NULL).
# @noRd
.qtlSumStatsCheckContainment <- function(
    ldReference,
    entries,
    df,
    minLdOverlapWarn
) {
    if (is.null(ldReference)) {
        return(invisible(NULL))
    }
    for (i in seq_len(nrow(df))) {
        .checkLdContainment(
            ldReference,
            entries[[i]],
            minLdOverlapWarn,
            str_c(
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
.qtlSumStatsArgs <- function(df, entries, genome, ldReference, nSampleCol) {
    args <- list(
        study = as.character(df$study),
        context = as.character(df$context),
        trait = as.character(df$trait),
        entry = entries,
        genome = genome,
        ldReference = ldReference
    )
    c(
        args,
        compact(list(
            nSample = nSampleCol,
            varY = if (is_in("varY", names(df))) as.numeric(df$varY)
        ))
    )
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
#' @param genome,ldReference,which,minLdOverlapWarn Passed to the QtlSumStats
#'   loader for \code{sumStatsManifest}.
#' @param columnMapping,sampleSelect,formatMapping Passed to the QtlSumStats
#'   loader for \code{sumStatsManifest}.
#' @param transposeCovariates Transpose covariate TSVs (QTLtools layout).
#' @return A \code{MultiStudyQtlDataset} object.
#' @examples
#' d <- system.file("extdata", "qtl_mini", package = "pecotmr")
#' manifest <- data.frame(study = c("s1", "s2"), context = "brain",
#'   phenotypePath = file.path(d, "example_geneexpr.bed.gz"),
#'   genotypePath = file.path(d, "example.chr22"))
#' loadMultiStudyQtlDatasetFromManifest(qtlDatasetsManifest = manifest)
#' @importFrom checkmate assertFlag assertNumber
#' @export
loadMultiStudyQtlDatasetFromManifest <- function(
    qtlDatasetsManifest,
    sumStatsManifest = NULL,
    genome = NULL,
    ldReference = NULL,
    which = NULL,
    minLdOverlapWarn = 0.5,
    columnMapping = NULL,
    sampleSelect = NULL,
    formatMapping = NULL,
    transposeCovariates = FALSE
) {
    assertFlag(transposeCovariates)
    assertNumber(minLdOverlapWarn, lower = 0, upper = 1)
    qtlDatasets <- .msqBuildDatasets(qtlDatasetsManifest, transposeCovariates)
    sumStats <- .msqLoadSumStats(
        sumStatsManifest,
        genome,
        ldReference,
        which,
        minLdOverlapWarn,
        columnMapping,
        sampleSelect,
        formatMapping
    )
    MultiStudyQtlDataset(qtlDatasets = qtlDatasets, sumStats = sumStats)
}

# Read the qtlDatasets manifest and build one QtlDataset per study.
# @noRd
.msqBuildDatasets <- function(qtlDatasetsManifest, transposeCovariates) {
    base <- .manifestBase(qtlDatasetsManifest)
    df <- .canonManifestCols(
        .readManifest(qtlDatasetsManifest),
        .qtlDatasetManifestAliases,
        required = c("study", "genotypePath", "context", "phenotypePath"),
        label = "MultiStudyQtlDataset (qtlDatasetsManifest)"
    )
    studies <- unique(as.character(df$study))
    set_names(
        map(
            studies,
            .msqDatasetForStudy,
            df = df,
            base = base,
            transposeCovariates = transposeCovariates
        ),
        studies
    )
}

# Load the optional companion sumStats manifest (NULL when none supplied).
# @noRd
.msqLoadSumStats <- function(
    sumStatsManifest,
    genome,
    ldReference,
    which,
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
        ldReference = ldReference,
        which = which,
        minLdOverlapWarn = minLdOverlapWarn,
        columnMapping = columnMapping,
        sampleSelect = sampleSelect,
        formatMapping = formatMapping
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One sumstats entry's chromosomes (empty for NULL / 0-length gr).
# @noRd
.mlEntryChroms <- function(gr) {
    if (is.null(gr) || length(gr) == 0L) {
        return(character(0))
    }
    as.character(GenomicRanges::seqnames(gr))
}

# The tabix seqname matching query seqname `w` (chr-prefix tolerant); NA if
# none.
# @noRd
.mlMatchSeqname <- function(w, available) {
    hit <- available[canonChrom(available) == canonChrom(w)]
    if (length(hit) >= 1L) hit[[1L]] else NA_character_
}

# TRUE when a resolved field value is NA or NULL (missing-key check).
# @noRd
.mlIsNaOrNull <- function(x) {
    is.na(x) || is.null(x)
}

# Each VCF row's first ALT allele, as character: one unlist() of the alt sets
# and an index into it, not an as.list() walk that pulls the rows out one at
# a time (about two seconds on three thousand variants). A row with no ALT
# allele gives NA.
# @noRd
.mlFirstAlts <- function(altList) {
    n <- lengths(altList)
    flat <- as.character(unlist(altList, use.names = FALSE))
    first <- cumsum(n) - n + 1L
    out <- rep(NA_character_, length(n))
    out[n > 0L] <- flat[first[n > 0L]]
    out
}

# One QtlSumStats entry for manifest row `i` (per-row columnMapping override).
# @noRd
.loadQtlSumStatsEntry <- function(
    i,
    df,
    base,
    windows,
    columnMapping,
    sampleSelect,
    formatMapping,
    allowNoN,
    traitColumn = NULL,
    variantIdAlleles = "none"
) {
    label <- str_c(
        "QtlSumStats[",
        df$study[[i]],
        "/",
        df$context[[i]],
        "/",
        df$trait[[i]],
        "]"
    )
    mapping <- if (
        is_in("columnMapping", names(df)) &&
            !is.na(df$columnMapping[[i]]) &&
            str_length(as.character(df$columnMapping[[i]])) > 0L
    ) {
        .resolveRel(as.character(df$columnMapping[[i]]), base)
    } else {
        columnMapping
    }
    .loadSumStatsEntry(
        .resolveRel(as.character(df$sumStatsPath[[i]]), base),
        .manifestWindowFor(i, windows),
        mapping,
        sampleSelect,
        formatMapping,
        label,
        allowNoN = !is.null(allowNoN) && isTRUE(allowNoN[[i]]),
        prefilters = list(
            traitColumn = traitColumn,
            trait = df$trait[[i]],
            variantIdAlleles = variantIdAlleles
        )
    )
}

# One study's QtlDataset from its manifest rows.
# @noRd
.msqDatasetForStudy <- function(s, df, base, transposeCovariates) {
    .buildQtlDatasetFromRows(
        df[as.character(df$study) == s, , drop = FALSE],
        s,
        base,
        transposeCovariates
    )
}
