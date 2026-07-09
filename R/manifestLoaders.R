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
    return(as.data.frame(manifest, stringsAsFactors = FALSE, check.names = FALSE))
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
  if (is.character(manifest) && length(manifest) == 1L && file.exists(manifest)) {
    dirname(normalizePath(manifest))
  } else {
    NULL
  }
}

# Resolve a possibly-relative path against a manifest base directory.
.resolveRel <- function(p, base) {
  if (is.null(p) || length(p) != 1L || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^(/|[A-Za-z]:|~)", p) || is.null(base) || file.exists(p)) return(p)
  file.path(base, p)
}

# Rename accepted aliases to the canonical camelCase name and validate that all
# `required` canonical columns are present. `aliases` is a named list keyed by
# canonical name; each value is the character vector of accepted source names
# (including the canonical name itself). The first alias present wins.
.canonManifestCols <- function(df, aliases, required, label) {
  for (canon in names(aliases)) {
    if (canon %in% names(df)) next
    hit <- intersect(aliases[[canon]], names(df))
    if (length(hit) >= 1L) {
      names(df)[match(hit[[1L]], names(df))] <- canon
    }
  }
  missingCols <- setdiff(required, names(df))
  if (length(missingCols) > 0L) {
    stop(label, " manifest is missing required column(s): ",
         paste(missingCols, collapse = ", "))
  }
  df
}

# Resolve a single collection-level scalar (e.g. genome) that may be supplied
# both as a constant manifest column and as a function argument. The column (if
# present) must be constant; a supplied argument must agree with it.
.reconcileScalar <- function(colValues, arg, name, required = TRUE) {
  colVal <- NULL
  if (!is.null(colValues)) {
    v <- unique(as.character(colValues[!is.na(colValues) &
                                       nzchar(as.character(colValues))]))
    if (length(v) > 1L) {
      stop("`", name, "` column must be constant across the manifest; got: ",
           paste(v, collapse = ", "))
    }
    if (length(v) == 1L) colVal <- v
  }
  if (!is.null(arg) && !is.null(colVal)) {
    if (!identical(as.character(arg), colVal)) {
      stop("`", name, "` argument ('", arg, "') disagrees with the manifest ",
           "column value ('", colVal, "').")
    }
    return(arg)
  }
  if (!is.null(arg)) return(arg)
  if (!is.null(colVal)) return(colVal)
  if (required) {
    stop("`", name, "` must be provided as an argument or a manifest column.")
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
  if (grepl("\\.vcf(\\.b?gz)?$", lower) || endsWith(lower, ".bcf") ||
      endsWith(lower, ".gds")) {
    return(GenotypeHandle(path = path))
  }
  if (endsWith(lower, ".bed")) {
    return(GenotypeHandle(plink1Prefix = sub("\\.bed$", "", path, ignore.case = TRUE)))
  }
  if (endsWith(lower, ".pgen")) {
    return(GenotypeHandle(plink2Prefix = sub("\\.pgen$", "", path, ignore.case = TRUE)))
  }
  if (file.exists(paste0(path, ".bed"))) return(GenotypeHandle(plink1Prefix = path))
  if (file.exists(paste0(path, ".pgen"))) return(GenotypeHandle(plink2Prefix = path))
  stop("Could not determine genotype format for '", path,
       "' (expected .bed/.pgen/.vcf[.gz]/.bcf/.gds or a PLINK prefix).")
}

# Turn an ldSketch specification into a GenotypeHandle. Accepts a prebuilt
# handle, a named chrom -> path vector (genoMeta), or a single string that is
# either a genotype file/prefix or a chrom-sharded genoMeta meta file.
.resolveLdSketch <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (methods::is(spec, "GenotypeHandle")) return(spec)
  if (is.character(spec) && !is.null(names(spec))) {
    return(GenotypeHandle(genoMeta = spec))
  }
  if (is.character(spec) && length(spec) == 1L) {
    lower <- tolower(spec)
    looksLikeGenotype <- grepl("\\.vcf(\\.b?gz)?$", lower) ||
      endsWith(lower, ".bcf") || endsWith(lower, ".gds") ||
      endsWith(lower, ".bed") || endsWith(lower, ".pgen") ||
      file.exists(paste0(spec, ".bed")) || file.exists(paste0(spec, ".pgen"))
    if (looksLikeGenotype) return(.detectGenotypeFormat(spec))
    # Otherwise treat it as a chrom-sharded genoMeta meta file.
    return(GenotypeHandle(genoMeta = spec))
  }
  stop("`ldSketch` must be a GenotypeHandle, a genotype path/prefix, or a ",
       "genoMeta spec (named chrom->path vector or meta-file path).")
}

# Resolve the LD sketch from the (argument, ldSketchPath column) pair.
.resolveLdSketchInput <- function(df, ldSketch, base) {
  if (methods::is(ldSketch, "GenotypeHandle")) return(ldSketch)
  colSpec <- NULL
  if (!is.null(df$ldSketchPath)) {
    v <- unique(as.character(df$ldSketchPath[!is.na(df$ldSketchPath) &
                                             nzchar(as.character(df$ldSketchPath))]))
    if (length(v) > 1L) stop("`ldSketchPath` column must be constant across the manifest.")
    if (length(v) == 1L) colSpec <- .resolveRel(v, base)
  }
  spec <- ldSketch
  if (!is.null(spec) && !is.null(colSpec) &&
      !identical(as.character(spec), as.character(colSpec))) {
    stop("`ldSketch` argument disagrees with the manifest `ldSketchPath` column.")
  }
  spec <- spec %||% colSpec
  if (is.null(spec)) {
    stop("`ldSketch` must be provided as an argument or an `ldSketchPath` column.")
  }
  .resolveLdSketch(spec)
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
  if (length(entryGr) == 0L) return(invisible(NULL))
  reqChrom <- unique(canonChrom(as.character(GenomicRanges::seqnames(entryGr))))
  availChrom <- .ldSketchChroms(ldSketch)
  missingChrom <- setdiff(reqChrom, availChrom)
  if (length(missingChrom) > 0L) {
    stop(label, ": chromosome(s) ", paste(missingChrom, collapse = ", "),
         " are absent from the LD sketch (available: ",
         paste(availChrom, collapse = ", "), ").")
  }
  snpInfo <- ldSketch@snpInfo
  if (is.null(snpInfo) || nrow(snpInfo) == 0L) {
    warning(label, ": LD sketch carries no variant metadata; ",
            "skipping the variant-overlap check.")
    return(invisible(NULL))
  }
  keepSketch <- canonChrom(as.character(snpInfo$CHR)) %in% reqChrom
  snpInfo <- snpInfo[keepSketch, , drop = FALSE]
  entryVid <- formatVariantId(
    chrom = canonChrom(as.character(GenomicRanges::seqnames(entryGr))),
    pos   = GenomicRanges::start(entryGr),
    A2    = as.character(entryGr$A2),
    A1    = as.character(entryGr$A1))
  sketchVid <- formatVariantId(
    chrom = canonChrom(as.character(snpInfo$CHR)),
    pos   = as.integer(snpInfo$BP),
    A2    = as.character(snpInfo$A2),
    A1    = as.character(snpInfo$A1))
  km <- matchVariants(entryVid, sketchVid)
  nOverlap <- length(km$idxA)
  if (nOverlap == 0L) {
    stop(label, ": none of the ", length(entryGr), " summary-statistic ",
         "variants are present in the LD sketch.")
  }
  frac <- nOverlap / length(entryGr)
  if (frac < minLdOverlapWarn) {
    warning(sprintf(
      "%s: only %.1f%% of summary-statistic variants (%d/%d) are present in the LD sketch.",
      label, 100 * frac, nOverlap, length(entryGr)))
  }
  invisible(NULL)
}

# =============================================================================
# Summary-statistic file readers
# =============================================================================

# Coerce a region argument (chr:start-end string, GRanges, or one-row
# data.frame with chrom/start/end) into a GRanges.
.asGRegion <- function(region) {
  if (methods::is(region, "GRanges")) return(region)
  if (is.data.frame(region)) {
    return(GenomicRanges::GRanges(
      seqnames = as.character(region$chrom),
      ranges = IRanges::IRanges(start = as.integer(region$start),
                                end = as.integer(region$end))))
  }
  if (is.character(region) && length(region) == 1L) {
    m <- regmatches(region, regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", region))[[1L]]
    if (length(m) != 4L) {
      stop("`region` must be a 'chr:start-end' string, a GRanges, or a ",
           "one-row data.frame with chrom/start/end.")
    }
    return(GenomicRanges::GRanges(
      seqnames = m[[2L]],
      ranges = IRanges::IRanges(start = as.integer(gsub(",", "", m[[3L]])),
                                end = as.integer(gsub(",", "", m[[4L]])))))
  }
  stop("`region` must be a 'chr:start-end' string, a GRanges, or a ",
       "one-row data.frame with chrom/start/end.")
}

# Remap a region's seqnames to the seqnames actually present in a tabix file
# (tolerant of chr-prefix differences). Returns NULL when no seqname matches.
.matchRegionSeqnames <- function(gr, available) {
  want <- as.character(GenomicRanges::seqnames(gr))
  mapped <- vapply(want, function(w) {
    hit <- available[canonChrom(available) == canonChrom(w)]
    if (length(hit) >= 1L) hit[[1L]] else NA_character_
  }, character(1))
  if (anyNA(mapped)) return(NULL)
  GenomicRanges::GRanges(
    seqnames = mapped,
    ranges = IRanges::IRanges(start = GenomicRanges::start(gr),
                              end = GenomicRanges::end(gr)))
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
    as.data.frame(readr::read_tsv(I(colLine), show_col_types = FALSE,
                                  progress = FALSE), check.names = FALSE)[0, , drop = FALSE]
  } else {
    data.frame()
  }
  if (is.null(gr)) return(emptyDf)
  lines <- unlist(Rsamtools::scanTabix(tf, param = gr), use.names = FALSE)
  if (length(lines) == 0L) return(emptyDf)
  txt <- paste(c(colLine, lines), collapse = "\n")
  # Read every column as character; .resolveSumstatCols() coerces each field to
  # its target type. This prevents readr from type-guessing an allele/ID column
  # that is uniformly T/F within the region as logical (the allele "T" -> TRUE),
  # which would corrupt the variant id and break LD-containment matching.
  as.data.frame(readr::read_tsv(I(txt), show_col_types = FALSE, progress = FALSE,
                                col_types = readr::cols(.default = readr::col_character())),
                check.names = FALSE)
}

# Default source-column aliases used when no explicit columnMapping is given
# (mirrors gwas_sumstats_construct.R).
.sumstatColumnAliases <- list(
  chrom      = c("chrom", "#chrom", "chr", "#chr"),
  pos        = c("pos", "position", "BP", "bp"),
  variant_id = c("variant_id", "SNP", "snp", "rsid"),
  A1         = c("A1", "a1"),
  A2         = c("A2", "a2"),
  Z          = c("z", "Z"),
  N          = c("n_sample", "N", "n"),
  BETA       = c("beta", "BETA"),
  SE         = c("se", "SE"),
  P          = c("p", "P", "pvalue", "pval"),
  MAF        = c("maf", "MAF", "effect_allele_frequency"),
  INFO       = c("info", "INFO"))

# Accepted columnMapping *key* spellings per canonical field. A mapping file
# uses the xqtl-protocol standard-key names (e.g. `z`, `n_sample`), which
# differ from the canonical mcol names (`Z`, `N`); this bridges the two.
.sumstatMappingKeys <- list(
  chrom      = c("chrom", "chr"),
  pos        = c("pos", "position"),
  variant_id = c("variant_id", "snp", "SNP"),
  A1         = c("A1", "a1"),
  A2         = c("A2", "a2"),
  Z          = c("z", "Z"),
  N          = c("n_sample", "N", "n"),
  BETA       = c("beta", "BETA"),
  SE         = c("se", "SE"),
  P          = c("p", "P", "pvalue"),
  MAF        = c("maf", "MAF", "effect_allele_frequency"),
  INFO       = c("info", "INFO"))

# Read a columnMapping specification into a named character vector mapping a
# standard key (chrom/pos/variant_id/...) to the source column name. Accepts a
# named list/vector, or a path to a YAML file of `standardName: sourceName`
# entries (the xqtl-protocol column-mapping format).
.readColumnMapping <- function(columnMapping) {
  if (is.null(columnMapping)) return(NULL)
  if (is.list(columnMapping) || (is.character(columnMapping) &&
                                 !is.null(names(columnMapping)))) {
    return(vapply(columnMapping, as.character, character(1)))
  }
  if (is.character(columnMapping) && length(columnMapping) == 1L) {
    if (!file.exists(columnMapping)) {
      stop("columnMapping file not found: ", columnMapping)
    }
    mapping <- yaml::read_yaml(columnMapping)
    if (!is.list(mapping) || is.null(names(mapping)) ||
        any(!nzchar(names(mapping)))) {
      stop("columnMapping file '", columnMapping,
           "' must be a YAML mapping of standardName: sourceName.")
    }
    vapply(mapping, as.character, character(1))
  } else {
    stop("`columnMapping` must be a named list/vector or a path to a YAML ",
         "mapping file.")
  }
}

# Standardise the columns of a raw sumstats data.frame into the canonical
# schema expected by .dfToEntryGranges: chrom, pos, SNP, A1, A2, Z, N (+
# optional BETA, SE, P, MAF, INFO). `mapping` (named std -> source) overrides
# the default aliases per key.
.resolveSumstatCols <- function(df, columnMapping, label) {
  # Strip a leading '#' from the first column name so a '#CHR'-style header
  # (common in GWAS TSVs) resolves the same way on the plain-text and tabix
  # read paths, and so explicit columnMapping values match the bare name.
  if (ncol(df) > 0L) names(df)[1L] <- sub("^#", "", names(df)[1L])
  mapping <- .readColumnMapping(columnMapping)
  resolveKey <- function(key) {
    # Look the field up in the mapping under any accepted standard-key spelling
    # (`z`/`Z`, `n_sample`/`N`, ...). `intersect` avoids the "subscript out of
    # bounds" a named character vector throws for an absent `[[` key.
    if (!is.null(mapping)) {
      cand <- intersect(.sumstatMappingKeys[[key]], names(mapping))
      if (length(cand) >= 1L) {
        src <- mapping[[cand[[1L]]]]
        if (!(src %in% names(df))) {
          stop(label, ": columnMapping['", cand[[1L]], "'] = '", src,
               "' is not a column in the sumstats file.")
        }
        return(src)
      }
    }
    intersect(.sumstatColumnAliases[[key]], names(df))[1L]
  }
  required <- c("chrom", "pos", "variant_id", "A1", "A2", "Z", "N")
  resolved <- setNames(lapply(required, resolveKey), required)
  missingKeys <- required[vapply(resolved, function(x) is.na(x) || is.null(x), logical(1))]
  if (length(missingKeys) > 0L) {
    stop(label, ": sumstats file is missing required field(s): ",
         paste(missingKeys, collapse = ", "),
         " (supply a columnMapping if the source names differ).")
  }
  out <- data.frame(
    chrom = as.character(df[[resolved$chrom]]),
    pos   = as.integer(df[[resolved$pos]]),
    SNP   = as.character(df[[resolved$variant_id]]),
    A1    = as.character(df[[resolved$A1]]),
    A2    = as.character(df[[resolved$A2]]),
    Z     = as.numeric(df[[resolved$Z]]),
    N     = as.numeric(df[[resolved$N]]),
    stringsAsFactors = FALSE)
  for (key in c("BETA", "SE", "P", "MAF", "INFO")) {
    src <- resolveKey(key)
    if (!is.na(src) && !is.null(src)) out[[key]] <- as.numeric(df[[src]])
  }
  out
}

# Default GWAS-VCF FORMAT tag mapping (canonical stat -> FORMAT field).
.gwasVcfFormatDefaults <- list(BETA = "ES", SE = "SE", P = "LP",
                               N = "SS", MAF = "EAF")

# Pick the FORMAT sample column of a GWAS-VCF to read. Defaults to the single
# sample when there is exactly one; otherwise `sampleSelect` must name it.
.pickVcfSample <- function(samples, sampleSelect, label) {
  if (!is.null(sampleSelect)) {
    if (!(sampleSelect %in% samples)) {
      stop(label, ": sampleSelect '", sampleSelect,
           "' is not a sample in the VCF (have: ",
           paste(samples, collapse = ", "), ").")
    }
    return(sampleSelect)
  }
  if (length(samples) == 1L) return(samples[[1L]])
  stop(label, ": the VCF has ", length(samples), " samples (",
       paste(samples, collapse = ", "),
       "); pass `sampleSelect` to choose the study column.")
}

# Convert a VCF (GWAS-VCF format) into the canonical sumstats data.frame. The
# effect allele (A1) is ALT; the other allele (A2) is REF. Stats come from the
# per-study FORMAT fields ES/SE/LP/SS/EAF, with Z = ES / SE.
.vcfToSumstatDf <- function(vcf, sampleSelect, formatMapping, label) {
  fmap <- modifyList(.gwasVcfFormatDefaults,
                     if (is.null(formatMapping)) list() else as.list(formatMapping))
  rr <- SummarizedExperiment::rowRanges(vcf)
  altList <- VariantAnnotation::alt(vcf)
  alt <- vapply(as.list(altList), function(a) as.character(a)[[1L]], character(1))
  ref <- as.character(VariantAnnotation::ref(vcf))
  geno <- VariantAnnotation::geno(vcf)
  sample <- .pickVcfSample(colnames(vcf), sampleSelect, label)
  getField <- function(tag) {
    if (!is.null(tag) && tag %in% names(geno)) {
      as.numeric(geno[[tag]][, sample])
    } else {
      NULL
    }
  }
  es  <- getField(fmap$BETA)
  se  <- getField(fmap$SE)
  lp  <- getField(fmap$P)
  ss  <- getField(fmap$N)
  eaf <- getField(fmap$MAF)
  if (is.null(es) || is.null(se)) {
    stop(label, ": GWAS-VCF must carry ES and SE FORMAT fields to derive Z ",
         "(configure via `formatMapping`).")
  }
  if (is.null(ss)) {
    stop(label, ": GWAS-VCF must carry an SS (sample size) FORMAT field.")
  }
  out <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(rr)),
    pos   = as.integer(GenomicRanges::start(rr)),
    SNP   = names(rr),
    A1    = alt,
    A2    = ref,
    Z     = es / se,
    N     = ss,
    BETA  = es,
    SE    = se,
    stringsAsFactors = FALSE)
  if (!is.null(lp))  out$P   <- 10^(-lp)
  if (!is.null(eaf)) out$MAF <- pmin(eaf, 1 - eaf)
  out
}

# Read one GWAS-VCF sumstats file (region-restricted only when bgzipped +
# tabix-indexed) into the canonical sumstats data.frame.
.readSumStatsVcf <- function(path, region, sampleSelect, formatMapping, label) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop(label, ": reading VCF sumstats requires the 'VariantAnnotation' ",
         "package; please install it.")
  }
  hasTbi <- file.exists(paste0(path, ".tbi"))
  if (!is.null(region) && hasTbi) {
    tf <- Rsamtools::TabixFile(path)
    gr <- .matchRegionSeqnames(.asGRegion(region), Rsamtools::seqnamesTabix(tf))
    vcf <- if (is.null(gr)) {
      VariantAnnotation::readVcf(path)[0, ]
    } else {
      VariantAnnotation::readVcf(tf, param = gr)
    }
  } else {
    if (!is.null(region) && !hasTbi) {
      warning(label, ": VCF '", path, "' is not bgzipped + tabix-indexed; ",
              "region ignored, reading the whole file.")
    }
    vcf <- VariantAnnotation::readVcf(path)
  }
  .vcfToSumstatDf(vcf, sampleSelect, formatMapping, label)
}

# Read one delimited-text sumstats file. Region reads only when a <path>.tbi
# sidecar exists; otherwise the whole file is read and any region is ignored.
.readSumStatsText <- function(path, region, columnMapping, label) {
  hasTbi <- file.exists(paste0(path, ".tbi"))
  raw <- if (!is.null(region) && hasTbi) {
    .readTabixRegion(path, region)
  } else {
    if (!is.null(region) && !hasTbi) {
      warning(label, ": '", path, "' is not tabix-indexed; ",
              "region ignored, reading the whole file.")
    }
    # Read every column as character; .resolveSumstatCols() coerces each field
    # to its target type. This prevents readr from type-guessing an allele/ID
    # column that is uniformly T/F as logical (the allele "T" -> TRUE).
    as.data.frame(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
                                  col_types = readr::cols(.default = readr::col_character())),
                  check.names = FALSE)
  }
  .resolveSumstatCols(raw, columnMapping, label)
}

# Dispatch a sumstats file to the text or GWAS-VCF reader. BCF is rejected.
.readSumStatsFile <- function(path, region, columnMapping, sampleSelect,
                              formatMapping, label) {
  lower <- tolower(path)
  if (endsWith(lower, ".bcf")) {
    stop(label, ": BCF sumstats are not supported; convert '", path,
         "' to a bgzipped, tabix-indexed VCF.")
  }
  if (grepl("\\.vcf(\\.b?gz)?$", lower)) {
    .readSumStatsVcf(path, region, sampleSelect, formatMapping, label)
  } else {
    .readSumStatsText(path, region, columnMapping, label)
  }
}

# Read one sumstats source into an entry GRanges (canonical mcols).
.loadSumStatsEntry <- function(path, region, columnMapping, sampleSelect,
                               formatMapping, label) {
  df <- .readSumStatsFile(path, region, columnMapping, sampleSelect,
                          formatMapping, label)
  .dfToEntryGranges(df)
}

# =============================================================================
# Phenotype SummarizedExperiment builder (individual-level QTL)
# =============================================================================

# Read a covariate TSV into a samples-as-rows numeric matrix. The first column
# holds row identifiers; `transpose = TRUE` handles QTLtools-format inputs
# (covariates as rows, samples as columns).
.readCovariateMatrix <- function(path, transpose = FALSE) {
  raw <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE,
                                       progress = FALSE), check.names = FALSE)
  rn <- as.character(raw[[1L]])
  m <- as.matrix(raw[, -1L, drop = FALSE])
  rownames(m) <- rn
  storage.mode(m) <- "double"
  if (transpose) m <- t(m)
  m
}

# Build one per-context SummarizedExperiment from a bgzipped BED phenotype
# file (QTLtools layout: #chr/start/end/ID + sample columns) and an optional
# per-context covariate TSV (attached as colData).
.buildContextSe <- function(bedPath, covPath, transposeCov) {
  bed <- as.data.frame(readr::read_tsv(bedPath, show_col_types = FALSE,
                                       progress = FALSE, comment = ""),
                       check.names = FALSE)
  chrCol   <- intersect(c("#chr", "#chrom", "chrom", "chr"), names(bed))[1L]
  startCol <- intersect(c("start", "Start"), names(bed))[1L]
  endCol   <- intersect(c("end", "End"), names(bed))[1L]
  geneCol  <- intersect(c("gene_id", "ID", "phenotype_id"), names(bed))[1L]
  if (any(is.na(c(chrCol, startCol, endCol, geneCol)))) {
    stop("phenotype BED is missing one of chrom/start/end/gene_id columns: ",
         bedPath)
  }
  meta <- bed[, c(chrCol, startCol, endCol, geneCol)]
  names(meta) <- c("chrom", "start", "end", "gene_id")
  sampleCols <- setdiff(names(bed),
    c("#chr", "#chrom", "chrom", "chr", "start", "Start", "end", "End",
      "gene_id", "ID", "phenotype_id", "strand", "Strand"))
  expr <- as.matrix(bed[, sampleCols, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- meta$gene_id
  rr <- GenomicRanges::GRanges(
    seqnames = meta$chrom,
    ranges = IRanges::IRanges(start = meta$start + 1L, end = meta$end))
  names(rr) <- meta$gene_id
  cd <- if (!is.null(covPath) && nzchar(covPath)) {
    pcov <- .readCovariateMatrix(covPath, transposeCov)
    common <- intersect(rownames(pcov), colnames(expr))
    if (length(common) == 0L) {
      stop("No shared samples between phenotype and covariate file: ", bedPath)
    }
    expr <- expr[, common, drop = FALSE]
    S4Vectors::DataFrame(pcov[common, , drop = FALSE], row.names = common)
  } else {
    S4Vectors::DataFrame(row.names = colnames(expr))
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr), rowRanges = rr, colData = cd)
}

# Aliases for the QtlDataset (phenotype) manifest.
.qtlDatasetManifestAliases <- list(
  context               = c("context", "cond", "condition"),
  phenotypePath         = c("phenotypePath", "phenotype_path", "path",
                            "phenotypeFile", "phenotype_file"),
  covariatePath         = c("covariatePath", "covariate_path", "cov_path",
                            "covPath"),
  study                 = c("study", "study_id", "studyId"),
  genotypePath          = c("genotypePath", "genotype_path", "genotypePrefix",
                            "genotype_prefix", "genotypeFile"),
  genotypeCovariatePath = c("genotypeCovariatePath", "genotype_covariate_path",
                            "genotype_covariates", "genotypeCovariates"))

# Build one QtlDataset from an already-canonicalised set of manifest rows for a
# single study. `genotypesOverride` / `genoCovOverride` short-circuit the
# genotypePath / genotypeCovariatePath columns when supplied as arguments.
.buildQtlDatasetFromRows <- function(rows, study, base, transposeCov, qc,
                                     genotypesOverride = NULL,
                                     genoCovOverride = NULL) {
  contexts <- unique(as.character(rows$context))
  phenotypes <- setNames(lapply(contexts, function(cx) {
    sub <- rows[as.character(rows$context) == cx, , drop = FALSE]
    pths <- unique(sub$phenotypePath[!is.na(sub$phenotypePath)])
    if (length(pths) != 1L) {
      stop("Context '", cx, "' (study '", study,
           "') must reference exactly one phenotypePath; got: ",
           paste(pths, collapse = ", "))
    }
    covPath <- NULL
    if ("covariatePath" %in% names(sub)) {
      covs <- unique(sub$covariatePath[!is.na(sub$covariatePath) &
                                       nzchar(as.character(sub$covariatePath))])
      if (length(covs) > 1L) {
        stop("Context '", cx, "' (study '", study,
             "') references multiple covariatePath values: ",
             paste(covs, collapse = ", "))
      }
      if (length(covs) == 1L) covPath <- .resolveRel(covs[[1L]], base)
    }
    .buildContextSe(.resolveRel(pths[[1L]], base), covPath, transposeCov)
  }), contexts)

  genoHandle <- if (methods::is(genotypesOverride, "GenotypeHandle")) {
    genotypesOverride
  } else {
    spec <- if (is.character(genotypesOverride)) {
      genotypesOverride
    } else {
      v <- unique(as.character(rows$genotypePath[!is.na(rows$genotypePath)]))
      if (length(v) != 1L) {
        stop("study '", study, "' must reference exactly one genotypePath; got: ",
             paste(v, collapse = ", "))
      }
      v[[1L]]
    }
    .detectGenotypeFormat(.resolveRel(spec, base))
  }

  genoCov <- if (is.matrix(genoCovOverride)) {
    genoCovOverride
  } else {
    spec <- if (is.character(genoCovOverride)) {
      genoCovOverride
    } else if ("genotypeCovariatePath" %in% names(rows)) {
      v <- unique(as.character(rows$genotypeCovariatePath[
        !is.na(rows$genotypeCovariatePath) &
        nzchar(as.character(rows$genotypeCovariatePath))]))
      if (length(v) > 1L) {
        stop("study '", study, "' references multiple genotypeCovariatePath values.")
      }
      if (length(v) == 1L) v[[1L]] else NULL
    } else {
      NULL
    }
    if (!is.null(spec)) {
      .readCovariateMatrix(.resolveRel(spec, base), transposeCov)
    } else {
      matrix(numeric(0), nrow = 0, ncol = 0)
    }
  }

  do.call(QtlDataset, c(list(
    study = study, genotypes = genoHandle, phenotypes = phenotypes,
    genotypeCovariates = genoCov), qc))
}

# =============================================================================
# Public loaders
# =============================================================================

#' @title Load a QtlDataset from a manifest
#' @description Build a single-study \code{\link{QtlDataset}} from a manifest
#'   describing one row per QTL context (or per trait x context). The manifest
#'   is a data.frame or a path to a delimited file (\code{.csv} -> CSV, else
#'   TSV). Column names may be snake_case or camelCase.
#' @param manifest A data.frame or a path to a manifest file. Recognised
#'   columns (canonical camelCase; snake_case aliases accepted):
#'   \code{context} (required), \code{phenotypePath} (required, bgzipped BED),
#'   \code{covariatePath} (optional, per-context covariates), and the
#'   single-valued \code{study} / \code{genotypePath} /
#'   \code{genotypeCovariatePath}.
#' @param study Study identifier; reconciled with a \code{study} column.
#' @param genotypes A \code{\link{GenotypeHandle}} or a genotype path/prefix;
#'   reconciled with a \code{genotypePath} column.
#' @param genotypeCovariates A numeric matrix (samples x covariates) or a path
#'   to a covariate TSV; reconciled with a \code{genotypeCovariatePath} column.
#' @param scaleResiduals,mafCutoff,macCutoff,xvarCutoff,imissCutoff,keepSamples,keepVariants,keepIndel
#'   Pass-through \code{\link{QtlDataset}} arguments (stored as lazy QC slots).
#' @param transposeCovariates Transpose covariate TSVs (QTLtools layout) before
#'   treating them as samples-as-rows.
#' @return A \code{QtlDataset} object.
#' @export
loadQtlDatasetFromManifest <- function(manifest, study = NULL,
                                       genotypes = NULL,
                                       genotypeCovariates = NULL,
                                       scaleResiduals = TRUE,
                                       mafCutoff = 0, macCutoff = 0,
                                       xvarCutoff = 0, imissCutoff = 0,
                                       keepSamples = character(0),
                                       keepVariants = character(0),
                                       keepIndel = TRUE,
                                       transposeCovariates = FALSE) {
  base <- .manifestBase(manifest)
  df <- .canonManifestCols(.readManifest(manifest),
                           .qtlDatasetManifestAliases,
                           required = c("context", "phenotypePath"),
                           label = "QtlDataset")
  study <- .reconcileScalar(df$study, study, "study")
  qc <- list(scaleResiduals = scaleResiduals, mafCutoff = mafCutoff,
             macCutoff = macCutoff, xvarCutoff = xvarCutoff,
             imissCutoff = imissCutoff, keepSamples = keepSamples,
             keepVariants = keepVariants, keepIndel = keepIndel)
  .buildQtlDatasetFromRows(df, study, base, transposeCovariates, qc,
                           genotypesOverride = genotypes,
                           genoCovOverride = genotypeCovariates)
}

# Aliases shared by the sumstats manifests.
.gwasSumStatsManifestAliases <- list(
  study         = c("study", "study_id", "studyId"),
  sumStatsPath  = c("sumStatsPath", "sumstats_path", "path", "file_path",
                    "gwas_tsv", "filePath"),
  columnMapping = c("columnMapping", "column_mapping", "column_mapping_file",
                    "columnMappingFile"),
  nCase         = c("nCase", "n_case"),
  nControl      = c("nControl", "n_control"),
  varY          = c("varY", "var_y"),
  genome        = c("genome"),
  ldSketchPath  = c("ldSketchPath", "ld_sketch_path", "ld_sketch"))

.qtlSumStatsManifestAliases <- list(
  study         = c("study", "study_id", "studyId"),
  context       = c("context", "cond", "condition"),
  trait         = c("trait", "gene_id", "geneId", "molecular_trait_id"),
  sumStatsPath  = c("sumStatsPath", "sumstats_path", "path", "file_path",
                    "filePath"),
  columnMapping = c("columnMapping", "column_mapping", "column_mapping_file",
                    "columnMappingFile"),
  varY          = c("varY", "var_y"),
  genome        = c("genome"),
  ldSketchPath  = c("ldSketchPath", "ld_sketch_path", "ld_sketch"))

#' @title Load a GwasSumStats collection from a manifest
#' @description Build a \code{\link{GwasSumStats}} from a manifest with one row
#'   per study. No QC is run (the result carries \code{qcInfo = list()}).
#' @param manifest A data.frame or path. Columns (snake_case aliases accepted):
#'   \code{study} (required, unique), \code{sumStatsPath} (required),
#'   \code{columnMapping} (optional), \code{nCase} / \code{nControl}
#'   (optional), \code{varY} (optional), and the single-valued \code{genome} /
#'   \code{ldSketchPath}.
#' @param genome Genome build; reconciled with a \code{genome} column.
#' @param ldSketch A \code{\link{GenotypeHandle}} or a spec
#'   (path/prefix/genoMeta); reconciled with an \code{ldSketchPath} column.
#' @param region Optional \code{chr:start-end} string, GRanges, or one-row
#'   data.frame restricting the variants read (honoured only for tabix-indexed
#'   text and bgzipped+tabixed VCF; ignored with a warning otherwise).
#' @param minLdOverlapWarn Warn when the fraction of sumstats variants present
#'   in the LD sketch falls below this (default 0.5).
#' @param columnMapping Optional default column mapping (a named list/vector,
#'   or a path to a YAML file of \code{standardName: sourceName} entries)
#'   applied when a row has no \code{columnMapping}.
#' @param sampleSelect Optional GWAS-VCF FORMAT sample (study) column to read.
#' @param formatMapping Optional GWAS-VCF FORMAT tag mapping (canonical stat ->
#'   FORMAT field), overriding ES/SE/LP/SS/EAF defaults.
#' @return A \code{GwasSumStats} object.
#' @export
loadGwasSumStatsFromManifest <- function(manifest, genome = NULL,
                                         ldSketch = NULL, region = NULL,
                                         minLdOverlapWarn = 0.5,
                                         columnMapping = NULL,
                                         sampleSelect = NULL,
                                         formatMapping = NULL) {
  base <- .manifestBase(manifest)
  df <- .canonManifestCols(.readManifest(manifest),
                           .gwasSumStatsManifestAliases,
                           required = c("study", "sumStatsPath"),
                           label = "GwasSumStats")
  if (anyDuplicated(as.character(df$study))) {
    stop("GwasSumStats manifest `study` values must be unique.")
  }
  genome <- .reconcileScalar(df$genome, genome, "genome")
  ldSketch <- .resolveLdSketchInput(df, ldSketch, base)

  entries <- lapply(seq_len(nrow(df)), function(i) {
    label <- paste0("GwasSumStats[study=", df$study[[i]], "]")
    mapping <- if ("columnMapping" %in% names(df) &&
                   !is.na(df$columnMapping[[i]]) &&
                   nzchar(as.character(df$columnMapping[[i]]))) {
      .resolveRel(as.character(df$columnMapping[[i]]), base)
    } else {
      columnMapping
    }
    gr <- .loadSumStatsEntry(.resolveRel(as.character(df$sumStatsPath[[i]]), base),
                             region, mapping, sampleSelect, formatMapping, label)
    .checkLdContainment(ldSketch, gr, minLdOverlapWarn, label)
    gr
  })

  args <- list(study = as.character(df$study), entry = entries,
               genome = genome, ldSketch = ldSketch)
  if ("nCase" %in% names(df))    args$nCase    <- as.numeric(df$nCase)
  if ("nControl" %in% names(df)) args$nControl <- as.numeric(df$nControl)
  if ("varY" %in% names(df))     args$varY     <- as.numeric(df$varY)
  do.call(GwasSumStats, args)
}

# Build the entry list + tuple vectors for a QtlSumStats manifest.
.loadQtlSumStatsEntries <- function(df, base, region, ldSketch,
                                    minLdOverlapWarn, columnMapping,
                                    sampleSelect, formatMapping) {
  lapply(seq_len(nrow(df)), function(i) {
    label <- paste0("QtlSumStats[", df$study[[i]], "/", df$context[[i]], "/",
                    df$trait[[i]], "]")
    mapping <- if ("columnMapping" %in% names(df) &&
                   !is.na(df$columnMapping[[i]]) &&
                   nzchar(as.character(df$columnMapping[[i]]))) {
      .resolveRel(as.character(df$columnMapping[[i]]), base)
    } else {
      columnMapping
    }
    gr <- .loadSumStatsEntry(.resolveRel(as.character(df$sumStatsPath[[i]]), base),
                             region, mapping, sampleSelect, formatMapping, label)
    if (!is.null(ldSketch)) {
      .checkLdContainment(ldSketch, gr, minLdOverlapWarn, label)
    }
    gr
  })
}

#' @title Load a QtlSumStats collection from a manifest
#' @description Build a \code{\link{QtlSumStats}} from a manifest with one row
#'   per \code{(study, context, trait)} tuple. No QC is run.
#' @param manifest A data.frame or path. Columns (snake_case aliases accepted):
#'   \code{study}, \code{context}, \code{trait} (required), \code{sumStatsPath}
#'   (required), \code{columnMapping} (optional), \code{varY} (optional), and
#'   the single-valued \code{genome} / \code{ldSketchPath}.
#' @param genome Genome build; reconciled with a \code{genome} column.
#' @param ldSketch A \code{\link{GenotypeHandle}} or spec; reconciled with an
#'   \code{ldSketchPath} column.
#' @param region,minLdOverlapWarn,columnMapping,sampleSelect,formatMapping As
#'   for \code{\link{loadGwasSumStatsFromManifest}}.
#' @return A \code{QtlSumStats} object.
#' @export
loadQtlSumStatsFromManifest <- function(manifest, genome = NULL,
                                        ldSketch = NULL, region = NULL,
                                        minLdOverlapWarn = 0.5,
                                        columnMapping = NULL,
                                        sampleSelect = NULL,
                                        formatMapping = NULL) {
  base <- .manifestBase(manifest)
  df <- .canonManifestCols(.readManifest(manifest),
                           .qtlSumStatsManifestAliases,
                           required = c("study", "context", "trait",
                                        "sumStatsPath"),
                           label = "QtlSumStats")
  genome <- .reconcileScalar(df$genome, genome, "genome")
  ldSketch <- .resolveLdSketchInput(df, ldSketch, base)
  entries <- .loadQtlSumStatsEntries(df, base, region, ldSketch,
                                     minLdOverlapWarn, columnMapping,
                                     sampleSelect, formatMapping)
  args <- list(study = as.character(df$study),
               context = as.character(df$context),
               trait = as.character(df$trait),
               entry = entries, genome = genome, ldSketch = ldSketch)
  if ("varY" %in% names(df)) args$varY <- as.numeric(df$varY)
  do.call(QtlSumStats, args)
}

#' @title Load a MultiStudyQtlDataset from manifests
#' @description Build a \code{\link{MultiStudyQtlDataset}} from a QtlDataset
#'   manifest (individual-level studies) and, optionally, a QtlSumStats
#'   manifest (summary-only studies).
#' @param qtlDatasetsManifest A data.frame or path. The QtlDataset schema with
#'   \code{study} and \code{genotypePath} \strong{mandatory} (grouping key =
#'   \code{study}); one \code{QtlDataset} is built per study group.
#' @param sumStatsManifest Optional QtlSumStats manifest (see
#'   \code{\link{loadQtlSumStatsFromManifest}}).
#' @param genome,ldSketch,region,minLdOverlapWarn,columnMapping,sampleSelect,formatMapping
#'   Passed to the QtlSumStats loader for \code{sumStatsManifest}.
#' @param transposeCovariates Transpose covariate TSVs (QTLtools layout).
#' @param scaleResiduals,mafCutoff,macCutoff,xvarCutoff,imissCutoff,keepSamples,keepVariants,keepIndel
#'   Pass-through \code{\link{QtlDataset}} arguments applied to every study.
#' @return A \code{MultiStudyQtlDataset} object.
#' @export
loadMultiStudyQtlDatasetFromManifest <- function(qtlDatasetsManifest,
                                                 sumStatsManifest = NULL,
                                                 genome = NULL, ldSketch = NULL,
                                                 region = NULL,
                                                 minLdOverlapWarn = 0.5,
                                                 columnMapping = NULL,
                                                 sampleSelect = NULL,
                                                 formatMapping = NULL,
                                                 transposeCovariates = FALSE,
                                                 scaleResiduals = TRUE,
                                                 mafCutoff = 0, macCutoff = 0,
                                                 xvarCutoff = 0, imissCutoff = 0,
                                                 keepSamples = character(0),
                                                 keepVariants = character(0),
                                                 keepIndel = TRUE) {
  base <- .manifestBase(qtlDatasetsManifest)
  df <- .canonManifestCols(.readManifest(qtlDatasetsManifest),
                           .qtlDatasetManifestAliases,
                           required = c("study", "genotypePath", "context",
                                        "phenotypePath"),
                           label = "MultiStudyQtlDataset (qtlDatasetsManifest)")
  qc <- list(scaleResiduals = scaleResiduals, mafCutoff = mafCutoff,
             macCutoff = macCutoff, xvarCutoff = xvarCutoff,
             imissCutoff = imissCutoff, keepSamples = keepSamples,
             keepVariants = keepVariants, keepIndel = keepIndel)
  studies <- unique(as.character(df$study))
  qtlDatasets <- setNames(lapply(studies, function(s) {
    .buildQtlDatasetFromRows(df[as.character(df$study) == s, , drop = FALSE],
                             s, base, transposeCovariates, qc)
  }), studies)

  sumStats <- NULL
  if (!is.null(sumStatsManifest)) {
    sumStats <- loadQtlSumStatsFromManifest(
      sumStatsManifest, genome = genome, ldSketch = ldSketch, region = region,
      minLdOverlapWarn = minLdOverlapWarn, columnMapping = columnMapping,
      sampleSelect = sampleSelect, formatMapping = formatMapping)
  }
  MultiStudyQtlDataset(qtlDatasets = qtlDatasets, sumStats = sumStats)
}
