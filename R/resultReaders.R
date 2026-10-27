# =============================================================================
# Published result readers
# -----------------------------------------------------------------------------
# The inverse of R/resultWriters.R: take a file that was published -- by
# `writeBed()` / `writeGwasVcf()`, or by another tool writing the same
# vocabulary -- and rebuild the S4 collection.
#
#   readFineMappingResult(path, ...)  -> Qtl / GwasFineMappingResult
#   readTwasWeights(path, ...)        -> TwasWeights
#
# Two rules shape everything here.
#
# PERMISSIVE about columns. Each consumer needs a different few: TWAS needs
# only the weight, Mendelian randomization beta/se/PIP, colocalization the LBF
# matrix on its non-enloc path. A file converted from another tool will carry
# some of those and not others, so a missing column is restored as absent, not
# as an error -- the consumers name what they lack when they are reached.
# What IS required is only what the class itself cannot exist without: a
# coordinate-encoding variant id, and for a fine-mapping result a PIP.
#
# STRICT about identity. A merged table holding several (study, context,
# trait, method) tuples reads back as a multi-row collection, one row per
# tuple. Where a per-row value is recorded per variant -- the analysis window,
# the provenance flags -- every variant of a tuple must agree, and a
# disagreement is an error: it means two different results were concatenated
# under one identity, and picking either value would silently publish one
# fit's window over another's.
# =============================================================================

#' @include AllGenerics.R AllClasses.R resultWriters.R
#' @importFrom dplyr slice mutate select bind_cols all_of
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map_chr map_lgl list_rbind compact set_names
#' @importFrom stringr str_c str_detect str_remove str_to_lower
#' @importFrom stringr str_flatten str_split_1
#' @importFrom magrittr is_in
#' @importFrom glue glue
#' @importFrom dplyr if_else
#' @importFrom rlang abort warn arg_match try_fetch .data
#' @importFrom checkmate assertString assertFlag
NULL

# ---- the inverse of the published vocabulary --------------------------------

# Published column -> the CANONICAL WIDE `topLoci` column it was projected
# from. Only the pure renames live here. The four quantities that are
# TRANSFORMED on the way out are inverted explicitly below, each next to a
# comment naming the forward transform: `start`/`end` from `pos`, `LP` from
# `p`, `EZ` with its supplied flag, and `CS<n>` from the credible-set label.
#
# `MAF` has no canonical slot -- the posterior view publishes `af`, the
# directional frequency -- so it is restored under its own name and travels as
# a passthrough column rather than being mistaken for `af`.
# @noRd
.PUB_RESTORE_VARIANT <- c(
    `#chr` = "chrom",
    ID = "variant_id",
    REF = "A2",
    ALT = "A1",
    AF = "af",
    MAF = "maf",
    NS = "N",
    ES = "marginal_beta",
    SE = "marginal_se",
    WGT = "posterior_mean",
    WGTSE = "posterior_sd",
    THETA = "posterior_theta",
    PIP = "pip",
    LBF = "logBF",
    LFSR = "lfsr"
)

# Published column -> the collection mcols column it was read from. These are
# per-ROW, so they are lifted out of the variant table and reconciled across
# each tuple's variants rather than stored per variant.
# @noRd
.PUB_RESTORE_IDENTITY <- c(
    study = "study",
    context = "context",
    trait = "trait",
    method = "method",
    block_id = "blockId",
    data_type = "dataType",
    data_source = "dataSource",
    weight_standardized = "standardized",
    joint_studies = "jointStudies",
    joint_contexts = "jointContexts",
    joint_traits = "jointTraits",
    region_given = "regionGiven"
)

# The published columns that are structural or per-row, and so are never
# carried through as a variant-level passthrough column.
# @noRd
.PUB_NON_PAYLOAD <- c(
    "#chr",
    "start",
    "end",
    "region_id",
    "region_start",
    "region_end",
    names(.PUB_RESTORE_IDENTITY)
)

# The identity axes a published table is grouped into collection rows by.
# @noRd
.PRD_GROUP_COLS <- c("study", "context", "trait", "method", "block_id")

# ---- reading the file -------------------------------------------------------

# One published table as a tibble of character columns, region-restricted when
# the file is bgzipped + tabix-indexed.
#
# Every column is read as character and coerced afterwards, for the same
# reason `.readTabixRegion()` does it: an allele or id column that happens to
# be uniformly "T" within a region would otherwise be type-guessed as logical,
# silently corrupting the variant identity.
# @noRd
.prdReadTable <- function(path, which, label) {
    if (!file.exists(path)) {
        abort(glue("{label}: file '{path}' does not exist."))
    }
    if (!is.null(which)) {
        if (file.exists(str_c(path, ".tbi"))) {
            return(.readTabixRegion(path, which))
        }
        warn(glue(
            "{label}: '{path}' is not bgzipped + tabix-indexed; `which` ",
            "ignored, reading the whole file."
        ))
    }
    readr::read_tsv(
        path,
        skip = length(.prdBedMetaLines(path)),
        show_col_types = FALSE,
        progress = FALSE,
        col_types = readr::cols(.default = readr::col_character())
    )
}

# Which reader a path asks for. The writers choose the extension, so `"auto"`
# is unambiguous on pecotmr's own output; `format` is for a file named
# something else.
# @noRd
.prdFormat <- function(path, format) {
    format <- arg_match(format, c("auto", "bed", "vcf"))
    if (format != "auto") {
        return(format)
    }
    stem <- str_remove(path, "\\.(gz|bgz|bz2)$")
    if (str_detect(stem, "\\.(vcf|bcf)$")) "vcf" else "bed"
}

# ---- reading a published VCF ------------------------------------------------

# One published GWAS-VCF as the SAME long table the BED reader produces: the
# structural columns, one column per FORMAT key, and the identity columns read
# off the `##SAMPLE` header block. Each sample is one fitted model, so each
# contributes its own block of rows.
# @noRd
.prdVcfTable <- function(path, which, sampleSelect, label) {
    vcf <- .prdReadVcf(path, which, label)
    samples <- colnames(vcf)
    if (!is.null(sampleSelect)) {
        absent <- setdiff(as.character(sampleSelect), samples)
        if (length(absent) > 0L) {
            abort(glue(
                "{label}: sampleSelect names {str_flatten(absent, ', ')}, ",
                "which the VCF does not have (have: ",
                "{str_flatten(samples, ', ')})."
            ))
        }
        samples <- as.character(sampleSelect)
    }
    if (length(samples) == 0L) {
        # A sites-only VCF -- no FORMAT, no sample columns -- carries no
        # per-variant statistics at all. Returning an empty table instead
        # only defers the failure to whichever required column is missed
        # first, which reported "no `PIP` column" for a file whose real
        # problem is that it holds no results.
        abort(glue(
            "{label}: the VCF has no samples, so it carries no published ",
            "results. A published file records each result as a sample ",
            "with its own FORMAT fields."
        ))
    }
    meta <- .prdVcfSampleMeta(vcf)
    pieces <- map(samples, .prdVcfSampleTable, vcf = vcf, meta = meta)
    list(
        build = .prdVcfBuild(vcf),
        table = list_rbind(pieces),
        fits = compact(list_flatten(unname(map2(
            samples,
            pieces,
            .prdVcfSampleFits,
            meta = meta
        ))))
    )
}

# The genome build a VCF records on its ##reference line, or NA.
# @noRd
.prdVcfBuild <- function(vcf) {
    .vcfMetaBuild(VariantAnnotation::header(vcf))
}

# @noRd
.vcfMetaBuild <- function(hdr) {
    ref <- VariantAnnotation::meta(hdr)$reference
    if (is.null(ref) || nrow(ref) == 0L) {
        return(NA_character_)
    }
    as.character(ref$Value[[1L]])
}

# A sample's partial fits, keyed by the group each belongs to: one for a
# single fit, or one per LD block of a sample that joined several, taken from
# the header's "|"-packed values in the order its block list gives. The keys
# are the ones .prdGroups() splits the stacked table by, so each partial fit
# is found again once every sample's rows are stacked together.
# @noRd
.prdVcfSampleFits <- function(sample, tbl, meta) {
    if (is.null(meta) || !is_in(sample, rownames(meta)) || nrow(tbl) == 0L) {
        return(list())
    }
    row <- meta[sample, , drop = FALSE]
    keys <- .pubIdentityKey(tbl, .PRD_GROUP_COLS)
    first <- which(!duplicated(keys))
    pos <- .prdBlockPositions(row, tbl$block_id[first]) %||%
        rep(1L, length(first))
    set_names(map(pos, .prdVcfPartialFit, row = row), keys[first])
}

# Where each block sits in a joined sample's header lists, or NULL for a
# sample holding one block.
# @noRd
.prdBlockPositions <- function(row, blocks) {
    listed <- .prdBlockList(row)
    if (length(listed) < 2L || is.null(blocks) || all(is.na(blocks))) {
        return(NULL)
    }
    # An NA block is a variant this sample does not have -- another sample of
    # the file contributed it -- and it stays NA.
    pos <- match(blocks, listed)
    unknown <- !is.na(blocks) & is.na(pos)
    if (any(unknown)) {
        abort(glue(
            "readFineMappingResult: the BID field names block(s) ",
            "{str_flatten(unique(blocks[unknown]), ', ')} that the ",
            "sample's header does not list."
        ))
    }
    pos
}

# The LD blocks a sample's `##SAMPLE` line lists, in packing order.
# @noRd
.prdBlockList <- function(row) {
    cell <- row[["block_id"]]
    if (is.null(cell) || is.na(cell[[1L]])) {
        return(character(0))
    }
    str_split_1(as.character(cell[[1L]]), fixed("|"))
}

# @noRd
.prdReadVcf <- function(path, which, label) {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        abort(glue(
            "{label}: reading a published VCF requires the ",
            "'VariantAnnotation' package."
        ))
    }
    if (!file.exists(path)) {
        abort(glue("{label}: file '{path}' does not exist."))
    }
    if (is.null(which) || !file.exists(str_c(path, ".tbi"))) {
        if (!is.null(which)) {
            warn(glue(
                "{label}: '{path}' is not bgzipped + tabix-indexed; ",
                "`which` ignored, reading the whole file."
            ))
        }
        return(VariantAnnotation::readVcf(path))
    }
    tf <- Rsamtools::TabixFile(path)
    gr <- .matchRegionSeqnames(
        which,
        Rsamtools::seqnamesTabix(tf)
    )
    if (is.null(gr)) {
        return(VariantAnnotation::readVcf(path)[0, ])
    }
    VariantAnnotation::readVcf(tf, param = gr)
}

# The per-fit `##SAMPLE` block as a plain data.frame keyed by sample id, or
# NULL when the file has none. This is where `writeGwasVcf()` records the
# identity, the analysis window and the two per-EFFECT fit quantities (`V`,
# `cs_index`) that no per-variant column could hold.
# @noRd
.prdVcfSampleMeta <- function(vcf) {
    # `meta()` is VariantAnnotation's own accessor for the header's structured
    # `##` blocks -- S4Vectors::metadata() is a different thing and returns
    # nothing here.
    md <- try_fetch(
        VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE,
        error = function(cnd) NULL
    )
    if (is.null(md) || nrow(md) == 0L) {
        return(NULL)
    }
    df <- as.data.frame(md, stringsAsFactors = FALSE)
    # Back to the shared vocabulary through the writer's own table, so the
    # rest of the reader sees the same names a published BED carries.
    known <- is_in(names(df), .PUB_SAMPLE_KEYS)
    names(df)[known] <- names(.PUB_SAMPLE_KEYS)[
        match(names(df)[known], .PUB_SAMPLE_KEYS)
    ]
    df
}

# One sample's rows: the fixed columns, its FORMAT values, and the identity
# the header recorded for it.
# @noRd
.prdVcfSampleTable <- function(sample, vcf, meta) {
    rr <- SummarizedExperiment::rowRanges(vcf)
    n <- length(rr)
    geno <- VariantAnnotation::geno(vcf)
    fixed <- tibble(
        `#chr` = as.character(GenomicRanges::seqnames(rr)),
        start = as.integer(GenomicRanges::start(rr)) - 1L,
        end = as.integer(GenomicRanges::start(rr)),
        ID = names(rr) %||% rep(NA_character_, n),
        REF = unname(as.character(VariantAnnotation::ref(vcf))),
        ALT = .mlFirstAlts(VariantAnnotation::alt(vcf))
    )
    fmt <- set_names(
        map(names(geno), .prdVcfFormatColumn, geno = geno, sample = sample),
        map_chr(names(geno), .prdVcfKeyName)
    )
    # BID is structure, not a statistic: it says which LD block of a joined
    # sample each variant belongs to, and so feeds the identity columns.
    ident <- .prdVcfIdentity(meta, sample, n, fmt[["BID"]])
    fmt[["BID"]] <- NULL
    # The records are the union of every sample's variants, and a sample
    # leaves the cells of a variant it does not have empty. Those rows are
    # not its results, so they are dropped rather than read back as variants
    # with no statistics.
    has <- reduce(map(fmt, .prdNotNa), `|`, .init = rep(FALSE, n))
    slice(
        bind_cols(fixed, as_tibble(fmt, .name_repair = "minimal"), ident),
        which(has)
    )
}

# @noRd
.prdNotNa <- function(v) !is.na(v)

# One FORMAT field's column for one sample. A field the sample does not report
# is all-NA, which is what an absent column means on the way back in.
# @noRd
.prdVcfFormatColumn <- function(key, geno, sample) {
    m <- geno[[key]]
    if (!is_in(sample, colnames(m))) {
        return(rep(NA_character_, nrow(m)))
    }
    v <- as.character(unname(m[, sample]))
    # A `Number=A` geno field comes back as a CHARACTER matrix whose missing
    # cells are the literal strings "NA" (or VCF's own "."), not NA. Left in
    # place they stop a numeric field from being recognised as numeric, which
    # is how a published Float came back as a character column.
    replace(v, !is.na(v) & is_in(v, c("NA", ".", "")), NA_character_)
}

# A FORMAT key as the published column name. The fixed vocabulary is already
# shared between the two formats, so only the per-credible-set `fullFit`
# fields need undoing: `writeGwasVcf()` uppercases those column names because a
# FORMAT key is conventionally uppercase, and nothing else in the schema is
# case-ambiguous.
# @noRd
.prdVcfKeyName <- function(key) {
    lower <- str_to_lower(key)
    if (str_detect(lower, "^(within_cs_pip|cs_logbf_|cs_effect_)")) {
        return(lower)
    }
    key
}

# The identity columns for one sample, recycled over its variants.
#
# A file with no `##SAMPLE` block still has the sample NAME. Splitting it back
# into a tuple would be guesswork, but taking it whole is not: it becomes the
# `study`, which keeps each sample a result of its own rather than letting
# several merge into one identity-less row, and the reader's own arguments
# supply whatever else is wanted.
# @noRd
.prdVcfIdentity <- function(meta, sample, n, bid = NULL) {
    if (is.null(meta) || !is_in(sample, rownames(meta))) {
        return(tibble(study = rep(sample, length.out = n)))
    }
    row <- meta[sample, , drop = FALSE]
    cols <- intersect(
        c(
            names(.PUB_RESTORE_IDENTITY),
            "region_id",
            "region_start",
            "region_end"
        ),
        colnames(row)
    )
    if (length(cols) == 0L) {
        return(tibble(study = rep(sample, length.out = n)))
    }
    pos <- .prdBlockPositions(row, bid)
    out <- as_tibble(
        set_names(
            map(cols, .prdMetaValue, row = row, n = n, pos = pos),
            cols
        ),
        .name_repair = "minimal"
    )
    .prdVcfRegionStart(out)
}

# The header fields a sample joining several LD blocks packs one value per
# block of, separated by "|" (see writeGwasVcf).
# @noRd
.PRD_BLOCK_FIELDS <- c(
    "block_id",
    "region_id",
    "region_start",
    "region_end",
    "region_given"
)

# A VCF header's `RegionStart` is 1-based, because VCF coordinates are; a
# published table's `region_start` is 0-based half-open, because BED
# coordinates are. Each format is right on its own terms, so the long form
# the readers share settles on the table's convention and the VCF value is
# shifted into it here -- once, rather than at every reader of the column.
# @noRd
.prdVcfRegionStart <- function(out) {
    if (!is_in("region_start", names(out))) {
        return(out)
    }
    mutate(
        out,
        region_start = as.character(
            suppressWarnings(as.integer(.data$region_start)) - 1L
        )
    )
}

# One header field's value, recycled over the sample's variants. "." is VCF's
# missing value and `.vcfHeaderValue()` writes it for a field one sample has
# and another does not, so it is restored as NA rather than as the string ".".
# @noRd
.prdMetaValue <- function(nm, row, n, pos = NULL) {
    v <- as.character(row[[nm]][[1L]])
    if (is.na(v) || v == ".") {
        return(rep(NA_character_, length.out = n))
    }
    if (is.null(pos) || !is_in(nm, .PRD_BLOCK_FIELDS)) {
        return(rep(v, length.out = n))
    }
    # A joined sample: each variant takes its own block's value.
    perBlock <- str_split_1(v, fixed("|"))[pos]
    replace(perBlock, perBlock == "", NA_character_)
}

# The per-effect fit quantities the `##SAMPLE` block carried, as the partial
# `susieFit` of one row -- `V` (the estimated prior variances) and the credible
# sets' `cs_index`. They are per EFFECT, not per variant, which is why they
# live in the header at all, and colocalization's credible-set filter reads
# them straight off the fit. Everything else a fit holds (alpha, lbf_variable,
# the posterior matrices) is not published, so this is deliberately partial:
# `NULL` when the header recorded neither.
# @noRd
.prdVcfPartialFit <- function(pos, row) {
    v <- .prdUnpackNumeric(.prdBlockPiece(row[["V"]], pos))
    csIndex <- .prdUnpackNumeric(.prdBlockPiece(row[["cs_index"]], pos))
    out <- compact(list(
        V = v,
        sets = if (!is.null(csIndex)) list(cs_index = as.integer(csIndex))
    ))
    if (length(out) == 0L) NULL else out
}

# The `pos`-th block's part of a header cell: a sample joining several LD
# blocks packs one value per block, separated by "|", and a single-block
# sample's cell is its one value whole.
# @noRd
.prdBlockPiece <- function(cell, pos) {
    if (is.null(cell) || length(cell) == 0L) {
        return(NULL)
    }
    chr <- as.character(cell[[1L]])
    if (is.na(chr) || chr == "." || is.na(pos)) {
        return(NULL)
    }
    pieces <- str_split_1(chr, fixed("|"))
    if (length(pieces) == 1L) pieces else pieces[[pos]]
}

# A ";"-packed header value as a numeric vector, or NULL when it was NA.
# @noRd
.prdUnpackNumeric <- function(cell) {
    if (is.null(cell) || length(cell) == 0L) {
        return(NULL)
    }
    chr <- as.character(cell[[1L]])
    if (length(chr) == 0L || is.na(chr) || chr == "") {
        return(NULL)
    }
    suppressWarnings(as.numeric(str_split_1(chr, ";")))
}

# ---- published table -> canonical per-variant table -------------------------

# The canonical wide `topLoci` table for one published table: the pure renames
# applied, the four transformed quantities inverted, and every other column
# carried through verbatim so a `columns = "all"` export and a third party's
# own annotations both survive the trip.
# @noRd
.prdCanonical <- function(df, method, priorVariance) {
    pos <- .prdPos(df)
    out <- tibble(
        variant_id = .prdVariantIds(df, pos),
        chrom = .prdChr(df),
        pos = pos
    )
    renamed <- intersect(names(.PUB_RESTORE_VARIANT), names(df))
    renamed <- setdiff(renamed, c("#chr", "ID"))
    out <- bind_cols(
        out,
        as_tibble(
            set_names(
                map(renamed, .prdRestoreColumn, df = df),
                .PUB_RESTORE_VARIANT[renamed]
            ),
            .name_repair = "minimal"
        )
    )
    out <- bind_cols(
        out,
        .prdMarginalZ(df),
        .prdMarginalP(df),
        .prdCsColumns(df, method),
        .prdPassthrough(df)
    )
    .prdWakefieldLbf(out, priorVariance)
}

# A published column restored to its canonical name, typed. Everything the
# renaming table covers is numeric except the four character columns, which
# are the alleles and the id-like fields.
# @noRd
.prdRestoreColumn <- function(nm, df) {
    v <- df[[nm]]
    if (is_in(nm, c("REF", "ALT"))) {
        return(as.character(v))
    }
    suppressWarnings(as.numeric(v))
}

# The 1-based position. A published range is 0-based half-open -- that is the
# BED format's requirement -- so `end` IS the 1-based position of a
# single-base record, and `start + 1` says the same thing. `end` is preferred
# because it is what the writer set from `pos` directly.
# @noRd
.prdPos <- function(df) {
    if (is_in("end", names(df))) {
        return(suppressWarnings(as.integer(df$end)))
    }
    if (is_in("start", names(df))) {
        return(suppressWarnings(as.integer(df$start)) + 1L)
    }
    rep(NA_integer_, nrow(df))
}

# @noRd
.prdChr <- function(df) {
    if (is_in("#chr", names(df))) {
        return(as.character(df[["#chr"]]))
    }
    rep(NA_character_, nrow(df))
}

# The variant ids, which are the one thing a result cannot be built without:
# an id renders a range and its alleles, so `.variantIdsToGRanges()` refuses
# one that does not encode coordinates. A published `ID` normally does, having
# been written from exactly that. When it does not -- an rsid-keyed file from
# another tool -- the id is rebuilt from the coordinates and alleles that are
# published alongside it, which is the same string the writer would have
# produced.
# @noRd
.prdVariantIds <- function(df, pos) {
    ids <- if (is_in("ID", names(df))) as.character(df$ID) else NULL
    if (!is.null(ids) && !any(is.na(parseVariantId(ids)$pos))) {
        return(ids)
    }
    parts <- map(c("#chr", "REF", "ALT"), .prdCol, df = df)
    if (any(map_lgl(parts, is.null)) || anyNA(pos)) {
        abort(str_c(
            "A published table must identify its variants: either an `ID` ",
            "column of chrom:pos:ref:alt ids, or `#chr`, `start`/`end`, ",
            "`REF` and `ALT` to build them from. An id without coordinates ",
            "has no variant identity to store."
        ))
    }
    formatVariantId(
        as.character(parts[[1L]]),
        pos,
        as.character(parts[[2L]]),
        as.character(parts[[3L]])
    )
}

# @noRd
.prdCol <- function(nm, df) {
    if (is_in(nm, names(df))) df[[nm]] else NULL
}

# `EZ` is the z the SOURCE supplied, and the spec means it that way: the z
# "provided if it was used to derive the ES and SE fields". So a published EZ
# restores `marginal_z_supplied = TRUE`, and its absence leaves z derived from
# ES / SE with the flag FALSE -- which is what keeps a re-export from
# publishing an EZ that pecotmr computed itself.
# @noRd
.prdMarginalZ <- function(df) {
    n <- nrow(df)
    ez <- .prdNumericOrNull(df, "EZ")
    if (!is.null(ez) && any(!is.na(ez))) {
        return(tibble(
            marginal_z = ez,
            marginal_z_supplied = !is.na(ez)
        ))
    }
    es <- .prdNumericOrNull(df, "ES")
    se <- .prdNumericOrNull(df, "SE")
    if (is.null(es) || is.null(se)) {
        return(tibble(.rows = n))
    }
    tibble(
        marginal_z = if_else(is.na(se) | se <= 0, NA_real_, es / se),
        marginal_z_supplied = rep(FALSE, n)
    )
}

# `LP` is -log10 p, which is the better published form (invertible, and no
# underflow at tiny p); the canonical table stores the p itself.
# @noRd
.prdMarginalP <- function(df) {
    lp <- .prdNumericOrNull(df, "LP")
    if (is.null(lp)) {
        return(tibble(.rows = nrow(df)))
    }
    tibble(marginal_p = 10^(-lp))
}

# @noRd
.prdNumericOrNull <- function(df, nm) {
    if (!is_in(nm, names(df))) {
        return(NULL)
    }
    suppressWarnings(as.numeric(df[[nm]]))
}

# The credible-set columns. Both formats publish an integer INDEX, because the
# label's prefix is the fitting method, which is already its own field. The
# canonical table stores the label, so it is rebuilt as `<method>_<index>`;
# NA stays NA, which is "in no credible set at this coverage".
# @noRd
.prdCsColumns <- function(df, method) {
    csNames <- names(df)[str_detect(names(df), "^CS[0-9.]+$")]
    if (length(csNames) == 0L) {
        return(tibble(.rows = nrow(df)))
    }
    levels <- str_remove(csNames, "^CS")
    cols <- c(
        set_names(
            map(csNames, .prdCsLabels, df = df, method = method),
            str_c("cs_", levels)
        ),
        compact(set_names(
            map(str_c("PUR", levels), .prdNumericOrNull, df = df),
            str_c("cs_", levels, "_purity")
        ))
    )
    as_tibble(cols, .name_repair = "minimal")
}

# @noRd
.prdCsLabels <- function(nm, df, method) {
    idx <- suppressWarnings(as.integer(df[[nm]]))
    prefix <- if (is.null(method) || is.na(method)) "cs" else method
    if_else(
        is.na(idx) | idx == 0L,
        NA_character_,
        str_c(prefix, "_", idx)
    )
}

# Every published column the schema does not already account for, carried
# through under its own name: the per-credible-set `fullFit` statistics, and
# whatever a `columns = "all"` export or another tool added.
#
# Each is numeric when every non-NA value parses as a number, and character
# otherwise. Guessing per column rather than per value is what keeps an allele
# or label column from being coerced to NA, and the file itself carries no
# types to consult.
# @noRd
.prdPassthrough <- function(df) {
    known <- c(
        names(.PUB_RESTORE_VARIANT),
        .PUB_NON_PAYLOAD,
        "ES",
        "SE",
        "EZ",
        "LP",
        names(df)[str_detect(names(df), "^(CS|PUR)[0-9.]+$")]
    )
    extra <- setdiff(names(df), known)
    if (length(extra) == 0L) {
        return(tibble(.rows = nrow(df)))
    }
    as_tibble(
        set_names(map(extra, .prdAutoType, df = df), extra),
        .name_repair = "minimal"
    )
}

# @noRd
.prdAutoType <- function(nm, df) {
    v <- as.character(df[[nm]])
    num <- suppressWarnings(as.numeric(v))
    if (all(is.na(num) == is.na(v))) num else v
}

# Wakefield's approximate Bayes factor, as an OPT-IN substitute for a `LBF`
# the file does not carry. Opt-in because it is an approximation with a prior
# the caller has to choose: `priorVariance` is Wakefield's W, and the
# resulting log Bayes factors are not what susie's `lbf_variable` would have
# been. Colocalization's non-enloc path needs a real LBF matrix and says so
# rather than accepting this.
# @noRd
.prdWakefieldLbf <- function(out, priorVariance) {
    if (is.null(priorVariance)) {
        return(out)
    }
    assertNumber(priorVariance, lower = 0)
    if (is_in("logBF", names(out)) && any(!is.na(out$logBF))) {
        abort(str_c(
            "`priorVariance` asks for Wakefield log Bayes factors, but this ",
            "file already carries LBF. Drop `priorVariance` to keep the ",
            "published values."
        ))
    }
    z <- .prdCol("marginal_z", out)
    se <- .prdCol("marginal_se", out)
    absent <- c("z (EZ, or ES/SE)", "standard error (SE)")[
        c(is.null(z), is.null(se))
    ]
    if (length(absent) > 0L) {
        abort(glue(
            "`priorVariance` asks for Wakefield log Bayes factors, which ",
            "need a z and a standard error. This file publishes no ",
            "{str_flatten(absent, ' and no ')}."
        ))
    }
    mutate(out, logBF = .approxBayesFactor(z, se, W = priorVariance)$lbf)
}

# ---- published table -> collection rows -------------------------------------

# One published table read into the long form both readers work from.
# @noRd
.prdRead <- function(path, which, format, sampleSelect, label) {
    assertString(path)
    if (!is.null(which)) {
        .assertRegionGRanges(which, "`which`")
    }
    if (.prdFormat(path, format) == "vcf") {
        return(.prdVcfTable(path, which, sampleSelect, label))
    }
    # A published table has no header to carry a fit in -- only the VCF's
    # `##SAMPLE` block does -- so a BED read restores no fit at all.
    list(
        table = .prdReadTable(path, which, label),
        fits = list(),
        build = .prdBedBuild(path)
    )
}

# The genome build a published BED records on its `##genome=` line, or NA.
# @noRd
.prdBedBuild <- function(path) {
    meta <- .prdBedMetaLines(path)
    hit <- meta[str_detect(meta, "^##genome=")]
    if (length(hit) == 0L) {
        NA_character_
    } else {
        str_remove(hit[[1L]], "^##genome=")
    }
}

# The `##` lines a published BED opens with, above its column header.
# @noRd
.prdBedMetaLines <- function(path) {
    head <- readr::read_lines(path, n_max = 20L, progress = FALSE)
    head[seq_len(.prdLeadingCount(str_detect(head, "^##")))]
}

# @noRd
.prdLeadingCount <- function(flags) {
    if (all(flags)) length(flags) else which(!flags)[[1L]] - 1L
}

# The row groups of a published table, in first-appearance order: one per
# (study, context, trait, method, block_id) tuple the file carries. A table
# with none of those axes is a single result, which is what a bare
# single-entry export writes.
# @noRd
.prdGroups <- function(df) {
    # A header-only file -- which is what an export of an empty result writes
    # -- has no groups at all, not one empty group: the identity axes are
    # column names there, with no values to key anything by.
    if (nrow(df) == 0L) {
        return(list())
    }
    cols <- intersect(.PRD_GROUP_COLS, names(df))
    if (length(cols) == 0L) {
        return(list(seq_len(nrow(df))))
    }
    key <- .pubIdentityKey(df, cols)
    split(seq_len(nrow(df)), factor(key, levels = unique(key)))
}

# The single value a per-ROW column takes over one group's variants.
#
# Several distinct values mean two different results were concatenated under
# one identity tuple, and picking either would silently publish one fit's
# window or provenance over another's -- so it is an error naming the column
# and what it held.
# @noRd
.prdGroupValue <- function(df, nm, idx, label) {
    if (!is_in(nm, names(df))) {
        return(NA_character_)
    }
    vals <- unique(as.character(df[[nm]][idx]))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) {
        return(NA_character_)
    }
    if (length(vals) > 1L) {
        abort(glue(
            "{label}: `{nm}` takes {length(vals)} different values ",
            "({str_flatten(vals, ', ')}) within one ",
            "(study, context, trait, method) group. It is a per-result ",
            "value, so two results were concatenated under one identity."
        ))
    }
    vals[[1L]]
}

# Reconcile what the file recorded with what the caller supplied. The
# argument FILLS IN a column the file does not carry -- which is the point of
# having it, for a file a user assembled themselves -- and a disagreement is
# an error rather than a silent override, because one of the two is wrong and
# the reader cannot tell which.
# `supplied` is never NULL here: both callers build it with `compact()`, so
# only the arguments the caller actually gave reach this function, and a
# NULL-guard would be dead by construction.
# @noRd
.prdResolve <- function(fileValue, supplied, nm, label) {
    sup <- as.character(supplied)[[1L]]
    if (!is.na(fileValue) && !identical(fileValue, sup)) {
        abort(glue(
            "{label}: `{nm}` was given as '{sup}' but the file records ",
            "'{fileValue}'. Drop the argument to use the file's value."
        ))
    }
    sup
}

# The per-row identity and provenance of one group.
# @noRd
.prdRowFields <- function(df, idx, supplied, label) {
    fromFile <- set_names(
        map(
            names(.PUB_RESTORE_IDENTITY),
            .prdGroupValue,
            df = df,
            idx = idx,
            label = label
        ),
        .PUB_RESTORE_IDENTITY
    )
    for (nm in names(supplied)) {
        fromFile[[nm]] <- .prdResolve(
            fromFile[[nm]],
            supplied[[nm]],
            nm,
            label
        )
    }
    fromFile
}

# The analysis WINDOW of one group as a one-range GRanges, or NULL when the
# file did not publish it. `region_start` is 0-based half-open like every
# other published coordinate, so it is shifted back by one; `region_id` stays
# the 1-based `chrom_start_end` identifier and is not a coordinate source.
# The chromosome is spelled in UCSC style, as the variants' ranges are.
# @noRd
.prdRowRegion <- function(df, idx, label) {
    start <- .prdGroupValue(df, "region_start", idx, label)
    end <- .prdGroupValue(df, "region_end", idx, label)
    chrom <- .prdRegionChrom(df, idx, label)
    if (is.na(start) || is.na(end) || is.na(chrom)) {
        return(NULL)
    }
    GenomicRanges::GRanges(
        withChrPrefix(chrom),
        IRanges::IRanges(
            start = as.integer(start) + 1L,
            end = as.integer(end)
        )
    )
}

# The window's chromosome: from `region_id` when it has one (that is what the
# identifier encodes), else the variants' own, which cannot disagree with it
# because a result is fine-mapped within one block.
# @noRd
.prdRegionChrom <- function(df, idx, label) {
    id <- .prdGroupValue(df, "region_id", idx, label)
    if (!is.na(id) && str_detect(id, "_")) {
        return(str_remove(id, "_[^_]*_[^_]*$"))
    }
    if (!is_in("#chr", names(df))) {
        return(NA_character_)
    }
    chrom <- unique(as.character(df[["#chr"]][idx]))
    if (length(chrom) == 1L) chrom else NA_character_
}

# Stack the per-group region ranges into the collections' `region` column:
# one range per row, or NULL when no row published one. A partial set would
# have to invent a window for the rest, so it is all or nothing.
# @noRd
.prdRegionColumn <- function(regions) {
    if (any(map_lgl(regions, is.null))) {
        return(NULL)
    }
    if (length(regions) == 0L) {
        return(NULL)
    }
    # Through a GRangesList, not `c()`: `c()` on a list of S4 objects falls
    # back to the list method and hands back a list, which the constructor
    # then rejects as "not a GRanges". Names come from the group keys and
    # would otherwise land on the stored ranges.
    unname(unlist(.rtlBindElements(regions)))
}

# A per-row character column for a constructor argument, or NULL when no row
# recorded it -- which is how the constructors are told to omit the column.
# @noRd
.prdRowColumn <- function(rows, nm) {
    v <- map_chr(rows, .prdRowField, nm = nm)
    if (all(is.na(v))) NULL else v
}

# @noRd
.prdRowField <- function(r, nm) {
    v <- r[[nm]]
    if (is.null(v) || length(v) == 0L) NA_character_ else as.character(v)[[1L]]
}

# The identity axes a class cannot be built without, checked once with every
# absent one named along with the argument that supplies it. The constructors
# reject NA on these axes, and a published file is free not to have recorded
# them -- a user who assembled the table themselves then says what they are.
# @noRd
.prdRequireIdentity <- function(rows, axes, label) {
    absent <- keep(axes, .prdAxisAbsent, rows = rows)
    if (length(absent) == 0L) {
        return(invisible(NULL))
    }
    args <- c(
        study = "studyName",
        context = "context",
        trait = "trait",
        method = "method"
    )
    abort(glue(
        "{label}: the file does not record {str_flatten(absent, ', ')}, ",
        "which this result is keyed by. Pass ",
        "{str_flatten(str_c('`', args[absent], '`'), ', ')} to supply ",
        "{if (length(absent) == 1L) 'it' else 'them'}."
    ))
}

# @noRd
.prdAxisAbsent <- function(axis, rows) {
    is.null(.prdRowColumn(rows, axis))
}

# `weight_standardized` describes the SCALE of the published weight, and it
# cannot be recovered from the numbers: the standardized-to-raw conversion is
# a per-variant factor 1/sqrt(2*af*(1-af)), so the two scales differ by an
# af-dependent amount at every variant and neither is recognisable on its own.
# So when a file publishes a weight without the flag, the caller has to say.
# @noRd
.prdRequireStandardized <- function(df, supplied, label) {
    if (!is.null(supplied)) {
        assertFlag(supplied)
        return(invisible(NULL))
    }
    if (!is_in("WGT", names(df)) || all(is.na(df$WGT))) {
        return(invisible(NULL))
    }
    if (is_in("weight_standardized", names(df))) {
        if (any(!is.na(df$weight_standardized))) {
            return(invisible(NULL))
        }
    }
    abort(glue(
        "{label}: this file publishes WGT but not `weight_standardized`, ",
        "and the scale cannot be recovered from the weights -- converting ",
        "between them rescales every variant by a different, allele-",
        "frequency-dependent factor. Pass `weightStandardized` to say ",
        "whether the genotypes behind these weights were standardized."
    ))
}

# @noRd
.prdFlag <- function(v) {
    if (is.null(v) || length(v) == 0L || is.na(v)) NA else as.logical(v)
}

# ---- entry point: fine-mapping results --------------------------------------

#' @title Read A Published Fine-Mapping Result
#' @description Rebuild a \code{QtlFineMappingResult} or
#'   \code{GwasFineMappingResult} from a file \code{\link{writeBed}} or
#'   \code{\link{writeGwasVcf}} wrote -- or from any file in
#'   the same vocabulary, which is what makes this the way into pecotmr for
#'   results produced elsewhere.
#'
#'   A merged table holding several \code{(study, context, trait, method)}
#'   tuples reads back as a multi-row collection, one row per tuple. Which
#'   class is returned follows the axes the file has: a result with a context
#'   or a trait is a QTL result, one without is a GWAS result.
#'
#' @section What is required, and what is not:
#'   Only what the class cannot exist without: variant ids that encode
#'   coordinates (or the \code{#chr} / \code{start} / \code{end} / \code{REF}
#'   / \code{ALT} columns to build them from) and a \code{PIP} column.
#'   Everything else is optional, because each consumer needs a different few
#'   columns -- Mendelian randomization the marginal effect and its standard
#'   error, colocalization the log Bayes factors -- and a file converted from
#'   another tool will have some and not others. A consumer reached without
#'   what it needs says which column is missing.
#'
#'   The fit itself is not published. A VCF's \code{##SAMPLE} header records
#'   the two per-effect quantities that no per-variant column could hold --
#'   the prior variances \code{V} and the credible sets' \code{cs_index} --
#'   and those are restored as a partial \code{susieFit}; a BED carries
#'   neither. Anything needing the posterior matrices (colocalization's LBF
#'   path) wants the \code{LBF} columns, not the fit.
#'
#' @section What a round trip changes:
#'   Every published quantity comes back identical, in either format. Three
#'   things are not identities:
#'   \itemize{
#'     \item \code{writeBed()} publishes the \strong{filtered} top-loci
#'       table, so a default export holds the variants that passed
#'       \code{signalCutoff} / \code{minPurity} rather than every variant of
#'       the fit. Export with \code{signalCutoff = 0} for the whole table.
#'     \item Credible-set membership comes back as \code{NA} where the fit
#'       recorded the \code{"<method>_0"} not-captured sentinel. Both formats
#'       publish an integer index, and \code{NA} is what every reader of the
#'       column already treats as "in no credible set".
#'     \item \code{data_source} and \code{weight_standardized} are
#'       populated on the way back -- the writer derives the first and the
#'       reader is told the second -- where a fit predating them left them
#'       empty.
#'   }
#'
#' @param path Path to the published file. Bgzipped and tabix-indexed files
#'   can be read in part, with \code{which}.
#' @param which Optional \code{GRanges}: read only the records that overlap
#'   it, as \code{which} does for \code{rtracklayer::import()} and
#'   \code{VariantAnnotation::ScanVcfParam()}. Ignored with a warning when
#'   the file is not tabix-indexed.
#' @param format \code{"auto"} (default) reads a \code{.vcf} / \code{.bcf}
#'   path as a VCF and anything else as a published table; \code{"bed"} /
#'   \code{"vcf"} say so explicitly for a file named otherwise.
#' @param studyName,context,trait,method Identity values for a file that did
#'   not record them. Each FILLS IN the corresponding column when the file
#'   lacks it; a value that disagrees with what the file does record is an
#'   error rather than a silent override.
#' @param dataType Optional assay / modality label, same fill-in rule.
#' @param weightStandardized Whether the published posterior effect
#'   (\code{WGT}) is on the standardized genotype scale. Required when the
#'   file publishes \code{WGT} without \code{weight_standardized}: the scale
#'   is not recoverable from the weights, since converting between the two
#'   rescales every variant by its own allele-frequency-dependent factor.
#' @param sampleSelect For a multi-sample VCF, which samples to read. Each
#'   sample is one fitted model. Defaults to all of them.
#' @param priorVariance Opt-in Wakefield prior variance (\code{W}) for
#'   deriving per-variant log Bayes factors from z when the file has no
#'   \code{LBF}. These are an approximation under a prior you chose, not
#'   susie's \code{lbf_variable}, so nothing derives them unasked.
#' @return A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @seealso \code{\link{writeBed}}, \code{\link{writeGwasVcf}},
#'   \code{\link{readTwasWeights}}
#' @examples
#' data(qtlFineMappingExample)
#' path <- file.path(tempdir(), "top_loci.bed")
#' writeBed(qtlFineMappingExample, path)
#' # This fit predates the provenance columns, so its published table does
#' # not say which genotype scale its posterior effects are on -- the reader
#' # will not guess.
#' readFineMappingResult(path, weightStandardized = TRUE)
#' @importFrom purrr keep imap
#' @importFrom checkmate assertNumber
#' @export
readFineMappingResult <- function(
    path,
    which = NULL,
    format = c("auto", "bed", "vcf"),
    studyName = NULL,
    context = NULL,
    trait = NULL,
    method = NULL,
    dataType = NULL,
    weightStandardized = NULL,
    sampleSelect = NULL,
    priorVariance = NULL
) {
    label <- "readFineMappingResult"
    src <- .prdRead(path, which, format, sampleSelect, label)
    df <- src$table
    .prdRequirePip(df, label)
    .prdRequireStandardized(df, weightStandardized, label)
    supplied <- compact(list(
        study = studyName,
        context = context,
        trait = trait,
        method = method,
        dataType = dataType,
        standardized = weightStandardized
    ))
    groups <- .prdGroups(df)
    rows <- imap(
        groups,
        .prdFmRow,
        df = df,
        supplied = supplied,
        label = label,
        priorVariance = priorVariance,
        fits = src$fits
    )
    .withGenome(
        .prdBuildFineMapping(rows, df, context, trait, label),
        src$build %||% NA_character_
    )
}

# A fine-mapping result is a posterior, and the PIP is the posterior. A table
# without one is some other kind of file -- a weights export, most likely --
# and building a result from it would produce rows whose every accessor is
# empty.
# @noRd
.prdRequirePip <- function(df, label) {
    if (is_in("PIP", names(df))) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: no `PIP` column. A fine-mapping result is a posterior ",
        "over variants, so the inclusion probability is the one statistic ",
        "it cannot be read without. For a weights file use readTwasWeights()."
    ))
}

# One collection row from one group of published rows.
# @noRd
.prdFmRow <- function(idx, key, df, supplied, label, priorVariance, fits) {
    fields <- .prdRowFields(df, idx, supplied, label)
    tl <- .prdCanonical(slice(df, idx), fields$method, priorVariance)
    list(
        fields = fields,
        region = .prdRowRegion(df, idx, label),
        entry = FineMappingRow(
            variantIds = tl$variant_id,
            susieFit = .prdFitFor(fits, key),
            topLoci = tl
        )
    )
}

# The partial fit the VCF header recorded for this group, if any. A table with
# no identity axes is one unnamed group, so the key is an index rather than a
# key -- which must not be allowed to pick the first fit by position.
# @noRd
.prdFitFor <- function(fits, key) {
    if (!is.character(key) || !is_in(key, names(fits))) {
        return(NULL)
    }
    fits[[key]]
}

# Which class the file describes, and build it. A QTL result is keyed by
# (study, context, trait, method) and a GWAS result by (study, method) plus
# its block, so the presence of a context or a trait is what separates them.
# An empty file has no values to look at, only its header, so there the
# COLUMN decides.
# @noRd
.prdBuildFineMapping <- function(rows, df, context, trait, label) {
    regions <- .prdRegionColumn(map(rows, "region"))
    fields <- map(rows, "fields")
    if (length(rows) == 0L) {
        return(.prdEmptyFineMapping(df, context, trait))
    }
    if (.prdIsQtl(df, fields, context, trait)) {
        return(.prdBuildQtl(rows, fields, regions, label))
    }
    .prdRequireIdentity(fields, c("study", "method"), label)
    out <- GwasFineMappingResult(
        studyName = .prdRowColumn(fields, "study"),
        method = .prdRowColumn(fields, "method"),
        entry = map(rows, "entry"),
        blockId = .prdRowColumn(fields, "blockId"),
        region = regions,
        dataType = .prdRowColumn(fields, "dataType"),
        dataSource = .prdRowColumn(fields, "dataSource"),
        weightStandardized = .prdStandardizedColumn(fields)
    )
    .prdRestoreRegionGiven(out, rows)
}

# Whether each block range was given is what the file says. The constructor
# counts every range it is handed as given, but a range read from a file is
# the analysis window, which says nothing of where it came from, so a file
# that does not record it reads back as computed. A group spanning several
# chromosomes is split into a row per chromosome, each keeping its group's.
# @noRd
.prdRestoreRegionGiven <- function(out, rows) {
    recorded <- as.logical(
        map_chr(map(rows, "fields"), .prdRowField, nm = "regionGiven")
    )
    fromIdx <- .rtlSplitBySeqname(map(map(rows, "entry"), variants))$fromIdx
    out$regionGiven <- !is.na(recorded[fromIdx]) & recorded[fromIdx]
    out
}


# A header-only file reads back as an EMPTY collection of the class its
# columns describe, rather than as an error about identity axes it has no
# rows to record. An empty export is a legitimate output -- a region where
# every entry was screened out still writes its schema.
# @noRd
.prdEmptyFineMapping <- function(df, context, trait) {
    if (.prdIsQtl(df, list(), context, trait)) {
        return(QtlFineMappingResult(
            studyName = character(0),
            context = character(0),
            trait = character(0),
            method = character(0),
            entry = list()
        ))
    }
    GwasFineMappingResult(
        studyName = character(0),
        method = character(0),
        entry = list()
    )
}

# @noRd
.prdBuildQtl <- function(rows, fields, regions, label) {
    .prdRequireIdentity(
        fields,
        c("study", "context", "trait", "method"),
        label
    )
    QtlFineMappingResult(
        studyName = .prdRowColumn(fields, "study"),
        context = .prdRowColumn(fields, "context"),
        trait = .prdRowColumn(fields, "trait"),
        method = .prdRowColumn(fields, "method"),
        entry = map(rows, "entry"),
        jointStudies = .prdRowColumn(fields, "jointStudies"),
        jointContexts = .prdRowColumn(fields, "jointContexts"),
        jointTraits = .prdRowColumn(fields, "jointTraits"),
        region = regions,
        dataType = .prdRowColumn(fields, "dataType"),
        dataSource = .prdRowColumn(fields, "dataSource"),
        weightStandardized = .prdStandardizedColumn(fields)
    )
}

# @noRd
.prdIsQtl <- function(df, fields, context, trait) {
    if (!is.null(context) || !is.null(trait)) {
        return(TRUE)
    }
    if (nrow(df) == 0L) {
        return(any(is_in(c("context", "trait"), names(df))))
    }
    !is.null(.prdRowColumn(fields, "context")) ||
        !is.null(.prdRowColumn(fields, "trait"))
}

# @noRd
.prdStandardizedColumn <- function(fields) {
    v <- .prdRowColumn(fields, "standardized")
    if (is.null(v)) NULL else as.logical(v)
}

# ---- entry point: TWAS weights ----------------------------------------------

#' @title Read Published TWAS Weights
#' @description Rebuild a \code{TwasWeights} from a file
#'   \code{\link{writeBed}} or \code{\link{writeGwasVcf}}
#'   wrote, or from any file in the same vocabulary.
#'
#'   A merged table reads back as a multi-row collection, one row per
#'   \code{(study, context, trait, method)} tuple. A fitted model whose
#'   weights were a MATRIX -- a joint fit, one column per member -- was
#'   published as one block per member with \code{trait} naming it, so it
#'   reads back as one row per member with \code{joint_traits} keeping the
#'   membership. The weights are identical either way; only the packing
#'   differs.
#'
#' @section What is required, and what is not:
#'   Variant ids that encode coordinates (or the columns to build them from)
#'   and a \code{WGT} column -- a weight is the whole content of a
#'   \code{TwasWeights}. The four identity axes must be recorded or supplied,
#'   because the collection is keyed by them. Everything else is optional and
#'   carried through if present.
#'
#'   Cross-validation metrics are never published (they describe a model
#'   scored on held-out folds, which is not per-variant and not
#'   data-secure to express per sample), so a read-back collection has no
#'   \code{cvResult} and no fitted models.
#'
#' @inheritParams readFineMappingResult
#' @param weightStandardized Whether the published weights are on the
#'   standardized genotype scale. Required when the file does not record it:
#'   the scale is not recoverable from the weights themselves.
#' @return A \code{TwasWeights}.
#' @seealso \code{\link{writeBed}}, \code{\link{readFineMappingResult}}
#' @examples
#' data(twasWeightsExample)
#' path <- file.path(tempdir(), "weights.bed")
#' writeBed(twasWeightsExample, path)
#' readTwasWeights(path)

#' @export
readTwasWeights <- function(
    path,
    which = NULL,
    format = c("auto", "bed", "vcf"),
    studyName = NULL,
    context = NULL,
    trait = NULL,
    method = NULL,
    dataType = NULL,
    weightStandardized = NULL,
    sampleSelect = NULL
) {
    label <- "readTwasWeights"
    src <- .prdRead(path, which, format, sampleSelect, label)
    df <- src$table
    .prdRequireWgt(df, label)
    .prdRequireStandardized(df, weightStandardized, label)
    supplied <- compact(list(
        study = studyName,
        context = context,
        trait = trait,
        method = method,
        dataType = dataType,
        standardized = weightStandardized
    ))
    rows <- imap(
        .prdGroups(df),
        .prdTwRow,
        df = df,
        supplied = supplied,
        label = label
    )
    .withGenome(.prdTwFromRows(rows, label), src$build %||% NA_character_)
}

# @noRd
.prdRequireWgt <- function(df, label) {
    if (!is_in("WGT", names(df))) {
        abort(glue(
            "{label}: no `WGT` column. The weight is the whole content of a ",
            "TwasWeights. For a fine-mapping table use ",
            "readFineMappingResult()."
        ))
    }
}

# Assemble the TwasWeights collection from the per-group `.prdTwRow()`
# records; an empty file reads back as an empty collection.
# @noRd
.prdTwFromRows <- function(rows, label) {
    if (length(rows) == 0L) {
        return(TwasWeights(
            studyName = character(0),
            context = character(0),
            trait = character(0),
            method = character(0),
            entry = list()
        ))
    }
    fields <- map(rows, "fields")
    .prdRequireIdentity(
        fields,
        c("study", "context", "trait", "method"),
        label
    )
    TwasWeights(
        studyName = .prdRowColumn(fields, "study"),
        context = .prdRowColumn(fields, "context"),
        trait = .prdRowColumn(fields, "trait"),
        method = .prdRowColumn(fields, "method"),
        entry = map(rows, "entry"),
        jointStudies = .prdRowColumn(fields, "jointStudies"),
        jointContexts = .prdRowColumn(fields, "jointContexts"),
        jointTraits = .prdRowColumn(fields, "jointTraits"),
        region = .prdRegionColumn(map(rows, "region"))
    )
}

# One TwasWeightsRow from one group of published rows.
#
# The ids go through `.prdCanonical()` so they are derived exactly as the
# fine-mapping reader derives them. Variant-level columns BEYOND the weight
# are dropped: a TwasWeightsRow holds alleles and a weight, and the published
# `AF` came from the collection's LD panel, which is not in the file either.
# @noRd
.prdTwRow <- function(idx, key, df, supplied, label) {
    fields <- .prdRowFields(df, idx, supplied, label)
    sub <- slice(df, idx)
    tl <- .prdCanonical(sub, fields$method, NULL)
    list(
        fields = fields,
        region = .prdRowRegion(df, idx, label),
        entry = TwasWeightsRow(
            variantIds = tl$variant_id,
            weights = suppressWarnings(as.numeric(sub$WGT)),
            weightStandardized = isTRUE(.prdFlag(fields$standardized)),
            weightsDataType = .prdNaNull(fields$dataType)
        )
    )
}

# @noRd
.prdNaNull <- function(v) {
    if (is.null(v) || length(v) == 0L || is.na(v)) NULL else as.character(v)
}
