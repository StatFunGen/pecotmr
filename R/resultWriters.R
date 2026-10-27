# =============================================================================
# Published result writers
# -----------------------------------------------------------------------------
# Every pecotmr output format that leaves the package as a FILE lives here,
# behind TWO generics rather than one entry point per class:
#
#   writeGwasVcf(obj, filename, ...) GWAS-VCF-style: one record per variant,
#                               one SAMPLE per study / fitted model, the
#                               statistics as FORMAT fields. GwasSumStats,
#                               every FineMappingResultBase, TwasWeights.
#                               VariantAnnotation (a Suggests) writes the file.
#   writeBed(obj, filename, ...) the published long-format tables, one row per
#                               published observation, bgzipped + tabix-indexed
#                               on request. FineMappingResultBase,
#                               TwasWeights.
#
# One file because they are one problem -- project an S4 collection onto a
# FIXED published schema, join the posterior view to the marginal one, sort
# genomically, compress and index -- and while the VCF and the BED writers were
# separate each kept its own copy of the column pick, the row narrowing, the
# identity-key join and the fullFit column scan, which had already drifted
# apart.
#
# Five layers, in order:
#   .pub*  the shared layer: column pick, identity keys, row narrowing, the
#          posterior/marginal join, genomic order, the table writer,
#          output-path decoration
#   .vcf*  the VariantAnnotation-backed VCF / BCF writer
#   .twv*  the TWAS-weight VCF projection (weights as per-model samples)
#   .fmb*  the fine-mapping BED schemas (topLoci, fsusie)
#   .twb*  the TWAS-weight BED schemas (weights, performance)
# =============================================================================

#' @include AllGenerics.R AllClasses.R gwasSumStats.R twasWeights.R
#' @importFrom S4Vectors DataFrame SimpleList mcols
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges DataFrameList CharacterList
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#' @importFrom Rsamtools asBcf bgzip indexTabix
#' @importFrom tools file_ext
#' @importFrom stringr str_sub
#' @importFrom purrr map map_chr keep set_names list_flatten
NULL

# ---- shared: schema assembly ------------------------------------------------

# The identity axes a view is split on when several collection rows are folded
# into one published table. `blockId` is deliberately absent: the fSuSiE views
# do not carry it, and both sides of a join must key on the same axes.
# @noRd
.PUB_ROW_COLS <- c("study", "context", "trait", "method")

# The identity axes a per-VARIANT join keys on -- the row identity plus the
# variant itself.
# @noRd
.PUB_VARIANT_COLS <- c(
    "study",
    "context",
    "trait",
    "blockId",
    "method",
    "variant_id"
)

# Column `col` of `df`, or NA when the view does not carry it. Published
# schemas are fixed, so a column a particular input happens not to have becomes
# an NA cell rather than a missing column. `n` makes that NA a full-length
# vector, which the VCF FORMAT matrices need; the BED writers let the table
# recycle a scalar instead.
# @noRd
.pubCol <- function(df, col, n = NULL) {
    if (!is.null(df) && is_in(col, names(df))) {
        return(df[[col]])
    }
    if (is.null(n)) NA else rep(NA, n)
}

# The row-identity key a published view joins on: the identity columns the
# frame actually carries, "\r"-joined. A bare single-entry view carries none of
# them, so the key degenerates to whatever is present (the variant alone),
# which is exactly the right join for one entry.
# @noRd
.pubIdentityKey <- function(df, cols) {
    present <- intersect(cols, names(df))
    if (length(present) == 0L) {
        return(rep("", nrow(df)))
    }
    exec(str_c, !!!map(as.list(df[present]), .pubKeyPart), sep = "\r")
}

# The same key read off a COLLECTION's tuple columns, so a long-format view's
# rows can be matched back to the collection row they came from.
# `.tupleColumn()` not `[[`: on a RangedTupleList `[[` extracts an ELEMENT,
# while the identity axes live in mcols.
# @noRd
.pubCollectionKey <- function(x, cols) {
    present <- intersect(cols, .tupleColumnNames(x))
    exec(
        str_c,
        !!!map(map(present, .pubTupleColumn, x = x), .pubKeyPart),
        sep = "\r"
    )
}

# @noRd
.pubTupleColumn <- function(nm, x) as.character(.tupleColumn(x, nm))

# One component of an identity key. str_c() propagates NA, so a single absent
# axis (a GWAS result has no trait) would make the WHOLE key NA, losing the
# axes that ARE set -- and match() pairs NA with NA, so two unrelated rows
# would silently join. Spelling the absent value keeps NA matching NA within a
# column, which is what an absent axis means, without erasing the rest.
# @noRd
.pubKeyPart <- function(v) {
    chr <- as.character(v)
    if_else(is.na(chr), "<NA>", chr)
}

# "chr22_10000000_18000000" -- the region_ID publication tables key on.
# @noRd
.pubRegionIds <- function(gr) {
    if (!methods::is(gr, "GRanges") || length(gr) == 0L) {
        return(character(0))
    }
    str_c(as.character(seqnames(gr)), start(gr), end(gr), sep = "_")
}

# A vector packed into one cell -- the published spelling for a column holding
# a curve, a grid or a set of names.
# @noRd
.pubPack <- function(z) str_flatten(as.character(z), collapse = ";")

# The per-CS variant-level columns a fullFit topLoci build adds: the
# `within_cs_pip` scalar plus the wide within_cs_pip_<lab> / cs_logbf_<lab> /
# cs_effect_<lab> sets. Data-driven in both the VCF and the BED writer, so the
# pattern lives once.
# @noRd
.pubFullFitColumns <- function(df) {
    names(df)[str_detect(names(df), "^(within_cs_pip|cs_logbf_|cs_effect_)")]
}

# Align a marginalEffects view to a posterior (topLoci) view, row for row, so
# the two can be read as one table. Both sides are keyed on (row identity,
# variant): across a whole collection that keeps two rows contributing the same
# variant apart, and for a single entry -- whose views carry no identity
# columns -- it degenerates to a variant_id join.
#
# `[` rather than slice(): a posterior variant with no marginal record must
# produce an NA row, and slice() rejects the NA index that expresses it.
# @noRd
.pubJoinMarginals <- function(base, marginals) {
    idx <- match(
        .pubIdentityKey(base, .PUB_VARIANT_COLS),
        .pubIdentityKey(marginals, .PUB_VARIANT_COLS)
    )
    marginals[idx, , drop = FALSE]
}

# Fold a list of collections into one, choosing the class's own combiner. The
# `list` methods exist so a genome's worth of per-region RDSes can be exported
# in one call -- the combiners also handle the single-element case, so there is
# no separate concatenate-and-sort step.
# @noRd
.pubCombineList <- function(x) {
    if (length(x) == 0L) {
        abort("`x` must contain at least one result to export.")
    }
    if (methods::is(x[[1L]], "TwasWeights")) {
        return(combineTwasWeights(x))
    }
    if (methods::is(x[[1L]], "FineMappingResultBase")) {
        return(combineFineMappingResults(x))
    }
    abort(glue(
        "`x` must be a TwasWeights, QtlFineMappingResult or ",
        "GwasFineMappingResult object, or a list of one of those; got a ",
        "list of {class(x[[1L]])[[1L]]}."
    ))
}

# Row indices whose identity axes equal every supplied selector. Axes the
# collection does not carry are ignored, so one call addresses a GWAS result
# (no context, no trait) and a QTL one alike.
# @noRd
.pubMatchingRows <- function(x, studyName, context, trait, method) {
    cols <- .tupleColumnNames(x)
    selectors <- compact(list(
        study = studyName,
        context = if (is_in("context", cols)) context,
        trait = if (is_in("trait", cols)) trait,
        method = method
    ))
    reduce(
        names(selectors),
        .pubNarrowRows,
        x = x,
        selectors = selectors,
        .init = seq_len(nrow(x))
    )
}

# Narrow a row-index vector to those rows whose `axis` column equals the
# caller's requested value. `.tupleColumn()` not `[[`: on a RangedTupleList
# `[[` extracts an ELEMENT, while the identity axes live in mcols.
# @noRd
.pubNarrowRows <- function(rows, axis, x, selectors) {
    rows[as.character(.tupleColumn(x, axis))[rows] == selectors[[axis]]]
}

# Which of a source view's OWN columns to append after the published schema.
# The published columns always lead -- that is what makes the format
# publishable -- so this only ever ADDS: NULL adds nothing, "all" adds every
# column the view carries that the schema did not already map, and a character
# vector adds exactly those.
#
# `mapped` names the source columns the schema already publishes (under their
# published spellings) plus the GRanges range columns, which duplicate the
# position columns. Without that, "all" would write the same value twice under
# two names, and an explicitly named one would collide outright.
# @noRd
.pubViewColumns <- function(view, columns, mapped) {
    if (is.null(columns)) {
        return(character(0))
    }
    available <- setdiff(names(view), mapped)
    if (identical(columns, "all")) {
        return(available)
    }
    assertCharacter(columns, any.missing = FALSE, .var.name = "columns")
    .pubCheckViewColumns(columns, view, available, mapped)
    columns
}

# @noRd
.pubCheckViewColumns <- function(columns, view, available, mapped) {
    clash <- intersect(columns, mapped)
    if (length(clash) > 0L) {
        abort(glue(
            "`columns` names {str_flatten(clash, ', ')}, which the published ",
            "schema already writes; drop them."
        ))
    }
    unknown <- setdiff(columns, names(view))
    if (length(unknown) > 0L) {
        abort(glue(
            "`columns` names {str_flatten(unknown, ', ')}, which this view ",
            "does not carry. Available: {str_flatten(available, ', ')}."
        ))
    }
    invisible(NULL)
}

# Refuse to export a collection whose rows lack a FULL-DATA fit.
#
# Cross-validation alone is not publishable. A CV-only run -- which is what
# twasWeightsPipeline() produces under `CrossValidationParam(fitFullData =
# FALSE)`, legal only WITH cross-validation -- stores each row as an empty
# element with no weights and nothing but its `cvResult`. All it holds is
# fold-level results, and those cannot go in a published file: fold
# membership is sample-level information and sample names can be protected.
# So an export needs the full-data fit, and its absence is an ERROR rather
# than a silently empty output file.
#
# The check is on the OBJECT, not on the flag: `fitFullData` is a run setting
# that may not have been the thing that emptied a row, and an export
# has to be safe whatever produced its input.
#
# It also keeps genomicRegion() honest. That derives a row's analysis window
# from range(x) when no `region` column was recorded, and a variant-less row
# contributes no range -- so the derived vector would be SHORTER than the
# collection and indexing it by row number would read a neighbour's window.
# @noRd
.pubRequireFullFit <- function(x, incomplete, what) {
    if (length(incomplete) == 0L) {
        return(invisible(x))
    }
    abort(glue(
        "{what}: {length(incomplete)} of {nrow(x)} row(s) carry no ",
        "full-data fit, only cross-validation results -- ",
        "{str_flatten(.pubRowLabels(x, incomplete), ', ')}. Fold-level ",
        "results cannot be published (fold membership is sample-level ",
        "information), so refit with ",
        "`CrossValidationParam(fitFullData = TRUE)` and export that. ",
        "The cross-validation record stays available through cvResult()."
    ))
}

# Rows of a TwasWeights with no full-data weights: a CV-only row has an empty
# element AND no `weight` metadata column (TwasWeightsRow() drops the column
# when `weights` is NULL).
# @noRd
.twbIncompleteRows <- function(x) {
    which(
        lengths(x) == 0L | !map_lgl(seq_len(nrow(x)), .twbRowHasWeights, x = x)
    )
}

# @noRd
.twbRowHasWeights <- function(i, x) !is.null(mcols(x[[i]])$weight)

# Rows of a fine-mapping result with no fitted variants: what a
# cross-validation only run (CrossValidationParam(fitFullData = FALSE)) leaves
# in each row.
# @noRd
.fmbIncompleteRows <- function(x) {
    which(lengths(x) == 0L)
}

# Identity labels for the named rows, for an error message.
# @noRd
.pubRowLabels <- function(x, rows) {
    keys <- .pubCollectionKey(x, .PUB_ROW_COLS)
    str_replace_all(keys[rows], "\r", "|")
}

# ---- shared: genomic order --------------------------------------------------

#' Order rows by genomic position, with chromosomes sorted naturally
#'
#' Returns the permutation that puts rows in genomic order: chromosomes
#' 1..22 numerically, then the non-numeric contigs (X, Y, MT, ...)
#' alphabetically after all of them, and positions numerically within a
#' chromosome. A plain character sort puts chr10 before chr2, which is why
#' published tables need this.
#'
#' Shared by every pecotmr writer and by the pipeline wrappers that merge
#' per-region tables, so one definition of "genomic order" covers both.
#'
#' @param chrom Chromosome labels. A leading \code{"chr"} is optional and is
#'   ignored for ordering, so \code{"chr2"} and \code{"2"} sort alike.
#' @param start Start positions; coerced to numeric, non-numeric entries sort
#'   last within their chromosome.
#' @param end Optional second position key, used only to break ties on
#'   \code{start} (a half-open BED's end column). Ignored when it is not
#'   numeric, so a table whose third column is not a position (an allele, say)
#'   can pass it harmlessly.
#' @return An integer permutation, as from \code{order()}.
#' @examples
#' genomicOrder(c("chr10", "chr2", "chrX", "chr2"), c(5, 90, 1, 10))
#' @export
genomicOrder <- function(chrom, start, end = NULL) {
    rank <- .chromRank(chrom)
    pos <- suppressWarnings(as.numeric(start))
    tie <- suppressWarnings(as.numeric(end))
    if (is.null(end) || all(is.na(tie))) {
        return(order(rank, as.character(chrom), pos))
    }
    order(rank, as.character(chrom), pos, tie)
}

# The sort key behind genomicOrder(): numeric contigs keep their number, the
# rest sort alphabetically among themselves after every numeric one. A key for
# ONE vector, not a global scale -- the non-numeric contigs are ranked among
# themselves, so values are comparable only within a single call (which is all
# order() needs).
# @noRd
.chromRank <- function(chrom) {
    bare <- str_remove(as.character(chrom), "^chr")
    num <- suppressWarnings(as.numeric(bare))
    alpha <- sort(unique(bare[is.na(num)]))
    if_else(is.na(num), 1000 + match(bare, alpha), num)
}

# ---- shared: the table writer -----------------------------------------------

# Append caller-supplied trailing columns. Publication tables sometimes carry
# a field no S4 object can know -- the legacy top-loci BED's `block_id` is
# parsed from the RDS FILENAME -- so the writer takes them rather than forcing
# the caller to re-open and rewrite the file it just produced.
# @noRd
.pubAppendExtra <- function(out, extraColumns) {
    if (is.null(extraColumns) || length(extraColumns) == 0L) {
        return(out)
    }
    if (!is.list(extraColumns) || !is_named(extraColumns)) {
        abort("`extraColumns` must be a NAMED list of columns.")
    }
    clash <- intersect(names(extraColumns), names(out))
    if (length(clash) > 0L) {
        abort(glue(
            "`extraColumns` would overwrite the published column(s) ",
            "{str_flatten(clash, ', ')}; rename them."
        ))
    }
    reduce(
        names(extraColumns),
        .pubAppendOneExtra,
        extra = extraColumns,
        .init = out
    )
}

# @noRd
.pubAppendOneExtra <- function(out, nm, extra) {
    v <- extra[[nm]]
    if (!is_in(length(v), c(1L, nrow(out)))) {
        abort(glue(
            "`extraColumns${nm}` has length {length(v)}; it must be 1 ",
            "(recycled) or {nrow(out)} (one value per row)."
        ))
    }
    out[[nm]] <- if (nrow(out) == 0L) v[0L] else v
    out
}

# Write a table as a published BED: plain when `path` has no .gz/.bgz
# extension, bgzipped + tabix-indexed when it has. Rsamtools rather than a
# shell bgzip/tabix, so no pipeline needs an htslib binary on PATH.
#
# The builders already order their columns to the schema, so this takes the
# frame as given; it supplies the header for an empty result and drops the
# columns that had nothing to say.
#
# utils::write.table(), not vroom_write(): the two differ in how they format a
# float (1e-06 vs 1e-6) and in how many digits they keep.
# `emptySchema` is the header an empty result still documents -- the schema of
# whichever table was being written, since the fsusie table's vocabulary is not
# the per-variant one.
# @noRd
.pubWriteTable <- function(
    out,
    path,
    naString,
    index,
    emptySchema = .pubSchema(character(0)),
    keepColumns = character(0),
    extraColumns = NULL,
    build = NA_character_
) {
    if (is.null(out) || nrow(out) == 0L) {
        # A header-only file documents the whole vocabulary rather than
        # collapsing to the required core, which is what dropping empty
        # columns would do when every column is vacuously absent.
        out <- as_tibble(
            set_names(
                rep(list(character()), length(emptySchema)),
                emptySchema
            ),
            .name_repair = "minimal"
        )
    } else {
        out <- .pubDropEmpty(out, keepColumns)
    }
    out <- .pubAppendExtra(out, extraColumns)
    if (nrow(out) > 0L) {
        out <- slice(out, genomicOrder(out[["#chr"]], out[["start"]]))
    }
    .pubEmitTable(out, path, naString, index, build)
}

# Materialise the ordered table at `path`, compressing and indexing it when the
# extension asks for it. A named genome build goes first, as a `##genome=`
# line above the column header; tabix skips both as comments.
# @noRd
.pubEmitTable <- function(out, path, naString, index, build = NA_character_) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    compress <- str_detect(path, "\\.(gz|bgz)$")
    plain <- if (compress) str_remove(path, "\\.(gz|bgz)$") else path
    con <- file(plain, "w")
    if (!is.na(build)) {
        writeLines(str_c("##genome=", build), con)
    }
    utils::write.table(
        out,
        con,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = naString
    )
    close(con)
    if (!compress) {
        return(invisible(path))
    }
    bgzip(plain, dest = path, overwrite = TRUE)
    file.remove(plain)
    # An empty table has no records to index, and indexTabix() errors on one.
    if (index && nrow(out) > 0L) {
        indexTabix(path, format = "bed")
    }
    invisible(path)
}

# ---- shared: output-path decoration -----------------------------------------

# Decorate `outputPath` with the spec's context / trait / block tags when split
# flags are set. Preserves the file extension. Examples: "out.vcf" +
# (context="brain") -> "out.brain.vcf" "out.vcf.bgz" + (context="brain",
# trait="ENSG1") -> "out.brain.ENSG1.vcf.bgz"
# @noRd
.decorateOutputPath <- function(
    outputPath,
    spec,
    splitByContext,
    splitByTrait,
    splitByBlock = FALSE
) {
    if (!isTRUE(splitByContext) && !isTRUE(splitByTrait) && !splitByBlock) {
        return(outputPath)
    }
    ext <- str_to_lower(tools::file_ext(outputPath))
    composite <- ext == "bgz" || ext == "gz"
    base <- if (composite) {
        str_remove(
            outputPath,
            regex("\\.[^.]+\\.(bgz|gz)$", ignore_case = TRUE)
        )
    } else {
        tools::file_path_sans_ext(outputPath)
    }
    ext_keep <- str_sub(
        outputPath,
        str_length(base) + 1L,
        str_length(outputPath)
    )
    tags <- c(
        if (isTRUE(splitByContext)) .pubPathTag(spec$context),
        if (isTRUE(splitByTrait)) .pubPathTag(spec$trait),
        if (splitByBlock) .pubPathTag(spec$block)
    )
    if (length(tags) == 0L) {
        return(outputPath)
    }
    str_c(base, ".", str_flatten(tags, "."), ext_keep)
}

# One path tag, or NULL when the axis carries no usable value (a context-less
# result split by context keeps its undecorated path).
# @noRd
.pubPathTag <- function(v) {
    if (is.null(v) || is.na(v) || str_length(v) == 0L) {
        return(NULL)
    }
    v
}

# ---- VCF: entry points ------------------------------------------------------

#' @rdname writeGwasVcf
#' @export
setMethod(
    "writeGwasVcf",
    signature("list"),
    function(obj, filename, sampleName = NULL, ...) {
        writeGwasVcf(
            .pubCombineList(obj),
            filename,
            sampleName = sampleName,
            ...
        )
    }
)

#' @rdname writeGwasVcf
#' @export
setMethod(
    "writeGwasVcf",
    signature("GwasSumStats"),
    function(obj, filename, sampleName = NULL, studyName = NULL) {
        .vcfRequireVariantAnnotation()
        studyName <- .vcfResolveStudy(obj, studyName)
        ss <- sumStats(obj, studyName = studyName)
        mc <- mcols(ss)
        .writeVcfImpl(
            chrom = as.character(seqnames(ss)),
            pos = start(ss),
            ref = mc$A2,
            alt = mc$A1,
            snpIds = mc$SNP,
            geno = .vcfSumstatsGeno(mc, length(ss)),
            genoHeader = .vcfSumstatsGenoHeader(),
            sampleName = sampleName %||% studyName,
            outputPath = filename,
            build = .oneBuild(obj)
        )
    }
)

# VariantAnnotation writes the file; it is a Suggests, so the writer checks.
# @noRd
.vcfRequireVariantAnnotation <- function() {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        abort("Package 'VariantAnnotation' is required for writeGwasVcf().")
    }
    invisible(NULL)
}

# Select which study to write (a GwasSumStats can hold many); a single-study
# collection defaults to its one study.
# @noRd
.vcfResolveStudy <- function(x, studyName) {
    if (!is.null(studyName)) {
        return(studyName)
    }
    # Studies, not rows: one study is several rows once split by chromosome
    # or LD block, and sumStats() stitches them back together.
    studies <- unique(as.character(x$study))
    if (length(studies) != 1L) {
        msg <- glue(
            "This GwasSumStats has {length(studies)} studies. ",
            "Pass `studyName = <name>` to select one."
        )
        abort(msg)
    }
    studies
}

# Per-sample geno matrices for the GWAS-VCF fields the mcols can fill.
#
# Three released bugs are fixed here, all of them the same mistake -- a
# reserved key carrying the wrong quantity:
#   * ES held the Z-SCORE. The spec defines ES as "effect size estimate
#     relative to the alternative allele", and defines EZ for a z "provided if
#     it was used to derive the ES and SE fields" -- which is exactly our case.
#   * SS is not a GWAS-VCF key at all. The spec's sample-size field is NS.
#   * AF held MAF. The spec's AF is the ALTERNATIVE allele frequency, and a
#     minor-allele frequency cannot be converted to one -- it has lost which
#     allele it refers to -- so it goes under our own MAF key instead of being
#     mislabelled.
# @noRd
.vcfSumstatsGeno <- function(mc, nSnps) {
    compact(list(
        ES = .vcfSumstatsField(mc, "BETA", nSnps),
        SE = .vcfSumstatsField(mc, "SE", nSnps),
        EZ = .vcfSuppliedZ(mc, nSnps),
        NS = .vcfSumstatsField(mc, "N", nSnps, as.integer),
        AF = .vcfSumstatsField(mc, "AF", nSnps),
        MAF = .vcfSumstatsField(mc, "MAF", nSnps)
    ))
}

# `EZ` only where the SOURCE supplied the z. The spec defines the field as the
# z "provided if it was used to derive the ES and SE fields", so publishing a
# Wald z that pecotmr computed from a recorded beta/se inverts its meaning --
# a reader would then believe ES and SE were derived from it. The loaders
# record `z_supplied` for exactly this; a collection built by hand without it
# is treated as not having supplied one, which is the safe direction.
# @noRd
.vcfSuppliedZ <- function(mc, nSnps) {
    z <- .vcfSumstatsField(mc, "Z", nSnps)
    if (is.null(z)) {
        return(NULL)
    }
    supplied <- as.logical(.pubCol(mc, "z_supplied", nSnps))
    if (!any(isTRUE(supplied) | (!is.na(supplied) & supplied))) {
        return(NULL)
    }
    matrix(replace(as.numeric(z), !supplied | is.na(supplied), NA_real_), nSnps)
}

# @noRd
.vcfSumstatsField <- function(mc, nm, nSnps, coerce = as.numeric) {
    if (!is_in(nm, colnames(mc))) {
        return(NULL)
    }
    matrix(coerce(mc[[nm]]), nSnps)
}

# The FORMAT header for the sumstats geno fields, with the spec's own wording
# for the keys it defines.
# @noRd
.vcfSumstatsGenoHeader <- function() {
    DataFrame(
        Number = rep("A", 6L),
        Type = c("Float", "Float", "Float", "Integer", "Float", "Float"),
        Description = c(
            "Effect size estimate relative to the alternative allele",
            "Standard error of effect size estimate",
            str_c(
                "Z-score provided if it was used to derive the ES and SE ",
                "fields"
            ),
            str_c(
                "Variant-specific number of samples/individuals with called ",
                "genotypes used to test association with specified trait"
            ),
            "Alternative allele frequency in trait subset",
            "Minor allele frequency (not the spec's AF: allele unknown)"
        ),
        row.names = c("ES", "SE", "EZ", "NS", "AF", "MAF")
    )
}

#' @rdname writeGwasVcf
#' @export
setMethod(
    "writeGwasVcf",
    signature("FineMappingResultBase"),
    function(
        obj,
        filename,
        sampleName = NULL,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        splitByContext = FALSE,
        splitByTrait = FALSE,
        splitByBlock = FALSE
    ) {
        .vcfRequireVariantAnnotation()
        rows <- .pubMatchingRows(obj, studyName, context, trait, method)
        .fitCoverageCheckedRows(obj[rows], "writeGwasVcf()")
        # Resolve what goes in each FILE. Without a split flag every matching
        # fit is a SAMPLE of one file; with one, each context, trait or LD
        # block gets its own file at a decorated path.
        groups <- .resolveFineMappingGroups(
            obj,
            outputPath = filename,
            rows = rows,
            splitByContext = splitByContext,
            splitByTrait = splitByTrait,
            splitByBlock = splitByBlock
        )
        invisible(map_chr(
            groups,
            .fmvWriteGroup,
            x = obj,
            sampleName = sampleName
        ))
    }
)

# Resolve what each OUTPUT FILE holds. Without a split flag every matching fit
# becomes a SAMPLE of one multi-sample VCF -- which is GWAS-VCF's own shape,
# where the samples are the studies -- so comparing two methods over a region
# means reading one file rather than two. A GWAS fit is every LD block of its
# (study, method), so its blocks share one sample. `splitByContext` /
# `splitByTrait` emit one file per context / trait and `splitByBlock` one per
# block, each at a decorated path; the fits that land on one path share it.
# @noRd
.resolveFineMappingGroups <- function(
    x,
    outputPath,
    rows,
    splitByContext,
    splitByTrait,
    splitByBlock
) {
    if (length(rows) == 0L) {
        abort("writeGwasVcf: no rows match the supplied selectors.")
    }
    .fmvCheckSplitByBlock(x, splitByBlock)
    specs <- map(
        .fmvFitRows(x, rows),
        .fitSpec,
        x = x,
        blockIds = .fmvBlockIds(x)
    )
    if (isTRUE(splitByBlock)) {
        specs <- list_flatten(map(specs, .fmvSpecPerBlock))
    }
    if (!isTRUE(splitByContext) && !isTRUE(splitByTrait) && !splitByBlock) {
        return(list(list(path = outputPath, specs = specs)))
    }
    paths <- map_chr(
        specs,
        .decorateOutputPath,
        outputPath = outputPath,
        splitByContext = splitByContext,
        splitByTrait = splitByTrait,
        splitByBlock = splitByBlock
    )
    unname(imap(
        split(specs, factor(paths, levels = unique(paths))),
        .fmvPathGroup
    ))
}

# @noRd
.fmvPathGroup <- function(specs, path) {
    list(path = path, specs = unname(specs))
}

# Blocks are a GWAS notion: a QTL fit is one row, so there is nothing to split.
# @noRd
.fmvCheckSplitByBlock <- function(x, splitByBlock) {
    assertFlag(splitByBlock)
    if (!splitByBlock || methods::is(x, "GwasFineMappingResult")) {
        return(invisible(NULL))
    }
    abort(glue(
        "writeGwasVcf: `splitByBlock` splits a GWAS fit into its LD blocks; a ",
        "{class(x)[[1L]]} fit is one row. Use `splitByContext` / ",
        "`splitByTrait` instead."
    ))
}

# The matching rows grouped into fits, in first-appearance order. A GWAS fit is
# every LD block of one (study, method); a QTL row is a fit of its own, its
# joint-membership columns being part of what identifies it.
# @noRd
.fmvFitRows <- function(x, rows) {
    if (!methods::is(x, "GwasFineMappingResult")) {
        return(as.list(rows))
    }
    key <- str_c(
        .fmrIdentityColumn(x, "study")[rows],
        .fmrIdentityColumn(x, "method")[rows],
        sep = "\r"
    )
    unname(split(rows, factor(key, levels = unique(key))))
}

# One spec per block of a fit, each carrying the block it is for.
# @noRd
.fmvSpecPerBlock <- function(spec) {
    map(unique(spec$blocks), .fmvBlockSpec, spec = spec)
}

# @noRd
.fmvBlockSpec <- function(block, spec) {
    on <- spec$blocks == block
    list_assign(
        spec,
        rows = spec$rows[on],
        blocks = spec$blocks[on],
        block = block
    )
}

# The LD block each row of a GWAS result covers: its `blockId`, which keys the
# block manifest, else its own range.
# @noRd
.fmvBlockIds <- function(x) {
    ids <- .fmrIdentityColumn(x, "blockId")
    if_else(is.na(ids), .rtlRangeKeys(x), ids)
}

# Build a (study, context, trait, method) spec for one fit: its rows -- one,
# or one per LD block of a GWAS fit -- their blocks, and the identity they
# share.
# @noRd
.fitSpec <- function(rows, x, blockIds) {
    r <- rows[[1L]]
    list(
        rows = rows,
        blocks = blockIds[rows],
        study = .fmrIdentityAt(x, "study", r),
        context = .fmrIdentityAt(x, "context", r),
        trait = .fmrIdentityAt(x, "trait", r),
        method = .fmrIdentityAt(x, "method", r)
    )
}

# ---- VCF: FORMAT field specs ------------------------------------------------

# One VCF FORMAT field spec: name, per-variant `values`, VCF `type`
# (Float/Integer), and header `desc`. Number is always "A" (per-allele).
# @noRd
.vcfSpec <- function(name, values, type, desc) {
    list(name = name, values = values, type = type, desc = desc)
}

# The per-variant FORMAT fields, named exactly as the BED columns are: the
# spec's keys for the quantities it defines, ours alongside them for the rest.
#
# `ES` is the MARGINAL effect and nothing else. It used to hold the posterior
# conditional effect when the fit reported one and the marginal beta otherwise,
# so a consumer could not tell which quantity it had. The posterior effect now
# has its own field -- `WGT`, which is also the prediction weight, because a
# susie posterior mean IS the TWAS coefficient.
# @noRd
.vcfCoreSpecs <- function(base, m, nSnps) {
    c(.vcfMarginalSpecs(m, nSnps), .vcfFitSpecs(base, nSnps))
}

# The marginal-association FORMAT fields, from the joined marginals `m`.
# @noRd
.vcfMarginalSpecs <- function(m, nSnps) {
    p <- suppressWarnings(as.numeric(.pubCol(m, "p", nSnps)))
    list(
        .vcfSpec(
            "ES",
            .pubCol(m, "beta", nSnps),
            "Float",
            "Effect size estimate relative to the alternative allele"
        ),
        .vcfSpec(
            "SE",
            .pubCol(m, "se", nSnps),
            "Float",
            "Standard error of effect size estimate"
        ),
        .vcfSpec(
            "EZ",
            .pubEz(m, nSnps),
            "Float",
            str_c(
                "Z-score provided if it was used to derive the ES and SE ",
                "fields"
            )
        ),
        .vcfSpec(
            "LP",
            if_else(is.na(p) | p <= 0, NA_real_, -log10(p)),
            "Float",
            "-log10 p-value for effect estimate"
        ),
        .vcfSpec(
            "NS",
            as.integer(.pubCol(m, "N", nSnps)),
            "Integer",
            str_c(
                "Variant-specific number of samples/individuals with called ",
                "genotypes used to test association with specified trait"
            )
        )
    )
}

# The per-variant FORMAT fields read off the fit itself (`base`): allele
# frequency and the posterior effect / prediction weight.
# @noRd
.vcfFitSpecs <- function(base, nSnps) {
    list(
        .vcfSpec(
            "AF",
            .pubCol(base, "af", nSnps),
            "Float",
            "Alternative allele frequency in trait subset"
        ),
        .vcfSpec(
            "WGT",
            .pubCol(base, "posterior_effect", nSnps),
            "Float",
            "Posterior effect estimate, i.e. the prediction weight"
        ),
        .vcfSpec(
            "WGTSE",
            .pubCol(base, "posterior_se", nSnps),
            "Float",
            "Standard error of the posterior effect estimate"
        ),
        .vcfSpec(
            "THETA",
            .pubCol(base, "theta", nSnps),
            "Float",
            str_c(
                "Infinitesimal component of WGT (susieInf / susieAsh); ",
                "already included in WGT, 0 for methods without one"
            )
        )
    )
}

# Credible-set FORMAT fields. pecotmr emits CS<coverage> (+ PUR<coverage> when
# a purity column exists) for EVERY cs_<coverage> column present, so coverage
# levels are data-driven rather than the 50/70/95 xqtl-protocol default.
# @noRd
.vcfCsSpecs <- function(base) {
    csCols <- names(base)[str_detect(names(base), "^cs_[0-9.]+$")]
    perCol <- map(csCols, .vcfCsColSpecs, base = base)
    list_flatten(perCol)
}

# Per-CS variant-level fullFit FORMAT fields (within_cs_pip scalar plus the
# wide within_cs_pip_<lab> / cs_logbf_<lab> / cs_effect_<lab> sets emitted when
# the topLoci was built with fullFit = TRUE). Each column becomes an uppercase
# field.
# @noRd
.vcfFullFitSpecs <- function(base) {
    map(.pubFullFitColumns(base), .vcfFullFitColSpec, base = base)
}

# Posterior FORMAT fields, only when a fine-mapping posterior table exists:
# PIP, per-variant log Bayes factor, local false sign rate, plus the dynamic
# credible-set and fullFit column sets.
# @noRd
.vcfPosteriorSpecs <- function(base, nSnps) {
    fixed <- list(
        .vcfSpec(
            "PIP",
            .pubCol(base, "pip", nSnps),
            "Float",
            "Posterior inclusion probability"
        ),
        .vcfSpec(
            "LBF",
            .pubCol(base, "logBF", nSnps),
            "Float",
            "Per-variant log Bayes factor (max single effect)"
        ),
        .vcfSpec(
            "LFSR",
            .pubCol(base, "lfsr", nSnps),
            "Float",
            "Local false sign rate (per-condition posterior)"
        )
    )
    c(fixed, .vcfCsSpecs(base), .vcfFullFitSpecs(base))
}

# FORMAT specs for one cs_<coverage> column: a CS index field plus a PUR purity
# field when a matching *_purity column exists.
# @noRd
.vcfCsColSpecs <- function(cc, base) {
    cov <- str_remove(cc, "^cs_")
    csSpec <- .vcfSpec(
        str_c("CS", cov),
        .pubCsIndex(base[[cc]]),
        "Integer",
        glue("Credible-set index at {cov}% coverage")
    )
    pc <- str_c(cc, "_purity")
    if (!is_in(pc, names(base))) {
        return(list(csSpec))
    }
    pur <- suppressWarnings(as.numeric(base[[pc]]))
    purSpec <- .vcfSpec(
        str_c("PUR", cov),
        pur,
        "Float",
        glue("Purity (min abs corr) of the {cov}% credible set")
    )
    list(csSpec, purSpec)
}

# The FORMAT spec for one fullFit posterior column (uppercased Float field).
# @noRd
.vcfFullFitColSpec <- function(cc, base) {
    v <- suppressWarnings(as.numeric(base[[cc]]))
    .vcfSpec(
        str_to_upper(cc),
        v,
        "Float",
        glue("Per-credible-set variant statistic ({cc})")
    )
}

# TRUE when a FORMAT spec has at least one non-NA value (drop all-NA specs).
# @noRd
.vcfSpecHasValues <- function(s) {
    any(!is.na(s$values))
}

# ---- VCF: fine-mapping assembly ---------------------------------------------

# Resolve the per-variant body: the fine-mapping POSTERIOR (topLoci at
# signalCutoff 0 = every fitted variant, carrying pip / CS membership /
# conditional effect) joined to the MARGINAL univariate sumstats
# (marginalEffects) where those exist. mvSuSiE / fSuSiE have no marginal
# sumstats, so the posterior drives the variant set; univariate susie adds
# marginal ES=beta / SE / LP / AF on top. Returns list(base, m, hasPost).
# @noRd
.vcfResolveBody <- function(entry, sn) {
    post <- try_fetch(
        as_tibble(as.data.frame(
            S4Vectors::mcols(topLoci(entry, signalCutoff = 0))
        )),
        error = function(cnd) NULL
    )
    marg <- try_fetch(
        as_tibble(as.data.frame(
            S4Vectors::mcols(marginalEffects(entry)),
            optional = TRUE
        )),
        error = function(cnd) NULL
    )
    hasPost <- !is.null(post) && nrow(post) > 0L
    hasMarg <- !is.null(marg) && nrow(marg) > 0L
    if (hasPost) {
        base <- post
        m <- if (hasMarg) .pubJoinMarginals(base, marg) else NULL
    } else if (hasMarg) {
        base <- marg
        m <- marg
    } else {
        msg <- glue("writeGwasVcf: entry [{sn}] has no variants to write")
        abort(msg)
    }
    list(base = base, m = m, hasPost = hasPost)
}

# The geno header DataFrame for a set of FORMAT specs.
# @noRd
.vcfGenoHeader <- function(specs) {
    DataFrame(
        Number = rep("A", length(specs)),
        Type = map_chr(specs, "type"),
        Description = map_chr(specs, "desc"),
        row.names = map_chr(specs, "name")
    )
}

# Sample names: the caller's, or the derived ones. A caller-supplied vector has
# to name EVERY sample -- the samples are what the file is keyed on, so a
# partial naming would silently mislabel the rest.
# @noRd
.vcfSampleNames <- function(derived, sampleName) {
    if (is.null(sampleName)) {
        return(derived)
    }
    if (length(sampleName) != length(derived)) {
        abort(glue(
            "`sampleName` has length {length(sampleName)} but ",
            "{length(derived)} model(s) will be written; pass one name per ",
            "model or NULL to derive them."
        ))
    }
    as.character(sampleName)
}

# Write one file's worth of fine-mapping fits: the records are the UNION of
# their variants, each fit is a sample, and a fit that did not include a
# variant leaves its cells empty.
# @noRd
.fmvWriteGroup <- function(group, x, sampleName) {
    # Every row's window at once: a genome-wide study's blocks join one
    # sample, and asking per block would rebuild them all for each.
    rows <- map(group$specs, .fmvRowSpecs, x = x, rg = genomicRegion(x))
    sn <- .vcfSampleNames(map_chr(rows, "sampleName"), sampleName)
    variants <- .fmvVariantUnion(rows)
    specs <- keep(
        .fmvFormatSpecs(rows, variants$variant_id),
        .vcfSpecHasValues
    )
    .writeVcfImpl(
        chrom = variants$chrom,
        pos = variants$pos,
        ref = variants$A2,
        alt = variants$A1,
        snpIds = variants$variant_id,
        geno = set_names(map(specs, "values"), map_chr(specs, "name")),
        genoHeader = .vcfGenoHeader(specs),
        sampleName = sn,
        outputPath = group$path,
        sampleMeta = .fmvSampleMeta(rows, sn),
        build = .oneBuild(x)
    )
    group$path
}

# One fit's body, its FORMAT specs, and the variants they are aligned to.
# @noRd
.fmvRowSpecs <- function(spec, x, rg = genomicRegion(x)) {
    sn <- .fmvSampleLabel(spec)
    blocks <- map2(
        spec$rows,
        spec$blocks,
        .fmvBlockBody,
        x = x,
        sn = sn,
        rg = rg
    )
    base <- bind_rows(map(blocks, "base"))
    nSnps <- nrow(base)
    list(
        sampleName = sn,
        spec = spec,
        perFit = .fmvJoinPerFit(map(blocks, "perFit")),
        base = base,
        nVariants = nSnps,
        variantIds = as.character(.pubCol(base, "variant_id", nSnps)),
        formats = c(
            .vcfCoreSpecs(base, .fmvJoinMarginals(blocks), nSnps),
            if (any(map_lgl(blocks, "hasPost"))) {
                .vcfPosteriorSpecs(base, nSnps)
            } else {
                list()
            },
            .fmvBlockIdSpec(blocks)
        )
    )
}

# One row's body and its per-fit header values. A GWAS fit fine-mapped block
# by block is several rows, and each contributes one of these to the sample.
# @noRd
.fmvBlockBody <- function(r, block, x, sn, rg) {
    entry <- x[r]
    list_assign(
        .vcfResolveBody(entry, sn),
        perFit = .pubPerFitFields(
            x,
            r,
            .fmrPartsSusieFit(.fmrRowParts(x, r)),
            rg
        ),
        block = block
    )
}

# The marginals of every block, row-aligned with the stacked bases. A block
# with none still holds its rows, as NA; with none anywhere the table has no
# columns, which the marginal fields read as missing, as they read NULL.
# @noRd
.fmvJoinMarginals <- function(blocks) {
    bind_rows(map2(map(blocks, "m"), blocks, .fmvMarginalsOrEmpty))
}

# @noRd
.fmvMarginalsOrEmpty <- function(m, block) {
    m %||% tibble(.rows = nrow(block$base))
}

# The BID field: each variant's LD block, written only when a sample joins
# several, since credible-set indices are numbered within a block.
# @noRd
.fmvBlockIdSpec <- function(blocks) {
    if (length(blocks) < 2L) {
        return(list())
    }
    list(.vcfSpec(
        "BID",
        list_c(map(blocks, .fmvBlockColumn)),
        "String",
        "LD block of the variant; credible sets are numbered within a block"
    ))
}

# @noRd
.fmvBlockColumn <- function(block) {
    rep(block$block, nrow(block$base))
}

# The per-fit header fields that differ between the blocks of one fit, packed
# in block order with "|" when a sample joins several.
# @noRd
.FMV_BLOCK_FIELDS <- c(
    "block_id",
    "region_id",
    "region_start",
    "region_end",
    "region_given",
    "V",
    "cs_index"
)

# The header values of a fit spread over several LD blocks: the per-block
# fields packed, the rest -- which describe the fit -- taken once, after
# checking every block agrees on them.
# @noRd
.fmvJoinPerFit <- function(perFits) {
    if (length(perFits) == 1L) {
        return(perFits[[1L]])
    }
    fields <- names(perFits[[1L]])
    set_names(map(fields, .fmvJoinField, perFits = perFits), fields)
}

# @noRd
.fmvJoinField <- function(nm, perFits) {
    vals <- map_chr(perFits, .fmvFieldText, nm = nm)
    if (all(is.na(vals))) {
        return(NA_character_)
    }
    if (is_in(nm, .FMV_BLOCK_FIELDS)) {
        return(str_flatten(replace_na(vals, ""), "|"))
    }
    if (n_distinct(vals) > 1L) {
        abort(glue(
            "writeGwasVcf: the LD blocks of one fit disagree on `{nm}` ",
            "({str_flatten(unique(vals), ', ')}), so they cannot be written ",
            "as one sample. Pass `splitByBlock = TRUE`."
        ))
    }
    perFits[[1L]][[nm]]
}

# @noRd
.fmvFieldText <- function(perFit, nm) {
    as.character(perFit[[nm]] %||% NA_character_)
}

# @noRd
.fmvSampleLabel <- function(spec) {
    as.character(glue(
        "{spec$study}|{spec$context %||% '_'}|",
        "{spec$trait %||% '_'}|{spec$method}"
    ))
}

# The union of every sample's variants, in genomic order, carrying the alleles
# the first fit to report each one recorded.
# @noRd
.fmvVariantUnion <- function(rows) {
    tbl <- bind_rows(map(rows, .fmvRowVariants))
    tbl <- slice(tbl, which(!duplicated(tbl$variant_id)))
    slice(tbl, genomicOrder(tbl$chrom, tbl$pos))
}

# @noRd
.fmvRowVariants <- function(r) {
    n <- r$nVariants
    tibble(
        variant_id = r$variantIds,
        chrom = as.character(.pubCol(r$base, "chrom", n)),
        pos = .pubCol(r$base, "pos", n),
        A1 = as.character(.pubCol(r$base, "A1", n)),
        A2 = as.character(.pubCol(r$base, "A2", n))
    )
}

# The multi-sample FORMAT specs: one per field name ANY sample reported, its
# values a nVariant x nSample matrix placed onto the union by variant id. The
# field set is the union because which posterior fields a fit reports depends
# on the method, and one file needs one FORMAT schema across its samples.
# @noRd
.fmvFormatSpecs <- function(rows, unionIds) {
    fields <- unique(list_c(c(
        list(character(0)),
        map(rows, .fmvRowFormatNames)
    )))
    map(fields, .fmvFormatSpec, rows = rows, unionIds = unionIds)
}

# @noRd
.fmvRowFormatNames <- function(r) map_chr(r$formats, "name")

# @noRd
.fmvFormatSpec <- function(nm, rows, unionIds) {
    template <- compact(map(rows, .fmvFindFormat, nm = nm))[[1L]]
    cols <- map(
        rows,
        .fmvFormatColumn,
        nm = nm,
        unionIds = unionIds,
        type = template$type
    )
    .vcfSpec(nm, exec(cbind, !!!cols), template$type, template$desc)
}

# @noRd
.fmvFindFormat <- function(r, nm) detect(r$formats, .fmvFormatNamed, nm = nm)

# @noRd
.fmvFormatNamed <- function(f, nm) identical(f$name, nm)

# One sample's column of a FORMAT field. The NA fill matches the field's
# declared VCF type, so an Integer field does not become a Float the moment one
# sample is missing a variant.
# @noRd
.fmvFormatColumn <- function(r, nm, unionIds, type) {
    fill <- switch(
        type,
        Integer = NA_integer_,
        String = NA_character_,
        NA_real_
    )
    out <- rep(fill, length(unionIds))
    f <- .fmvFindFormat(r, nm)
    if (is.null(f)) {
        return(out)
    }
    out[match(r$variantIds, unionIds)] <- f$values
    out
}

# The ##SAMPLE header block for a fine-mapping VCF: each fit's identity tuple
# and how many variants it contributed, so the file says which fits it holds
# without reparsing the sample names.
# @noRd
.fmvSampleMeta <- function(rows, sampleNames) {
    tbl <- bind_rows(map(rows, .fmvSampleMetaRow))
    .pubSampleMetaFrame(tbl, sampleNames)
}

# The `##SAMPLE` key for each per-fit field.
#
# The DATA vocabulary is shared between the two formats, but a structured VCF
# header key cannot be: Rsamtools' header parser splits the `<...>` on
# `,(?=[[:alnum:]]+=)`, so a key containing an underscore is not recognised as
# a key at all -- it glues onto the previous field's value, and everything
# after it on the line is lost. Six of these fields were unreadable for that
# reason. CamelCase is also what the GWAS-VCF spec's own SAMPLE keys use
# (`TotalVariants`, `StudyType`), so this follows the format's convention
# rather than departing from ours. This table is the only place the two
# spellings meet, and `.prdVcfSampleMeta()` reads them back through it.
# @noRd
.PUB_SAMPLE_KEYS <- c(
    total_variants = "TotalVariants",
    study = "Study",
    context = "Context",
    trait = "Trait",
    method = "Method",
    block_id = "BlockId",
    data_type = "DataType",
    data_source = "DataSource",
    weight_standardized = "WeightStandardized",
    joint_studies = "JointStudies",
    joint_contexts = "JointContexts",
    joint_traits = "JointTraits",
    region_id = "RegionId",
    region_start = "RegionStart",
    region_end = "RegionEnd",
    region_given = "RegionGiven",
    V = "V",
    cs_index = "CsIndex"
)

# The ##SAMPLE block as a DataFrame, with the fields that are empty for EVERY
# sample dropped -- the same rule the BED columns follow, so a header is not
# padded with `BlockId=.,JointStudies=.,V=.` for an input that has none.
# @noRd
.pubSampleMetaFrame <- function(tbl, sampleNames) {
    kept <- names(tbl)[!map_lgl(tbl, .pubAllAbsent)]
    df <- as.data.frame(mutate(
        select(tbl, all_of(kept)),
        across(everything(), .vcfHeaderValue)
    ))
    names(df) <- .pubSampleKeys(names(df))
    `rownames<-`(as(df, "DataFrame"), .vcfHeaderValue(sampleNames))
}

# Field names as `##SAMPLE` keys. An unmapped name is a new per-fit field that
# has not been given a key: letting it through would write a line that reads
# back short by that field and every field after it, so it stops here.
# @noRd
.pubSampleKeys <- function(nms) {
    absent <- nms[!is_in(nms, names(.PUB_SAMPLE_KEYS))]
    if (length(absent) > 0L) {
        abort(glue(
            "Per-fit field(s) {str_flatten(absent, ', ')} have no ",
            "##SAMPLE key. Add them to .PUB_SAMPLE_KEYS: an underscored key ",
            "is silently unreadable."
        ))
    }
    unname(.PUB_SAMPLE_KEYS[nms])
}

# @noRd
.fmvSampleMetaRow <- function(r) {
    # As text: a sample joining several LD blocks packs its per-block fields
    # into one string, and one header holds both kinds of sample.
    bind_cols(
        tibble(total_variants = r$nVariants),
        as_tibble(map(r$perFit, as.character))
    )
}

# The per-fit fields a `##SAMPLE` line carries, named exactly as the BED's
# per-fit COLUMNS are -- one vocabulary across both formats. The header line
# spells them in CamelCase, because a VCF structured-header key must be
# alphanumeric to be parseable; `.PUB_SAMPLE_KEYS` is that one translation.
#
# `V` and `cs_index` are the reason this block is load-bearing rather than
# decoration: they are per-EFFECT, so no per-variant column can hold them, and
# colocPipeline's row filter reads them off the fit -- `sets$cs_index` or `V`.
# Without them a round-tripped LBF matrix could be read but not filtered.
# @noRd
.pubPerFitFields <- function(x, row, fit, rg = genomicRegion(x)) {
    inRange <- !is.null(row) && length(rg) >= row
    md <- mcols(x)
    list(
        study = .pubRowValue(md, "study", row),
        context = .pubRowValue(md, "context", row),
        trait = .pubRowValue(md, "trait", row),
        method = .pubRowValue(md, "method", row),
        block_id = .pubRowValue(md, "blockId", row),
        data_type = .pubRowValue(md, "dataType", row),
        data_source = .pubDataSource(x),
        joint_studies = .pubRowValue(md, "jointStudies", row),
        joint_contexts = .pubRowValue(md, "jointContexts", row),
        joint_traits = .pubRowValue(md, "jointTraits", row),
        weight_standardized = .pubRowValue(md, "standardized", row),
        region_id = if (inRange) .pubRegionIds(rg)[[row]] else NA_character_,
        region_start = if (inRange) start(rg)[[row]] else NA_integer_,
        region_end = if (inRange) end(rg)[[row]] else NA_integer_,
        region_given = .pubRowValue(md, "regionGiven", row),
        V = .pubPackOrNa(.fitField(fit, "V")),
        cs_index = .pubPackOrNa(.fitField(fit, "sets", "cs_index"))
    )
}

# One mcols cell as a string. `[[row]]` then pack, rather than
# as.character() on the column: several of these are List columns (a
# TwasWeights `dataType` is one) and as.character() on a CompressedList
# errors.
# @noRd
.pubRowValue <- function(md, nm, row) {
    if (is.null(row) || !is_in(nm, colnames(md))) {
        return(NA_character_)
    }
    cell <- md[[nm]][[row]]
    if (length(cell) == 0L) NA_character_ else .pubPack(cell)
}

# @noRd
.pubPackOrNa <- function(v) {
    if (length(v) == 0L) NA_character_ else .pubPack(v)
}

# ---- VCF: TWAS weights ------------------------------------------------------

#' @rdname writeGwasVcf
#' @export
setMethod(
    "writeGwasVcf",
    signature("TwasWeights"),
    function(
        obj,
        filename,
        sampleName = NULL,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        .vcfRequireVariantAnnotation()
        .pubRequireFullFit(obj, .twbIncompleteRows(obj), "writeGwasVcf")
        rows <- .pubMatchingRows(obj, studyName, context, trait, method)
        if (length(rows) == 0L) {
            abort("writeGwasVcf: no rows match the supplied selectors.")
        }
        # The same rule the BED path applies, and for the same reason: a
        # fine-mapping method's "weights" are the posterior re-expressed, so
        # the fine-mapping result says strictly more. This path used to publish
        # them, which made the two formats disagree about what a weights file
        # contains.
        rows <- intersect(rows, .twbExportableRows(obj, "writeGwasVcf"))
        if (length(rows) == 0L) {
            abort(str_c(
                "writeGwasVcf: every selected row is fine-mapping-derived, so ",
                "there is nothing to write as weights. Write the ",
                "fine-mapping result."
            ))
        }
        .twvWriteSamples(
            obj,
            .twvSampleSpecs(obj, rows),
            sampleName,
            filename
        )
    }
)

# Project the selected models onto one GWAS-VCF-style file: the records are the
# UNION of their variants, each model is a sample, and a model that did not
# weight a variant leaves that cell empty.
#
# Multi-sample rather than one file per model because that is what GWAS-VCF
# does with several studies, and because a region's TwasWeights routinely holds
# a dozen methods for one gene -- a dozen near-identical files would be the
# wrong shape for the same information.
# @noRd
.twvWriteSamples <- function(x, samples, sampleName, outputPath) {
    gr <- .twvVariantUnion(x, samples)
    ids <- .grVariantIds(gr)
    mc <- mcols(gr)
    wgt <- .twvWeightMatrix(x, samples, ids)
    meta <- .twbRowMetadata(x)
    sn <- .twvSampleNames(samples, meta, sampleName)
    .writeVcfImpl(
        chrom = as.character(seqnames(gr)),
        pos = start(gr),
        ref = as.character(mc$A2),
        alt = as.character(mc$A1),
        snpIds = ids,
        geno = list(WGT = wgt),
        genoHeader = .twvGenoHeader(),
        sampleName = sn,
        outputPath = outputPath,
        sampleMeta = .twvSampleMeta(
            samples,
            sn,
            meta,
            colSums(!is.na(wgt)),
            x
        ),
        build = .oneBuild(x)
    )
}

# The ##SAMPLE header block: what GWAS-VCF records per study, here per fitted
# MODEL -- its identity tuple, how many variants it actually weighted, and the
# scale its weights are on.
#
# Cross-validation metrics are deliberately NOT here. They describe how a model
# scored on held-out FOLDS, and a published file has no data-secure way to say
# anything fold-level: fold membership is sample-level information and sample
# names can be protected. The export carries the full-data fit and nothing
# else; callers who want the CV record read cvResult() and tabulate it
# themselves.
# @noRd
.twvSampleMeta <- function(samples, sampleNames, meta, nWeighted, x) {
    tbl <- bind_rows(map(
        seq_along(samples),
        .twvSampleMetaRow,
        samples = samples,
        meta = meta,
        nWeighted = nWeighted,
        x = x
    ))
    .pubSampleMetaFrame(tbl, sampleNames)
}

# @noRd
.twvSampleMetaRow <- function(j, samples, meta, nWeighted, x) {
    s <- samples[[j]]
    ident <- slice(meta, s$row)
    fields <- .pubPerFitFields(x, s$row, NULL)
    # A weight-matrix column is a member of a joint fit, so `trait` names the
    # member while `joint_traits` keeps the membership. Assigned rather than
    # bound, so there is one `trait` field and not a repaired duplicate.
    fields$trait <- .twbTraitOf(s$label, ident$trait)
    bind_cols(tibble(total_variants = nWeighted[[j]]), as_tibble(fields))
}

# One ##SAMPLE field value. These are written unquoted, so a comma, an equals
# sign, an angle bracket or whitespace would break the line's key=value
# structure; an absent value is spelled "." as elsewhere in VCF.
# @noRd
.vcfHeaderValue <- function(v) {
    chr <- str_replace_all(as.character(v), "[,=<>[:space:]]+", "_")
    if_else(is.na(chr), ".", chr)
}

# @noRd
.twvGenoHeader <- function() {
    DataFrame(
        Number = "A",
        Type = "Float",
        Description = str_c(
            "TWAS prediction weight for the effect (ALT) allele; empty where ",
            "the model did not weight the variant"
        ),
        row.names = "WGT"
    )
}

# One sample spec per (collection row, weight column): the row it reads, which
# weight column, and that column's label. A plain weight vector contributes one
# sample; a weight MATRIX (fSuSiE collapses a functional fit to variants x
# features) contributes one per column.
# @noRd
.twvSampleSpecs <- function(x, rows) {
    list_flatten(map(rows, .twvRowSampleSpecs, x = x))
}

# @noRd
.twvRowSampleSpecs <- function(i, x) {
    # .pubRequireFullFit() has already established that every row has weights.
    labels <- .twbWeightLabels(mcols(x[[i]])$weight)
    map(seq_along(labels), .twvOneSampleSpec, row = i, labels = labels)
}

# @noRd
.twvOneSampleSpec <- function(k, row, labels) {
    list(row = row, column = k, label = labels[[k]])
}

# The union of every selected model's variants, in genomic order. The weight
# mcols is dropped first: it is numeric on one row and a matrix on another, and
# c() on GRanges requires one mcols schema.
# @noRd
.twvVariantUnion <- function(x, samples) {
    rows <- unique(map_int(samples, "row"))
    all <- exec(c, !!!map(rows, .twvRowAlleles, x = x))
    gr <- all[!duplicated(.grVariantIds(all))]
    gr[genomicOrder(as.character(seqnames(gr)), start(gr))]
}

# @noRd
.twvRowAlleles <- function(i, x) {
    gr <- x[[i]]
    S4Vectors::`mcols<-`(
        gr,
        value = mcols(gr)[, c("A1", "A2"), drop = FALSE]
    )
}

# The nVariant x nSample weight matrix: one column per model, placed onto the
# union by variant id.
# @noRd
.twvWeightMatrix <- function(x, samples, unionIds) {
    exec(
        cbind,
        !!!map(samples, .twvSampleColumn, x = x, unionIds = unionIds)
    )
}

# @noRd
.twvSampleColumn <- function(s, x, unionIds) {
    gr <- x[[s$row]]
    out <- rep(NA_real_, length(unionIds))
    out[match(.grVariantIds(gr), unionIds)] <-
        .twbWeightColumn(s$column, mcols(gr)$weight)
    out
}

# Sample names: the caller's, or the model's identity tuple -- with the weight
# column's label appended when the row carries a weight matrix, since those
# columns are separate samples and must have separate names.
# @noRd
.twvSampleNames <- function(samples, meta, sampleName) {
    .vcfSampleNames(
        map_chr(samples, .twvSampleLabel, meta = meta),
        sampleName
    )
}

# @noRd
.twvSampleLabel <- function(s, meta) {
    ident <- slice(meta, s$row)
    tags <- c(ident$study, ident$context, ident$trait, ident$method)
    if (!is.na(s$label)) {
        tags <- c(tags, s$label)
    }
    str_flatten(if_else(is.na(tags), "_", tags), "|")
}

# ---- VCF: the VariantAnnotation writer --------------------------------------

# Internal implementation shared by all methods. `sampleName` may name several
# samples, in which case every geno matrix is nVariant x nSample.
# `sampleMeta` is an optional per-sample DataFrame written as ##SAMPLE header
# lines (GWAS-VCF's per-study metadata block).
# @noRd
.writeVcfImpl <- function(
    chrom,
    pos,
    ref,
    alt,
    snpIds,
    geno,
    genoHeader,
    sampleName,
    outputPath,
    sampleMeta = NULL,
    build = NA_character_
) {
    # UCSC chromosome names (chr1, chrX, chrM), as the package's ranges use.
    chrom <- withChrPrefix(chrom)
    gr <- .vcfBuildRanges(chrom, pos, ref, alt, snpIds)
    hdr <- .vcfBuildHeader(
        sampleName,
        genoHeader,
        geno,
        sampleMeta,
        build,
        unique(chrom)
    )
    vcf <- .vcfBuildObject(gr, sampleName, hdr, geno, ref, alt)
    .vcfWriteByFormat(vcf, outputPath, chrom)
    invisible(outputPath)
}

# Per-variant rowRanges (end extended to the wider of ref / alt for indels).
# @noRd
.vcfBuildRanges <- function(chrom, pos, ref, alt, snpIds) {
    GRanges(
        chrom,
        IRanges(
            start = as.integer(pos),
            end = as.integer(pos) + pmax(str_length(ref), str_length(alt)) - 1L,
            names = snpIds
        )
    )
}

# VCF header (VCFv4.2, the genome build when the object names one, the
# optional ##SAMPLE block, and the geno header subset to the fields actually
# present). The build is written as ##reference and as the assembly of each
# ##contig, the two places a VCF reader looks for it.
# @noRd
.vcfBuildHeader <- function(
    sampleName,
    genoHeader,
    geno,
    sampleMeta = NULL,
    build = NA_character_,
    contigs = character(0)
) {
    named <- !is.na(build)
    sections <- compact(list(
        fileformat = DataFrame(Value = "VCFv4.2", row.names = "fileformat"),
        reference = if (named) {
            DataFrame(Value = build, row.names = "reference")
        },
        contig = if (named && length(contigs) > 0L) {
            DataFrame(
                assembly = rep(build, length(contigs)),
                row.names = contigs
            )
        },
        SAMPLE = sampleMeta
    ))
    hdr <- VariantAnnotation::VCFHeader(
        header = exec(DataFrameList, !!!sections),
        sample = sampleName
    )
    VariantAnnotation::`geno<-`(
        hdr,
        value = genoHeader[
            is_in(rownames(genoHeader), names(geno)),
            ,
            drop = FALSE
        ]
    )
}

# Assemble + finalize the VCF object (ref / alt / FILTER, sorted).
# @noRd
.vcfBuildObject <- function(gr, sampleName, hdr, geno, ref, alt) {
    vcf <- VariantAnnotation::VCF(
        rowRanges = gr,
        colData = DataFrame(Samples = sampleName, row.names = sampleName),
        exptData = list(header = hdr),
        geno = SimpleList(geno)
    )
    alleled <- VariantAnnotation::`alt<-`(
        VariantAnnotation::`ref<-`(vcf, value = DNAStringSet(ref)),
        value = .vcfAltList(alt)
    )
    filtered <- VariantAnnotation::`fixed<-`(
        alleled,
        value = `[[<-`(
            VariantAnnotation::fixed(alleled),
            "FILTER",
            value = "PASS"
        )
    )
    sort(filtered)
}

# The ALT column's container. A DNAStringSetList is the normal one, but ALT is
# the one allele field the VCF spec lets be non-DNA: "*" for a variant spanned
# by an upstream deletion, or a symbolic <...> allele. DNAStringSet rejects
# those outright ("key 42 (char '*') not in lookup table"), so such a file
# carries its ALT as a CharacterList. REF is DNA by the spec and unaffected.
# @noRd
.vcfAltList <- function(alt) {
    if (isTRUE(all(str_detect(alt, "^[ACGTNacgtn]+$")))) {
        return(DNAStringSetList(as.list(alt)))
    }
    CharacterList(as.list(alt))
}

# Write the VCF in the format implied by the output extension.
# VariantAnnotation::writeVcf() appends
# ".bgz" when index = TRUE, so paths are passed WITHOUT the .bgz/.gz suffix.
# @noRd
.vcfWriteByFormat <- function(vcf, outputPath, chrom) {
    ext <- file_ext(outputPath)
    if (ext == "bcf") {
        .vcfWriteBcf(vcf, outputPath, chrom)
    } else if (ext == "gz" || ext == "bgz") {
        .vcfWriteBgz(vcf, outputPath)
    } else {
        VariantAnnotation::writeVcf(vcf, outputPath)
    }
    invisible(NULL)
}

# BCF path: write a temporary bgzipped VCF, then convert to BCF via asBcf.
# @noRd
#' @importFrom rlang try_fetch
.vcfWriteBcf <- function(vcf, outputPath, chrom) {
    tmpVcfStem <- tempfile(fileext = ".vcf")
    tmpVcfBgz <- str_c(tmpVcfStem, ".bgz")
    on.exit(
        unlink(c(tmpVcfBgz, str_c(tmpVcfBgz, ".tbi")), force = TRUE),
        add = TRUE
    )
    VariantAnnotation::writeVcf(vcf, tmpVcfStem, index = TRUE)
    # asBcf appends ".bcf" to destination, so strip the extension.
    # Rsamtools disabled asBcf() (>= 2.26 raises "temporarily disabled"), so
    # the bare upstream error is translated into something actionable rather
    # than surfacing as an opaque failure from a documented output format.
    try_fetch(
        asBcf(
            tmpVcfBgz,
            dictionary = unique(chrom),
            destination = str_remove(outputPath, "\\.bcf$")
        ),
        error = function(cnd) {
            msg <- glue(
                "writeGwasVcf: BCF output needs a working ",
                "Rsamtools::asBcf(), which the installed Rsamtools does not ",
                "provide. Write a bgzipped VCF instead by giving the output ",
                "path a .vcf.bgz extension."
            )
            abort(msg, parent = cnd)
        }
    )
}

# bgzip / gzip path: VariantAnnotation::writeVcf() always creates .bgz;
# renamed to .gz on request.
# @noRd
.vcfWriteBgz <- function(vcf, outputPath) {
    vcfStem <- str_remove(outputPath, "\\.(bgz|gz)$")
    VariantAnnotation::writeVcf(vcf, vcfStem, index = TRUE)
    actualPath <- str_c(vcfStem, ".bgz")
    if (actualPath != outputPath && file.exists(actualPath)) {
        file.rename(actualPath, outputPath)
        tbiActual <- str_c(actualPath, ".tbi")
        if (file.exists(tbiActual)) {
            file.rename(tbiActual, str_c(outputPath, ".tbi"))
        }
    }
}

# ---- BED: the fSuSiE schema -------------------------------------------------

# @noRd
.FMB_FSUSIE_COLUMNS <- c(
    "#chr",
    "start",
    "end",
    "ID",
    "REF",
    "ALT",
    "AF",
    "PIP",
    ".cs.",
    "cs_id",
    "cs_root",
    "grid_resolution",
    "grid_positions",
    "grid_effects",
    "grid_band_halfwidth",
    "trait_positions",
    "trait_effects",
    "joint_traits",
    "study",
    "context",
    "trait",
    "method",
    "region_id",
    "region_start",
    "region_end"
)

# fSuSiE carries the SAME effect twice, at two resolutions, and the published
# table holds both: the fitted curve on the wavelet grid, and that curve
# reinterpolated onto the positions of the traits the fit was run over
# (linearly between grid points, constant past the ends). They are
# complementary -- one is the model, the other is what it predicts where the
# data actually sit -- so both are written by default and `levels` drops a
# group when only one is wanted.
#
# Dropping a group REMOVES its columns from the published order rather than
# re-ordering what is left, so a single-level table is a strict subset of the
# full one.
# @noRd
.FMB_FSUSIE_WAVELET_COLUMNS <- c(
    "grid_resolution",
    "grid_positions",
    "grid_effects",
    "grid_band_halfwidth"
)

# The trait-level group. `joint_traits` rather than a new `trait_names`: these
# ARE the joint members (fsusieTraitEffects() reads them from the row's
# `jointTraits` column), so they use the vocabulary the rest of the schema
# already has for joint membership. "epi_mark" also presumed epigenomic data,
# which a functional fit does not require.
# @noRd
.FMB_FSUSIE_TRAIT_COLUMNS <- c(
    "trait_positions",
    "trait_effects",
    "joint_traits"
)

# @noRd
.fmbFsusieColumns <- function(levels, csCols = character(0)) {
    at <- match(".cs.", .FMB_FSUSIE_COLUMNS)
    full <- c(
        .FMB_FSUSIE_COLUMNS[seq_len(at - 1L)],
        csCols,
        .FMB_FSUSIE_COLUMNS[seq(at + 1L, length(.FMB_FSUSIE_COLUMNS))]
    )
    setdiff(
        full,
        c(
            if (!is_in("wavelet", levels)) .FMB_FSUSIE_WAVELET_COLUMNS,
            if (!is_in("trait", levels)) .FMB_FSUSIE_TRAIT_COLUMNS
        )
    )
}

# The coverage that KEYS the rows. A functional fit's table is one row per
# (context, credible set), so one level has to define the row; the caller's
# first requested level wins, else the first the fit carries. Membership at
# every other level rides along as its own CS<n>/PUR<n> pair, so a secondary
# coverage is represented rather than silently dropped.
# @noRd
.fmbPrimaryCsLevel <- function(view, coverage) {
    levels <- .pubCsLevels(view, coverage)
    if (length(levels) == 0L) {
        abort(str_c(
            "writeBed(schema = 'fsusie'): the result carries no ",
            "credible-set column, so there are no rows to key on."
        ))
    }
    levels[[1L]]
}

# The source columns .FMB_FSUSIE_COLUMNS already publishes, plus the range
# columns an as.data.frame() adds. The credibleSets view carries a good deal
# more per-variant detail (beta / se / N / logBF / the other coverages and
# their purities), which `columns = "all"` appends.
# @noRd
.FMB_FSUSIE_SOURCE <- c(
    "seqnames",
    "start",
    "end",
    "width",
    "strand",
    "chrom",
    "pos",
    "variant_id",
    "A1",
    "A2",
    "af",
    "pip",
    "cs_95",
    "study",
    "context",
    "trait",
    "method"
)

# One (row, credible set) block of the fSuSiE table. `v` is that CS's
# credibleSets() rows; `band` / `marks` are the matching curve views, already
# narrowed to the same row.
#
# The "chr" prefix is kept, as in the other two schemas and in the protocol's
# own phenotype BEDs; stripping it here was unique to this table and made one
# reader rule impossible. AF, not MAF, is published: the same `af` the other
# table publishes. Writing pmin(af, 1 - af) here lost which allele the
# frequency referred to, so one export of an object was recoverable and the
# other was not. Membership is written at EVERY coverage the fit carries, not
# just the one keying the rows -- a secondary coverage would otherwise be
# silently dropped.
# @noRd
.fmbFsusieBlock <- function(
    v,
    band,
    marks,
    rowRegionIds,
    rowRegions,
    added,
    levels,
    primary,
    csCols
) {
    key <- .pubIdentityKey(v, .PUB_ROW_COLS)[[1L]]
    cs <- v[[str_c("cs_", primary)]][[1L]]
    context <- v$context[[1L]]
    b <- .fmbCurveRows(band, key, cs)
    if (nrow(b) == 0L) {
        return(NULL)
    }
    m <- .fmbCurveRows(marks, key, cs)
    i <- match(key, names(rowRegionIds))
    regionId <- rowRegionIds[[i]]
    rg <- rowRegions[i]
    csId <- str_c(context, regionId, str_remove(cs, "^fsusie_"), sep = ":")
    out <- tibble(
        `#chr` = as.character(v$chrom),
        start = v$pos - 1L,
        end = v$pos,
        ID = v$variant_id,
        REF = .pubCol(v, "A2"),
        ALT = .pubCol(v, "A1"),
        AF = .pubCol(v, "af"),
        PIP = v$pip,
        cs_id = csId,
        cs_root = csId,
        study = .pubCol(v, "study"),
        context = context,
        trait = .pubCol(v, "trait"),
        method = .pubCol(v, "method"),
        region_id = regionId,
        region_start = start(rg) - 1L,
        region_end = end(rg)
    )
    out <- bind_cols(out, .fmbFsusieGrid(b, m), .pubFrameCs(v, csCols))
    bind_cols(
        select(out, all_of(.fmbFsusieColumns(levels, csCols))),
        select(v, all_of(added))
    )
}

# The rows of a curve view (credible band or trait marks) belonging to one
# (row identity, credible set) pair.
# @noRd
.fmbCurveRows <- function(view, key, cs) {
    keep <- .pubIdentityKey(view, .PUB_ROW_COLS) == key & view$cs == cs
    view[keep, , drop = FALSE]
}

# The packed functional-grid columns of one fSuSiE block: the credible band
# `b` over the grid and the per-trait marks `m`.
# @noRd
.fmbFsusieGrid <- function(b, m) {
    tibble(
        grid_resolution = nrow(b),
        grid_positions = .pubPack(b$pos),
        grid_effects = .pubPack(b$effect),
        grid_band_halfwidth = .pubPack(b$halfwidth),
        trait_positions = .pubPack(m$traitPos),
        trait_effects = .pubPack(m$effect),
        joint_traits = .pubPack(m$jointTrait)
    )
}

# The fSuSiE BED table: one block per (row, credible set) of the fsusie rows.
#
# genomicRegion() is the window each row was ANALYSED over (the requested
# region, or the trait's cis window), which is what a publication table means
# by its region -- not the narrower span of whichever variants survived QC.
#
# Splitting the views on (context, cs) alone is only safe for a single
# region: once several regions are combined into one table, two regions can
# each contribute context "c1" / "fsusie_1", and they must stay separate rows
# with their own windows.
# @noRd
.fmbFsusieTable <- function(x, columns, levels, coverage) {
    fs <- x[as.character(.tupleColumn(x, "method")) == "fsusie"]
    if (nrow(fs) == 0L) {
        return(NULL)
    }
    band <- as.data.frame(
        S4Vectors::mcols(fsusieCredibleBand(fs)),
        optional = TRUE
    )
    if (nrow(band) == 0L) {
        abort(str_c(
            "writeBed(schema = 'fsusie'): the functional fields ",
            "are absent from every row. Refit with a retention level that ",
            "keeps the fSuSiE curves."
        ))
    }
    csTbl <- as.data.frame(S4Vectors::mcols(credibleSets(fs)), optional = TRUE)
    primary <- .fmbPrimaryCsLevel(csTbl, coverage)
    csCols <- .pubCsColumns(csTbl, coverage)
    marks <- as.data.frame(.fsusieTraitEffectsTable(fs))
    rowRegions <- genomicRegion(fs)
    rowRegionIds <- set_names(
        .pubRegionIds(rowRegions),
        .pubCollectionKey(fs, .PUB_ROW_COLS)
    )
    pieces <- compact(map(
        split(
            csTbl,
            list(
                .pubIdentityKey(csTbl, .PUB_ROW_COLS),
                csTbl[[str_c("cs_", primary)]]
            ),
            drop = TRUE
        ),
        .fmbFsusieBlock,
        band = band,
        marks = marks,
        rowRegionIds = rowRegionIds,
        rowRegions = rowRegions,
        added = .pubViewColumns(csTbl, columns, .FMB_FSUSIE_SOURCE),
        levels = levels,
        primary = primary,
        csCols = csCols
    ))
    if (length(pieces) == 0L) NULL else list_rbind(pieces)
}

# ---- BED: the topLoci schema ------------------------------------------------

# The ONE published vocabulary, shared by both formats.
#
# Fine-mapping results and TWAS weights are the same kind of record -- the
# fine-mapping posterior effect IS the prediction weight -- so they are one
# schema, not two. A file carries the subset of these columns it has anything
# to say about (see .pubDropEmpty), so a weights table is not padded with
# twenty empty posterior columns.
#
# Names are the GWAS-VCF FORMAT keys for every per-variant quantity the spec
# defines, our own uppercase keys alongside them for the quantities it does
# not, and snake_case for the per-fit fields -- which are BED columns here and
# `##SAMPLE` keys in the VCF, so they cannot inherit a spec name. One
# vocabulary across both formats; no translation layer to drift.
#
# The CS columns are NOT here: which coverages exist is fixed at FIT time by
# CredibleSetParam, so they are discovered per input (see .pubCsColumns) and
# spliced in at `.cs.` below.
# @noRd
.PUB_COLUMNS <- c(
    # Structural. `ID`/`REF`/`ALT` mirror the VCF's fixed columns -- REF/ALT
    # says which allele is which, where `a1`/`a2` left it to convention.
    "#chr",
    "start",
    "end",
    "ID",
    "REF",
    "ALT",
    # Frequency and sample size. `AF` is the ALT-allele frequency the spec
    # defines; `MAF` is ours, for a source that only recorded the minor
    # frequency -- which is not derivable to AF, since it has lost which
    # allele it refers to. `NS` is the spec's sample-size key (there is no
    # `SS` in GWAS-VCF).
    "AF",
    "MAF",
    "NS",
    # Marginal association. `ES`/`SE` are the MARGINAL effect and its error --
    # never a z-score, never silently the posterior. `EZ` carries z, and the
    # spec defines it as exactly that: the z "provided if it was used to
    # derive the ES and SE fields". `LP` is -log10 p, which is also the better
    # representation: invertible, and no underflow at tiny p.
    "ES",
    "SE",
    "EZ",
    "LP",
    # Posterior effect and its error. This IS the prediction weight: a susie
    # posterior mean is sum_l alpha_l * mu_l, which is the TWAS coefficient.
    # `conditional_effect` is gone entirely: it was `WGT / PIP` once the
    # posterior mean moved onto the scale `coef()` uses, so it duplicated two
    # columns that are both here.
    "WGT",
    "WGTSE",
    # The infinitesimal half of WGT. susieInf / susieAsh split the posterior
    # into a sparse and an infinitesimal term and coef.susie() sums them, so
    # WGT is right either way -- this is what makes the split visible, which
    # is otherwise unrecoverable from the published file. 0 elsewhere.
    "THETA",
    # Fine-mapping posterior. CS<n>/PUR<n> are spliced in after these.
    "PIP",
    "LBF",
    "LFSR",
    ".cs.",
    # Identity: the tuple the collections are keyed by.
    "study",
    "context",
    "trait",
    "method",
    "block_id",
    # Provenance. `data_type` is the assay/modality; `data_source` is
    # individual-level vs summary-statistics, which `data_type` was being
    # mistaken for and which is derivable from whether an ldReference is
    # present.
    "data_type",
    "data_source",
    "weight_standardized",
    # Joint-fit membership. The identity column carries the "joint" sentinel
    # and these list the members; the writers used to drop them entirely.
    "joint_studies",
    "joint_contexts",
    "joint_traits",
    # The window ANALYSED, without which a read-back result's genomicRegion()
    # silently degrades to the span of the variants that survived QC, and,
    # for a GWAS block, whether its range was given or is that span.
    "region_id",
    "region_start",
    "region_end",
    "region_given"
)

# Always written, so the output is always a valid BED even when every optional
# column is empty.
# @noRd
.PUB_REQUIRED_COLUMNS <- c("#chr", "start", "end", "ID")

# The credible-set column pair for every coverage the input actually carries.
# Which coverages exist is a FIT-time decision (CredibleSetParam), and
# .fmrRowCs() READS a stored `cs_<n>` column rather than recomputing, so a
# writer can restrict what it publishes but cannot invent a level the fit never
# computed. The VCF path already discovered these; the BED path hardcoded
# 95/70/50, which is the inconsistency this removes.
# @noRd
.pubCsColumns <- function(view, coverage = NULL) {
    levels <- .pubCsLevels(view, coverage)
    list_c(c(
        list(character(0)),
        map(levels, .pubCsPair)
    ))
}

# @noRd
.pubCsPair <- function(lv) c(str_c("CS", lv), str_c("PUR", lv))

# The available-levels phrase. `str_flatten(character(0))` is "", not NULL, so
# a `%||%` fallback never fires and the message read "Available: ."
# @noRd
.pubLevelList <- function(have) {
    if (length(have) == 0L) "<none>" else str_flatten(have, ", ")
}

# The coverage levels to publish: every `cs_<n>` column the view carries, or
# the requested subset -- with a requested level the fit never computed an
# error that names what IS available, rather than a silently empty column.
# @noRd
.pubCsLevels <- function(view, coverage = NULL) {
    have <- str_remove(
        names(view)[str_detect(names(view), "^cs_[0-9.]+$")],
        "^cs_"
    )
    if (is.null(coverage)) {
        return(have)
    }
    want <- as.character(as.numeric(coverage) * 100)
    absent <- setdiff(want, have)
    if (length(absent) > 0L) {
        abort(glue(
            "`coverage` asks for {str_flatten(absent, ', ')}, which this fit ",
            "did not compute; credible sets are fixed at fit time by ",
            "CredibleSetParam. Available: {.pubLevelList(have)}."
        ))
    }
    want
}

# Splice the discovered CS pairs into the schema at the `.cs.` placeholder.
# @noRd
.pubSchema <- function(csCols) {
    at <- match(".cs.", .PUB_COLUMNS)
    c(
        .PUB_COLUMNS[seq_len(at - 1L)],
        csCols,
        .PUB_COLUMNS[seq(at + 1L, length(.PUB_COLUMNS))]
    )
}

# Build one record kind's rows with every schema column present -- the ones it
# carries, and NA for the ones its class cannot know. .pubDropEmpty() then
# removes whatever stayed entirely NA, so the file shows a user only columns
# that say something.
# @noRd
.pubFrame <- function(cols, schema) {
    absent <- setdiff(schema, names(cols))
    full <- c(cols, set_names(rep(list(NA), length(absent)), absent))
    select(exec(tibble, !!!full), all_of(schema))
}

# Drop the columns that are entirely absent. A fixed schema padded with twenty
# empty columns is worse for a reader than a varying one, and nothing
# downstream breaks: `.pubCol()` already returns NA for a column a table does
# not have, so absent and all-NA are the same thing on the way back in.
# @noRd
.pubDropEmpty <- function(out, keepColumns = character(0)) {
    keep <- names(out)[
        is_in(names(out), c(.PUB_REQUIRED_COLUMNS, keepColumns)) |
            !map_lgl(out, .pubAllAbsent)
    ]
    select(out, all_of(keep))
}

# @noRd
.pubAllAbsent <- function(v) all(is.na(v))

# The columns the caller named explicitly. `"all"` is a request for whatever
# happens to be there, so those stay subject to the drop-empty rule; a named
# column is a request for that column, so it is written even if it is empty.
# @noRd
.pubNamedColumns <- function(columns) {
    if (is.null(columns) || identical(columns, "all")) {
        return(character(0))
    }
    as.character(columns)
}

# @noRd
.FMB_TOPLOCI_SOURCE <- c(
    "seqnames",
    "start",
    "end",
    "width",
    "strand",
    "chrom",
    "pos",
    "A1",
    "A2",
    "variant_id",
    "af",
    "beta",
    "se",
    "z",
    "trait",
    "context",
    "cs_95",
    "cs_70",
    "cs_50",
    "cs_95_purity",
    "cs_70_purity",
    "cs_50_purity",
    "posterior_effect",
    "posterior_se",
    "pip",
    "logBF",
    "lfsr",
    "blockId",
    "study",
    "method",
    "N",
    "p"
)

# The marginal estimates of `x`'s variants joined onto its posterior table.
# @noRd
.fmbMarginalsFor <- function(tl, x) {
    .pubJoinMarginals(
        tl,
        as.data.frame(S4Vectors::mcols(marginalEffects(x)), optional = TRUE)
    )
}

# The topLoci BED table.
#
# `tl$beta`/`tl$se` are the POSTERIOR mean and sd -- topLoci() is the
# posterior view -- while `mm$beta`/`mm$se` are the MARGINAL estimates. Same
# column name, different quantity, which is why they land in different
# published fields. Per-CS variant-level columns from a fullFit build ride
# after the fixed schema; `columns` adds whatever else the caller asked for
# after those.
# @noRd
.fmbTopLociTable <- function(x, signalCutoff, minPurity, columns, coverage) {
    tl <- as.data.frame(
        topLoci(x, signalCutoff = signalCutoff, minPurity = minPurity)
    )
    if (nrow(tl) == 0L) {
        return(NULL)
    }
    mm <- .fmbMarginalsFor(tl, x)
    win <- .fmbRowWindows(x, tl)
    pos <- .pubCol(tl, "pos")
    p <- suppressWarnings(as.numeric(.pubCol(mm, "p")))
    csCols <- .pubCsColumns(tl, coverage)
    cols <- c(
        list(
            `#chr` = .pubCol(tl, "chrom"),
            start = pos - 1L,
            end = pos,
            ID = .pubCol(tl, "variant_id"),
            REF = .pubCol(tl, "A2"),
            ALT = .pubCol(tl, "A1"),
            AF = .pubCol(tl, "af"),
            NS = .pubCol(tl, "N"),
            ES = .pubCol(mm, "beta"),
            SE = .pubCol(mm, "se"),
            EZ = .pubEz(mm),
            LP = if_else(is.na(p) | p <= 0, NA_real_, -log10(p)),
            WGT = .pubCol(tl, "posterior_effect"),
            WGTSE = .pubCol(tl, "posterior_se"),
            THETA = .pubCol(tl, "theta"),
            PIP = .pubCol(tl, "pip"),
            LBF = .pubCol(tl, "logBF"),
            LFSR = .pubCol(tl, "lfsr")
        ),
        .fmbCsValues(tl, csCols),
        .pubIdentityValues(x, tl),
        list(
            region_id = win$id,
            region_start = win$start,
            region_end = win$end,
            region_given = win$given
        )
    )
    out <- .pubFrame(cols, .pubSchema(csCols))
    added <- union(
        .pubFullFitColumns(tl),
        .pubViewColumns(tl, columns, .FMB_TOPLOCI_SOURCE)
    )
    bind_cols(out, select(tl, all_of(added)))
}

# The CS membership / purity values as a tibble, for binding into a block.
#
# `csCols` is never empty here. The only caller runs downstream of
# `.fmbPrimaryCsLevel()`, which aborts outright for a table carrying no
# `cs_<n>` column, and both columns of every level that IS present come back
# from `.pubCsPair()` unconditionally -- so an empty-`csCols` branch would be
# dead by construction.
# @noRd
.pubFrameCs <- function(view, csCols) {
    as_tibble(
        .fmbCsValues(view, csCols),
        .name_repair = "minimal"
    )
}

# @noRd
.fmbCsValues <- function(tl, csCols) {
    set_names(map(csCols, .fmbCsValue, tl = tl), csCols)
}

# @noRd
.fmbCsValue <- function(nm, tl) {
    lv <- str_remove(nm, "^(CS|PUR)")
    if (str_starts(nm, "PUR")) {
        return(.pubCol(tl, str_c("cs_", lv, "_purity")))
    }
    .pubCsIndex(.pubCol(tl, str_c("cs_", lv)))
}

# A credible-set label ("susie_1", "fsusie_2") as its integer index, which is
# what both formats publish: the prefix is the fitting method, already its own
# field, and an index is filterable in a table and representable in a VCF
# Integer FORMAT field. NA -- not 0 -- for "in no credible set at this
# coverage": NA is missing in both formats, where 0 is a sentinel that could
# collide with a real index.
# @noRd
.pubCsIndex <- function(v) {
    idx <- suppressWarnings(as.integer(str_remove(as.character(v), ".*_")))
    replace(idx, is.na(idx) | idx == 0L, NA_integer_)
}

# The identity and provenance values shared by both record kinds. `data_source`
# is derived rather than stored: a collection's ldReference is the LD reference
# the fit was run against, so its absence means the fit saw individual-level
# data.
# @noRd
.pubIdentityValues <- function(x, view) {
    list(
        study = .pubCol(view, "study"),
        context = .pubCol(view, "context"),
        trait = .pubCol(view, "trait"),
        method = .pubCol(view, "method"),
        block_id = .pubCol(view, "blockId"),
        data_type = .pubCol(view, "dataType"),
        data_source = .pubRowOrDerived(view, "dataSource", .pubDataSource(x)),
        weight_standardized = .pubCol(view, "standardized"),
        joint_studies = .pubCol(view, "jointStudies"),
        joint_contexts = .pubCol(view, "jointContexts"),
        joint_traits = .pubCol(view, "jointTraits")
    )
}

# The fallback for `data_source`: a row's LD reference is the one its fit ran
# against, so a row with none is individual-level data. Only a fallback -- a
# row that recorded its own wins. Per row in .pubDataSourceOf(); for the
# collection as a whole it is NA when its rows are of both kinds.
# @noRd
.pubDataSource <- function(x) {
    kinds <- unique(.pubDataSourceOf(x))
    if (length(kinds) == 1L) kinds else NA_character_
}

# @noRd
.pubDataSourceOf <- function(x) {
    if_else(is.na(.ldRefKeys(x)), "individual", "sumstat")
}

# `EZ` only when the SOURCE supplied z. The spec defines the field as the z
# "provided if it was used to derive the ES and SE fields", so publishing a z
# that pecotmr itself computed as beta/se would invert its meaning -- and a
# reader would then believe ES/SE were derived when they were recorded.
# @noRd
.pubEz <- function(view, n = NULL) {
    # `isTRUE()` is length-1 only: `isTRUE(c(TRUE, TRUE))` is FALSE, so the
    # all-or-nothing test this used never fired for a view with more than one
    # variant and EZ was never published at all. Per ROW anyway, matching
    # `.vcfSuppliedZ()`: a file can hold a mix, and a derived z must not
    # borrow a neighbour's provenance.
    supplied <- as.logical(.pubCol(view, "z_supplied", n))
    z <- as.numeric(.pubCol(view, "z", n))
    if (!any(!is.na(supplied) & supplied)) {
        return(if (is.null(n)) NA_real_ else rep(NA_real_, n))
    }
    replace(z, is.na(supplied) | !supplied, NA_real_)
}

# Absent OR entirely empty both mean "the rows did not record this", because
# .fmrRowMetadata() NA-fills a column the collection lacks -- so testing only
# for absence would return the NA column and never reach the fallback.
# @noRd
.pubRowOrDerived <- function(view, nm, fallback) {
    v <- .pubCol(view, nm)
    if (all(is.na(v))) {
        return(fallback)
    }
    v
}

# The analysis window of the collection row each long-view row came from,
# matched on the row identity so a folded multi-row table keeps each row's own
# window rather than the first one's. Both keys are built from the SAME
# columns -- those the view and the collection both carry -- since a GWAS view
# NA-fills a context and trait its collection does not have; and `blockId` is
# among them, as it is what tells the blocks of one GWAS fit apart.
# @noRd
.fmbRowWindows <- function(x, df) {
    rg <- genomicRegion(x)
    cols <- intersect(
        intersect(c(.PUB_ROW_COLS, "blockId"), names(df)),
        .tupleColumnNames(x)
    )
    idx <- match(.pubIdentityKey(df, cols), .pubCollectionKey(x, cols))
    # 0-based half-open, like every other coordinate in a BED -- that is the
    # format's requirement, not a pecotmr choice. `region_id` stays the
    # canonical 1-based `chrom_start_end`, because it is an IDENTIFIER matched
    # against upstream region lists rather than a coordinate pair.
    list(
        id = .pubRegionIds(rg)[idx],
        start = start(rg)[idx] - 1L,
        end = end(rg)[idx],
        given = .pubRegionGiven(x)[idx]
    )
}

# Whether each row's block range was given, or NA where the collection does
# not record it (only a GWAS result does).
# @noRd
.pubRegionGiven <- function(x) {
    if (!is_in("regionGiven", .tupleColumnNames(x))) {
        return(rep(NA, nrow(x)))
    }
    as.logical(.tupleColumn(x, "regionGiven"))
}

# ---- BED: fine-mapping entry point ------------------------------------------

#' @rdname writeBed
#' @export
setMethod(
    "writeBed",
    signature("list"),
    function(obj, filename, ...) {
        writeBed(.pubCombineList(obj), filename, ...)
    }
)

# `region_ID` / `region_start` / `region_end` come from genomicRegion() -- the
# window each row was ANALYSED over. A nominal TAD boundary is what a
# publication table means by its region, not the narrower span of whichever
# variants survived QC.
#' @rdname writeBed
#' @export
setMethod(
    "writeBed",
    signature("FineMappingResultBase"),
    function(
        obj,
        filename,
        schema = c("topLoci", "fsusie"),
        signalCutoff = 0.025,
        minPurity = NULL,
        coverage = NULL,
        levels = c("wavelet", "trait"),
        columns = NULL,
        extraColumns = NULL
    ) {
        schema <- arg_match(schema)
        levels <- arg_match(levels, multiple = TRUE)
        assertString(filename)
        .pubRequireFullFit(obj, .fmbIncompleteRows(obj), "writeBed")
        if (identical(schema, "fsusie")) {
            return(.pubWriteTable(
                .fmbFsusieTable(obj, columns, levels, coverage),
                filename,
                naString = "NA",
                index = TRUE,
                emptySchema = .fmbFsusieColumns(levels),
                extraColumns = extraColumns,
                build = .oneBuild(obj)
            ))
        }
        .pubWriteTable(
            .fmbTopLociTable(obj, signalCutoff, minPurity, columns, coverage),
            filename,
            naString = "NA",
            index = TRUE,
            keepColumns = .pubNamedColumns(columns),
            extraColumns = extraColumns,
            build = .oneBuild(obj)
        )
    }
)

# ---- BED: the TWAS-weight schemas -------------------------------------------

# The source columns the weights schema already publishes. Unlike the
# fine-mapping views, a TwasWeights element carries nothing else today -- its
# variants hold exactly A1 / A2 / weight -- so `columns = "all"` is a no-op on
# a pipeline-built collection and exists for the custom-built ones that
# annotate their variants further.
# @noRd
.TWB_WEIGHT_SOURCE <- c("A1", "A2", "weight")

# A zero-row frame naming every variant-level column ANY row's variants carry.
# `columns` resolves against this rather than against one row, because the
# published schema must not vary row by row: a column one row has and another
# lacks has to become an NA cell, not a column that appears halfway down.
# @noRd
.twbVariantColumns <- function(x) {
    nms <- unique(list_c(c(
        list(character(0)),
        map(seq_len(nrow(x)), .twbRowVariantColumns, x = x)
    )))
    as_tibble(
        set_names(rep(list(character()), length(nms)), nms),
        .name_repair = "minimal"
    )
}

# @noRd
.twbRowVariantColumns <- function(i, x) colnames(mcols(x[[i]]))

# Per-row identity and payload provenance, as a tibble parallel to the
# collection's rows. The identity axes are mcols columns on a TwasWeights, not
# element fields, so they are read through `.tupleColumn()`.
# @noRd
.twbRowMetadata <- function(x) {
    tibble(
        study = .pubTupleColumn("study", x),
        trait = .pubTupleColumn("trait", x),
        context = .pubTupleColumn("context", x),
        method = .pubTupleColumn("method", x),
        standardized = .twbColumnOr(x, "standardized", as.logical, NA),
        dataType = .twbListColumn(x, "dataType"),
        dataSource = .pubDataSourceOf(x),
        jointStudies = .twbColumnOr(
            x,
            "jointStudies",
            as.character,
            NA_character_
        ),
        jointContexts = .twbColumnOr(
            x,
            "jointContexts",
            as.character,
            NA_character_
        ),
        jointTraits = .twbColumnOr(
            x,
            "jointTraits",
            as.character,
            NA_character_
        )
    )
}

# An mcols column read through `coerce`, or a recycled default when the
# collection does not carry it.
# @noRd
.twbColumnOr <- function(x, nm, coerce, default) {
    if (!is_in(nm, .tupleColumnNames(x))) {
        return(rep(default, nrow(x)))
    }
    coerce(.tupleColumn(x, nm))
}

# A per-row List column flattened to one string per row.
#
# No absent-column branch, unlike `.twbColumnOr()` above: this is only ever
# asked for `dataType`, which `TwasWeights()` materializes for every row even
# when no row declares one -- and that constructor is the only way to build
# the class, so the column is always there.
# @noRd
.twbListColumn <- function(x, nm) {
    map_chr(as.list(.tupleColumn(x, nm)), .twbScalarChr)
}

# One List cell as a single string: an absent value becomes NA, and a
# multi-value cell is ";"-packed rather than silently truncated.
# @noRd
.twbScalarChr <- function(v) {
    if (length(v) == 0L) NA_character_ else .pubPack(v)
}

# The weight-matrix column labels of one row. A plain numeric vector -- the
# univariate case, and what a multivariate fit is split into one row per
# context -- has no feature axis, so it gets a single NA label. A matrix keeps
# one label per column: its colname when it has them (fSuSiE's functional grid
# positions), else its index.
# @noRd
.twbWeightLabels <- function(w) {
    if (is.null(dim(w))) {
        return(NA_character_)
    }
    colnames(w) %||% as.character(seq_len(ncol(w)))
}

# The k-th weight column of a row's weights, as a plain numeric vector.
# @noRd
.twbWeightColumn <- function(k, w) {
    if (is.null(dim(w))) as.numeric(w) else as.numeric(w[, k])
}

# Fine-mapping methods have no business in a weights export: their "weights"
# are the fine-mapping posterior re-expressed, so the fine-mapping result says
# strictly more. Those rows are skipped with a message rather than written.
# @noRd
.TWB_FINE_MAPPING_METHODS <- c(
    "susie",
    "susie_inf",
    "susie_ash",
    "mvsusie",
    "fsusie"
)

# @noRd
.twbExportableRows <- function(x, label = "writeBed") {
    methods <- .pubTupleColumn("method", x)
    skip <- which(is_in(methods, .TWB_FINE_MAPPING_METHODS))
    if (length(skip) > 0L) {
        inform(glue(
            "{label}: skipping {length(skip)} fine-mapping-derived ",
            "row(s) ({str_flatten(unique(methods[skip]), ', ')}). Those ",
            "weights are the fine-mapping posterior re-expressed -- export ",
            "the fine-mapping result instead, which carries the posterior, ",
            "its standard error and the credible sets."
        ))
    }
    setdiff(seq_len(nrow(x)), skip)
}

# One (row, weight column) block of the weights table: every variant of the
# row, carrying that column's weight and the row's identity. A weight MATRIX
# column is a member of a joint fit, so `trait` names the member and
# `joint_traits` keeps the full membership.
# @noRd
.twbWeightBlock <- function(k, gr, w, labels, ident, win, af, schema, added) {
    mc <- mcols(gr)
    n <- length(gr)
    cols <- c(
        list(
            `#chr` = as.character(seqnames(gr)),
            start = start(gr) - 1L,
            end = start(gr),
            ID = .grVariantIds(gr),
            REF = as.character(.pubCol(mc, "A2", n)),
            ALT = as.character(.pubCol(mc, "A1", n)),
            AF = af,
            WGT = .twbWeightColumn(k, w),
            study = ident$study,
            context = ident$context,
            trait = .twbTraitOf(labels[[k]], ident$trait),
            method = ident$method,
            data_type = ident$dataType,
            data_source = ident$dataSource,
            weight_standardized = ident$standardized,
            joint_studies = ident$jointStudies,
            joint_contexts = ident$jointContexts,
            joint_traits = ident$jointTraits,
            region_id = win$id,
            region_start = win$start,
            region_end = win$end
        )
    )
    bind_cols(.pubFrame(cols, schema), .twbAddedColumns(mc, added, n))
}

# The caller-requested variant-level columns of one block, NA-filled for a row
# whose variants do not carry them.
# @noRd
.twbAddedColumns <- function(mc, added, n) {
    # A zero-column tibble still has to declare its row count, or bind_cols
    # refuses to recycle the block it is being joined to.
    if (length(added) == 0L) {
        return(tibble(.rows = n))
    }
    as_tibble(
        set_names(map(added, .pubCol, df = mc, n = n), added),
        .name_repair = "minimal"
    )
}

# The trait a weight column belongs to: the matrix column's label when it has
# one (a joint fit's member), else the row's own trait. `%||%` is no use here --
# a plain weight vector's label is NA, not NULL.
# @noRd
.twbTraitOf <- function(label, trait) {
    if (is.na(label)) trait else label
}

# @noRd
.twbWeightRows <- function(i, x, meta, regions, regionIds, schema, added) {
    gr <- x[[i]]
    w <- mcols(gr)$weight
    labels <- .twbWeightLabels(w)
    list_rbind(map(
        seq_along(labels),
        .twbWeightBlock,
        gr = gr,
        w = w,
        labels = labels,
        ident = slice(meta, i),
        win = list(
            id = .pubRegionIds(regions)[[i]],
            start = start(regions)[[i]] - 1L,
            end = end(regions)[[i]]
        ),
        af = .twbVariantAf(.rowLdRef(x, i), gr),
        schema = schema,
        added = added
    ))
}

# Per-variant allele frequency, from the row's LD reference when it records
# one (an LdMatrix with an AF column, or a stochastic sketch, whose real
# frequencies are recorded from its sidecar). A TwasWeights element records
# only alleles and a weight, so otherwise there is no frequency to publish --
# and `AF` is what the standardized-to-raw weight conversion needs, so its
# absence is worth seeing rather than faking.
# @noRd
.twbVariantAf <- function(ref, gr) {
    af <- if (is.null(ref)) NULL else .ldRefAf(ref)
    if (is.null(af)) {
        return(NA_real_)
    }
    rr <- .ldRefRanges(ref)
    as.numeric(af)[match(.grVariantIds(gr), .grVariantIds(rr))]
}

# @noRd
.twbWeightsTable <- function(x, columns, coverage) {
    rows <- .twbExportableRows(x)
    if (length(rows) == 0L) {
        abort(str_c(
            "writeBed: every row is fine-mapping-derived, so there is ",
            "nothing to write as weights. Export the fine-mapping result."
        ))
    }
    schema <- .pubSchema(.pubCsColumns(tibble(), coverage))
    regions <- genomicRegion(x)
    pieces <- compact(map(
        rows,
        .twbWeightRows,
        x = x,
        meta = .twbRowMetadata(x),
        regions = regions,
        regionIds = .pubRegionIds(regions),
        schema = schema,
        added = .pubViewColumns(
            .twbVariantColumns(x),
            columns,
            .TWB_WEIGHT_SOURCE
        )
    ))
    if (length(pieces) == 0L) NULL else list_rbind(pieces)
}

# ---- BED: TWAS-weight entry point -------------------------------------------

# One table, one row per (variant, weight column). There is no companion
# performance table: cross-validation metrics describe a MODEL scored on
# held-out FOLDS, nothing about them is per-variant, and a published file has
# no data-secure way to express anything fold-level. They stay on the object,
# where cvResult() returns them and the caller tabulates whatever they need.
#' @rdname writeBed
#' @export
setMethod(
    "writeBed",
    signature("TwasWeights"),
    function(
        obj,
        filename,
        coverage = NULL,
        columns = NULL,
        extraColumns = NULL
    ) {
        assertString(filename)
        .pubRequireFullFit(obj, .twbIncompleteRows(obj), "writeBed")
        .pubWriteTable(
            .twbWeightsTable(obj, columns, coverage),
            filename,
            naString = "NA",
            index = TRUE,
            keepColumns = .pubNamedColumns(columns),
            extraColumns = extraColumns,
            build = .oneBuild(obj)
        )
    }
)
