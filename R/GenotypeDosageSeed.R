# =============================================================================
# GenotypeDosageSeed -- a DelayedArray backend over a VCF, GDS or PLINK1 file
# -----------------------------------------------------------------------------
# A VCF, GDS or PLINK1 panel reads its dosages through this seed, as a PLINK2
# or BGEN panel reads through a GenotypeHandle's GhSeed. Each row is one
# variant: one ALT allele of one record of the file, so a record with several
# ALT alleles is read as one variant per ALT allele. The seed knows where each
# row's record is in the file and which allele the row counts, and turns what
# the file stores into that allele's dosage. What it reads is its `field`:
#
#   vcfGt               a VCF's GT calls, through VariantAnnotation
#   vcfDs               a VCF's FORMAT/DS values, through VariantAnnotation
#   snpGdsSnpOrder      a SNPRelate GDS's genotype codes, stored SNP by SNP
#   snpGdsSampleOrder   or sample by sample, through GDSArray
#   seqArrayGt          a SeqArray GDS's allele indices, through GDSArray
#   seqArrayDs          a SeqArray GDS's annotation/format/DS, through GDSArray
#   plinkBed            a PLINK1 .bed's genotype codes, through PlinkMatrix
#
# VCFArray's seed is not built on. It reads a row through every record that
# overlaps the row's range, so in a VCF whose records overlap -- a split
# multi-allelic site, a deletion spanning a SNP -- it returns more rows than
# were asked for and DelayedArray rejects the block. The VCF fields are read
# through VariantAnnotation, as VCFArray reads them, and each record is
# picked out of what its range returns.
#
# Orientation: variants x samples, as for GhSeed.
# =============================================================================

#' @include AllGenerics.R
NULL

#' DelayedArray seed over a VCF, GDS or PLINK1 file
#'
#' Lets the dosages of a VCF, GDS or PLINK1 file back a \code{DelayedArray}
#' that reads the file a block at a time, as \code{GhSeed} does for the
#' formats a \code{GenotypeHandle} reads.
#'
#' @slot backend What the file is read through: a
#'   \code{VariantAnnotation::VcfFile}, the \code{GDSArray} seed of the
#'   GDS node holding the field, or the stem of a PLINK1 fileset, which
#'   \code{PlinkMatrix::read_bed_subset()} reads.
#' @slot field Which field is read, and so how it becomes a dosage:
#'   \code{"vcfGt"}, \code{"vcfDs"}, \code{"snpGdsSnpOrder"},
#'   \code{"snpGdsSampleOrder"}, \code{"seqArrayGt"},
#'   \code{"seqArrayDs"} or \code{"plinkBed"}.
#' @slot rowRecord Integer; for each row, the element of \code{records} it
#'   reads.
#' @slot rowAllele Integer; for each row, the ALT allele it counts, 1 for a
#'   record's first ALT allele.
#' @slot records Where each record is in the file. For a VCF, a
#'   \code{GRanges} of each record's REF allele, with \code{posRank}, the
#'   record's rank among those starting at its position; for a GDS, an
#'   \code{IRanges} of the columns (or SeqArray's 2-bit planes) of the node
#'   holding it; for a PLINK1 fileset, an \code{IRanges} of each variant's
#'   position in the \code{.bed}.
#' @slot dn List of dimnames (variant ids, sample ids).
#' @name GenotypeDosageSeed-class
#' @keywords internal
#' @importClassesFrom S4Vectors Vector
#' @importFrom BiocGenerics width
setClass(
    "GenotypeDosageSeed",
    representation(
        backend = "ANY",
        field = "character",
        rowRecord = "integer",
        rowAllele = "integer",
        records = "Vector",
        dn = "list"
    )
)

#' @rdname GenotypeDosageSeed-class
setMethod("dim", "GenotypeDosageSeed", function(x) {
    c(length(x@rowRecord), length(x@dn[[2L]]))
})

#' @rdname GenotypeDosageSeed-class
setMethod("dimnames", "GenotypeDosageSeed", function(x) x@dn)

#' @rdname GenotypeDosageSeed-class
#' @importFrom S4Arrays extract_array
setMethod("extract_array", "GenotypeDosageSeed", function(x, index) {
    # A NULL index means everything along that dimension, as for GhSeed.
    dm <- dim(x)
    i <- index[[1L]] %||% seq_len(dm[[1L]])
    j <- index[[2L]] %||% seq_len(dm[[2L]])
    if (length(i) == 0L || length(j) == 0L) {
        return(matrix(numeric(0), length(i), length(j)))
    }
    # Each distinct row and sample is read once, in file order, and the block
    # then laid out as asked: DelayedArray may ask in any order, and repeat.
    rows <- sort(unique(i))
    samples <- sort(unique(j))
    block <- .gdSeedRead(x, rows, samples)
    block[match(i, rows), match(j, samples), drop = FALSE]
})

# Internal accessors: GenotypeDosageSeed is not exported, so these add no
# public surface, and they keep the `@` next to the class.
# @noRd
.gdSeedBackend <- function(x) x@backend

# @noRd
.gdSeedField <- function(x) x@field

# @noRd
.gdSeedRowRecord <- function(x) x@rowRecord

# @noRd
.gdSeedRowAllele <- function(x) x@rowAllele

# @noRd
.gdSeedRecords <- function(x) x@records

# The dosages of rows `rows` for samples `samples`, both ascending and
# distinct: a length(rows) x length(samples) matrix. A GDS seed's backend is
# GDSArray's, whose methods are there only once its namespace is loaded,
# which reading a saved panel back does not do.
# @noRd
.gdSeedRead <- function(x, rows, samples) {
    if (str_starts(.gdSeedField(x), "snpGds|seqArray")) {
        loadNamespace("GDSArray")
    }
    switch(
        .gdSeedField(x),
        vcfGt = ,
        vcfDs = .vcfSeedRead(x, rows, samples),
        snpGdsSnpOrder = ,
        snpGdsSampleOrder = .snpGdsSeedRead(x, rows, samples),
        seqArrayGt = .seqArrayGtSeedRead(x, rows, samples),
        seqArrayDs = .seqArrayDsSeedRead(x, rows, samples),
        plinkBed = .plinkBedSeedRead(x, rows, samples)
    )
}

# The seed cut to rows `keep` (ascending positions), carrying only the
# records those rows read, so a panel cut to a few variants does not carry
# every record of the file.
# @noRd
.gdSeedRows <- function(seed, keep) {
    rowRecord <- .gdSeedRowRecord(seed)[keep]
    used <- sort(unique(rowRecord))
    methods::initialize(
        seed,
        rowRecord = match(rowRecord, used),
        rowAllele = .gdSeedRowAllele(seed)[keep],
        records = .gdSeedRecords(seed)[used],
        dn = list(dimnames(seed)[[1L]][keep], dimnames(seed)[[2L]])
    )
}

# ---- VCF -------------------------------------------------------------------

# The VCF FORMAT field a seed reads.
# @noRd
.vcfSeedFieldName <- function(x) {
    if (.gdSeedField(x) == "vcfGt") "GT" else "DS"
}

# A VCF seed's rows: the records they read are fetched by their ranges, each
# picked out of what its range returns, and the field turned into the
# dosage of each row's ALT allele.
# @noRd
.vcfSeedRead <- function(x, rows, samples) {
    recordOfRow <- .gdSeedRowRecord(x)[rows]
    wanted <- sort(unique(recordOfRow))
    spans <- .gdSeedRecords(x)[wanted]
    queryRanges <- IRanges::reduce(spans)
    queryRanges <- `names<-`(
        queryRanges,
        as.character(seq_along(queryRanges))
    )
    sampleIds <- dimnames(x)[[2L]][samples]
    vcf <- VariantAnnotation::readVcf(
        .vcfSeedFile(x),
        param = VariantAnnotation::ScanVcfParam(
            fixed = NA,
            info = NA,
            geno = .vcfSeedFieldName(x),
            samples = sampleIds,
            which = queryRanges
        )
    )
    returned <- .vcfMatchRecords(spans, queryRanges, vcf, x)
    values <- VariantAnnotation::geno(vcf)[[.vcfSeedFieldName(x)]]
    picked <- values[
        returned[match(recordOfRow, wanted)],
        match(sampleIds, colnames(values)),
        drop = FALSE
    ]
    alleles <- .gdSeedRowAllele(x)[rows]
    if (.gdSeedField(x) == "vcfGt") {
        return(vcfGtDosageCpp(picked, alleles))
    }
    .vcfDsDosage(picked, alleles)
}

# A fresh VcfFile over the seed's file and index. The one the seed holds is
# not read through directly: one saved with the panel and loaded again would
# carry a stale connection.
# @noRd
.vcfSeedFile <- function(x) {
    file <- .gdSeedBackend(x)
    VariantAnnotation::VcfFile(path(file), index = Rsamtools::index(file))
}

# Which row of `vcf` holds each record of `spans`.
#
# A range returns every record overlapping it, so a record is not simply the
# n-th row returned. Each record lies in exactly one of the disjoint
# `queryRanges`, and that range returns every record starting at the
# record's position, in file order; so the record is the row from its range
# that starts at its position with its rank there, which `posRank` recorded
# when the panel was read. Matching on integers keeps the seed free of a
# string per record.
# @noRd
.vcfMatchRecords <- function(spans, queryRanges, vcf, x) {
    got <- SummarizedExperiment::rowRanges(vcf)
    have <- tibble(
        range = as.integer(as.character(got$paramRangeID)),
        start = as.integer(start(got))
    ) |>
        mutate(
            rank = row_number(),
            .by = all_of(c("range", "start"))
        )
    want <- tibble(
        range = findOverlaps(spans, queryRanges, select = "first"),
        start = as.integer(start(spans)),
        rank = as.integer(S4Vectors::mcols(spans)$posRank)
    )
    hit <- vctrs::vec_match(want, have)
    if (anyNA(hit)) {
        abort(c(
            glue(
                "{path(.gdSeedBackend(x))} no longer holds the records the ",
                "panel read from it."
            ),
            i = "Read it again with readGenotypes()."
        ))
    }
    hit
}

# The dosage of each row's ALT allele from a block of FORMAT/DS values: one
# value per ALT allele (Number=A, list cells), or one per record (Number=1),
# which is then its first ALT allele's. A row asking for an allele the
# record holds no value for is missing.
# @noRd
.vcfDsDosage <- function(values, alleles) {
    if (!is.list(values)) {
        return(unname(`[<-`(values, alleles != 1L, , value = NA_real_)))
    }
    # A list keeps its dim through as.vector(); drop it to walk the cells.
    cells <- `dim<-`(values, NULL)
    nValues <- lengths(cells)
    allele <- rep(alleles, times = ncol(values))
    held <- allele <= nValues
    dosage <- rep(NA_real_, length(cells))
    if (any(held)) {
        flat <- as.numeric(list_c(cells))
        at <- cumsum(nValues) - nValues + allele
        dosage <- replace(dosage, held, flat[at[held]])
    }
    matrix(dosage, nrow = nrow(values))
}

# ---- SNPRelate GDS ---------------------------------------------------------

# A SNPRelate GDS stores genotype codes 0-2 counting the first of a SNP's two
# alleles, which is therefore A1, and 3 for a missing call.
# @noRd
.snpGdsSeedRead <- function(x, rows, samples) {
    variant <- start(.gdSeedRecords(x))[.gdSeedRowRecord(x)[rows]]
    codes <- if (.gdSeedField(x) == "snpGdsSnpOrder") {
        extract_array(.gdSeedBackend(x), list(variant, samples))
    } else {
        t(extract_array(.gdSeedBackend(x), list(samples, variant)))
    }
    replace(`storage.mode<-`(codes, "double"), codes == 3L, NA_real_)
}

# ---- SeqArray GDS ----------------------------------------------------------

# A SeqArray GDS stores each call's allele indices (0 for REF, k for the
# k-th ALT allele) in 2-bit planes, ploidy x samples x planes: one plane per
# record, or more for a record with too many alleles for one, its index then
# the planes' sum weighted by 1, 4, 16, .... A missing allele sets every bit
# of the record's planes. A call with a missing allele is missing, as
# SeqArray's own dosages read it.
# @noRd
.seqArrayGtSeedRead <- function(x, rows, samples) {
    recordOfRow <- .gdSeedRowRecord(x)[rows]
    wanted <- sort(unique(recordOfRow))
    planes <- .gdSeedRecords(x)[wanted]
    bits <- extract_array(
        .gdSeedBackend(x),
        list(NULL, samples, sequence(width(planes), start(planes)))
    )
    last <- cumsum(width(planes))
    calls <- map2(
        last - width(planes) + 1L,
        width(planes),
        .seqArrayRecordCalls,
        bits = bits
    )
    copies <- map2(
        calls[match(recordOfRow, wanted)],
        .gdSeedRowAllele(x)[rows],
        .seqArrayCopies
    )
    matrix(list_c(copies), nrow = length(rows), byrow = TRUE)
}

# One record's allele indices, ploidy x samples, from its `n` planes starting
# at plane `first` of `bits`; NA where an allele is missing.
# @noRd
.seqArrayRecordCalls <- function(first, n, bits) {
    weighted <- map(
        seq_len(n),
        .seqArrayPlaneValue,
        first = first,
        bits = bits
    )
    value <- reduce(weighted, `+`)
    replace(value, value == 4^n - 1, NA)
}

# Plane `k` of a record whose planes start at `first`, weighted by its place.
# @noRd
.seqArrayPlaneValue <- function(k, first, bits) {
    plane <- matrix(bits[,, first + k - 1L], nrow = dim(bits)[[1L]])
    plane * 4^(k - 1L)
}

# The copies of ALT allele `allele` in each sample's call; NA for a call with
# a missing allele.
# @noRd
.seqArrayCopies <- function(calls, allele) {
    colSums(calls == allele)
}

# A SeqArray GDS stores FORMAT/DS as samples x values, each record's values
# in consecutive columns: one per ALT allele (Number=A), or one (Number=1),
# which is then its first ALT allele's. A missing value is stored as NaN.
# @noRd
.seqArrayDsSeedRead <- function(x, rows, samples) {
    spans <- .gdSeedRecords(x)[.gdSeedRowRecord(x)[rows]]
    alleles <- .gdSeedRowAllele(x)[rows]
    column <- start(spans) + alleles - 1L
    held <- alleles <= width(spans)
    missing <- matrix(NA_real_, length(rows), length(samples))
    if (!any(held)) {
        return(missing)
    }
    columns <- sort(unique(column[held]))
    values <- t(extract_array(.gdSeedBackend(x), list(samples, columns)))
    dosage <- `[<-`(
        missing,
        held,
        ,
        value = values[match(column[held], columns), , drop = FALSE]
    )
    replace(dosage, is.nan(dosage), NA_real_)
}

# ---- PLINK1 ----------------------------------------------------------------

# A PLINK1 .bed stores 2-bit genotype codes SNP by SNP. PlinkMatrix reads them
# as the copies of each variant's A1 allele -- the .bim's fifth column, as
# pecotmr counts it -- with NA for a missing call. It is not told the sample
# count: it counts the .fam itself, and a wrong count would misread silently.
# The stem may be a bundled example's portable reference, resolved here.
# PlinkMatrix is checked for here too, not only when the panel is built: a
# panel loaded from a saved object, such as the bundled examples, was built
# elsewhere.
# @noRd
.plinkBedSeedRead <- function(x, rows, samples) {
    .requireReaderPackage("PlinkMatrix", "PLINK1")
    variant <- start(.gdSeedRecords(x))[.gdSeedRowRecord(x)[rows]]
    PlinkMatrix::read_bed_subset(
        .resolveGenotypeResourcePath(.gdSeedBackend(x)),
        variant,
        samples
    )
}
