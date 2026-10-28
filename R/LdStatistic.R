# =============================================================================
# LdStatistic S4 virtual class
# -----------------------------------------------------------------------------
# Abstract container for pre-computed LD statistics. Subclasses
# (LdEigen, LdScore) provide method-specific representations: LdEigen
# for eigendecomposition-based methods (LDER/HDL/sHDL), LdScore for
# LD-score-based methods (S-LDSC/g-LDSC).
# =============================================================================

#' @include AllGenerics.R AllClasses.R LdBlocks.R
NULL

#' @title LD Statistic (Virtual Base Class)
#' @description Abstract container for pre-computed LD statistics. Subclasses
#'   provide method-specific representations: eigendecompositions (for
#'   LDER/HDL/sHDL) and LD score matrices (for S-LDSC/g-LDSC).
#'
#'   An \code{LdStatistic} \strong{is} a \code{GRanges} of the reference
#'   variants: one range per SNP, with \code{SNP} / \code{A1} / \code{A2}
#'   (and optionally \code{MAF}) in \code{mcols}. Per-variant statistics
#'   live in \code{mcols} too, so subsetting narrows the ranges and the
#'   statistics together. The genome build comes from \code{seqinfo()} rather
#'   than a slot of its own.
#'
#'   A statistic built from a genotype panel (or a
#'   \code{\link{GenotypeMixture}}) keeps that reference rather than its
#'   blocks' LD: each block's LD is read from the panel when an estimator
#'   reaches the block, so only one block's is ever in memory, and
#'   \code{mcols()} records each variant's row in the panel as
#'   \code{panelRow}. One built from \code{\link{LdMatrix}} blocks stores
#'   what it needs of their LD instead.
#' @slot ldBlocks A \code{GRanges} of LD block intervals.
#' @slot nSamples Integer, sample size of the LD reference panel.
#' @slot inSample Logical, whether the LD reference is from the same cohort as
#'   the GWAS (affects bias correction).
#' @slot ldReference The genotype panel or \code{GenotypeMixture} the
#'   statistic reads its blocks' LD from, or \code{NULL} when it was built
#'   from \code{LdMatrix} blocks.
#' @export
setClass(
    "LdStatistic",
    contains = c("VIRTUAL", "GRanges"),
    representation(
        ldBlocks = "GRanges",
        nSamples = "integer",
        inSample = "logical",
        ldReference = "LdReference_OR_NULL"
    ),
    prototype(ldReference = NULL),
    validity = function(object) .validateLdStatistic(object)
)

# @noRd
#' @importFrom checkmate makeAssertCollection assertCount assertLogical
.validateLdStatistic <- function(object) {
    coll <- makeAssertCollection()
    assertCount(
        object@nSamples,
        positive = TRUE,
        .var.name = "nSamples",
        add = coll
    )
    # assertLogical(len = 1) rather than assertFlag: the slot's declared type
    # already excludes non-logicals, and NA is tolerated here as it was before.
    assertLogical(object@inSample, len = 1L, .var.name = "inSample", add = coll)
    if (length(object) == 0L) {
        coll$push("an LdStatistic must carry at least one variant")
    }
    referenceCheck <- .ldStatReferenceCheck(object)
    if (!isTRUE(referenceCheck)) {
        coll$push(referenceCheck)
    }
    coll$getMessages()
}

# A statistic read from a genotype reference finds each variant there by its
# `panelRow`, so every variant needs one, and it must be a row of the panel.
# @noRd
.ldStatReferenceCheck <- function(object) {
    ref <- object@ldReference
    if (is.null(ref)) {
        return(TRUE)
    }
    if (!.ldIsGenotypeRef(ref)) {
        return(str_c(
            "`ldReference` must be a genotype panel or a GenotypeMixture: ",
            "LdMatrix blocks are stored, not read on demand"
        ))
    }
    rows <- S4Vectors::mcols(object, use.names = FALSE)$panelRow
    nRows <- length(.ldRefRanges(ref))
    valid <- is.integer(rows) &&
        !anyNA(rows) &&
        all(rows >= 1L & rows <= nRows)
    if (valid) {
        return(TRUE)
    }
    str_c(
        "a statistic read from `ldReference` needs an integer `panelRow` ",
        "column naming each variant's row in it"
    )
}

# The variants of an LD reference, as the GRanges every subclass is built on.
# Shared by the LdScore() and LdEigen() constructors so the two agree on what
# a reference panel's variant table looks like.
# @noRd
.ldStatRanges <- function(snpInfo, genome) {
    required <- c("SNP", "CHR", "BP", "A1", "A2")
    missingCols <- setdiff(required, colnames(snpInfo))
    if (length(missingCols) > 0L) {
        abort(glue(
            "`snpInfo` is missing column(s): ",
            "{str_flatten(missingCols, ', ')}."
        ))
    }
    gr <- GenomicRanges::GRanges(
        seqnames = withChrPrefix(as.character(snpInfo$CHR)),
        ranges = IRanges::IRanges(as.integer(snpInfo$BP), width = 1L)
    )
    labelled <- S4Vectors::`mcols<-`(
        `names<-`(gr, as.character(snpInfo$SNP)),
        value = S4Vectors::DataFrame(
            select(snpInfo, -any_of(c("CHR", "BP"))),
            row.names = NULL
        )
    )
    .withGenomeBuild(
        labelled,
        !is.null(genome) && length(genome) == 1L && nzchar(genome),
        genome
    )
}

# The ranges carrying `genome` as their seqinfo build, when the caller named
# one. GRanges has no constructor argument for it, so it is written back.
# @noRd
.withGenomeBuild <- function(gr, named, genome) {
    if (!named) {
        return(gr)
    }
    GenomeInfoDb::`genome<-`(gr, value = genome)
}

# =============================================================================
# Building an LD statistic from an LD reference
# -----------------------------------------------------------------------------
# buildLdEigen() and buildLdScore() both start from an LD reference -- one
# LdMatrix per LD block, or a genotype panel cut at the blocks' spans -- and
# differ only in what they compute per block. Everything before that
# (splitting into blocks, assembling the variant table, resolving nRef and
# the genome build) is shared, so it lives here next to .ldStatRanges().
# =============================================================================

# Shared front half of buildLdEigen()/buildLdScore(): validate the input,
# split it into blocks, and resolve what the LdStatistic base needs. A
# genotype reference is kept, with each variant's row in it, for the
# statistic to read its blocks from.
# @noRd
.h2RefPrepare <- function(ldReference, ldBlocks, nSamples, genome) {
    blocks <- .h2RefBlocks(ldReference, ldBlocks)
    resolvedGenome <- .h2RefResolveGenome(ldReference, genome)
    lazy <- .ldIsGenotypeRef(ldReference)
    list(
        blocks = blocks,
        snpInfo = .h2RefSnpInfo(blocks),
        snpIdx = .h2RefSnpIdx(blocks),
        ldBlocks = .h2RefBlockRanges(blocks, resolvedGenome),
        nRef = .h2RefResolveNRef(ldReference, nSamples),
        genome = resolvedGenome,
        ldReference = if (lazy) ldReference,
        panelRow = if (lazy) list_c(map(blocks, "rows"))
    )
}

# The reference's LD blocks: list(R, gr) per LdMatrix, or list(rows, gr) per
# block a genotype panel is cut into at `ldBlocks` -- the panel rows whose LD
# is computed when it is needed, not here.
# @noRd
.h2RefBlocks <- function(ldReference, ldBlocks) {
    if (.ldIsGenotypeRef(ldReference)) {
        return(.h2RefPanelBlocks(ldReference, ldBlocks))
    }
    if (!is.null(ldBlocks)) {
        abort(c(
            "`ldBlocks` cuts a genotype panel into LD blocks.",
            i = "A list of LdMatrix objects is already one block per element."
        ))
    }
    map(.h2RefMatrixList(ldReference), .h2RefMatrixBlock)
}

# A single LdMatrix is one block; a list must hold LdMatrix objects only.
# @noRd
.h2RefMatrixList <- function(ldReference) {
    if (.ldIsMatrixRef(ldReference)) {
        return(list(ldReference))
    }
    isList <- is.list(ldReference) &&
        length(ldReference) > 0L &&
        all(map_lgl(ldReference, .ldIsMatrixRef))
    if (!isList) {
        abort(c(
            str_c(
                "`ldReference` must be a genotype panel, an LdMatrix, or a ",
                "non-empty list of LdMatrix objects (one per LD block)."
            ),
            i = glue("Got {class(ldReference)[[1L]]}.")
        ))
    }
    ldReference
}

# @noRd
.h2RefMatrixBlock <- function(x) {
    list(R = ldMatrix(x), gr = as(x, "GRanges"))
}

# A genotype panel cut into the LD blocks `ldBlocks` spans. Each block holds
# the panel variants inside it, whose LD is computed from their dosages;
# variants outside every block are not part of the reference, and a variant
# inside two blocks would be scored twice, so overlapping blocks are refused.
# @noRd
.h2RefPanelBlocks <- function(panel, ldBlocks) {
    if (is.null(ldBlocks)) {
        abort(c(
            "`ldBlocks` is required with a genotype panel.",
            i = "It says where each LD block starts and ends."
        ))
    }
    spans <- .asLdBlockRanges(ldBlocks)
    gr <- .ldRefRanges(panel)
    hits <- GenomicRanges::findOverlaps(gr, spans)
    if (anyDuplicated(S4Vectors::queryHits(hits)) > 0L) {
        abort("`ldBlocks` overlap: a variant may sit in only one LD block.")
    }
    idx <- compact(map(
        seq_along(spans),
        .h2RefSpanIdx,
        hits = hits
    ))
    if (length(idx) == 0L) {
        abort("None of the panel's variants falls inside `ldBlocks`.")
    }
    map(idx, .h2RefPanelBlockAt, panel = panel, gr = gr)
}

# The panel rows inside block `k`, in panel order, or NULL for an empty block.
# @noRd
.h2RefSpanIdx <- function(k, hits) {
    idx <- sort(S4Vectors::queryHits(hits)[S4Vectors::subjectHits(hits) == k])
    if (length(idx) == 0L) NULL else idx
}

# @noRd
.h2RefPanelBlockAt <- function(idx, panel, gr) {
    list(rows = idx, gr = gr[idx])
}

# Global (concatenated) variant positions of each block.
# @noRd
.h2RefSnpIdx <- function(blocks) {
    sizes <- map_int(blocks, .h2RefBlockSize)
    ends <- cumsum(sizes)
    map2(as.integer(ends - sizes + 1L), as.integer(ends), seq.int)
}

# @noRd
.h2RefBlockSize <- function(block) length(block$gr)

# The variant table .ldStatRanges() expects, stacked over the blocks in
# order. Both kinds of reference name their variants `SNP`; frequency comes
# from `AF` (an LdMatrix's, or a stochastic sketch's recorded one) as MAF.
# @noRd
.h2RefSnpInfo <- function(blocks) {
    bind_rows(map(blocks, .h2RefSnpInfoOne))
}

# @noRd
.h2RefSnpInfoOne <- function(block) {
    gr <- block$gr
    md <- S4Vectors::mcols(gr, use.names = FALSE)
    info <- tibble(
        SNP = .h2RefVariantIds(gr, md),
        CHR = as.character(seqnames(gr)),
        BP = as.integer(start(gr)),
        A1 = .h2RefAllele(md, "A1"),
        A2 = .h2RefAllele(md, "A2")
    )
    maf <- .h2RefMaf(md)
    if (is.null(maf)) info else mutate(info, MAF = maf)
}

# @noRd
.h2RefVariantIds <- function(gr, md) {
    if (is_in("SNP", colnames(md))) {
        return(as.character(md$SNP))
    }
    if (!is.null(names(gr))) {
        return(as.character(names(gr)))
    }
    str_c(as.character(seqnames(gr)), ":", start(gr))
}

# @noRd
.h2RefAllele <- function(md, which) {
    if (!is_in(which, colnames(md))) {
        abort(glue(
            "the LD reference's variants carry no `{which}` column; an ",
            "LdStatistic needs both alleles."
        ))
    }
    as.character(md[[which]])
}

# @noRd
.h2RefMaf <- function(md) {
    if (!is_in("AF", colnames(md))) {
        return(NULL)
    }
    af <- as.numeric(md$AF)
    pmin(af, 1 - af)
}

# One range per block, spanning the variants it covers. These are the blocks
# the jackknife leaves out one at a time, so they must line up with
# eigenList / ldMatrixList element for element; a statistic read from a panel
# finds each block's variants by placing them in these ranges.
# @noRd
.h2RefBlockRanges <- function(blocks, genome) {
    gr <- GRanges(
        seqnames = map_chr(blocks, .h2RefBlockChrom),
        ranges = IRanges::IRanges(
            start = map_int(blocks, .h2RefBlockStart),
            end = map_int(blocks, .h2RefBlockEnd)
        )
    )
    .withGenomeBuild(gr, .h2RefNamedGenome(genome), genome)
}

# An LD block is within-chromosome by construction, and .h2RefBlockRanges()
# collapses each to a single range, so a block spanning chromosomes would
# silently get a range covering everything between them.
# @noRd
.h2RefBlockChrom <- function(block) {
    chrom <- unique(as.character(seqnames(block$gr)))
    if (length(chrom) != 1L) {
        abort(glue(
            "an LD block spans {length(chrom)} chromosomes ",
            "({str_flatten(chrom, ', ')}); each block must lie on one."
        ))
    }
    chrom
}

# @noRd
.h2RefBlockStart <- function(block) as.integer(min(start(block$gr)))

# @noRd
.h2RefBlockEnd <- function(block) as.integer(max(start(block$gr)))

# @noRd
.h2RefNamedGenome <- function(genome) {
    !is.null(genome) &&
        length(genome) == 1L &&
        !is.na(genome) &&
        nzchar(genome)
}

# nRef is the number of independent samples the reference's LD is worth --
# a panel's sample count, a sketch's effective size, or what each LdMatrix
# records -- unless the caller names one with `nSamples`. Differing sizes
# across LdMatrix blocks are an error rather than a silent pick: nRef drives
# HDL's finite-reference correction and the n of ldsc's unbiased r^2 in the
# LD scores, so the wrong one shifts h2.
# @noRd
.h2RefResolveNRef <- function(ldReference, nSamples) {
    if (!is.null(nSamples)) {
        return(as.integer(nSamples))
    }
    if (.ldIsGenotypeRef(ldReference)) {
        return(as.integer(.ldRefNoiseN(ldReference)))
    }
    found <- discard(
        unique(map_int(.h2RefMatrixList(ldReference), .h2RefNSamplesOf)),
        is.na
    )
    if (length(found) == 0L) {
        abort(c(
            "`nSamples` is required: no LdMatrix records its panel's size.",
            i = "Pass `nSamples`, or build them with LdMatrix(nSamples = )."
        ))
    }
    if (length(found) > 1L) {
        abort(glue(
            "the LdMatrix blocks record differing reference panel sizes ",
            "({str_flatten(found, ', ')}); pass `nSamples` to choose one."
        ))
    }
    found[[1L]]
}

# A named helper rather than `nSamples` itself: inside .h2RefResolveNRef()
# that name is the formal, which would shadow the accessor in FUN position.
# @noRd
.h2RefNSamplesOf <- function(x) {
    nSamples(x)
}

# @noRd
.h2RefResolveGenome <- function(ldReference, genome) {
    if (.h2RefNamedGenome(genome)) {
        return(genome)
    }
    refs <- if (.ldIsGenotypeRef(ldReference)) {
        list(.ldRefRanges(ldReference))
    } else {
        .h2RefMatrixList(ldReference)
    }
    found <- discard(unique(map_chr(refs, .h2RefGenomeOf)), is.na)
    if (length(found) == 1L) found[[1L]] else NA_character_
}

# @noRd
.h2RefGenomeOf <- function(x) {
    g <- discard(unique(GenomeInfoDb::genome(as(x, "GRanges"))), is.na)
    if (length(g) == 1L) g[[1L]] else NA_character_
}

#' @rdname ldPanelNSamples
#' @export
setMethod("ldPanelNSamples", "LdStatistic", function(x) x@nSamples)

#' @rdname inSample
#' @export
setMethod("inSample", "LdStatistic", function(x) x@inSample)

#' @rdname ldBlocks
#' @export
setMethod("ldBlocks", "LdStatistic", function(x) x@ldBlocks)

#' @rdname ldReference
#' @export
setMethod("ldReference", "LdStatistic", function(x) x@ldReference)

# An LdStatistic saved before the reference's sample size was named
# `nSamples` holds it as `nRef`, and one saved before statistics kept their LD
# reference has no `ldReference` slot: it stored what it needs of its blocks'
# LD, which a NULL reference now says.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "LdStatistic",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message(str_c("updateObject(object = '", class(object)[[1L]], "')"))
        }
        object <- .renameSavedSlots(object, c(nRef = "nSamples"))
        if (!methods::.hasSlot(object, "ldReference")) {
            object <- `slot<-`(
                object,
                "ldReference",
                check = FALSE,
                value = NULL
            )
        }
        methods::callNextMethod(object, ..., verbose = verbose)
    }
)

# =============================================================================
# Reading an LD statistic's blocks
# -----------------------------------------------------------------------------
# The estimators walk a statistic block by block. Each block is a record:
# `idx`, the block's variant positions in the statistic, and where its LD
# comes from -- `R`, the matrix an LdScore stores; `eigen`, the decomposition
# an LdEigen stores; or `rows`, the block's rows in the genotype reference it
# is read from. A scores-only LdScore has `idx` alone. The readers below turn
# a record into the LD form a computation needs, reading a genotype
# reference's block only then.
# =============================================================================

# The records of a statistic that stores nothing per block: its variants
# placed in its LD blocks, with their rows in its genotype reference when it
# has one.
# @noRd
.ldStatLazyBlocks <- function(x) {
    idx <- .ldStatPlacedIdx(x)
    if (is.null(x@ldReference)) {
        return(map(idx, .ldStatIdxBlock))
    }
    rows <- S4Vectors::mcols(x, use.names = FALSE)$panelRow
    map(idx, .ldStatPanelBlock, rows = rows)
}

# Each LD block's variant positions, in block order, found by placing the
# variants in `ldBlocks(x)`. A block no variant falls in has no entry: the
# placement only groups the hits there are.
# @noRd
.ldStatPlacedIdx <- function(x) {
    snpInfo <- tibble(CHR = as.character(seqnames(x)), BP = start(x))
    idx <- unname(snpsPerBlock(snpInfo, x@ldBlocks))
    if (length(idx) == 0L) {
        abort("no LD block covers any reference variant.")
    }
    idx
}

# @noRd
.ldStatIdxBlock <- function(idx) {
    list(idx = idx)
}

# @noRd
.ldStatPanelBlock <- function(idx, rows) {
    list(idx = idx, rows = rows[idx])
}

# A block's LD source (see .ldSourceRows()): the stored matrix, the stored
# decomposition, or a factor of it read from the genotype reference.
# @noRd
.ldStatBlockLd <- function(block, ldReference) {
    if (!is.null(block$R)) {
        return(list(R = block$R))
    }
    if (!is.null(block$eigen)) {
        return(block$eigen[c("vectors", "values")])
    }
    list(F = .ldRefFactor(ldReference, block$rows))
}

# A block's eigendecomposition, keeping the leading components that carry
# `truncation` of its eigenvalue mass: the stored one, or one computed from
# the factor of the block read from the genotype reference. A computed one
# lists only the min(m, p) directions the factor's SVD returns; `nNull` counts
# the block's other directions, whose eigenvalue is exactly zero, when
# nothing is truncated.
# @noRd
.ldStatBlockEigen <- function(block, ldReference, truncation) {
    if (!is.null(block$eigen)) {
        return(list_assign(block$eigen, nNull = 0L))
    }
    dec <- .ldFactorEigen(.ldRefFactor(ldReference, block$rows))
    keep <- .ldEigenKeep(dec$values, truncation)
    list(
        values = dec$values[keep],
        vectors = dec$vectors[, keep, drop = FALSE],
        snpIdx = block$idx,
        nNull = if (truncation >= 1) dec$nNull else 0L
    )
}

# `obj` reading its blocks from the genotype reference it was built from: the
# reference in its slot and each variant's row in it as `panelRow`. A
# statistic built from LdMatrix blocks is returned as it is.
# @noRd
.ldStatAttachReference <- function(obj, prep) {
    if (is.null(prep$ldReference)) {
        return(obj)
    }
    S4Vectors::mcols(obj)$panelRow <- prep$panelRow
    obj@ldReference <- prep$ldReference
    validObject(obj)
    obj
}
