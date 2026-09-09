# =============================================================================
# LdStatistic S4 virtual class
# -----------------------------------------------------------------------------
# Abstract container for pre-computed LD statistics. Subclasses
# (LdEigen, LdScore) provide method-specific representations: LdEigen
# for eigendecomposition-based methods (LDER/HDL/sHDL), LdScore for
# LD-score-based methods (S-LDSC/g-LDSC).
# =============================================================================

#' @include AllGenerics.R LdBlocks.R
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
#' @slot ldBlocks A \code{GRanges} of LD block intervals.
#' @slot nRef Integer, sample size of the LD reference panel.
#' @slot inSample Logical, whether the LD reference is from the same cohort as
#'   the GWAS (affects bias correction).
#' @export
setClass(
    "LdStatistic",
    contains = c("VIRTUAL", "GRanges"),
    representation(
        ldBlocks = "GRanges",
        nRef = "integer",
        inSample = "logical"
    ),
    validity = function(object) .validateLdStatistic(object)
)

# @noRd
.validateLdStatistic <- function(object) {
    errors <- character()
    if (length(object@nRef) != 1L || object@nRef <= 0L) {
        errors <- c(errors, "'nRef' must be a single positive integer")
    }
    if (length(object@inSample) != 1L) {
        errors <- c(errors, "'inSample' must be a single logical value")
    }
    if (length(object) == 0L) {
        errors <- c(errors, "an LdStatistic must carry at least one variant")
    }
    if (length(errors) == 0) TRUE else errors
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
    names(gr) <- as.character(snpInfo$SNP)
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        select(snpInfo, -any_of(c("CHR", "BP"))),
        row.names = NULL
    )
    if (!is.null(genome) && length(genome) == 1L && nzchar(genome)) {
        GenomeInfoDb::genome(gr) <- genome
    }
    gr
}

# =============================================================================
# Building an LD statistic from already-loaded LD
# -----------------------------------------------------------------------------
# buildLdEigen() and buildLdScore() both start from LD that has already been
# read -- one or more LdData, as loadLdMatrix() returns them --
# and differ only in what they compute per block. Everything before that
# (splitting into blocks, assembling the variant table, resolving nRef and the
# genome build) is shared, so it lives here next to .ldStatRanges().
# =============================================================================

# Shared front half of buildLdEigen()/buildLdScore(): validate the input,
# split it into blocks, and resolve what the LdStatistic base needs.
# @noRd
.ldRefPrepare <- function(ldBlockData, nRef, genome) {
    dataList <- .ldRefAsDataList(ldBlockData)
    blocks <- list_flatten(map(dataList, .ldRefBlocksOf))
    if (length(blocks) == 0L) {
        abort("`ldBlockData` yielded no LD blocks.")
    }
    resolvedGenome <- .ldRefResolveGenome(dataList, genome)
    list(
        blocks = blocks,
        snpInfo = .ldRefSnpInfo(blocks),
        snpIdx = .ldRefSnpIdx(blocks),
        ldBlocks = .ldRefBlockRanges(blocks, resolvedGenome),
        nRef = .ldRefResolveNRef(dataList, nRef),
        genome = resolvedGenome
    )
}

# Accept a single LdData as well as a list of them.
# @noRd
.ldRefAsDataList <- function(ldBlockData) {
    if (is(ldBlockData, "LdData")) {
        return(list(ldBlockData))
    }
    if (!is.list(ldBlockData) || length(ldBlockData) == 0L) {
        abort(glue(
            "`ldBlockData` must be an LdData or a non-empty list of them ",
            "(got {class(ldBlockData)[[1L]]})."
        ))
    }
    bad <- which(!map_lgl(ldBlockData, is, class2 = "LdData"))
    if (length(bad) > 0L) {
        abort(glue(
            "`ldBlockData` element(s) {str_flatten(bad, ', ')} are not ",
            "LdData objects."
        ))
    }
    ldBlockData
}

# The per-block LD matrices of one LdData, each paired with the variants it
# describes. An LdData holds either a single matrix (one block) or a list of
# them (block-diagonal LD across a multi-block region).
# @noRd
.ldRefBlocksOf <- function(x) {
    R <- getCorrelation(x)
    gr <- getVariantInfo(x)
    if (!is.list(R)) {
        return(list(.ldRefOneBlock(R, gr)))
    }
    idx <- .ldRefBlockIndices(x, R, length(gr))
    map2(R, idx, .ldRefBlockFromIdx, gr = gr)
}

# @noRd
.ldRefBlockFromIdx <- function(R, idx, gr) {
    .ldRefOneBlock(R, gr[idx])
}

# @noRd
.ldRefOneBlock <- function(R, gr) {
    R <- as.matrix(R)
    if (nrow(R) != length(gr)) {
        abort(glue(
            "an LD block's correlation matrix is {nrow(R)}x{ncol(R)} but ",
            "covers {length(gr)} variant(s)."
        ))
    }
    list(R = R, gr = gr)
}

# Which variants each block matrix covers. blockMetadata's startIdx/endIdx is
# authoritative when it lines up with the matrices; otherwise the blocks are
# taken in order, each covering as many variants as its matrix has rows.
# @noRd
.ldRefBlockIndices <- function(x, R, nVariants) {
    sizes <- map_int(R, nrow)
    if (sum(sizes) != nVariants) {
        abort(glue(
            "the LD block matrices cover {sum(sizes)} variant(s) but the ",
            "LdData carries {nVariants}."
        ))
    }
    md <- getBlockMetadata(x)
    hasIdx <- is.data.frame(md) &&
        all(is_in(c("startIdx", "endIdx"), names(md))) &&
        nrow(md) == length(R)
    if (hasIdx) {
        return(map2(as.integer(md$startIdx), as.integer(md$endIdx), seq.int))
    }
    ends <- cumsum(sizes)
    map2(as.integer(ends - sizes + 1L), as.integer(ends), seq.int)
}

# Global (concatenated) variant positions of each block.
# @noRd
.ldRefSnpIdx <- function(blocks) {
    sizes <- map_int(blocks, .ldRefBlockSize)
    ends <- cumsum(sizes)
    map2(as.integer(ends - sizes + 1L), as.integer(ends), seq.int)
}

# @noRd
.ldRefBlockSize <- function(block) length(block$gr)

# The variant table .ldStatRanges() expects, stacked over the blocks in
# order. An LdData names its variants `variant_id`, and reports frequency as
# allele_freq rather than MAF.
# @noRd
.ldRefSnpInfo <- function(blocks) {
    bind_rows(map(blocks, .ldRefSnpInfoOne))
}

# @noRd
.ldRefSnpInfoOne <- function(block) {
    gr <- block$gr
    md <- S4Vectors::mcols(gr, use.names = FALSE)
    info <- tibble(
        SNP = .ldRefVariantIds(gr, md),
        CHR = as.character(seqnames(gr)),
        BP = as.integer(start(gr)),
        A1 = .ldRefAllele(md, "A1"),
        A2 = .ldRefAllele(md, "A2")
    )
    maf <- .ldRefMaf(md)
    if (is.null(maf)) info else mutate(info, MAF = maf)
}

# @noRd
.ldRefVariantIds <- function(gr, md) {
    if (is_in("variant_id", colnames(md))) {
        return(as.character(md$variant_id))
    }
    if (is_in("SNP", colnames(md))) {
        return(as.character(md$SNP))
    }
    if (!is.null(names(gr))) {
        return(as.character(names(gr)))
    }
    str_c(as.character(seqnames(gr)), ":", start(gr))
}

# @noRd
.ldRefAllele <- function(md, which) {
    if (!is_in(which, colnames(md))) {
        abort(glue(
            "the LD reference's variants carry no `{which}` column; an ",
            "LdStatistic needs both alleles."
        ))
    }
    as.character(md[[which]])
}

# @noRd
.ldRefMaf <- function(md) {
    if (is_in("MAF", colnames(md))) {
        return(as.numeric(md$MAF))
    }
    if (is_in("allele_freq", colnames(md))) {
        af <- as.numeric(md$allele_freq)
        return(pmin(af, 1 - af))
    }
    NULL
}

# One range per block, spanning the variants it covers. These are the blocks
# the jackknife leaves out one at a time, so they must line up with
# eigenList / ldMatrixList element for element.
# @noRd
.ldRefBlockRanges <- function(blocks, genome) {
    gr <- GRanges(
        seqnames = map_chr(blocks, .ldRefBlockChrom),
        ranges = IRanges::IRanges(
            start = map_int(blocks, .ldRefBlockStart),
            end = map_int(blocks, .ldRefBlockEnd)
        )
    )
    if (.ldRefNamedGenome(genome)) {
        GenomeInfoDb::genome(gr) <- genome
    }
    gr
}

# An LD block is within-chromosome by construction, and .ldRefBlockRanges()
# collapses each to a single range, so a block spanning chromosomes would
# silently get a range covering everything between them.
# @noRd
.ldRefBlockChrom <- function(block) {
    chrom <- unique(withChrPrefix(as.character(seqnames(block$gr))))
    if (length(chrom) != 1L) {
        abort(glue(
            "an LD block spans {length(chrom)} chromosomes ",
            "({str_flatten(chrom, ', ')}); each block must lie on one."
        ))
    }
    chrom
}

# @noRd
.ldRefBlockStart <- function(block) as.integer(min(start(block$gr)))

# @noRd
.ldRefBlockEnd <- function(block) as.integer(max(start(block$gr)))

# @noRd
.ldRefNamedGenome <- function(genome) {
    !is.null(genome) &&
        length(genome) == 1L &&
        !is.na(genome) &&
        nzchar(genome)
}

# nRef comes from the LdData unless the caller names one. Differing panel
# sizes are an error rather than a silent pick: nRef drives HDL's
# finite-reference correction, so the wrong one shifts h2.
# @noRd
.ldRefResolveNRef <- function(dataList, nRef) {
    if (!is.null(nRef)) {
        return(as.integer(nRef))
    }
    found <- unique(map_int(dataList, .ldRefNRefOf))
    found <- found[!is.na(found)]
    if (length(found) == 0L) {
        abort(glue(
            "`nRef` is required: none of the supplied LdData records a ",
            "reference panel sample size."
        ))
    }
    if (length(found) > 1L) {
        abort(glue(
            "the supplied LdData record differing reference panel sizes ",
            "({str_flatten(found, ', ')}); pass `nRef` to choose one."
        ))
    }
    found[[1L]]
}

# @noRd
.ldRefNRefOf <- function(x) {
    n <- getNRef(x)
    if (length(n) == 0L) NA_integer_ else as.integer(n[[1L]])
}

# @noRd
.ldRefResolveGenome <- function(dataList, genome) {
    if (.ldRefNamedGenome(genome)) {
        return(genome)
    }
    found <- unique(map_chr(dataList, .ldRefGenomeOf))
    found <- found[!is.na(found)]
    if (length(found) == 1L) found[[1L]] else NA_character_
}

# @noRd
.ldRefGenomeOf <- function(x) {
    g <- unique(GenomeInfoDb::genome(as(x, "GRanges")))
    g <- g[!is.na(g)]
    if (length(g) == 1L) g[[1L]] else NA_character_
}

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdStatistic", function(x) x@nRef)

#' @rdname getInSample
#' @export
setMethod("getInSample", "LdStatistic", function(x) x@inSample)

#' @rdname getLdBlocks
#' @export
setMethod("getLdBlocks", "LdStatistic", function(x) x@ldBlocks)

#' @rdname getGenome
#' @export
setMethod("getGenome", "LdStatistic", function(x, ...) {
    # The build lives in seqinfo, not a slot: a GRanges already has somewhere
    # to keep it, and storing it twice is what the retired LdBlocks class did.
    build <- unique(GenomeInfoDb::genome(x))
    build <- build[!is.na(build)]
    if (length(build) == 0L) NA_character_ else build[[1L]]
})
