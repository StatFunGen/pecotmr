# =============================================================================
# LdData S4 class
# -----------------------------------------------------------------------------
# Container for LD information used by fine-mapping / colocalization. Holds
# a pre-computed correlation matrix (or a list of per-block matrices) and/or
# a GenotypeHandle (or list of handles for mixture panels) for on-demand
# correlation computation.
# =============================================================================

#' @include GenotypeHandle.R AllClasses.R
NULL

#' @title LD Data Container
#' @description S4 container for LD information. Stores either a pre-computed
#'   correlation matrix or a \code{GenotypeHandle} (or list of handles for
#'   mixture panels) for lazy genotype/correlation access.
#'
#'   An \code{LdData} \strong{is} a \code{GRanges} over the variants it
#'   covers -- carrying A1, A2, variant_id and optionally allele_freq,
#'   variance and n_nomiss as metadata columns -- so \code{length()},
#'   \code{seqnames()} and \code{start()} answer directly, as they do for
#'   \code{\link{LdScore}} and \code{\link{LdEigen}}. The correlation may
#'   be NULL, in which case those ranges are the object's only record of which
#'   variants it describes.
#'
#' @slot correlation A correlation matrix, a list of per-block matrices
#'   (block-diagonal LD), or NULL if genotypes are available and R should be
#'   computed on demand.
#' @slot genotypeHandle Where genotypes are read from: a genotype handle, a
#'   list of them (for mixture panels), a matrix of dosages already extracted
#'   and filtered, or NULL when only pre-computed R is available. Pass a
#'   genotype panel (see \code{\link{readGenotypes}}) to the constructor and
#'   it is unwrapped to its handle.
#' @slot snpIdx Integer vector of 1-based SNP indices into the handle's
#'   \code{snpInfo}. NULL when correlation is pre-computed, or when the
#'   source is a matrix (which is already the subset).
#' @slot blockMetadata A \code{GRanges} with one range per block: its
#'   \code{seqnames} and range are the block's chromosome and genomic span,
#'   and \code{mcols} carry \code{blockId}, \code{size}, \code{startIdx}
#'   and \code{endIdx} -- the index range into the correlation matrix the
#'   block addresses. Extra columns the caller supplied are kept as mcols.
#'   A block whose index range is out of bounds has no span to derive and
#'   gets a width-0 range.
#' @slot nRef Integer, reference panel sample size.
#' @slot mixtureWeights NULL when \code{genotypeHandle} is a single
#'   \code{GenotypeHandle}; a numeric vector of mixing proportions (one per
#'   panel, summing to 1) when \code{genotypeHandle} is a list of
#'   \code{GenotypeHandle}s. Used by \code{getCorrelation()} to compute a
#'   weighted-average mixture LD matrix; required whenever \code{genotypeHandle}
#'   is a list and \code{getCorrelation()} will be called.
#' @export
setClass(
    "LdData",
    contains = "GRanges",
    representation(
        correlation = "LdCorrelation",
        genotypeHandle = "LdGenotypeSource",
        snpIdx = "LdSnpIndex",
        blockMetadata = "GRanges",
        nRef = "integer",
        mixtureWeights = "LdMixtureWeights"
    ),
    validity = function(object) {
        errors <- c(
            if (is.null(object@correlation) && is.null(object@genotypeHandle)) {
                str_c(
                    "At least one of 'correlation' or ",
                    "'genotypeHandle' must be non-NULL"
                )
            },
            if (length(object) == 0) "an LdData must cover >= 1 variant",
            .ldCheckGenotypeSource(object@genotypeHandle),
            .ldCheckCorrelation(object@correlation),
            .ldCheckMixtureWeights(object),
            .ldCheckBlockMetadata(object)
        )
        if (length(errors) == 0) TRUE else errors
    }
)

# Mixture weights are only meaningful over a LIST of panels, and must then be
# a proper simplex over them. NULL weights are always valid.
# @noRd
.ldCheckMixtureWeights <- function(object) {
    if (is.null(object@mixtureWeights)) {
        return(NULL)
    }
    if (!is.list(object@genotypeHandle)) {
        return(str_c(
            "'mixtureWeights' may only be set when ",
            "'genotypeHandle' is a list of panels"
        ))
    }
    w <- object@mixtureWeights
    if (!is.numeric(w) || length(w) != length(object@genotypeHandle)) {
        return(str_c(
            "'mixtureWeights' must be numeric of length ",
            "equal to the genotypeHandle list"
        ))
    }
    if (any(w < 0) || abs(sum(w) - 1) > 1e-6) {
        return("'mixtureWeights' must be non-negative and sum to 1")
    }
    NULL
}

#' @describeIn LdData-class Refused. Subsetting would narrow the variants
#'   while \code{correlation} -- variant-by-variant, and \code{snpIdx}, which
#'   indexes the reference panel -- stayed as they were, leaving an LD matrix
#'   describing variants the object no longer has. Build the LdData over the
#'   variant set you want instead.
#' @param x An \code{LdData}.
#' @param i,j,... Subscripts; any use is an error.
#' @param drop Ignored.
#' @return Nothing: this method always signals an error.
#' @export
setMethod("[", "LdData", function(x, i, j, ..., drop = TRUE) {
    abort(glue(
        "an LdData cannot be subset: its correlation is variant-by-variant ",
        "and its snpIdx addresses the reference panel, so narrowing the ",
        "variants would leave both describing variants that are no longer ",
        "present. Build the LdData over the variant set you want instead."
    ))
})

#' @rdname show-methods
#' @export
setMethod("show", "LdData", function(object) {
    n_var <- length(object)
    has_R <- !is.null(object@correlation)
    has_geno <- !is.null(object@genotypeHandle)
    r_type <- if (has_R && is.list(object@correlation)) {
        "block-diagonal"
    } else {
        "single"
    }
    cat(glue("LdData: {n_var} variants\n", .trim = FALSE))
    cat(glue(
        "  Correlation: {if (has_R) r_type else 'NULL'}, ",
        "Genotype handle: {if (has_geno) 'available' else 'NULL'}\n",
        .trim = FALSE
    ))
    cat(glue("  Reference N: {object@nRef}\n", .trim = FALSE))
})

# One element of a mixture list: unwrap a panel, pass anything else through
# for validity to judge.
# @noRd
.ldDataSourceElement <- function(x) {
    open <- .openGenotypeHandle(x)
    if (is.null(open)) x else open
}

# Normalise whatever the caller passed for `genotypeHandle` into one of the
# shapes the slot admits. A panel is the public shape (readGenotypes()); the
# slot stores the handle behind it, because `snpIdx` selects into the handle's
# whole snpInfo and a panel adds nothing the handle does not already carry.
# Anything else is passed through for the slot's class union to accept or
# reject -- this is a normaliser, not a second type check.
# @noRd
.ldDataGenotypeSource <- function(x) {
    if (is.null(x) || is.matrix(x)) {
        return(x)
    }
    if (is.list(x)) {
        # Unwrap the panels but do not reject here: validity reports a bad
        # mixture element by its index, which a purrr-wrapped abort would bury.
        return(map(x, .ldDataSourceElement))
    }
    open <- .openGenotypeHandle(x)
    if (!is.null(open)) {
        return(open)
    }
    # The slot's class union would reject this anyway, but its message names
    # the union rather than what to pass instead.
    hint <- if (methods::is(x, "DelayedArray")) {
        str_c(
            " -- pass the panel it came from rather than its assay, so the ",
            "variant selection in `snpIdx` still means something"
        )
    } else {
        ""
    }
    abort(glue(
        "`genotypeHandle` must be a genotype panel from readGenotypes(), a ",
        "list of them, a matrix of dosages, or NULL (got ",
        "{class(x)[[1L]]}){hint}."
    ))
}

# --- block metadata --------------------------------------------------------
#
# Five paths used to build this table -- the two region loaders, the two
# bare-matrix wrappers, and callers handing over a bare genomic span -- and
# each emitted a different set of columns, so consumers read `startIdx`,
# `size` and `chrom` that were there by luck rather than by contract. One
# normaliser at the only door (LdData()) fills whatever is absent and stores
# a single shape, the way .ldDataGenotypeSource() does for the handle.
#
# These are the columns consumers actually read. `blockStart` / `blockEnd`
# are deliberately NOT required: the two region loaders record a block's
# genomic span, nothing reads it, and demanding it would reject every
# multi-block table built without it. Supplied columns are carried through.
# @noRd
.ldBlockColumns <- c("blockId", "chrom", "size", "startIdx", "endIdx")

# What survives as mcols once chrom / blockStart / blockEnd become the
# GRanges itself.
# @noRd
.ldBlockMcols <- c("blockId", "size", "startIdx", "endIdx")

# Supplied values are never replaced: a caller handing over deliberately
# out-of-range indices keeps them, because partitionLdMatrix() still has to
# reject those and substituting valid ones would hide the error.
# @noRd
.ldBlockMetadata <- function(blockMetadata, variants) {
    if (length(variants) == 0L) {
        # Nothing to derive a block from. Validity rejects the object for
        # covering no variant; crashing here would pre-empt that message.
        return(GenomicRanges::GRanges())
    }
    tbl <- .ldBlockMetadataTibble(blockMetadata)
    if (nrow(tbl) == 0L) {
        # A placeholder table means "one block", not "no blocks": an LdData
        # always covers at least one variant.
        tbl <- tibble(blockId = 1L)
    }
    tbl <- .ldBlockMetadataIndices(tbl, length(variants))
    tbl <- .ldBlockMetadataOrder(.ldBlockMetadataDerive(tbl, variants))
    .ldBlockMetadataGRanges(tbl, variants)
}

# The stored shape: one range per block, so the invariants the merge logic
# depends on are the ones GRanges already answers -- isDisjoint() for
# overlap, sort() for order, seqnames() for the same-chromosome guard. The
# index payload the matrix is addressed by rides along in mcols.
# @noRd
.ldBlockMetadataGRanges <- function(tbl, variants) {
    span <- .ldBlockSpans(tbl, variants)
    gr <- GenomicRanges::GRanges(
        seqnames = tbl$chrom,
        ranges = IRanges::IRanges(start = span$start, end = span$end)
    )
    keep <- c(
        intersect(.ldBlockMcols, names(tbl)),
        setdiff(names(tbl), c(.ldBlockColumns, "blockStart", "blockEnd"))
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        as.data.frame(select(tbl, all_of(keep)))
    )
    gr
}

# Each block's genomic span. Supplied blockStart / blockEnd win; otherwise
# the span is read off the variants the index range covers.
#
# A block whose index range is out of bounds has no span to read, and a
# GRanges cannot hold an NA one (seqnames reject NAs outright). It gets a
# WIDTH-0 range at the first variant instead: legal, and visibly not an
# interval, so nothing downstream mistakes it for a real span.
# partitionLdMatrix() rejects the block on its indices regardless.
# @noRd
.ldBlockSpans <- function(tbl, variants) {
    if (is_in("blockStart", names(tbl)) && is_in("blockEnd", names(tbl))) {
        return(list(start = tbl$blockStart, end = tbl$blockEnd))
    }
    parts <- map(
        seq_len(nrow(tbl)),
        .ldBlockSpanAt,
        tbl = tbl,
        variants = variants
    )
    list(
        start = map_dbl(parts, "start"),
        end = map_dbl(parts, "end")
    )
}

# @noRd
.ldBlockSpanAt <- function(i, tbl, variants) {
    s <- tbl$startIdx[[i]]
    e <- tbl$endIdx[[i]]
    if (.ldBlockRangeUsable(s, e, length(variants))) {
        part <- variants[seq.int(s, e)]
        return(list(
            start = min(GenomicRanges::start(part)),
            end = max(GenomicRanges::end(part))
        ))
    }
    anchor <- GenomicRanges::start(variants)[[1L]]
    list(start = anchor, end = anchor - 1L)
}

# @noRd
.ldBlockRangeUsable <- function(s, e, n) {
    !is.na(s) && !is.na(e) && s >= 1L && e >= s && e <= n
}

# Any accepted shape as a tibble. `start`/`end` is how a bare genomic span
# arrives; they name the BLOCK's span, so they become blockStart/blockEnd.
# @noRd
.ldBlockMetadataTibble <- function(x) {
    if (!is(x, "GRanges") && !is.data.frame(x) && !is(x, "DataFrame")) {
        # as.data.frame() would happily turn a string into a one-column
        # table, so the shapes have to be named rather than coerced blindly.
        abort(glue(
            "`blockMetadata` must be a GRanges, a data.frame or a ",
            "DataFrame; got {class(x)[[1L]]}"
        ))
    }
    tbl <- if (is(x, "GRanges")) {
        .ldBlockMetadataFromGranges(x)
    } else {
        as_tibble(as.data.frame(x))
    }
    .ldBlockMetadataRenameSpan(tbl)
}

# @noRd
.ldBlockMetadataFromGranges <- function(x) {
    md <- as_tibble(as.data.frame(S4Vectors::mcols(x)))
    # An LdBlocks-style GRanges carries chrom / blockStart / blockEnd in its
    # mcols as well as in its ranges. The mcols win: binding both would
    # duplicate the names and leave no column called `chrom` at all.
    fromRanges <- tibble(
        chrom = as.character(GenomicRanges::seqnames(x)),
        blockStart = GenomicRanges::start(x),
        blockEnd = GenomicRanges::end(x)
    )
    bind_cols(
        md,
        select(
            fromRanges,
            all_of(setdiff(names(fromRanges), names(md)))
        )
    )
}

# @noRd
.ldBlockMetadataRenameSpan <- function(tbl) {
    if (!is_in("blockStart", names(tbl)) && is_in("start", names(tbl))) {
        tbl <- rename(tbl, blockStart = "start")
    }
    if (!is_in("blockEnd", names(tbl)) && is_in("end", names(tbl))) {
        tbl <- rename(tbl, blockEnd = "end")
    }
    tbl
}

# The index range is the one thing that cannot be derived for several
# blocks: only the caller knows which variants each covers. A single block
# covers all of them.
# @noRd
.ldBlockMetadataIndices <- function(tbl, n) {
    absent <- setdiff(c("startIdx", "endIdx"), names(tbl))
    if (length(absent) == 0L) {
        return(tbl)
    }
    if (length(absent) == 1L) {
        abort(glue(
            "`blockMetadata` gives one of startIdx / endIdx but not the ",
            "other; supply both or neither."
        ))
    }
    if (is_in("size", names(tbl))) {
        # Blocks are taken in order, so the sizes give the ranges. This is
        # the LdBlocks shape: a block table that says how many variants each
        # block holds but not where they sit in the matrix.
        sizes <- as.integer(tbl$size)
        ends <- cumsum(sizes)
        return(mutate(tbl, startIdx = ends - sizes + 1L, endIdx = ends))
    }
    if (nrow(tbl) > 1L) {
        abort(glue(
            "`blockMetadata` describes {nrow(tbl)} blocks but gives neither ",
            "an index range nor `size`, so which variants each block covers ",
            "cannot be determined."
        ))
    }
    mutate(tbl, startIdx = 1L, endIdx = n)
}

# `blockId`, `chrom` and `size` all follow from the index range, for one
# block or many.
# @noRd
.ldBlockMetadataDerive <- function(tbl, variants) {
    absent <- setdiff(.ldBlockColumns, names(tbl))
    if (length(absent) == 0L) {
        return(tbl)
    }
    derived <- tibble(
        blockId = seq_len(nrow(tbl)),
        chrom = .ldBlockChroms(tbl, variants),
        size = as.integer(tbl$endIdx - tbl$startIdx + 1L),
        startIdx = tbl$startIdx,
        endIdx = tbl$endIdx
    )
    bind_cols(tbl, select(derived, all_of(absent)))
}

# @noRd
.ldBlockChroms <- function(tbl, variants) {
    map_chr(
        seq_len(nrow(tbl)),
        .ldBlockChromAt,
        tbl = tbl,
        variants = variants
    )
}

# An out-of-range block falls back to the data's own chromosome rather than
# NA: seqnames cannot hold NA, and partitionLdMatrix() rejects the block on
# its indices anyway. Its range is width 0 (see .ldBlockSpanAt).
# @noRd
.ldBlockChromAt <- function(i, tbl, variants) {
    s <- tbl$startIdx[[i]]
    e <- tbl$endIdx[[i]]
    part <- if (.ldBlockRangeUsable(s, e, length(variants))) {
        variants[seq.int(s, e)]
    } else {
        variants
    }
    as.character(GenomicRanges::seqnames(part))[[1L]]
}

# Canonical columns first, anything the caller added after.
# @noRd
.ldBlockMetadataOrder <- function(tbl) {
    select(tbl, all_of(.ldBlockColumns), everything())
}

# Every LdData describes at least one block, with the columns its consumers
# read. The normaliser guarantees both; this is the backstop that keeps a
# future construction path from quietly dropping one.
# @noRd
.ldCheckBlockMetadata <- function(object) {
    bm <- object@blockMetadata
    absent <- setdiff(.ldBlockMcols, names(S4Vectors::mcols(bm)))
    if (length(absent) > 0L) {
        return(str_c(
            "'blockMetadata' is missing required column(s): ",
            str_flatten(absent, ", ")
        ))
    }
    if (length(bm) == 0L) {
        return("'blockMetadata' must describe at least one block")
    }
    NULL
}

#' @title Create an LdData Object
#' @description Construct an \code{LdData} from a correlation matrix and/or
#'   genotype handle, plus variant metadata as a GRanges.
#' @param correlation A correlation matrix, list of matrices, or NULL.
#' @param genotypeHandle A genotype panel (see
#'   \code{\link{readGenotypes}}), a list of panels for a mixture reference,
#'   a matrix of already-extracted dosages, or NULL.
#' @param snpIdx Vector of 1-based SNP indices, coerced to integer, or NULL.
#' @param variants A GRanges with variant metadata (must have variant_id in
#'   mcols, plus A1, A2).
#' @param blockMetadata Block boundaries: a \code{GRanges}, a
#'   \code{data.frame} or a \code{DataFrame}. Whatever is absent for a
#'   single block is derived from \code{variants}, so a bare genomic span
#'   -- or nothing at all -- is enough; with several blocks either the index
#'   columns or \code{size} must be supplied. Normalised on the way in and
#'   stored as a \code{GRanges} (see the \code{blockMetadata} slot).
#' @param nRef Integer, reference panel sample size.
#' @param mixtureWeights Optional numeric vector of mixing proportions, one per
#'   panel in \code{genotypeHandle} when it is a list. Must be non-negative and
#'   sum to 1. Required whenever \code{genotypeHandle} is a list and downstream
#'   code will call \code{getCorrelation()}.
#' @return An \code{LdData} object.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' ld
#' @importFrom checkmate assert checkMatrix checkList checkNull
#' @importFrom checkmate assertNumeric
#' @export
LdData <- function(
    correlation = NULL,
    genotypeHandle = NULL,
    snpIdx = NULL,
    variants,
    blockMetadata,
    nRef = 0L,
    mixtureWeights = NULL
) {
    # correlation is documented as a matrix OR a list of matrices OR NULL,
    # so this must be an or-combination, not assertMatrix.
    assert(
        checkMatrix(correlation),
        checkList(correlation),
        checkNull(correlation),
        .var.name = "correlation"
    )
    assertNumeric(mixtureWeights, null.ok = TRUE)
    obj <- new(
        "LdData",
        variants,
        correlation = correlation,
        genotypeHandle = .ldDataGenotypeSource(genotypeHandle),
        snpIdx = if (is.null(snpIdx)) NULL else as.integer(snpIdx),
        blockMetadata = .ldBlockMetadata(blockMetadata, variants),
        nRef = as.integer(nRef),
        mixtureWeights = mixtureWeights
    )
    validObject(obj)
    obj
}

# Internal: convert a refPanel data.frame (chrom/pos/A1/A2/variant_id, with
# optional allele_freq/variance/n_nomiss) into the GRanges an LdData is.
.refPanelToGranges <- function(refPanel) {
    chr <- withChrPrefix(refPanel$chrom)
    pos <- as.integer(refPanel$pos)

    gr <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = pos, width = 1L)
    )

    mcolsData <- DataFrame(
        variant_id = refPanel$variant_id,
        A1 = refPanel$A1,
        A2 = refPanel$A2
    )

    optional <- intersect(
        c("allele_freq", "variance", "n_nomiss"),
        names(refPanel)
    )
    `mcols<-`(
        gr,
        value = cbind(
            mcolsData,
            DataFrame(refPanel[optional], check.names = FALSE)
        )
    )
}

# One panel scaled by its mixture weight.
# @noRd
.ldScalePanel <- function(w, panel) {
    w * panel
}

#' @rdname getCorrelation
#' @export
setMethod("getCorrelation", "LdData", function(x) {
    if (!is.null(x@correlation)) {
        return(x@correlation)
    }
    if (is.null(x@genotypeHandle)) {
        abort("No correlation matrix or genotype handle available")
    }
    if (is.list(x@genotypeHandle)) {
        if (is.null(x@mixtureWeights)) {
            msg <- glue(
                "Cannot compute mixture LD: `mixtureWeights` is NULL. ",
                "Construct LdData with mixtureWeights = <numeric vector> ",
                "when supplying a list of GenotypeHandles."
            )
            abort(msg)
        }
        perPanel <- map(x@genotypeHandle, .ldPanelLd, snpIdx = x@snpIdx)
        dims <- map_int(perPanel, nrow)
        if (length(unique(dims)) != 1L) {
            msg <- glue(
                "Mixture panels yielded LD matrices of differing ",
                "dimensions: {str_flatten(dims, ', ')}. All panels ",
                "must be aligned on the same variant subset."
            )
            abort(msg)
        }
        # The mixture is a weighted sum over the panels, so it is a fold.
        weighted <- map2(x@mixtureWeights, perPanel, .ldScalePanel)
        return(`dimnames<-`(
            reduce(
                weighted,
                `+`,
                .init = matrix(0, nrow = dims[[1L]], ncol = dims[[1L]])
            ),
            dimnames(perPanel[[1L]])
        ))
    }
    computeLd(.ldSourceDosages(x@genotypeHandle, x@snpIdx), method = "sample")
})

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "LdData", function(x) {
    if (is.null(x@genotypeHandle)) {
        return(NULL)
    }
    if (is.list(x@genotypeHandle)) {
        # A mixture list may hold matrices as well as handles.
        map(x@genotypeHandle, .ldSourceDosages, x@snpIdx)
    } else {
        .ldSourceDosages(x@genotypeHandle, x@snpIdx)
    }
})

#' @rdname hasGenotypes
#' @export
setMethod("hasGenotypes", "LdData", function(x) {
    !is.null(getGenotypeHandle(x))
})

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "LdData", function(x) {
    mcols(x)$variant_id
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "LdData", function(x) {
    # A plain GRanges view: callers subset and re-wrap it, and carrying the
    # LD payload along would make those copies quietly expensive.
    as(x, "GRanges")
})

#' @rdname getBlockMetadata
#' @export
setMethod("getBlockMetadata", "LdData", function(x) {
    x@blockMetadata
})

#' @rdname getRefPanel
#' @export
setMethod("getRefPanel", "LdData", function(x) {
    mutate(
        as_tibble(as.data.frame(mcols(x))),
        chrom = as.character(seqnames(x)),
        pos = start(x)
    )
})

#' @rdname getGenotypeHandle
#' @keywords internal
setMethod("getGenotypeHandle", "LdData", function(x) x@genotypeHandle)

#' @rdname getMixtureWeights
#' @export
setMethod("getMixtureWeights", "LdData", function(x) x@mixtureWeights)

#' @rdname getSnpIdx
#' @export
setMethod("getSnpIdx", "LdData", function(x) x@snpIdx)

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdData", function(x) x@nRef)

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Sample-LD matrix for one mixture panel's genotype handle (over `snpIdx`).
# @noRd
.ldPanelLd <- function(h, snpIdx) {
    computeLd(.ldSourceDosages(h, snpIdx), method = "sample")
}

# Dosages from a genotype source. A matrix IS the dosages already -- extracted
# and filtered upstream, with `snpIdx` NULL because the matrix is the subset --
# so it is returned as-is rather than fed to the file readers, which would
# dispatch on it and fail.
# @noRd
.ldSourceDosages <- function(x, snpIdx) {
    if (is.matrix(x)) {
        return(x)
    }
    .dosageMatrix(x, snpIdx)
}

# Every element of a mixture list has to be something LD can be computed from.
# The class union admits any list, so the elements are checked here.
# @noRd
.ldCheckGenotypeSource <- function(x) {
    if (!is.list(x)) {
        return(character())
    }
    ok <- map_lgl(x, .ldIsGenotypeSource)
    if (all(ok)) {
        return(character())
    }
    str_c(
        "'genotypeHandle' list elements must each be a genotype panel or a ",
        "dosage matrix; element(s) ",
        str_flatten(which(!ok), ", "),
        " are not."
    )
}

# @noRd
.ldIsGenotypeSource <- function(x) {
    methods::is(x, "GenotypeHandle") || is.matrix(x)
}

# Block-diagonal LD is one matrix per block. The union admits any list, so the
# elements are checked here -- a list of something else would otherwise sail
# through and fail later inside the LD arithmetic.
# @noRd
.ldCheckCorrelation <- function(x) {
    if (!is.list(x)) {
        return(character())
    }
    ok <- map_lgl(x, is.matrix)
    if (all(ok)) {
        return(character())
    }
    str_c(
        "'correlation' list elements must each be a matrix (block-diagonal ",
        "LD is one matrix per block); element(s) ",
        str_flatten(which(!ok), ", "),
        " are not."
    )
}
