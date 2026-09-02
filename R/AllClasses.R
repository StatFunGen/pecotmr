# =============================================================================
# AllClasses.R
# -----------------------------------------------------------------------------
# Virtual base classes shared across the package. Concrete subclasses live
# in their own per-class files (QtlSumStats.R, GwasSumStats.R, QtlDataset.R,
# QtlFineMappingResult.R, GwasFineMappingResult.R, etc.).
#
# Per Bioconductor convention this file is loaded first in the Collate
# ordering (the "AllClasses.R" filename sorts to the top of the alphabet),
# and every method-bearing file uses `@include AllClasses.R` so roxygen
# topologically orders the Collate field for us.
# =============================================================================

#' @include AllGenerics.R GenotypeHandle.R RangedTupleList.R
#' @importFrom methods setClass setMethod new is validObject
NULL

# The LD reference panel a collection was harmonized against, as the same
# RangedSummarizedExperiment shape QtlDataset uses for its genotypes: variants
# in rowRanges, dosages in a DelayedArray assay that reads lazily through a
# GenotypeHandle. NULL when a collection carries no panel.
#
# A union rather than "ANY": the slot has always been documented as holding a
# panel, and an untyped slot with a typed docstring is the kind of thing that
# is only true until someone puts something else in it.
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setClassUnion(
    "LdSketchOrNULL",
    c("RangedSummarizedExperiment", "NULL")
)

# What an LdData can read genotypes from. Four shapes, all real:
#
#   GenotypeHandle  a lazy file-backed panel; `snpIdx` selects into its whole
#                   snpInfo, so the handle must not be pre-narrowed
#   list            one handle per panel for a mixture reference, averaged by
#                   `mixtureWeights`
#   matrix          dosages already extracted and filtered, kept so
#                   getGenotypes() answers without reopening the file (see
#                   .loadLdFromBlocks); `snpIdx` is NULL in this case because
#                   the matrix is already the subset
#   NULL            no genotypes -- the object carries a pre-computed R
#
# A panel (RangedSummarizedExperiment) is what callers pass; the constructor
# unwraps it to its handle, so the slot itself never holds one.
#
# A union rather than "ANY", for the same reason as above: the slot's
# docstring named two of these four shapes and nothing enforced even that.
setClassUnion(
    "LdGenotypeSource",
    c("GenotypeHandle", "matrix", "list", "NULL")
)

# The rest of LdData's payload, typed for the same reason. Each union is the
# set of shapes the class is actually built with, confirmed by recording slot
# classes at validity across the LD test files (validity runs for every
# object, however it was constructed).
#
# `correlation` is a single matrix, or one matrix per block for
# block-diagonal LD, or NULL when it has to be computed from genotypes.
setClassUnion("LdCorrelation", c("matrix", "list", "NULL"))

# `snpIdx` selects into the genotype source's snpInfo; NULL when the
# correlation is pre-computed, or when the source is already the subset. The
# constructor coerces to integer, so a caller may pass doubles.
setClassUnion("LdSnpIndex", c("integer", "NULL"))

# Block boundaries, as ranges or as a table. Not nullable: every LdData is
# built for some block, and the constructor has always required it.
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
setClassUnion("LdBlockMetadata", c("GRanges", "data.frame", "DataFrame"))

# Mixing proportions, one per panel, when `genotypeHandle` is a list.
setClassUnion("LdMixtureWeights", c("numeric", "NULL"))

# =============================================================================
# SumStatsBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS summary statistics collections.
# Concrete subclasses (QtlSumStats, GwasSumStats) inherit from
# RangedTupleList and share the ldSketch / qcInfo slots (the genome build
# lives in seqinfo, not a slot). Each element
# is one tuple's per-variant GRanges: x[[i]], formerly x$entry[[i]].
#
# getZ / getN / getMaf / nSnps are
# defined once on SumStatsBase (they only delegate to getSumStats); subsetChr /
# getVarY / getSumStats / getSumstatDf stay on the concrete subclass because
# they rely on the tuple shape (3-tuple QtlSumStats, 1-tuple GwasSumStats).
# =============================================================================

#' @title Summary Statistics Base Class
#' @description Virtual base class for QTL and GWAS summary statistics
#'   collections. Concrete subclasses (\code{QtlSumStats}, \code{GwasSumStats})
#'   inherit from \code{\linkS4class{RangedTupleList}} and share the
#'   \code{ldSketch} / \code{qcInfo} slots, and the genome build in
#'   \code{seqinfo()}.
#'
#'   Each element is the per-variant \code{GRanges} of one tuple, so
#'   \code{x[[i]]} is that tuple's summary statistics and the identity columns
#'   live in \code{mcols(x)}. There is no \code{entry} column: what used to be
#'   \code{x$entry[[i]]} is now simply \code{x[[i]]}.
#' @slot ldSketch The \code{GenotypeHandle} the QC pipeline harmonized against,
#'   or \code{NULL}. Optional: LD-free workflows (e.g. mash, which operates
#'   across conditions per variant) carry \code{NULL}; pipelines that need LD
#'   validate its presence when they consume the collection.
#' @slot qcInfo A \code{list} recording which QC steps ran. Empty \code{list()}
#'   on construction; populated by \code{summaryStatsQc()} with a per-step audit
#'   record (filter names, drop counts, liftover target, RAISS settings, etc.).
#'   Fine-mapping and TWAS-weights pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0L} -- the slot serves as both the gating
#'   flag and the audit trail.
#' @export
setClass(
    "SumStatsBase",
    contains = c("VIRTUAL", "RangedTupleList"),
    representation(
        ldSketch = "LdSketchOrNULL",
        qcInfo = "list"
    )
)

# TRUE when two elements share both their identity tuple AND their span.
# Splitting a multi-chromosome study makes the tuple alone non-unique, so the
# key is the tuple plus the element's range (spec 4.4: uniqueness is enforced
# on `(identity tuple..., range)`).
# @noRd
.ssHasDuplicateKeys <- function(object, tupleCols) {
    if (nrow(object) == 0L) {
        return(FALSE)
    }
    md <- mcols(object)
    tupleKey <- map(tupleCols, .ssKeyColumn, md = md)
    spans <- range(object)
    rangeKey <- map_chr(as.list(spans), .ssSpanLabel)
    keys <- exec(str_c, !!!c(tupleKey, list(rangeKey)), sep = "|")
    anyDuplicated(keys) > 0L
}

# @noRd
.ssKeyColumn <- function(cn, md) {
    as.character(md[[cn]])
}

# A stable label for one element's span; empty elements collapse to "".
# @noRd
.ssSpanLabel <- function(g) {
    if (length(g) == 0L) {
        return("")
    }
    str_c(
        as.character(seqnames(g))[[1L]],
        ":",
        min(start(g)),
        "-",
        max(end(g))
    )
}

# Stitch a tuple's elements back into one GRanges, optionally restricted to a
# region. The seqname split means one tuple can own several elements, but the
# accessor contract is still "one GRanges per tuple"; `ranges` lets a caller
# pull just the part they want instead of materialising the whole span.
# @noRd
.ssStitchElements <- function(x, idx, ranges = NULL) {
    gr <- .rtlGatherElements(x, idx)
    if (is.null(ranges)) {
        return(gr)
    }
    win <- .asGRegion(ranges)
    onWindowChrom <- as.character(seqnames(gr)) %in%
        as.character(seqnames(win))
    if (!any(onWindowChrom)) {
        return(gr[0L])
    }
    gr[onWindowChrom & IRanges::overlapsAny(gr, win)]
}

# The build recorded in seqinfo must be exactly one non-NA value. A missing
# build is an error rather than a default: every downstream liftover / LD
# join keys on it, and silently guessing hg38 is how a mismatched panel gets
# through. Mixed builds mean the parts were never comparable.
# @noRd
.sumStatsCheckGenome <- function(object) {
    # seqinfo records the build per SEQLEVEL, so a collection spanning none
    # has nowhere to keep one -- whether it has no elements at all or only
    # empty ones (a PIP-screened region emptied by summaryStatsQc). That is
    # the one case where a missing build is not a defect: there is nothing
    # for it to describe. A subset that empties an existing collection keeps
    # its seqinfo (see .rtlRebuild), so this exempts only what was built with
    # no ranges in the first place.
    if (length(GenomeInfoDb::seqlevels(object)) == 0L) {
        return(NULL)
    }
    build <- unique(GenomeInfoDb::genome(object))
    build <- build[!is.na(build)]
    if (length(build) == 1L && str_length(build) > 0L) {
        return(NULL)
    }
    if (length(build) == 0L) {
        return("no genome build in seqinfo(); set one with genome(x) <- ...")
    }
    str_c(
        "seqinfo() names more than one genome build (",
        str_flatten(build, ", "),
        ")"
    )
}

#' @rdname getGenome
#' @examples
#' data(qtlSumStatsExample)
#' getGenome(qtlSumStatsExample)
#' @export
setMethod("getGenome", "SumStatsBase", function(x, ...) {
    # The build lives in seqinfo, exactly as it does on LdStatistic: a
    # GRangesList already has somewhere to keep it, and a parallel `genome`
    # slot went stale against it (getGenome() said hg19 while genome(x) said
    # NA, so every Bioconductor path that reads genome(x) saw nothing).
    build <- unique(GenomeInfoDb::genome(x))
    build <- build[!is.na(build)]
    if (length(build) == 0L) NA_character_ else build[[1L]]
})

#' @rdname getQcInfo
#' @examples
#' data(qtlSumStatsExample)
#' getQcInfo(qtlSumStatsExample)
#' @export
setMethod("getQcInfo", "SumStatsBase", function(x, ...) x@qcInfo)

#' @rdname getQcDiagnostics
#' @examples
#' data(qtlSumStatsExample)
#' getQcDiagnostics(qtlSumStatsExample)
#' @export
setMethod("getQcDiagnostics", "SumStatsBase", function(x, entry = 1L, ...) {
    qc <- x@qcInfo
    if (length(qc) == 0L) {
        return(NULL)
    }
    audits <- qc$entryAudit
    if (is.null(audits)) {
        return(NULL)
    }
    if (is.null(entry)) {
        out <- map(audits, "ldMismatchDiagnostics")
        keep <- !map_lgl(out, is.null)
        if (!any(keep)) {
            return(NULL)
        }
        set_names(out[keep], seq_along(audits)[keep])
    } else {
        if (
            !is.numeric(entry) ||
                length(entry) != 1L ||
                entry < 1L ||
                entry > length(audits)
        ) {
            msg <- glue(
                "`entry` must be a single integer in 1:{length(audits)}."
            )
            abort(msg)
        }
        audits[[as.integer(entry)]]$ldMismatchDiagnostics
    }
})

#' @rdname getLdSketch
#' @examples
#' data(qtlSumStatsExample)
#' getLdSketch(qtlSumStatsExample)
#' @export
setMethod("getLdSketch", "SumStatsBase", function(x, ...) x@ldSketch)

#' @rdname getStudy
#' @examples
#' data(qtlDatasetExample)
#' getStudy(qtlDatasetExample)
#' @export
setMethod("getStudy", "SumStatsBase", function(x) unique(as.character(x$study)))

# getZ / getN / getMaf / nSnps delegate purely to getSumStats (which the
# concrete subclass dispatches with its own tuple shape), so they live once on
# the base and forward `...` through.
#' @rdname getZ
#' @examples
#' data(qtlSumStatsExample)
#' getZ(qtlSumStatsExample)
#' @export
setMethod("getZ", "SumStatsBase", function(x, ...) mcols(getSumStats(x, ...))$Z)

#' @rdname getN
#' @examples
#' data(qtlSumStatsExample)
#' getN(qtlSumStatsExample)
#' @export
setMethod("getN", "SumStatsBase", function(x, ...) mcols(getSumStats(x, ...))$N)

# getP / getBeta / getSe are first-class alongside getZ / getN: they read the
# optional P / BETA / SE mcols and return NULL when the entry does not carry
# them (DataFrame `$` semantics), so a p-value-primary sumstats (e.g. TensorQTL
# cis output) is an equal citizen to a Z-primary GWAS sumstats.
#' @rdname getP
#' @examples
#' data(qtlSumStatsExample)
#' getP(qtlSumStatsExample)
#' @export
setMethod("getP", "SumStatsBase", function(x, ...) mcols(getSumStats(x, ...))$P)

#' @rdname getBeta
#' @examples
#' data(qtlSumStatsExample)
#' getBeta(qtlSumStatsExample)
#' @export
setMethod("getBeta", "SumStatsBase", function(x, ...) {
    mcols(getSumStats(x, ...))$BETA
})

#' @rdname getSe
#' @examples
#' data(qtlSumStatsExample)
#' getSe(qtlSumStatsExample)
#' @export
setMethod("getSe", "SumStatsBase", function(x, ...) {
    mcols(getSumStats(x, ...))$SE
})

#' @rdname getMaf
#' @examples
#' data(qtlDatasetExample)
#' getMaf(qtlDatasetExample)
#' @export
setMethod("getMaf", "SumStatsBase", function(x, ...) {
    mc <- mcols(getSumStats(x, ...))
    if (is_in("MAF", colnames(mc))) mc$MAF else NULL
})

#' @rdname nSnps
#' @examples
#' data(qtlSumStatsExample)
#' nSnps(qtlSumStatsExample)
#' @export
setMethod("nSnps", "SumStatsBase", function(x, ...) length(getSumStats(x, ...)))

# =============================================================================
# FineMappingResultBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS fine-mapping result collections.
# Concrete subclasses (QtlFineMappingResult, GwasFineMappingResult) carry
# a DFrame of per-fit rows plus a shared ldSketch slot. Downstream
# pipelines dispatch on FineMappingResultBase for behaviors that apply to
# either flavour, and on the concrete subclass when the tuple shape
# matters.
# =============================================================================

#' @title Fine-Mapping Result Base Class
#' @description Virtual base class for fine-mapping result collections. Concrete
#'   subclasses (\code{QtlFineMappingResult}, \code{GwasFineMappingResult})
#'   carry a \code{DFrame} of per-fit rows and a shared \code{ldSketch} slot.
#'   Downstream pipelines should dispatch on \code{FineMappingResultBase} for
#'   behaviors that apply to either flavour, and on the concrete subclass when
#'   the tuple shape matters.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the fits were computed
#'   against, or \code{NULL} when the fits were derived from individual-level
#'   data (no LD reference). Used downstream for cross-pipeline LD-sketch
#'   identity validation.
#' @export
setClass(
    "FineMappingResultBase",
    contains = c("VIRTUAL", "RangedTupleList"),
    representation(ldSketch = "LdSketchOrNULL")
)

#' @rdname getStudy
#' @export
setMethod("getStudy", "FineMappingResultBase", function(x) {
    unique(as.character(x$study))
})

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "FineMappingResultBase", function(x, ...) x@ldSketch)

#' @rdname getMethodNames
#' @examples
#' data(qtlFineMappingExample)
#' getMethodNames(qtlFineMappingExample)
#' @export
setMethod("getMethodNames", "FineMappingResultBase", function(x) {
    unique(as.character(x$method))
})

#' @noRd
setMethod(
    "adjustPips",
    "FineMappingResultBase",
    function(x, keepVariants, ...) {
        if (nrow(x) == 0L) {
            return(x)
        }
        # Entries sharing no variant with `keepVariants` have nothing to
        # renormalize, so they are DROPPED -- matching subsetRegion's rule that
        # elements trimming to zero variants go away. Every other failure (a
        # fit that cannot be honestly subset, a slot whose width disagrees with
        # the variant count) PROPAGATES: silently leaving those entries
        # unadjusted mixes adjusted and unadjusted fits in one object, which is
        # the bug this replaces.
        overlaps <- map_lgl(
            .collectionEntries(x),
            .fmrEntryOverlaps,
            keepVariants = keepVariants
        )
        if (!any(overlaps)) {
            msg <- glue(
                "adjustPips: no entry shares a variant with `keepVariants`; ",
                "the two variant sets are disjoint."
            )
            abort(msg)
        }
        if (!all(overlaps)) {
            msg <- glue(
                "adjustPips: dropping {sum(!overlaps)} of {length(overlaps)} ",
                "entries that share no variant with `keepVariants`."
            )
            inform(msg)
        }
        out <- x[overlaps, ]
        # The entry is a derived view now, so adjusting means rebuilding each
        # element (and its fit payload) from the adjusted entry. The @listData
        # write this replaced is gone with the DFrame representation.
        adjusted <- map(
            .collectionEntries(out),
            adjustPips,
            keepVariants = keepVariants,
            ...
        )
        .fmrFromEntries(out, adjusted)
    }
)

# Rebuild a fine-mapping collection from adjusted entries, keeping every
# identity column and collection-level slot.
# @noRd
.fmrFromEntries <- function(x, entries) {
    grl <- GenomicRanges::GRangesList(map(entries, rowVariants))
    md <- mcols(x, use.names = FALSE)
    md$susieFit <- S4Vectors::SimpleList(map(entries, getSusieFit))
    md$cvResult <- S4Vectors::SimpleList(map(entries, getCvResult))
    mcols(grl) <- md
    new(class(x), grl, ldSketch = .asLdSketch(getLdSketch(x)))
}

# TRUE when an entry shares at least one variant with `keepVariants`, matched
# the same (chrom, pos, allele) way adjustPips() itself matches them.
# @noRd
.fmrEntryOverlaps <- function(entry, keepVariants) {
    matched <- matchVariants(
        .fmrPartsVariantIds(entry),
        as.character(keepVariants)
    )
    length(matched$idxA) > 0L
}

# =============================================================================
# Variant reconciliation between two fine-mapping collections
# -----------------------------------------------------------------------------
# Coloc needs both sides scored on the SAME variant set; TWAS / MR / cTWAS need
# only the QTL side adjusted to the GWAS one. Both go through adjustPips()
# above, which renormalizes each retained single-effect posterior over the
# variants that survive.
#
# What reconciliation does NOT do is recover the information the dropped
# variants carried. Coverage falls as overlap shrinks -- equally for a
# renormalization and for a full refit -- so that loss is irreducible, not an
# approximation error. getRetainedMass() is what makes it visible instead of
# letting a heavily-trimmed fit look as confident as an untouched one.
# =============================================================================

#' @rdname intersectVariants
#' @export
setMethod(
    "intersectVariants",
    signature(x = "FineMappingResultBase", y = "FineMappingResultBase"),
    function(x, y, oneSided = FALSE, ...) {
        shared <- .rcShared(x, y)
        if (isTRUE(oneSided)) {
            return(adjustPips(x, shared, ...))
        }
        list(
            x = adjustPips(x, shared, ...),
            y = adjustPips(y, shared, ...)
        )
    }
)

# The variants two collections have in common, matched allele-aware so a
# chr-prefix or separator difference does not read as no-overlap. Erroring on
# an empty intersection is deliberate: returning two empty collections would
# be indistinguishable from "reconciled fine, nothing colocalizes", which is a
# materially different scientific conclusion.
# @noRd
.rcShared <- function(x, y) {
    xv <- .rcAllVariants(x)
    yv <- .rcAllVariants(y)
    matched <- matchVariants(xv, yv)
    if (length(matched$idxA) == 0L) {
        msg <- glue(
            "intersectVariants: the two collections share no variants ",
            "({length(xv)} vs {length(yv)} distinct). Check that they were ",
            "built against the same genome and variant-id convention."
        )
        abort(msg)
    }
    xv[matched$idxA]
}

# Every distinct variant in a collection, in element order.
# @noRd
.rcAllVariants <- function(x) {
    if (nrow(x) == 0L) {
        return(character(0))
    }
    unique(.grVariantIds(unlist(x, use.names = FALSE)))
}

# =============================================================================
# Retained-mass diagnostic
# =============================================================================

#' @rdname getRetainedMass
#' @export
setMethod("getRetainedMass", "FineMappingResultBase", function(x, ...) {
    if (nrow(x) == 0L) {
        return(.rcEmptyMass(x))
    }
    parts <- map(seq_len(nrow(x)), .rcMassForRow, x = x)
    parts <- compact(parts)
    if (length(parts) == 0L) {
        return(.rcEmptyMass(x))
    }
    bind_rows(parts)
})

# The zero-row result, carrying the SAME identity columns the populated one
# would. Returning a bare (effect, retainedMass, nVariants) tibble instead
# would make the empty case a different shape from the non-empty case, so a
# caller that selects `study` breaks only when there is nothing to report --
# the worst time to find out. Which identity columns exist depends on the
# concrete class (a GWAS collection has no context / trait), so they are read
# off `x` rather than hard-coded.
# @noRd
.rcEmptyMass <- function(x) {
    cols <- names(.rcIdentityCols(x, integer(0)))
    ident <- set_names(rep(list(character(0)), length(cols)), cols)
    as_tibble(c(
        ident,
        list(
            effect = integer(0),
            retainedMass = numeric(0),
            nVariants = integer(0)
        )
    ))
}

# Per-effect retained mass for one element, with its identity columns
# attached.
#
# The mass is read off the STORED alpha: after a reconciliation each effect's
# row has been renormalized over the retained variants, so what is recoverable
# here is how concentrated that effect now is, not the pre-subset share. The
# pre-subset share is recorded at adjustment time (see .adjustPipsRetainedMass)
# and carried on the fit.
# @noRd
.rcMassForRow <- function(i, x) {
    fit <- mcols(x)$susieFit[[i]]
    mass <- .rcFitRetainedMass(fit)
    if (is.null(mass)) {
        return(NULL)
    }
    ident <- .rcIdentityCols(x, i)
    bind_cols(
        as_tibble(ident),
        as_tibble(list(
            effect = seq_along(mass),
            retainedMass = as.numeric(mass),
            nVariants = rep(length(x[[i]]), length(mass))
        ))
    )
}

# The retained mass a fit recorded when it was last adjusted, or NULL when the
# fit has never been through a reconciliation (nothing was dropped, so there is
# no mass to report rather than a vector of 1s implying a subset happened).
# @noRd
.rcFitRetainedMass <- function(fit) {
    if (is.null(fit) || !is.list(fit)) {
        return(NULL)
    }
    fit[["retained_mass"]]
}

# The element's identity columns, as a one-row list for recycling.
# @noRd
.rcIdentityCols <- function(x, i) {
    md <- mcols(x)
    cols <- intersect(
        c("study", "context", "trait", "method", "blockId"),
        colnames(md)
    )
    set_names(map(cols, .rcIdentityValue, md = md, i = i), cols)
}

# @noRd
.rcIdentityValue <- function(cn, md, i) {
    if (length(i) == 0L) {
        return(character(0))
    }
    as.character(md[[cn]])[[i]]
}

# Select the single FineMappingRow addressed by a (tuple / region) key. Each
# concrete subclass implements this with its own row selector; the delegating
# accessors below then live once on the base and route through it.
setGeneric(".fmrSelectEntry", function(x, ...) {
    standardGeneric(".fmrSelectEntry")
})

#' @rdname getCs
#' @export
setMethod(
    "getCs",
    "FineMappingResultBase",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        coverage = 0.95,
        minPurity = NULL,
        ...
    ) {
        # Selectors pinning one entry -> that entry's bare credible-set table;
        # no /
        # partial selectors -> aggregate every matching entry's credible sets,
        # prefixed with the row identity (study/context/trait/blockId/method).
        # `minPurity` is an independent CS-quality filter, orthogonal to
        # coverage.
        .fmrAggregateView(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = .fmrRowCs,
            coverage = coverage,
            minPurity = minPurity
        )
    }
)

#' @rdname getLbf
#' @export
setMethod("getLbf", "FineMappingResultBase", function(x, ...) {
    .fmrAggregateView(x, perEntry = .fmrRowLbf)
})

#' @rdname getCredibleSetSummary
#' @export
setMethod(
    "getCredibleSetSummary",
    "FineMappingResultBase",
    function(x, coverage = 0.95, ...) {
        .fmrAggregateView(
            x,
            perEntry = .fmrRowCredibleSetSummary,
            coverage = coverage
        )
    }
)

#' @rdname fsusieCredibleBand
#' @export
setMethod("fsusieCredibleBand", "FineMappingResultBase", function(x, ...) {
    .fmrAggregateView(x, perEntry = .fmrRowFsusieCredibleBand)
})

#' @rdname fsusieAffectedRegions
#' @export
setMethod("fsusieAffectedRegions", "FineMappingResultBase", function(x, ...) {
    grs <- map(seq_len(nrow(x)), .fsusieEntryAffectedRegions, x = x)
    grs <- grs[lengths(grs) > 0L]
    if (length(grs) == 0L) {
        return(GenomicRanges::GRanges())
    }
    exec(c, !!!grs)
})

#' @rdname getTopLoci
#' @export
setMethod(
    "getTopLoci",
    "FineMappingResultBase",
    function(
        x,
        type = c("data.frame", "GRanges"),
        signalCutoff = 0.025,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        minPurity = NULL,
        ...
    ) {
        type <- arg_match(type)
        # type = "GRanges" is honored only for a single pinned entry; the
        # aggregate (identity-prefixed) form is data.frame-only.
        if (type == "GRanges") {
            return(.fmrbTopLociGranges(
                x,
                study,
                context,
                trait,
                method,
                region,
                signalCutoff,
                minPurity
            ))
        }
        # data.frame: bare per-variant table when selectors pin one entry, else
        # the matching rows' per-variant tables stacked with row-identity cols.
        .fmrAggregateView(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = .fmrRowTopLoci,
            type = "data.frame",
            signalCutoff = signalCutoff,
            minPurity = minPurity
        )
    }
)

# Resolve the single pinned entry and return its GRanges view (aggregating
# across multiple entries requires type = "data.frame").
# @noRd
.fmrbTopLociGranges <- function(
    x,
    study,
    context,
    trait,
    method,
    region,
    signalCutoff,
    minPurity
) {
    sel <- tryCatch(
        .fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ),
        error = function(e) e
    )
    if (inherits(sel, "error")) {
        msg <- glue(
            "getTopLoci: aggregating across multiple entries requires ",
            "type = 'data.frame'."
        )
        abort(msg)
    }
    .fmrRowTopLoci(
        sel,
        type = "GRanges",
        signalCutoff = signalCutoff,
        minPurity = minPurity
    )
}

#' @rdname getMarginalEffects
#' @export
setMethod(
    "getMarginalEffects",
    "FineMappingResultBase",
    function(
        x,
        maxPval = NULL,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        # Selectors pinning one entry -> that entry's bare marginal table; no /
        # partial selectors -> aggregate every matching entry's marginals,
        # prefixed with the row identity (study/context/trait/blockId/method).
        .fmrAggregateView(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = .fmrRowMarginalEffects,
            maxPval = maxPval
        )
    }
)

#' @rdname getRegion
#' @examples
#' data(qtlFineMappingExample)
#' getRegion(qtlFineMappingExample)
#' @export
setMethod("getRegion", "FineMappingResultBase", function(x, ...) {
    .getRegionColumn(x)
})

#' @rdname getTraitPosition
#' @examples
#' data(qtlDatasetExample)
#' getTraitPosition(qtlDatasetExample)
#' @export
setMethod("getTraitPosition", "FineMappingResultBase", function(x, ...) {
    .getTraitPosColumn(x)
})

#' @rdname getSusieFit
#' @export
setMethod(
    "getSusieFit",
    "FineMappingResultBase",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        .fmrPartsSusieFit(.fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)

#' @rdname resolveWeights
#' @export
setMethod("resolveWeights", "FineMappingResultBase", function(x, ...) {
    # The per-variant weight of the row a selector pins. Defined on the
    # collection because that is what getFineMappingResult() now returns; the
    # body is the per-row primitive, so the two cannot drift.
    .fmrRowResolveWeights(.fmrSelectEntry(x, ...), ...)
})

#' @rdname getVariantIds
#' @export
setMethod(
    "getVariantIds",
    "FineMappingResultBase",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        .fmrPartsVariantIds(.fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "SumStatsBase", function(x, ...) {
    # One method covers GwasSumStats and QtlSumStats: each class's own
    # getSumStats() knows its selectors (study, or study/context/trait) and
    # raises the ambiguity error when a multi-row collection is addressed
    # without one, so `...` carries them through untouched.
    #
    # Rendered with the same .grVariantIds() the row classes use, so an id
    # means the same string whichever object produced it.
    .grVariantIds(getSumStats(x, ...))
})

#' @rdname subsetChr
#' @export
setMethod("subsetChr", "SumStatsBase", function(x, chr) {
    # The whole-seqname special case of subsetRegion(): one verb, one set of
    # semantics. Elements on other chromosomes are dropped rather than kept as
    # empties, which is what the seqname split makes natural anyway -- after
    # splitting, an element belongs to exactly one chromosome.
    chrName <- withChrPrefix(chr)
    subsetRegion(
        x,
        GenomicRanges::GRanges(
            chrName,
            IRanges::IRanges(1L, .Machine$integer.max)
        )
    )
})
