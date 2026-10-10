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
    "LdSketch_OR_NULL",
    c("RangedSummarizedExperiment", "NULL")
)

# What an LdData can read genotypes from. Four shapes, all real:
#
#   GenotypeHandle  a lazy file-backed panel; `snpIdx` selects into its whole
#                   snpInfo, so the handle must not be pre-narrowed
#   list            one handle per panel for a mixture reference, averaged by
#                   `mixtureWeights`
#   matrix          dosages already extracted and filtered, kept so
#                   genotypes() answers without reopening the file (see
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

# Mixing proportions, one per panel, when `genotypeHandle` is a list.
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
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
# z / nSamples / maf / nSnps are
# defined once on SumStatsBase (they only delegate to sumStats); subsetChr /
# varY / sumStats / as.data.frame stay on the concrete subclass because
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
#'   \code{length(qcInfo(x)) == 0L} -- the slot serves as both the gating
#'   flag and the audit trail.
#' @export
setClass(
    "SumStatsBase",
    contains = c("VIRTUAL", "RangedTupleList"),
    representation(
        ldSketch = "LdSketch_OR_NULL",
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
    build <- discard(unique(GenomeInfoDb::genome(object)), is.na)
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

#' @rdname qcInfo
#' @examples
#' data(qtlSumStatsExample)
#' qcInfo(qtlSumStatsExample)
#' @export
setMethod("qcInfo", "SumStatsBase", function(x) x@qcInfo)

#' @rdname qcDiagnostics
#' @examples
#' data(qtlSumStatsExample)
#' qcDiagnostics(qtlSumStatsExample)
#' @export
setMethod("qcDiagnostics", "SumStatsBase", function(x, entry = 1L) {
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

#' @rdname ldSketch
#' @examples
#' data(qtlSumStatsExample)
#' ldSketch(qtlSumStatsExample)
#' @export
setMethod("ldSketch", "SumStatsBase", function(x) x@ldSketch)

#' @rdname studyName
#' @examples
#' data(qtlDatasetExample)
#' studyName(qtlDatasetExample)
#' @export
setMethod("studyName", "SumStatsBase", function(x) {
    unique(as.character(x$study))
})

# z / nSamples / maf / nSnps delegate purely to sumStats (which the
# concrete subclass dispatches with its own tuple shape), so they live once on
# the base and forward `...` through.
#' @rdname z
#' @examples
#' data(qtlSumStatsExample)
#' z(qtlSumStatsExample)
#' @export
setMethod(
    "z",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mcols(exec(sumStats, x, !!!sel))$Z
    }
)

#' @rdname nSamples
#' @examples
#' data(qtlSumStatsExample)
#' nSamples(qtlSumStatsExample)
#' @export
setMethod(
    "nSamples",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mcols(exec(sumStats, x, !!!sel))$N
    }
)

# pval / coef / se are first-class alongside z / nSamples: they read the
# optional P / BETA / SE mcols and return NULL when the entry does not carry
# them (DataFrame `$` semantics), so a p-value-primary sumstats (e.g. TensorQTL
# cis output) is an equal citizen to a Z-primary GWAS sumstats.
#' @rdname pval
#' @examples
#' data(qtlSumStatsExample)
#' pval(qtlSumStatsExample)
#' @export
setMethod(
    "pval",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mcols(exec(sumStats, x, !!!sel))$P
    }
)

#' @rdname coef-methods
#' @examples
#' data(qtlSumStatsExample)
#' coef(qtlSumStatsExample)
#' @export
setMethod(
    "coef",
    "SumStatsBase",
    function(
        object,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL,
        ...
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mcols(exec(sumStats, object, !!!sel))$BETA
    }
)

#' @rdname se
#' @examples
#' data(qtlSumStatsExample)
#' se(qtlSumStatsExample)
#' @export
setMethod(
    "se",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mcols(exec(sumStats, x, !!!sel))$SE
    }
)

#' @rdname maf
#' @examples
#' data(qtlDatasetExample)
#' maf(qtlDatasetExample)
#' @export
setMethod(
    "maf",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mc <- mcols(exec(sumStats, x, !!!sel))
        if (is_in("MAF", colnames(mc))) mc$MAF else NULL
    }
)

#' @rdname af
#' @examples
#' data(qtlSumStatsExample)
#' af(qtlSumStatsExample)
#' @export
setMethod(
    "af",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        mc <- mcols(exec(sumStats, x, !!!sel))
        if (is_in("AF", colnames(mc))) mc$AF else NULL
    }
)

# =============================================================================
# Variant-level replacement
# =============================================================================

# The element index the selectors pin. GwasSumStats is keyed by study alone
# and QtlSumStats by the (study, context, trait) tuple; branching here rather
# than adding an internal generic keeps one copy of each replacement body.
# @noRd
.ssRowIndexFor <- function(x, studyName, context, trait) {
    if (methods::is(x, "QtlSumStats")) {
        return(.qtlSumStatsSelectRow(x, studyName, context, trait))
    }
    .gwasRefuseQtlSelectors(context, trait, NULL)
    .gwasSelectStudy(x, studyName)
}

# Write one per-variant column onto the pinned element and put it back. The
# `[[<-` route rebuilds the collection and re-runs validity, so a write that
# would break the identity invariants is refused rather than stored.
# @noRd
.ssWriteVariantColumn <- function(x, idx, column, values) {
    el <- x[[idx]]
    mcols(el)[[column]] <- values
    x[[idx]] <- el
    x
}

# Positional replacement. The vector is taken to be in the element's own
# variant order, so length is the only thing checkable -- and it is checked
# hard: a recycled or truncated frequency vector is silently wrong, which is
# exactly the failure the GRanges form exists to avoid.
# @noRd
.ssVariantColumnPositional <- function(x, idx, column, value) {
    n <- length(x[[idx]])
    if (length(value) != n) {
        abort(glue(
            "{column}: a positional replacement needs one value per variant ",
            "({n}); got {length(value)}. Pass a GRanges carrying A1/A2 to ",
            "match on (chrom, pos, A1, A2) instead."
        ))
    }
    .ssWriteVariantColumn(x, idx, column, as.numeric(value))
}

# Variant ids for either side of a match, rebuilt from the coordinates and
# alleles rather than read off a stored id column, so the two sides agree on
# format. NOTE formatVariantId() takes A2 BEFORE A1.
# @noRd
.ssVariantKeys <- function(gr) {
    mc <- mcols(gr)
    formatVariantId(
        chrom = as.character(seqnames(gr)),
        pos = start(gr),
        A2 = as.character(mc$A2),
        A1 = as.character(mc$A1)
    )
}

# The payload column of a replacement GRanges. A1/A2 are the match key, so
# exactly one other mcols column is the values; requiring that rather than a
# fixed name keeps the caller from having to know the stored column names.
# @noRd
.ssReplacementPayload <- function(value, column) {
    mc <- mcols(value)
    nms <- colnames(mc)
    if (!all(is_in(c("A1", "A2"), nms))) {
        abort(glue(
            "{column}: a GRanges replacement must carry A1 and A2 in ",
            "mcols() so an allele swap can be detected."
        ))
    }
    payload <- setdiff(nms, c("A1", "A2"))
    if (length(payload) != 1L) {
        abort(glue(
            "{column}: a GRanges replacement must carry exactly one value ",
            "column besides A1/A2 (got {length(payload)}",
            "{if (length(payload)) str_c(': ', str_flatten(payload, ', '))",
            " else ''})."
        ))
    }
    as.numeric(mc[[payload]])
}

# Matched replacement. Ids are rebuilt from both sides and paired with
# matchVariants, so a chr-prefix or separator difference does not read as
# no-overlap. A variant the replacement does not name becomes NA -- its value
# is unknown, not unchanged. `directional` says whether an allele swap
# inverts the quantity: AF goes to 1 - AF, while MAF (min(af, 1 - af)) and N
# are swap-invariant, and complementing those would corrupt them.
# @noRd
.ssVariantColumnMatched <- function(x, idx, column, value, directional) {
    vals <- .ssReplacementPayload(value, column)
    el <- x[[idx]]
    m <- matchVariants(.ssVariantKeys(el), .ssVariantKeys(value))
    if (length(m$idxA) == 0L) {
        abort(glue(
            "{column}: no variant in the replacement matches this entry on ",
            "(chrom, pos, A1, A2)."
        ))
    }
    picked <- vals[m$idxB]
    if (directional) {
        picked <- ifelse(m$sign < 0, 1 - picked, picked)
    }
    out <- rep(NA_real_, length(el))
    out[m$idxA] <- picked
    unused <- length(value) - length(m$idxB)
    if (unused > 0L) {
        warn(glue(
            "{column}: {unused} of {length(value)} replacement variant(s) ",
            "are not in this entry and were ignored."
        ))
    }
    missed <- length(el) - length(m$idxA)
    if (missed > 0L) {
        warn(glue(
            "{column}: {missed} of {length(el)} entry variant(s) are not in ",
            "the replacement and were set NA."
        ))
    }
    .ssWriteVariantColumn(x, idx, column, out)
}

# @noRd
.ssReplacePositional <- function(x, column, value, studyName, context, trait) {
    .ssVariantColumnPositional(
        x,
        .ssRowIndexFor(x, studyName, context, trait),
        column,
        value
    )
}

# @noRd
.ssReplaceMatched <- function(
    x,
    column,
    value,
    directional,
    studyName,
    context,
    trait
) {
    .ssVariantColumnMatched(
        x,
        .ssRowIndexFor(x, studyName, context, trait),
        column,
        value,
        directional
    )
}

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "af",
    signature(x = "SumStatsBase", value = "numeric"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplacePositional(x, "AF", value, studyName, context, trait)
    }
)

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "af",
    signature(x = "SumStatsBase", value = "GRanges"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplaceMatched(
            x,
            "AF",
            value,
            directional = TRUE,
            studyName,
            context,
            trait
        )
    }
)

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "maf",
    signature(x = "SumStatsBase", value = "numeric"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplacePositional(x, "MAF", value, studyName, context, trait)
    }
)

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "maf",
    signature(x = "SumStatsBase", value = "GRanges"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplaceMatched(
            x,
            "MAF",
            value,
            directional = FALSE,
            studyName,
            context,
            trait
        )
    }
)

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "nSamples",
    signature(x = "SumStatsBase", value = "numeric"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplacePositional(x, "N", value, studyName, context, trait)
    }
)

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "nSamples",
    signature(x = "SumStatsBase", value = "GRanges"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplaceMatched(
            x,
            "N",
            value,
            directional = FALSE,
            studyName,
            context,
            trait
        )
    }
)

#' @rdname nSnps
#' @examples
#' data(qtlSumStatsExample)
#' nSnps(qtlSumStatsExample)
#' @export
setMethod(
    "nSnps",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        length(exec(sumStats, x, !!!sel))
    }
)

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
    representation(ldSketch = "LdSketch_OR_NULL")
)

#' @rdname studyName
#' @export
setMethod("studyName", "FineMappingResultBase", function(x) {
    unique(as.character(x$study))
})

#' @rdname ldSketch
#' @export
setMethod("ldSketch", "FineMappingResultBase", function(x) x@ldSketch)

#' @rdname methodNames
#' @examples
#' data(qtlFineMappingExample)
#' methodNames(qtlFineMappingExample)
#' @export
setMethod("methodNames", "FineMappingResultBase", function(x) {
    unique(as.character(x$method))
})

#' @noRd
setMethod(
    "adjustPips",
    "FineMappingResultBase",
    function(x, keepVariants) {
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
            keepVariants = keepVariants
        )
        .fmrFromEntries(out, adjusted)
    }
)

# The collection's identity columns plus the per-entry payload columns.
# @noRd
.withEntryPayload <- function(md, entries) {
    # `[[<-` not cbind(): the collection may already carry these columns, and
    # they must be REPLACED -- cbind would append a second copy of each, which
    # every reader then shadows with the stale one.
    withFit <- `[[<-`(
        md,
        "susieFit",
        value = S4Vectors::SimpleList(map(entries, susieFit))
    )
    `[[<-`(
        withFit,
        "cvResult",
        value = S4Vectors::SimpleList(map(entries, cvResult))
    )
}

# Rebuild a fine-mapping collection from adjusted entries, keeping every
# identity column and collection-level slot.
# @noRd
.fmrFromEntries <- function(x, entries) {
    grl <- `mcols<-`(
        GenomicRanges::GRangesList(map(entries, variants)),
        value = .withEntryPayload(mcols(x, use.names = FALSE), entries)
    )
    new(class(x), grl, ldSketch = .asLdSketch(ldSketch(x)))
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
# approximation error. retainedMass() is what makes it visible instead of
# letting a heavily-trimmed fit look as confident as an untouched one.
# =============================================================================

#' @rdname intersectVariants
#' @export
setMethod(
    "intersectVariants",
    signature(x = "FineMappingResultBase", y = "FineMappingResultBase"),
    function(x, y, oneSided = FALSE) {
        shared <- .rcShared(x, y)
        if (isTRUE(oneSided)) {
            return(adjustPips(x, shared))
        }
        list(
            x = adjustPips(x, shared),
            y = adjustPips(y, shared)
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

#' @rdname retainedMass
#' @export
setMethod("retainedMass", "FineMappingResultBase", function(x) {
    if (nrow(x) == 0L) {
        return(.rcEmptyMass(x))
    }
    parts <- compact(map(seq_len(nrow(x)), .rcMassForRow, x = x))
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

#' @rdname credibleSets
#' @export
setMethod(
    "credibleSets",
    "FineMappingResultBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        coverage = 0.95,
        minPurity = NULL
    ) {
        # Selectors pinning one entry -> that entry's bare credible-set table;
        # no /
        # partial selectors -> aggregate every matching entry's credible sets,
        # prefixed with the row identity (study/context/trait/blockId/method).
        # `minPurity` is an independent CS-quality filter, orthogonal to
        # coverage.
        .fmrAggregateView(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = .fmrRowCs,
            viewArgs = list(
                coverage = coverage,
                minPurity = minPurity
            )
        )
    }
)

#' @rdname lbf
#' @export
setMethod("lbf", "FineMappingResultBase", function(x) {
    .fmrAggregateView(x, perEntry = .fmrRowLbf)
})

#' @rdname credibleSetSummary
#' @export
setMethod(
    "credibleSetSummary",
    "FineMappingResultBase",
    function(x, coverage = 0.95) {
        .fmrAggregateView(
            x,
            perEntry = .fmrRowCredibleSetSummary,
            viewArgs = list(coverage = coverage)
        )
    }
)

#' @rdname fsusieCredibleBand
#' @export
setMethod("fsusieCredibleBand", "FineMappingResultBase", function(x) {
    .fmrAggregateView(x, perEntry = .fmrRowFsusieCredibleBand)
})

#' @rdname fsusieAffectedRegions
#' @export
setMethod("fsusieAffectedRegions", "FineMappingResultBase", function(x) {
    perRow <- map(seq_len(nrow(x)), .fsusieEntryAffectedRegions, x = x)
    grs <- perRow[lengths(perRow) > 0L]
    if (length(grs) == 0L) {
        return(GenomicRanges::GRanges())
    }
    exec(c, !!!grs)
})

#' @rdname topLoci
#' @export
setMethod(
    "topLoci",
    "FineMappingResultBase",
    function(
        x,
        signalCutoff = 0.025,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        minPurity = NULL,
        raw = FALSE
    ) {
        # The five selectors travel together to both branches, so they are
        # named once and spliced rather than relisted twice.
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            method = method,
            region = region
        )
        # raw = TRUE hands back the stored canonical table verbatim, which is
        # a data.frame by definition -- there is no posterior projection to
        # put on a range.
        if (isTRUE(raw)) {
            return(exec(.fmrbTopLociTable, x, !!!sel, list(raw = TRUE)))
        }
        # The per-variant table first -- bare when the selectors pin one
        # entry, else the matching rows' tables stacked with row-identity
        # columns -- then one conversion to ranges. Converting last is what
        # lets the aggregate form be a GRanges too: the identity columns just
        # become mcols, where the old type = "GRanges" branch had to refuse
        # anything but a single pinned entry.
        tbl <- exec(
            .fmrbTopLociTable,
            x,
            !!!sel,
            list(signalCutoff = signalCutoff, minPurity = minPurity)
        )
        .fmeTopLociGRanges(tbl)
    }
)


# The per-variant top-loci table: the stored table verbatim (raw) or the
# projected posterior view, selected by `viewArgs`. Both stack with
# row-identity columns when the selectors match more than one entry. The
# caller converts to ranges; this stays a table so the aggregate form can be
# stacked before conversion.
# @noRd
.fmrbTopLociTable <- function(
    x,
    studyName,
    context,
    trait,
    method,
    region,
    viewArgs
) {
    .fmrAggregateView(
        x,
        studyName = studyName,
        context = context,
        trait = trait,
        method = method,
        region = region,
        perEntry = .fmrRowTopLoci,
        viewArgs = viewArgs
    )
}
# raw = TRUE hands back the stored canonical table verbatim, so the
# posterior-view projection does not apply -- and marginal effects are
# not a posterior view, so they stay a table rather than ranges.
#' @rdname marginalEffects
#' @export
setMethod(
    "marginalEffects",
    "FineMappingResultBase",
    function(
        x,
        maxPval = NULL,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL
    ) {
        # Selectors pinning one entry -> that entry's bare marginal table; no /
        # partial selectors -> aggregate every matching entry's marginals,
        # prefixed with the row identity (study/context/trait/blockId/method).
        .fmrAggregateView(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = .fmrRowMarginalEffects,
            viewArgs = list(maxPval = maxPval)
        )
    }
)

#' @rdname genomicRegion
#' @examples
#' data(qtlFineMappingExample)
#' genomicRegion(qtlFineMappingExample)
#' @export
setMethod("genomicRegion", "FineMappingResultBase", function(x) {
    .getRegionColumn(x)
})

#' @rdname traitPosition
#' @examples
#' data(qtlDatasetExample)
#' traitPosition(qtlDatasetExample)
#' @export
setMethod("traitPosition", "FineMappingResultBase", function(x) {
    .getTraitPosColumn(x)
})

#' @rdname susieFit
#' @export
setMethod(
    "susieFit",
    "FineMappingResultBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL
    ) {
        .fmrPartsSusieFit(.fmrSelectEntry(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)

#' @rdname resolveWeights
#' @export
setMethod(
    "resolveWeights",
    "FineMappingResultBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL
    ) {
        # The per-variant weight of the row a selector pins. Defined on the
        # collection because that is what fineMappingResult() now
        # returns; the body is the per-row primitive, so they cannot drift.
        #
        # The selectors are named, not absorbed by `...`: .fmrSelectEntry is the
        # only consumer, and each concrete class refuses the ones it does not
        # index by (a QtlFineMappingResult has no `region` axis).
        .fmrRowResolveWeights(.fmrSelectEntry(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)

#' @rdname variantIds
#' @export
setMethod(
    "variantIds",
    "FineMappingResultBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL
    ) {
        .fmrPartsVariantIds(.fmrSelectEntry(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)

#' @rdname variantIds
#' @export
setMethod(
    "variantIds",
    "SumStatsBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL
    ) {
        # One method covers GwasSumStats and QtlSumStats. The selectors are the
        # union: GwasSumStats reads only `study` / `ranges` and refuses the
        # QtlSumStats-only ones rather than accepting and ignoring them.
        #
        # Rendered with the same .grVariantIds() the row classes use, so an id
        # means the same string whichever object produced it.
        sel <- list(
            study = studyName,
            context = context,
            trait = trait,
            annotateSignificance = annotateSignificance,
            ranges = ranges
        )
        .grVariantIds(exec(sumStats, x, !!!sel))
    }
)

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
