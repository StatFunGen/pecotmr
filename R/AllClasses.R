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

#' @include AllGenerics.R GenotypeHandle.R RangedTupleList.R LdMatrix.R
#' @include GenotypeMixture.R
#' @importFrom methods setClass setMethod new is validObject
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
NULL

# The LD reference a collection was harmonized against -- the one place every
# pipeline reads LD from. It is one of two things:
#
#   RangedSummarizedExperiment  a genotype panel, the same shape QtlDataset
#                               uses for its genotypes: variants in rowRanges,
#                               dosages in a DelayedArray assay that reads
#                               lazily through a GenotypeHandle
#   GenotypeMixture             several such panels, their LD averaged
#   LdMatrix                    a correlation matrix the caller built
#
# NULL when a collection carries no LD reference.
#
# A union rather than "ANY": the slot has always been documented as holding
# one of these, and an untyped slot with a typed docstring is the kind of
# thing that is only true until someone puts something else in it.
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setClassUnion(
    "LdReference_OR_NULL",
    c("RangedSummarizedExperiment", "GenotypeMixture", "LdMatrix", "NULL")
)

# =============================================================================
# SumStatsBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS summary statistics collections.
# Concrete subclasses (QtlSumStats, GwasSumStats) inherit from
# RangedTupleList and share the ldReference / qcOptions slots (the genome build
# lives in seqinfo, not a slot). Each element
# is one tuple's per-variant GRanges: x[[i]], formerly x$entry[[i]].
#
# z / nSamples / maf / nSnps are
# defined once on SumStatsBase (they only delegate to sumStats);
# varY / sumStats / as.data.frame stay on the concrete subclass because
# they rely on the tuple shape (3-tuple QtlSumStats, 1-tuple GwasSumStats).
# =============================================================================

#' @title Summary Statistics Base Class
#' @description Virtual base class for QTL and GWAS summary statistics
#'   collections. Concrete subclasses (\code{QtlSumStats}, \code{GwasSumStats})
#'   inherit from \code{\linkS4class{RangedTupleList}} and share the
#'   \code{ldReference} / \code{qcOptions} slots, and the genome build in
#'   \code{seqinfo()}.
#'
#'   Each element is the per-variant \code{GRanges} of one tuple, so
#'   \code{x[[i]]} is that tuple's summary statistics and the identity columns
#'   live in \code{mcols(x)}. There is no \code{entry} column: what used to be
#'   \code{x$entry[[i]]} is now simply \code{x[[i]]}.
#' @slot ldReference A named list of the LD references the rows were
#'   harmonized against and every pipeline reads LD from, each a genotype
#'   panel (see \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}}
#'   or an \code{\link{LdMatrix}}. Each distinct reference is stored once, and
#'   the \code{ldReference} column of \code{mcols} names the entry each row
#'   uses, so rows can sit on different references. Optional: an LD-free row
#'   (e.g. for mash, which operates across conditions per variant) has
#'   \code{NA} there; pipelines that need LD validate its presence when they
#'   consume the collection.
#' @slot qcOptions A \code{list} recording the QC that ran over the
#'   collection: the options \code{summaryStatsQc()} was called with and when.
#'   Empty \code{list()} on construction. Pipelines that need QC'd summary
#'   statistics refuse inputs where \code{length(qcOptions(x)) == 0L}. What
#'   QC did to each row is in the \code{qcAudit} column of \code{mcols()},
#'   read as a table with \code{\link{qcInfo}}.
#' @export
setClass(
    "SumStatsBase",
    contains = c("VIRTUAL", "RangedTupleList"),
    representation(
        ldReference = "list",
        qcOptions = "list"
    ),
    prototype(ldReference = list())
)

methods::setValidity("SumStatsBase", function(object) {
    .ldRefDictProblems(object@ldReference, .ldRefKeys(object)) %||% TRUE
})

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
    rangeKey <- .keyEmptyElements(object, .rtlRangeKeys(object))
    keys <- exec(str_c, !!!c(tupleKey, list(rangeKey)), sep = "|")
    anyDuplicated(keys) > 0L
}

# Span keys with every EMPTY element given a key of its own. An empty element
# has no span to key on -- it is a region QC screened out, or one a filter()
# emptied -- and several of them in one study are the normal outcome of a PIP
# screen over a block-split study, not duplicates. Its block names it when
# the collection carries `blockId`; otherwise its position does, so two empty
# elements only collide when they claim the same block.
# @noRd
.keyEmptyElements <- function(object, spanKeys) {
    empty <- lengths(object) == 0L
    if (!any(empty)) {
        return(spanKeys)
    }
    md <- mcols(object, use.names = FALSE)
    block <- if (is_in("blockId", colnames(md))) {
        as.character(md$blockId)
    } else {
        rep(NA_character_, length(object))
    }
    emptyKey <- if_else(
        is.na(block),
        str_c("\u0001row:", seq_along(spanKeys)),
        str_c("\u0001block:", block)
    )
    if_else(empty, emptyKey, spanKeys)
}

# @noRd
.ssKeyColumn <- function(cn, md) {
    as.character(md[[cn]])
}

# Stitch a tuple's elements back into one GRanges, optionally restricted to a
# window. The seqname split means one tuple can own several elements, but the
# accessor contract is still "one GRanges per tuple". `window` is for the
# pipelines, which score a pair over the span its weights reach; the public
# accessors return the whole tuple and leave narrowing to subsetByOverlaps().
# @noRd
.ssStitchElements <- function(x, idx, window = NULL) {
    gr <- .rtlGatherElements(x, idx)
    if (is.null(window)) {
        return(gr)
    }
    win <- .assertRegionGRanges(window, "`window`")
    onWindowChrom <- is_in(
        as.character(seqnames(gr)),
        as.character(seqnames(win))
    )
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

# The genome build as one string for a show() line. genome() answers per
# SEQLEVEL, so interpolating it directly prints the line once per chromosome
# of a multi-chromosome collection, and not at all for one spanning none.
# @noRd
.showGenomeBuild <- function(x) {
    build <- unique(discard(GenomeInfoDb::genome(x), is.na))
    if (length(build) == 0L) {
        return("NA")
    }
    str_flatten(build, ", ")
}

# ---- The genome build a result carries --------------------------------------

# The builds `x` names: a ranged object's seqinfo genome, a QtlDataset's
# experiments', a MultiStudyQtlDataset's studies' and summary statistics',
# a plain list's elements'. An unset (NA) build names nothing, and neither
# does an object without ranges (NULL, an unranged SummarizedExperiment).
# @noRd
.namedBuilds <- function(x) {
    builds <- if (is.list(x) && !is.object(x)) {
        list_c(map(x, .namedBuilds), ptype = character())
    } else if (is(x, "MultiStudyQtlDataset")) {
        c(
            list_c(map(as.list(x), .namedBuilds), ptype = character()),
            .namedBuilds(sumStats(x))
        )
    } else if (is(x, "MultiAssayExperiment")) {
        list_c(
            map(as.list(MultiAssayExperiment::experiments(x)), .namedBuilds),
            ptype = character()
        )
    } else if (.isRangedInput(x)) {
        GenomeInfoDb::genome(x)
    } else {
        character(0)
    }
    unique(discard(unname(as.character(builds)), is.na))
}

# @noRd
.isRangedInput <- function(x) {
    is(x, "GenomicRanges") ||
        is(x, "GRangesList") ||
        is(x, "RangedSummarizedExperiment")
}

# The one build `x` names, or NA when it names none -- or several, which a
# collection assembled from mismatched parts can.
# @noRd
.oneBuild <- function(x) {
    builds <- .namedBuilds(x)
    if (length(builds) == 1L) builds else NA_character_
}

# The build a result computed from `inputs` (a list) carries: the one they
# name, or NA when none names one. Inputs naming different builds are
# refused, as GenomeInfoDb refuses to merge ranges whose genomes differ:
# positions on one build do not place variants on another.
# @noRd
.inputsBuild <- function(inputs, caller) {
    builds <- unique(list_c(map(inputs, .namedBuilds), ptype = character()))
    if (length(builds) > 1L) {
        abort(c(
            glue(
                "{caller}: the inputs name different genome builds ",
                "({str_flatten(builds, ', ')})."
            ),
            i = str_c(
                "Positions on one build do not place variants on another: ",
                "lift the inputs onto one build, or, when the names only ",
                "spell the same build differently, give them one name with ",
                "genome(x) <- \"<build>\"."
            )
        ))
    }
    if (length(builds) == 0L) NA_character_ else builds
}

# `x` with its seqinfo naming `build` for every sequence, or unchanged when
# the build is unknown (NA).
# @noRd
.withGenome <- function(x, build) {
    if (is.null(x) || is.na(build)) {
        return(x)
    }
    GenomeInfoDb::`genome<-`(x, value = build)
}

#' @rdname qcInfo
#' @examples
#' data(qtlSumStatsExample)
#' qcInfo(qtlSumStatsExample)
#' @export
setMethod("qcInfo", "SumStatsBase", function(x) {
    .asViewDataFrame(.ssQcInfoTable(x))
})

#' @rdname qcOptions
#' @examples
#' data(qtlSumStatsExample)
#' qcOptions(qtlSumStatsExample)
#' @export
setMethod("qcOptions", "SumStatsBase", function(x) x@qcOptions)

# The per-row QC audit as a table: one row per tuple, its identity columns,
# then each scalar count the audits record, a nested count named by its path
# (contentFilters.nDropped). A row a step did not run on is NA there; the
# per-variant tables an audit can also hold are read with qcDiagnostics().
# @noRd
.ssQcInfoTable <- function(x) {
    ids <- intersect(
        c("study", "context", "trait", "blockId"),
        .tupleColumnNames(x)
    )
    base <- as_tibble(set_names(map(ids, .ssIdentityColumn, x = x), ids))
    audits <- mcols(x, use.names = FALSE)$qcAudit
    if (is.null(audits)) {
        return(base)
    }
    bind_cols(base, .qcAuditRows(as.list(audits)))
}

# @noRd
.ssIdentityColumn <- function(cn, x) {
    as.character(.tupleColumn(x, cn))
}

# One column per scalar count any audit records, NA for the rows without it.
# @noRd
.qcAuditRows <- function(audits) {
    flat <- map(audits, .qcAuditScalars)
    cols <- unique(list_c(map(flat, names), ptype = character()))
    as_tibble(
        set_names(map(cols, .qcAuditColumn, flat = flat), cols),
        .rows = length(audits)
    )
}

# @noRd
.qcAuditColumn <- function(col, flat) {
    list_c(map(flat, col, .default = NA))
}

# One audit's scalar counts as a flat named list, a nested count named by its
# path. Tables and other non-scalar values are left out.
# @noRd
.qcAuditScalars <- function(a, prefix = "") {
    if (!is.list(a) || is.data.frame(a) || length(a) == 0L) {
        return(list())
    }
    list_flatten(unname(map2(a, str_c(prefix, names(a)), .qcAuditLeaf)))
}

# @noRd
.qcAuditLeaf <- function(v, name) {
    if (is.list(v) && !is.data.frame(v)) {
        return(.qcAuditScalars(v, str_c(name, ".")))
    }
    if (is.atomic(v) && length(v) == 1L) {
        return(set_names(list(v), name))
    }
    list()
}

#' @rdname qcDiagnostics
#' @examples
#' data(qtlSumStatsExample)
#' qcDiagnostics(qtlSumStatsExample)
#' @export
setMethod(
    "qcDiagnostics",
    "SumStatsBase",
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        audits <- mcols(x, use.names = FALSE)$qcAudit
        if (is.null(audits)) {
            return(NULL)
        }
        idx <- .ssRowIndexFor(x, studyName, context, trait)
        tables <- compact(map(audits[idx], "ldMismatchDiagnostics"))
        if (length(tables) == 0L) {
            return(NULL)
        }
        # A GWAS study split across rows checked each row on its own, so its
        # tables stack, one row per variant, as sumStats() stitches them.
        DataFrame(bind_rows(tables), check.names = FALSE)
    }
)

# Row `i`'s QC audit, or NULL when QC recorded none for it.
# @noRd
.rowQcAudit <- function(x, i) {
    audits <- mcols(x, use.names = FALSE)$qcAudit
    if (is.null(audits) || i > length(audits)) {
        return(NULL)
    }
    audits[[i]]
}

# Before each row's QC audit moved into mcols(), a summary-statistics
# collection held its whole QC record in a `qcInfo` slot: the options, and
# under `entryAudit` one audit per row. An object saved then still carries
# that slot as an attribute; it is split into the `qcOptions` slot and the
# `qcAudit` column. An object older still held its genome build in a
# `genome` slot (.ssMigrateGenome()). An object that old needs the updates
# every collection shares too, so the inherited update runs either way.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "SumStatsBase",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message(str_c("updateObject(object = '", class(object)[[1L]], "')"))
        }
        old <- attr(object, "qcInfo", exact = TRUE)
        if (!is.null(old)) {
            object <- .ssMigrateQcInfo(object, old)
        }
        methods::callNextMethod(
            .ssMigrateGenome(object),
            ...,
            verbose = verbose
        )
    }
)

# Before the genome build moved into seqinfo(), summary statistics held it
# in a `genome` slot; an object saved then still carries it as an
# attribute, and it becomes the build seqinfo() names (none, when the slot
# was left empty).
# @noRd
.ssMigrateGenome <- function(object) {
    build <- attr(object, "genome", exact = TRUE)
    if (is.null(build)) {
        return(object)
    }
    attr(object, "genome") <- NULL
    .withGenome(object, build[1L])
}

# @noRd
.ssMigrateQcInfo <- function(object, old) {
    attr(object, "qcInfo") <- NULL
    attr(object, "qcOptions") <- old[setdiff(names(old), "entryAudit")]
    audits <- old$entryAudit
    if (length(audits) == 0L) {
        return(object)
    }
    if (length(audits) != length(object)) {
        warn(glue(
            "updateObject(): the saved QC record holds {length(audits)} row ",
            "audit(s) for {length(object)} row(s), so they cannot be matched ",
            "to rows and are dropped; qcOptions() keeps the settings."
        ))
        return(object)
    }
    `$<-`(object, "qcAudit", S4Vectors::SimpleList(audits))
}

#' @rdname ldReference
#' @examples
#' data(qtlSumStatsExample)
#' # One entry per distinct reference; each row names its own.
#' ldReference(qtlSumStatsExample)
#' qtlSumStatsExample$ldReference
#' @export
setMethod("ldReference", "SumStatsBase", function(x) {
    .ldRefsInUse(x)
})

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
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        mcols(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))$Z
    }
)

#' @rdname sampleSize
#' @examples
#' data(qtlSumStatsExample)
#' head(sampleSize(qtlSumStatsExample))
#' @export
setMethod(
    "sampleSize",
    "SumStatsBase",
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        mcols(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))$N
    }
)

# pval / coef / se are first-class alongside z / sampleSize: they read the
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
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        mcols(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))$P
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
    function(object, studyName = NULL, context = NULL, trait = NULL, ...) {
        mcols(sumStats(
            object,
            studyName = studyName,
            context = context,
            trait = trait
        ))$BETA
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
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        mcols(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))$SE
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
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        mc <- mcols(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))
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
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        mc <- mcols(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))
        if (is_in("AF", colnames(mc))) mc$AF else NULL
    }
)

# =============================================================================
# Variant-level replacement
# =============================================================================

# The element indices the selectors pin: every row of the study for a
# GwasSumStats, which is keyed by study alone and split by chromosome or LD
# block, and the one (study, context, trait) row for a QtlSumStats. Branching
# here rather than adding an internal generic keeps one copy of each
# replacement body.
# @noRd
.ssRowIndexFor <- function(x, studyName, context, trait) {
    if (methods::is(x, "QtlSumStats")) {
        return(.qtlSumStatsSelectRow(x, studyName, context, trait))
    }
    .gwasRefuseQtlSelectors(context, trait, NULL)
    .gwasSelectStudy(x, studyName)
}

# Write one per-variant column onto the pinned elements and put them back.
# `idx` can name several elements -- a GWAS study split by chromosome or LD
# block -- and `values` runs over them in `idx` order, the order sumStats()
# stitches them in. The rebuild re-runs validity, so a write that would break
# the identity invariants is refused rather than stored.
# @noRd
.ssWriteVariantColumn <- function(x, idx, column, values) {
    pieces <- split(
        values,
        factor(rep(seq_along(idx), lengths(x)[idx]), levels = seq_along(idx))
    )
    elements <- as.list(x)
    elements[idx] <- map2(
        elements[idx],
        pieces,
        .ssSetVariantColumn,
        column = column
    )
    .rtlRebuild(x, elements, seq_along(x))
}

# @noRd
.ssSetVariantColumn <- function(el, values, column) {
    mcols(el)[[column]] <- unname(values)
    el
}

# Positional replacement. The vector is taken to be in sumStats()'s variant
# order for the selection, so length is the only thing checkable -- and it is
# checked hard: a recycled or truncated frequency vector is silently wrong,
# which is exactly the failure fillFromReference() exists to avoid.
# @noRd
.ssVariantColumnPositional <- function(x, idx, column, value) {
    n <- sum(lengths(x)[idx])
    if (length(value) != n) {
        abort(glue(
            "{column}: a positional replacement needs one value per variant ",
            "({n}, in sumStats() order); got {length(value)}. Use ",
            "fillFromReference() to match on (chrom, pos, A1, A2) instead."
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

# The values column of a fillFromReference() reference. A1/A2 are the match
# key; the values are the column named like the target when the reference
# has one (a panel GRanges carries several), else its one other column.
# @noRd
.ssReplacementPayload <- function(value, column) {
    mc <- mcols(value)
    nms <- colnames(mc)
    if (!all(is_in(c("A1", "A2"), nms))) {
        abort(glue(
            "fillFromReference: `reference` must carry A1 and A2 in ",
            "mcols() so an allele swap can be detected."
        ))
    }
    if (is_in(column, nms)) {
        return(as.numeric(mc[[column]]))
    }
    payload <- setdiff(nms, c("A1", "A2"))
    if (length(payload) != 1L) {
        abort(glue(
            "fillFromReference: `reference` has no `{column}` column, so ",
            "its values must be the one column besides A1/A2 (got ",
            "{length(payload)}",
            "{if (length(payload)) str_c(': ', str_flatten(payload, ', '))",
            " else ''})."
        ))
    }
    as.numeric(mc[[payload]])
}

# Matched fill. Ids are rebuilt from both sides and paired with matchVariants,
# so a chr-prefix or separator difference does not read as no-overlap. A
# variant the reference does not name becomes NA -- its value is unknown, not
# unchanged. `directional` says whether an allele swap inverts the quantity:
# AF goes to 1 - AF, while MAF (min(af, 1 - af)) and N are swap-invariant,
# and complementing those would corrupt them.
# @noRd
.ssVariantColumnMatched <- function(x, idx, column, value, directional) {
    vals <- .ssReplacementPayload(value, column)
    # The whole selection at once, so a variant the reference names is
    # "unused" only when it is in none of the study's rows.
    el <- .ssStitchElements(x, idx)
    m <- matchVariants(.ssVariantKeys(el), .ssVariantKeys(value))
    if (length(m$idxA) == 0L) {
        abort(glue(
            "fillFromReference: no variant in `reference` matches the ",
            "selection on (chrom, pos, A1, A2)."
        ))
    }
    picked <- vals[m$idxB]
    if (directional) {
        picked <- ifelse(m$sign < 0, 1 - picked, picked)
    }
    out <- rep(NA_real_, length(el))
    out[m$idxA] <- picked
    .ssWarnUnmatched(
        column,
        length(value) - length(m$idxB),
        length(value),
        length(el) - length(m$idxA),
        length(el)
    )
    .ssWriteVariantColumn(x, idx, column, out)
}

# The two ways a matched fill falls short, each said with its count: reference
# variants the selection does not hold, and selected variants left NA.
# @noRd
.ssWarnUnmatched <- function(column, unused, nReference, missed, nSelected) {
    if (unused > 0L) {
        warn(glue(
            "fillFromReference({column}): {unused} of {nReference} reference ",
            "variant(s) are not in the selection and were ignored."
        ))
    }
    if (missed > 0L) {
        warn(glue(
            "fillFromReference({column}): {missed} of {nSelected} selected ",
            "variant(s) are not in `reference` and were set NA."
        ))
    }
    invisible(NULL)
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

#' @title Fill A Per-Variant Column From A Reference
#' @description Fill the \code{AF}, \code{MAF} or \code{N} column of a
#'   \code{GwasSumStats} / \code{QtlSumStats} selection from a reference
#'   \code{GRanges}, pairing variants on \code{(chrom, pos, A1, A2)} -- the
#'   usual case being effect-allele frequencies the source file did not
#'   carry, taken from a reference panel. Order does not matter.
#'
#'   An allele swap sends \code{AF} to \code{1 - AF}, so a flipped match is
#'   complemented; \code{MAF} is \code{min(af, 1 - af)} and \code{N} is a
#'   count, both swap-invariant, so neither is. A selected variant the
#'   reference does not name is set \code{NA} -- its value is unknown, not
#'   unchanged -- with a warning, and so is a reference variant the selection
#'   does not hold; no overlap at all is an error.
#'
#'   Values already in \code{\link{sumStats}} order are assigned positionally
#'   instead, with \code{af(x) <- value} (see
#'   \code{\link{variantColumnReplace}}).
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param reference A \code{GRanges} carrying \code{A1} and \code{A2} in
#'   \code{mcols()} plus the values: a column named like \code{column}, or
#'   else its one column besides \code{A1}/\code{A2}.
#' @param column Which column to fill: \code{"AF"} (default), \code{"MAF"} or
#'   \code{"N"}.
#' @param studyName,context,trait Selectors pinning one study (every row of a
#'   GWAS study) or one \code{QtlSumStats} entry; \code{context} and
#'   \code{trait} are \code{QtlSumStats}-only.
#' @return \code{x}, with the selection's column filled.
#' @examples
#' data(qtlSumStatsExample)
#' ref <- sumStats(qtlSumStatsExample)
#' S4Vectors::mcols(ref) <- S4Vectors::DataFrame(
#'     A1 = ref$A1,
#'     A2 = ref$A2,
#'     AF = runif(length(ref))
#' )
#' x <- fillFromReference(qtlSumStatsExample, ref, "AF")
#' head(af(x))
#' @importFrom checkmate assertClass
#' @export
fillFromReference <- function(
    x,
    reference,
    column = c("AF", "MAF", "N"),
    studyName = NULL,
    context = NULL,
    trait = NULL
) {
    assertClass(x, "SumStatsBase")
    assertClass(reference, "GRanges")
    column <- arg_match(column)
    .ssVariantColumnMatched(
        x,
        .ssRowIndexFor(x, studyName, context, trait),
        column,
        reference,
        directional = identical(column, "AF")
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
    "maf",
    signature(x = "SumStatsBase", value = "numeric"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplacePositional(x, "MAF", value, studyName, context, trait)
    }
)

#' @rdname variantColumnReplace
#' @export
setReplaceMethod(
    "sampleSize",
    signature(x = "SumStatsBase", value = "numeric"),
    function(x, studyName = NULL, context = NULL, trait = NULL, ..., value) {
        .ssReplacePositional(x, "N", value, studyName, context, trait)
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
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        length(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))
    }
)

# =============================================================================
# FineMappingResultBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS fine-mapping result collections.
# Concrete subclasses (QtlFineMappingResult, GwasFineMappingResult) carry
# a DFrame of per-fit rows plus a shared ldReference slot. Downstream
# pipelines dispatch on FineMappingResultBase for behaviors that apply to
# either flavour, and on the concrete subclass when the tuple shape
# matters.
# =============================================================================

#' @title Fine-Mapping Result Base Class
#' @description Virtual base class for fine-mapping result collections. Concrete
#'   subclasses (\code{QtlFineMappingResult}, \code{GwasFineMappingResult})
#'   carry a \code{DFrame} of per-fit rows and an \code{ldReference} slot
#'   listing the LD references those fits used.
#'   Downstream pipelines should dispatch on \code{FineMappingResultBase} for
#'   behaviors that apply to either flavour, and on the concrete subclass when
#'   the tuple shape matters.
#' @slot ldReference A named list of the LD references the fits were
#'   computed against, each a genotype panel, a \code{GenotypeMixture} or an
#'   \code{LdMatrix}, stored once apiece; the \code{ldReference} column of
#'   \code{mcols} names the entry each row used (\code{NA} for a fit derived
#'   from individual-level data). Used downstream to check that the rows a
#'   pipeline combines share an LD reference.
#' @section Narrowed rows: A row is built from its fit, holding the fit's
#'   variables in the fit's order, so its fit is read off it by position. A
#'   row narrowed since -- by \code{endoapply()} or \code{[[<-} -- keeps its
#'   fit as it was, and is found in the fit by the names the fit gives its
#'   variables (a narrowed row whose fit names none, such as a trimmed fit,
#'   cannot be read). Every reader of the fit warns, with condition class
#'   \code{pecotmrFitCoverage}, when a row holds fewer variants than its fit
#'   describes. A row replaced by one of the same length is read by
#'   position. The readers are the accessors
#'   (\code{susieFit()}, \code{lbf()}, \code{cvResult()},
#'   \code{credibleSetSummary()} and the fSuSiE views),
#'   \code{computeCsCorrelation()}, \code{susieResult()},
#'   \code{extractCsInfo()}, \code{extractTopPipInfo()},
#'   \code{writeGwasVcf()}, and the pipelines that score the fits:
#'   \code{colocPipeline(adjustPips = FALSE)}, \code{qtlEnrichmentPipeline()}
#'   and \code{twasWeightsPipeline()}. Values held per variant, such as
#'   \code{topLoci()}, are narrowed with the row and do not warn. Nothing is
#'   rescaled: \code{\link{intersectVariants}} restricts a fit explicitly, as
#'   \code{colocPipeline()} does by default.
#' @export
setClass(
    "FineMappingResultBase",
    contains = c("VIRTUAL", "RangedTupleList"),
    representation(ldReference = "list"),
    prototype(ldReference = list())
)

methods::setValidity("FineMappingResultBase", function(object) {
    .ldRefDictProblems(object@ldReference, .ldRefKeys(object)) %||% TRUE
})

#' @rdname studyName
#' @export
setMethod("studyName", "FineMappingResultBase", function(x) {
    unique(as.character(x$study))
})

#' @rdname ldReference
#' @export
setMethod("ldReference", "FineMappingResultBase", function(x) {
    .ldRefsInUse(x)
})

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
        # renormalize, so they are DROPPED, leaving no entry with zero
        # variants. Every other failure (a
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
    withCv <- `[[<-`(
        withFit,
        "cvResult",
        value = S4Vectors::SimpleList(map(entries, cvResult))
    )
    `[[<-`(
        withCv,
        "derivedColumns",
        value = IRanges::CharacterList(map(entries, .fmRowDerivedColumns))
    )
}

# Rebuild a fine-mapping collection from adjusted entries, keeping every
# identity column and collection-level slot.
# @noRd
.fmrFromEntries <- function(x, entries) {
    # .rtlRebuild() carries the names, seqinfo (the genome build) and every
    # slot -- the rows keep their keys in mcols, so the references go across
    # as is -- and only the payload columns are new.
    rebuilt <- .rtlRebuild(x, map(entries, variants), seq_len(length(x)))
    withPayload <- `mcols<-`(
        rebuilt,
        value = .withEntryPayload(mcols(rebuilt, use.names = FALSE), entries)
    )
    S4Vectors::`metadata<-`(withPayload, value = S4Vectors::metadata(x))
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

# ---- Trimming a collection to what its role reads ---------------------------

#' @title Trim A Fine-Mapping Collection To What Its Role Reads
#' @description Keep, of each row's fit, only what a fine-mapping collection
#'   is read for in its \code{role}, so a genome-wide collection fits in
#'   memory. Nothing is sampled or summarized ahead of time: a trimmed
#'   collection gives the same enrichment estimates as the untrimmed one, and
#'   can still be resampled.
#'
#'   \describe{
#'     \item{\code{"outcome"}}{The side whose PIPs are scanned, as in
#'       \code{\link{qtlEnrichmentPipeline}}'s \code{gwasFineMappingResult}:
#'       each fit keeps its \code{pip} and its credible sets (\code{sets},
#'       \code{sets_secondary}).}
#'     \item{\code{"annotation"}}{The side whose single effects are sampled,
#'       as in \code{qtlFineMappingResult}: each fit also keeps \code{alpha}
#'       and the prior variances \code{V} of its effects of positive prior
#'       variance -- the only ones the enrichment samples -- with rows and
#'       entries named by effect id (\code{"L<k>"}). A fit with no such
#'       effect keeps a zero-row \code{alpha}.}
#'   }
#'
#'   Each row's per-variant columns are cut to the identity columns, the
#'   PIPs and the credible-set columns, and those two are stored only where
#'   the trimmed fit cannot give them back. The posterior summaries, marginal
#'   statistics and per-fit metadata columns go, so \code{\link{topLoci}} of a
#'   trimmed collection reports PIPs and credible sets only. The
#'   cross-validation results go as well, so a trimmed collection cannot feed
#'   \code{\link{twasWeightsPipeline}}. Identity columns, names, the genome
#'   build, LD references and \code{metadata()} are kept. Trimming a trimmed
#'   collection changes nothing.
#' @param x A \code{\link{GwasFineMappingResult}} or
#'   \code{\link{QtlFineMappingResult}}.
#' @param role \code{"outcome"} or \code{"annotation"}. \code{NULL} (default)
#'   picks by class: \code{"outcome"} for a \code{GwasFineMappingResult},
#'   \code{"annotation"} for a \code{QtlFineMappingResult}.
#' @return \code{x}, trimmed.
#' @examples
#' data(qtlFineMappingExample)
#' trimmed <- trimFineMappingResult(qtlFineMappingExample)
#' names(susieFit(trimmed))
#' c(
#'     before = utils::object.size(qtlFineMappingExample),
#'     after = utils::object.size(trimmed)
#' )
#' @importFrom checkmate assertMultiClass
#' @importFrom rlang arg_match0
#' @export
trimFineMappingResult <- function(x, role = NULL) {
    assertMultiClass(x, c("GwasFineMappingResult", "QtlFineMappingResult"))
    role <- role %||%
        if (methods::is(x, "GwasFineMappingResult")) "outcome" else "annotation"
    role <- arg_match0(role, c("outcome", "annotation"))
    entries <- map(seq_len(nrow(x)), .trimRow, x = x, role = role)
    .fmrFromEntries(x, entries)
}

# Row i of `x`, trimmed for `role`: its fit trimmed, its cross-validation
# results dropped, and its table cut to the PIPs and credible-set columns,
# which the row builder stores only where the trimmed fit cannot give them
# back. Each variant id is held once, by the row's ranges and alleles: the
# trimmed fit's vectors go unnamed, read by position, except for a row
# narrowed since its fit, which is found in the fit by name, so its fit keeps
# the names of its variables (.partsFitIds()).
# @noRd
.trimRow <- function(i, x, role) {
    parts <- .fmrRowParts(x, i)
    gr <- variants(parts)
    if (length(gr) == 0L) {
        return(parts)
    }
    tl <- .fmrPartsTopLoci(parts)
    csCols <- names(tl)[str_detect(names(tl), "^cs_[0-9]+(_purity)?$")]
    fit <- susieFit(parts)
    fitIds <- if (.rowHoldsFit(length(gr), .fitWidth(fit))) {
        NULL
    } else {
        .partsFitIds(parts)
    }
    .fmRowBuild(
        `mcols<-`(
            GenomicRanges::granges(gr, use.mcols = FALSE),
            value = mcols(gr, use.names = FALSE)[, c("A1", "A2"), drop = FALSE]
        ),
        .trimFit(fit, role, fitIds),
        tl[intersect(c("variant_id", "pip", csCols), names(tl))],
        NULL
    )
}

# A fit cut to what `role` reads (see trimFineMappingResult()), its PIPs and
# effect columns named by `fitIds` (a narrowed row's) or unnamed. A fit that
# is not a list -- none at all -- is left as it is.
# @noRd
.trimFit <- function(fit, role, fitIds = NULL) {
    if (!is.list(fit) || length(fit) == 0L) {
        return(fit)
    }
    kept <- compact(list(
        pip = .trimNamed(fit[["pip"]], fitIds),
        sets = fit[["sets"]],
        sets_secondary = .trimSecondarySets(fit[["sets_secondary"]])
    ))
    if (role == "annotation") {
        kept <- c(.trimEffects(fit, fitIds), kept)
    }
    `class<-`(kept, class(fit))
}

# `pip` named by `fitIds`, or unnamed when there are none to name it by.
# @noRd
.trimNamed <- function(pip, fitIds) {
    if (is.null(pip)) {
        return(NULL)
    }
    `names<-`(pip, if (length(fitIds) == length(pip)) fitIds)
}

# The credible sets of each other coverage, without the per-coverage copies
# of the PIPs and the rest that a fit's sets_secondary may carry: only the
# sets are read (.csSetsForCoverage()).
# @noRd
.trimSecondarySets <- function(sec) {
    if (length(sec) == 0L) {
        return(NULL)
    }
    map(sec, function(s) list(sets = s$sets))
}

# The alpha rows and prior variances of the fit's effects of positive prior
# variance, named by effect id, and alpha's columns by `fitIds` when given.
# An already trimmed fit's names are its ids.
# @noRd
.trimEffects <- function(fit, fitIds = NULL) {
    alpha <- fit[["alpha"]]
    priorVariance <- fit[["V"]] %||% fit[["prior_variance"]]
    if (is.null(alpha) || is.null(priorVariance)) {
        return(list())
    }
    alpha <- as.matrix(alpha)
    ids <- if (.namedByEffect(alpha)) {
        rownames(alpha)
    } else {
        str_c("L", seq_len(nrow(alpha)))
    }
    keep <- which(as.numeric(priorVariance) > 0)
    list(
        alpha = `dimnames<-`(
            alpha[keep, , drop = FALSE],
            list(
                ids[keep],
                if (length(fitIds) == ncol(alpha)) fitIds
            )
        ),
        V = set_names(as.numeric(priorVariance)[keep], ids[keep])
    )
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
setMethod(
    "retainedMass",
    "FineMappingResultBase",
    function(x) {
        .asViewDataFrame(.rcMassTable(x = x))
    }
)

# retainedMass() as a tibble.
# @noRd
.rcMassTable <- function(x) {
    if (nrow(x) == 0L) {
        return(.rcEmptyMass(x))
    }
    parts <- compact(map(seq_len(nrow(x)), .rcMassForRow, x = x))
    if (length(parts) == 0L) {
        return(.rcEmptyMass(x))
    }
    bind_rows(parts)
}

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

# Select the single FineMappingRow addressed by a tuple key or a GWAS block's
# `blockId`. Each concrete subclass
# implements this with its own row selector; the delegating accessors below
# then live once on the base and route through it.
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
        blockId = NULL,
        coverage = 0.95,
        minPurity = NULL
    ) {
        # Selectors pinning one entry -> that entry's bare credible-set table;
        # no /
        # partial selectors -> aggregate every matching entry's credible sets,
        # prefixed with the row identity (study/context/trait/blockId/method).
        # `minPurity` is an independent CS-quality filter, orthogonal to
        # coverage. Ranged either way, so it narrows with subsetByOverlaps().
        tbl <- .fmrAggregateView(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            blockId = blockId,
            perEntry = .fmrRowCs,
            viewArgs = list(
                coverage = coverage,
                minPurity = minPurity
            )
        )
        .withGenome(.variantTableGRanges(tbl), .oneBuild(x))
    }
)

#' @rdname lbf
#' @export
setMethod("lbf", "FineMappingResultBase", function(x) {
    .fitCoverageCheckedRows(x, "lbf()")
    .withGenome(
        .variantTableGRanges(.fmrAggregateView(x, perEntry = .fmrRowLbf)),
        .oneBuild(x)
    )
})

#' @rdname credibleSetSummary
#' @export
setMethod(
    "credibleSetSummary",
    "FineMappingResultBase",
    function(x, coverage = 0.95) {
        .fitCoverageCheckedRows(x, "credibleSetSummary()")
        .asViewDataFrame(.fmrCredibleSetSummaryTable(
            x = x,
            coverage = coverage
        ))
    }
)

# credibleSetSummary() as a tibble.
# @noRd
.fmrCredibleSetSummaryTable <- function(x, coverage = 0.95) {
    .fmrAggregateView(
        x,
        perEntry = .fmrRowCredibleSetSummary,
        viewArgs = list(coverage = coverage)
    )
}

#' @rdname fsusieCredibleBand
#' @export
setMethod("fsusieCredibleBand", "FineMappingResultBase", function(x) {
    .fitCoverageCheckedRows(x, "fsusieCredibleBand()")
    .withGenome(
        .bandTableGRanges(
            .fmrAggregateView(x, perEntry = .fmrRowFsusieCredibleBand)
        ),
        .oneBuild(x)
    )
})

#' @rdname fsusieTraitEffects
#' @export
setMethod(
    "fsusieTraitEffects",
    "FineMappingResultBase",
    function(x) {
        .fitCoverageCheckedRows(x, "fsusieTraitEffects()")
        .asViewDataFrame(.fsusieTraitEffectsTable(x = x))
    }
)

# fsusieTraitEffects() as a tibble, for the writers.
# @noRd
.fsusieTraitEffectsTable <- function(x) {
    # Not .fmrAggregateView(): the per-entry payload carries only the fit, and
    # this view also needs the row's `jointTraits` column for the trait names.
    meta <- .fmrRowMetadata(x)
    parts <- compact(map(
        seq_len(nrow(x)),
        .fsusieEntryTraitEffects,
        x = x,
        meta = meta
    ))
    if (length(parts) == 0L) {
        return(bind_cols(slice(meta, integer(0)), .emptyTraitEffects()))
    }
    .rbindAligned(parts)
}

#' @rdname fsusieAffectedRegions
#' @export
setMethod("fsusieAffectedRegions", "FineMappingResultBase", function(x) {
    .fitCoverageCheckedRows(x, "fsusieAffectedRegions()")
    perRow <- map(seq_len(nrow(x)), .fsusieEntryAffectedRegions, x = x)
    grs <- perRow[lengths(perRow) > 0L]
    if (length(grs) == 0L) {
        return(GenomicRanges::GRanges())
    }
    .withGenome(exec(c, !!!grs), .oneBuild(x))
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
        blockId = NULL,
        minPurity = NULL
    ) {
        sel <- list(
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            blockId = blockId
        )
        # The per-variant table first -- bare when the selectors pin one
        # entry, else the matching rows' tables stacked with row-identity
        # columns -- then one conversion to ranges. Converting last is what
        # lets the aggregate form be a GRanges too: the identity columns just
        # become mcols, where the old type = "GRanges" branch had to refuse
        # anything but a single pinned entry.
        tbl <- .fmrbTopLociTable(
            x,
            sel,
            list(signalCutoff = signalCutoff, minPurity = minPurity)
        )
        .withGenome(.variantTableGRanges(tbl), .oneBuild(x))
    }
)

# The per-variant top-loci table: the projected posterior view `viewArgs`
# configures, stacked with row-identity columns when the selectors match more
# than one entry. The
# caller converts to ranges; this stays a table so the aggregate form can be
# stacked before conversion.
# @noRd
.fmrbTopLociTable <- function(x, selectors, viewArgs) {
    exec(
        .fmrAggregateView,
        x,
        !!!selectors,
        perEntry = .fmrRowTopLoci,
        viewArgs = viewArgs
    )
}
# Marginal effects are per variant, so they come back ranged like the
# posterior view.
#' @rdname marginalEffects
#' @export
setMethod(
    "marginalEffects",
    "FineMappingResultBase",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        # Selectors pinning one entry -> that entry's bare marginal table; no /
        # partial selectors -> aggregate every matching entry's marginals,
        # prefixed with the row identity (study/context/trait/blockId/method).
        tbl <- .fmrAggregateView(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            blockId = blockId,
            perEntry = .fmrRowMarginalEffects
        )
        .withGenome(.variantTableGRanges(tbl), .oneBuild(x))
    }
)

# Collections saved while rows recorded the variants their fits were made
# with carry a `fitVariants` column. A row's variants are its fit's
# variables, read by position, so the column is dropped. The inherited update
# runs either way, for the updates every collection shares.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "FineMappingResultBase",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message(str_c("updateObject(object = '", class(object)[[1L]], "')"))
        }
        md <- mcols(object, use.names = FALSE)
        if (is_in("fitVariants", colnames(md))) {
            md$fitVariants <- NULL
            object <- `mcols<-`(object, value = md)
        }
        methods::callNextMethod(object, ..., verbose = verbose)
    }
)

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
        blockId = NULL
    ) {
        parts <- .fmrSelectEntry(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            blockId = blockId
        )
        .fmrPartsSusieFit(.fitCoverageChecked(parts, "susieFit()"))
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
        blockId = NULL
    ) {
        # The per-variant weight of the row a selector pins. Defined on the
        # collection because that is what fineMappingResult() now
        # returns; the body is the per-row primitive, so they cannot drift.
        #
        # The selectors are named, not absorbed by `...`: .fmrSelectEntry is the
        # only consumer, and each concrete class refuses the ones it does not
        # index by (a QtlFineMappingResult has no `blockId` axis).
        .fmrRowResolveWeights(.fmrSelectEntry(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            blockId = blockId
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
        blockId = NULL
    ) {
        .fmrPartsVariantIds(.fmrSelectEntry(
            x,
            studyName = studyName,
            context = context,
            trait = trait,
            method = method,
            blockId = blockId
        ))
    }
)

#' @rdname variantIds
#' @export
setMethod(
    "variantIds",
    "SumStatsBase",
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        # One method covers GwasSumStats and QtlSumStats. The selectors are the
        # union: GwasSumStats reads only `studyName` and refuses the
        # QtlSumStats-only ones rather than accepting and ignoring them.
        #
        # Rendered with the same .grVariantIds() the row classes use, so an id
        # means the same string whichever object produced it.
        .grVariantIds(sumStats(
            x,
            studyName = studyName,
            context = context,
            trait = trait
        ))
    }
)
