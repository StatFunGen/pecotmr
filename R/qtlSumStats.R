# =============================================================================
# QtlSumStats S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study, context,
# trait). Each row holds a per-tuple GRanges of summary statistics
# (variant_id + per-variant Z/N/MAF mcols). Class-level slots ldSketch
# (a GenotypeHandle for the LD reference) + genome (the genome build)
# apply uniformly across rows. Built-in qcInfo slot tracks which
# summaryStatsQc() passes have been run.
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

setClass(
    "QtlSumStats",
    contains = "SumStatsBase",
    validity = function(object) .validateQtlSumStats(object)
)

# ---- QtlSumStats validity helpers ------------------------------------------

# Collect all contract violations (empty vector = valid). The per-entry checks
# run only once the basic slot/column checks pass (they assume those columns).
# @noRd
.validateQtlSumStats <- function(object) {
    errors <- c(
        .qssCheckLdSketch(object),
        .qssCheckRequiredCols(object),
        .qssCheckGenome(object),
        .qssCheckQcInfo(object),
        .validateTraitPosColumn(object)
    )
    if (length(errors) == 0L) {
        errors <- .qssCheckEntries(object)
    }
    if (length(errors) == 0L) TRUE else errors
}

# ldSketch must be a GenotypeHandle or NULL.
# @noRd
.qssCheckLdSketch <- function(object) {
    if (
        !is.null(object@ldSketch) &&
            !methods::is(object@ldSketch, "GenotypeHandle")
    ) {
        return("'ldSketch' must be a GenotypeHandle or NULL")
    }
    NULL
}

# The study/context/trait metadata columns must be present. `names(object)` is
# element names on a RangedTupleList, so the check reads mcols directly.
# @noRd
.qssCheckRequiredCols <- function(object) {
    missingCols <- setdiff(
        c("study", "context", "trait"),
        colnames(mcols(object))
    )
    if (length(missingCols) > 0L) {
        return(str_c("missing columns: ", str_flatten(missingCols, ", ")))
    }
    NULL
}

# genome slot: a single non-empty string.
# @noRd
.qssCheckGenome <- function(object) {
    if (length(object@genome) != 1L || str_length(object@genome) == 0L) {
        return("'genome' slot must be a single non-empty character string")
    }
    NULL
}

# qcInfo slot must be a list.
# @noRd
.qssCheckQcInfo <- function(object) {
    if (!is.list(object@qcInfo)) {
        return("'qcInfo' slot must be a list")
    }
    NULL
}

# Element contract. The elements ARE GRanges by construction now -- the
# container is a GRangesList -- and RangedTupleList's validity enforces the
# one-seqname/one-strand invariant, so only tuple uniqueness is left to check.
# @noRd
.qssCheckEntries <- function(object) {
    .qssCheckTupleUniqueness(object)
}

# (study, context, trait) tuple uniqueness.
# @noRd
.qssCheckTupleUniqueness <- function(object) {
    if (.ssHasDuplicateKeys(object, c("study", "context", "trait"))) {
        return("(study, context, trait, range) uniqueness violated")
    }
    NULL
}

#' @title QTL Summary Statistics Handling
#' @description Constructor and accessor methods for \code{QtlSumStats}, the
#'   DFrame-subclass collection keyed by \code{(study, context, trait)}.
#' @name pecotmr-qtl-sumstats
#' @keywords internal
#' @importFrom GenomicRanges GRanges seqnames start
#' @importFrom S4Vectors DataFrame SimpleList mcols
#' @include AllGenerics.R
NULL

# =============================================================================
# Constructor
# =============================================================================

#' @title Create a QtlSumStats Collection Object
#' @description Construct a \code{QtlSumStats} S4 DFrame-subclass collection
#'   from per-tuple vectors and a list of \code{GRanges} entries (one per
#'   tuple), plus a single LD sketch handle and a single genome build that apply
#'   to the whole collection. Each \code{GRanges} entry must carry per-variant
#'   statistics in its mcols (\code{SNP}, \code{A1}, \code{A2}, \code{Z},
#'   \code{N}; plus optional \code{MAF}, \code{INFO}, \code{BETA}, \code{SE},
#'   \code{P}).
#' @param study Character vector of study identifiers (per tuple).
#' @param context Character vector of context labels (per tuple).
#' @param trait Character vector of trait identifiers (per tuple).
#' @param entry A list / \code{SimpleList} of \code{GRanges}, one per tuple.
#'   Same length as \code{study}, \code{context}, and \code{trait}.
#' @param genome Single character string giving the genome build (e.g.,
#'   \code{"hg19"}, \code{"hg38"}). Uniform across the collection because all
#'   entries share the same LD sketch.
#' @param ldSketch A \code{GenotypeHandle} carrying the LD reference.
#' @param varY Optional numeric vector of per-tuple phenotype variances
#'   (\code{NA_real_} entries allowed).
#' @param nSample Optional per-tuple total sample size (numeric; default
#'   \code{NULL}). Attached only when supplied (length 1 or length(study)). Used
#'   as the study-level fallback for the per-variant \code{N} when a tuple has
#'   no per-variant \code{N} column. Named \code{nSample} to avoid clashing with
#'   \code{getNSamples()} (the LD-panel sample size). Unlike GWAS, QTL
#'   collections carry no case/control counts (molecular traits are
#'   quantitative), so only this total-N fallback is exposed.
#' @param ... Additional per-tuple columns to attach to the collection.
#' @param qcInfo A \code{list} recording which QC steps ran. Empty \code{list()}
#'   on construction; populated by \code{summaryStatsQc()} with a per-step audit
#'   record. Fine-mapping / TWAS pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0}.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @return A \code{QtlSumStats} object.
#' @examples
#' gh <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100 * 1:3, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(SNP = paste0("rs", 1:3),
#'   A1 = "A", A2 = "G", Z = rnorm(3), N = 100L)
#' QtlSumStats(study = "s1", context = "brain", trait = "g1", entry = list(gr),
#'   genome = "hg38", ldSketch = gh)
#' @export
QtlSumStats <- function(
    study,
    context,
    trait,
    entry,
    genome,
    ldSketch = NULL,
    varY = NA_real_,
    nSample = NULL,
    qcInfo = list(),
    traitPos = NULL,
    ...
) {
    if (
        missing(study) ||
            missing(context) ||
            missing(trait) ||
            missing(entry) ||
            missing(genome)
    ) {
        msg <- glue(
            "`study`, `context`, `trait`, `entry`, and `genome` ",
            "are all required."
        )
        abort(msg)
    }
    n <- length(study)
    varY <- .qssValidateArgs(context, trait, entry, genome, varY, n)
    entry <- .qssAddTraitDistances(entry, traitPos, n)
    cols <- .qssBaseCols(study, context, trait, varY)
    cols <- .qssAppendNSample(cols, nSample, n)
    cols <- .appendTraitPosCol(cols, traitPos, n)
    cols <- .qssAppendExtras(cols, list(...))
    dfArgs <- c(cols, list(check.names = FALSE))
    # The per-tuple GRanges become the collection's ELEMENTS; the tuple keys
    # and per-tuple scalars go in mcols. There is no `entry` column.
    # mcols are attached to the GRangesList BEFORE new(), because new()
    # validates during initialize() and the validity method needs the identity
    # columns to already be there.
    # A multi-seqname entry (e.g. a genome-wide GWAS) is split into one
    # element per chromosome, with its metadata row replicated alongside.
    split <- .rtlSplitBySeqname(entry)
    grl <- GenomicRanges::GRangesList(split$entry)
    md <- exec(S4Vectors::DataFrame, !!!dfArgs)
    mcols(grl) <- md[split$fromIdx, , drop = FALSE]
    obj <- methods::new(
        "QtlSumStats",
        grl,
        ldSketch = ldSketch,
        genome = as.character(genome),
        qcInfo = as.list(qcInfo)
    )
    methods::validObject(obj)
    obj
}

# Validate genome / entry / length consistency and recycle varY. Returns the
# (possibly recycled) varY.
# @noRd
.qssValidateArgs <- function(context, trait, entry, genome, varY, n) {
    if (length(genome) != 1L) {
        msg <- glue(
            "`genome` must be a single character string (one build per ",
            "collection, because all entries share the LD sketch)."
        )
        abort(msg)
    }
    if (!is.list(entry)) {
        abort(
            "`entry` must be a list (or SimpleList) of GRanges, one per tuple."
        )
    }
    if (length(context) != n || length(trait) != n || length(entry) != n) {
        msg <- glue(
            "`study`, `context`, `trait`, and `entry` must all have the ",
            "same length."
        )
        abort(msg)
    }
    .qssRecycleTo(varY, n, "varY")
}

# Recycle a length-1 vector to n, or require length n exactly.
# @noRd
.qssRecycleTo <- function(x, n, name) {
    if (length(x) == 1L && n > 1L) {
        x <- rep(x, n)
    }
    if (length(x) != n) {
        msg <- glue("`{name}` must have length 1 or length(study).")
        abort(msg)
    }
    x
}

# Base column list (with tss/tes distances auto-filled onto each entry from the
# trait position when supplied; the authoritative shape check is in
# .appendTraitPosCol).
# @noRd
.qssBaseCols <- function(study, context, trait, varY) {
    list(
        study = as.character(study),
        context = as.character(context),
        trait = as.character(trait),
        varY = as.numeric(varY)
    )
}

# Add the per-variant distance-to-trait-position column to each element before
# it becomes a collection element. Split out of .qssBaseCols because the
# elements are no longer one of the metadata columns -- they are the container
# itself.
# @noRd
.qssAddTraitDistances <- function(entry, traitPos, n) {
    computable <- !is.null(traitPos) &&
        methods::is(traitPos, "GRanges") &&
        length(traitPos) == n
    if (!computable) {
        return(entry)
    }
    .appendTraitDistances(entry, traitPos)
}

# Attach the OPTIONAL per-tuple total-sample-size column (default NULL leaves
# the
# original schema; summaryStatsQc() fills a missing per-variant N from it).
# @noRd
.qssAppendNSample <- function(cols, nSample, n) {
    if (is.null(nSample)) {
        return(cols)
    }
    cols$nSample <- as.numeric(.qssRecycleTo(nSample, n, "nSample"))
    cols
}

# Append any user-supplied extra columns (from `...`).
# @noRd
.qssAppendExtras <- function(cols, extras) {
    for (nm in names(extras)) {
        cols[[nm]] <- extras[[nm]]
    }
    cols
}

# Annotate each entry's variants with tss_distance / tes_distance from the
# trait-position GRanges (one range per tuple, aligned to `entry`): variant
# position minus the trait's TSS (= start of the trait-position range) and TES
# (= end). For a single-base trait position start() == end(), so the two are
# identical (as expected for point positions). Existing tss_distance /
# tes_distance mcols are preserved (not clobbered). `@noRd`
.appendTraitDistances <- function(entry, traitPos) {
    if (is.null(traitPos)) {
        return(entry)
    }
    map(
        seq_along(entry),
        .qssAppendTraitDist,
        entry = entry,
        traitPos = traitPos
    )
}

# =============================================================================
# Accessors
# =============================================================================

# Internal: resolve a (study, context, trait) tuple to its element indices.
# Returns a VECTOR: a tuple whose entry spanned several chromosomes was split
# into one element per seqname at construction.
.qtlSumStatsSelectRow <- function(x, study, context, trait) {
    if (nrow(x) == 0L) {
        abort("QtlSumStats has no rows.")
    }
    anyUnset <- missing(study) ||
        is.null(study) ||
        missing(context) ||
        is.null(context) ||
        missing(trait) ||
        is.null(trait)
    if (anyUnset) {
        if (.qssSingleTuple(x)) {
            return(seq_len(nrow(x)))
        }
        msg <- glue(
            "This QtlSumStats has {nrow(x)} entries. Pass `study`, ",
            "`context`, and `trait` to select one."
        )
        abort(msg)
    }
    if (length(study) != 1L || length(context) != 1L || length(trait) != 1L) {
        abort("`study`, `context`, and `trait` must each be length 1.")
    }
    .qssMatchTuple(x, study, context, trait)
}

# TRUE when every element belongs to the same (study, context, trait) tuple,
# so an unqualified accessor call is unambiguous even across a seqname split.
# @noRd
.qssSingleTuple <- function(x) {
    md <- mcols(x)
    keys <- str_c(
        as.character(md$study),
        as.character(md$context),
        as.character(md$trait),
        sep = "|"
    )
    n_distinct(keys) == 1L
}

# Resolve the element indices for a (study, context, trait) tuple.
# @noRd
.qssMatchTuple <- function(x, study, context, trait) {
    idx <- .matchTupleRows(
        x,
        list(study = study, context = context, trait = trait)
    )
    if (length(idx) == 0L) {
        msg <- glue(
            "No entry for (study='{study}', context='{context}', ",
            "trait='{trait}')."
        )
        abort(msg)
    }
    # More than one match is expected, not a violation: uniqueness is keyed on
    # (tuple, range), so a tuple split across chromosomes owns several
    # elements. The caller stitches them.
    idx
}

#' @rdname getSumStats
#' @param annotateSignificance Optional correction-method name
#'   (\code{"permutation"} / \code{"bonferroni_original"} /
#'   \code{"bonferroni_filtered"} / \code{"qvalue"}). When set on a QtlSumStats
#'   enriched by \code{\link{qtlAssociationPostprocess}}, a logical
#'   \code{significant} mcol for that method is added to the returned entry (the
#'   significance is derived on the fly, not stored). Flat export flattens this
#'   full entry GRanges (all mcols) directly; note \code{\link{getSumstatDf}} is
#'   a fixed GWAS-schema view and does not carry the association columns.
#' @export
setMethod(
    "getSumStats",
    signature(x = "QtlSumStats"),
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL,
        ranges = NULL,
        ...
    ) {
        idx <- .qtlSumStatsSelectRow(x, study, context, trait)
        gr <- .ssStitchElements(x, idx, ranges)
        if (!is.null(annotateSignificance)) {
            m <- arg_match(
                annotateSignificance,
                c(
                    "permutation",
                    "bonferroni_original",
                    "bonferroni_filtered",
                    "qvalue"
                )
            )
            S4Vectors::mcols(gr)[["significant"]] <- .qapSignificanceMask(
                x,
                m
            )[[idx]]
        }
        gr
    }
)

# getZ / getN / getMaf / nSnps are provided once by SumStatsBase (AllClasses.R);
# they only delegate to getSumStats().

#' @rdname getSumstatDf
#' @export
setMethod(
    "getSumstatDf",
    "QtlSumStats",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        require = character(0),
        derive = c("none", "zFromBetaSe"),
        keepChrPrefix = TRUE
    ) {
        derive <- arg_match(derive)
        gr <- getSumStats(x, study = study, context = context, trait = trait)
        .entryToSumstatDf(
            gr,
            require = require,
            derive = derive,
            keepChrPrefix = keepChrPrefix,
            label = glue(
                "QtlSumStats[{if (is.null(study)) '<auto>' else study}/",
                "{if (is.null(context)) '<auto>' else context}/",
                "{if (is.null(trait)) '<auto>' else trait}]"
            )
        )
    }
)


#' @rdname getVarY
#' @export
setMethod(
    "getVarY",
    "QtlSumStats",
    function(x, study = NULL, context = NULL, trait = NULL) {
        idx <- .qtlSumStatsSelectRow(x, study, context, trait)
        val <- x$varY[[idx]]
        if (is.na(val)) NULL else val
    }
)

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlSumStats", function(x) {
    unique(as.character(x$context))
})

#' @rdname getTraits
#' @export
setMethod("getTraits", "QtlSumStats", function(x) unique(as.character(x$trait)))

# Trait position provenance. Optional for a QtlSumStats -- it cannot be inferred
# from summary statistics, so it is only present when the caller supplied it.
# Returns the traitPos GRanges (whole column, or the rows matching `traitId`),
# or a scalar NA when no trait position was supplied.
#' @rdname getTraitPosition
#' @export
setMethod("getTraitPosition", "QtlSumStats", function(x, traitId = NULL, ...) {
    tp <- .getTraitPosColumn(x)
    if (!methods::is(tp, "GRanges") || is.null(traitId)) {
        return(tp)
    }
    idx <- which(is_in(as.character(x$trait), as.character(traitId)))
    if (length(idx) == 0L) {
        return(NA)
    }
    # The trait position is per-TUPLE, but a tuple split across chromosomes
    # owns several elements and the column is replicated across them, so the
    # matches collapse back to the one position.
    unique(tp[idx])
})

# =============================================================================
# Show method
# =============================================================================

#' @rdname show-methods
#' @export
setMethod("show", "QtlSumStats", function(object) {
    cat(glue(
        "QtlSumStats: {nrow(object)} entries, ",
        "genome build {getGenome(object)}\n",
        .trim = FALSE
    ))
    if (nrow(object) > 0L) {
        cat(glue(
            "  {n_distinct(object$study)} studies, ",
            "{n_distinct(object$context)} contexts, ",
            "{n_distinct(object$trait)} traits\n",
            .trim = FALSE
        ))
    }
    ld <- getLdSketch(object)
    ldSrc <- if (is.null(ld)) {
        "none (LD-free)"
    } else {
        glue("{getFormat(ld)} @ {getPath(ld)}")
    }
    cat(glue("  LD sketch: {ldSrc}\n", .trim = FALSE))
})

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Column `cn` of `object`, via `[[` (preserves the required `entry` column that
# `object[, cn]` would drop).
# @noRd

# Append tss_distance / tes_distance mcols to entry `i` from traitPos[i] (no-op
# when the entry is empty or the trait position is unset). Existing distance
# mcols are preserved.
# @noRd
.qssAppendTraitDist <- function(i, entry, traitPos) {
    gr <- entry[[i]]
    if (length(gr) == 0L) {
        return(gr)
    }
    tp <- traitPos[i]
    tssPos <- suppressWarnings(GenomicRanges::start(tp))
    tesPos <- suppressWarnings(GenomicRanges::end(tp))
    if (length(tssPos) == 0L || is.na(tssPos)) {
        return(gr)
    }
    pos <- GenomicRanges::start(gr)
    mc <- S4Vectors::mcols(gr)
    if (is.null(mc) || ncol(mc) == 0L) {
        mc <- S4Vectors::DataFrame(row.names = NULL)
    }
    if (is.null(mc[["tss_distance"]])) {
        mc[["tss_distance"]] <- pos - tssPos
    }
    if (is.null(mc[["tes_distance"]])) {
        mc[["tes_distance"]] <- pos - tesPos
    }
    S4Vectors::mcols(gr) <- mc
    gr
}

# Entry `i` of a SumStats collection restricted to chromosome `chrName`. Shared
# by the QtlSumStats and GwasSumStats subsetChr methods.
# @noRd
