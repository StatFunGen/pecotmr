# =============================================================================
# GwasSumStats S4 class
# -----------------------------------------------------------------------------
# GRangesList-subclass collection keyed by the identity tuple (study). A study
# is split into one row per chromosome, or per window when the constructor is
# given a `region`, each holding that part of its GWAS summary statistics;
# the accessors stitch a study's rows back together. Class-level slots genome
# + qcOptions apply uniformly across rows; ldReference holds each distinct LD
# reference once, and the ldReference mcols column names the one every row
# uses.
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

#' @title GWAS Summary-Statistic Collection
#' @description S4 collection of GWAS summary statistics keyed by the identity
#'   tuple \code{(study)}. A study is split into one element per chromosome,
#'   or per window when \code{\link{GwasSumStats}} is given a \code{region},
#'   each a per-variant \code{GRanges}; \code{\link{sumStats}} and the other
#'   accessors stitch a study's elements back together.
#' @details Required column: \code{study}; the constructor also adds
#'   \code{blockId}, naming each element's chromosome or window, and the
#'   element's analysis window: \code{region}, a \code{GRanges} with one
#'   range per element, and \code{regionGiven}. \code{regionGiven} is
#'   \code{TRUE} where the window was given as a \code{region}, and
#'   \code{FALSE} where it is the span of the element's variants when the
#'   collection was first built, before any QC. A QC pass carries both
#'   unchanged; the span of the variants left is \code{range(x)}.
#'   \code{\link{colocPipeline}},
#'   \code{\link{colocboostPipeline}} and
#'   \code{\link{causalInferencePipeline}} check the QTL variants against
#'   these windows. Variants past a given window lie in one that was not
#'   loaded, which is an error. Variants past a computed window may only be
#'   ones the GWAS lacks, so more than one in a row is a warning.
#'   A study repeats across its elements, so uniqueness is on
#'   \code{(study, range)}. The class-level
#'   slots inherited from \code{\linkS4class{SumStatsBase}} -- \code{genome}
#'   and \code{qcOptions} -- apply uniformly to every row; \code{ldReference}
#'   holds each distinct LD reference once, and the \code{ldReference} column
#'   names the one each row uses.
#' @seealso \code{\link{GwasSumStats}} for the constructor and
#'   \code{\linkS4class{QtlSumStats}} for the QTL counterpart.
#' @importFrom checkmate makeAssertCollection assertNames assertList
#' @export
setClass(
    "GwasSumStats",
    contains = "SumStatsBase",
    validity = function(object) {
        # The ldReference slot's class union enforces its type.
        coll <- makeAssertCollection()
        assertNames(
            colnames(mcols(object)) %||% character(0),
            must.include = "study",
            what = "colnames",
            .var.name = "mcols",
            add = coll
        )
        coll$push(.sumStatsCheckGenome(object))
        assertList(object@qcOptions, .var.name = "qcOptions", add = coll)
        slotErrors <- coll$getMessages()
        if (length(slotErrors) > 0L) {
            return(slotErrors)
        }
        # The elements ARE GRanges by construction now -- the container is
        # a GRangesList -- so the old per-element type and length checks
        # are gone. The one-seqname/one-strand invariant is enforced by
        # RangedTupleList's own validity.
        # Keyed on (study, range), not study alone: a study split across
        # chromosomes contributes one element per seqname, so the study
        # label legitimately repeats.
        if (.ssHasDuplicateKeys(object, "study")) {
            return("(study, range) must be unique")
        }
        TRUE
    }
)


#' @rdname show-methods
setMethod("show", "GwasSumStats", function(object) {
    cat(glue(
        "GwasSumStats: {n_distinct(object$study)} studies in ",
        "{nrow(object)} rows, genome build {.showGenomeBuild(object)}\n",
        .trim = FALSE
    ))
    ldSrc <- .ldRefsLabel(ldReference(object), "none (LD-free)")
    cat(glue("  LD reference: {ldSrc}\n", .trim = FALSE))
})


#' @title GWAS Summary Statistics Handling
#' @description Constructor, accessors, and converters for
#'   \code{GwasSumStats}, the \code{RangedTupleList} collection of GWAS
#'   summary statistics keyed by \code{study}.
#' @name pecotmr-gwas-sumstats
#' @keywords internal
#' @importFrom GenomicRanges GRanges seqnames start
#' @importFrom S4Vectors DataFrame mcols mcols<- SimpleList
#' @importFrom IRanges IRanges
#' @include AllGenerics.R
NULL

# =============================================================================
# Constructor
# =============================================================================

# Recycle a length-1 per-study scalar to one value per study (or validate an
# already-per-study vector), coerced to numeric.
# @noRd
.recyclePerStudy <- function(v, nm, studyName) {
    if (length(v) == 1L && length(studyName) > 1L) {
        v <- rep(v, length(studyName))
    }
    if (length(v) != length(studyName)) {
        msg <- glue("`{nm}` must have length 1 or length(study).")
        abort(msg)
    }
    as.numeric(v)
}

#' @title Create a GwasSumStats Collection Object
#' @description Construct a \code{GwasSumStats} S4 collection from per-study
#'   tuple vectors and a list of \code{GRanges} entries (one per study), plus
#'   the LD reference(s) the rows were harmonized against and a single genome
#'   build that applies to the whole collection.
#'
#' Each \code{GRanges} entry must carry per-variant statistics in its mcols (at
#' minimum \code{SNP}, \code{A1}, \code{A2}, \code{Z}, \code{N}; optionally
#' \code{MAF}, \code{INFO}, \code{BETA}, \code{SE}, \code{P}).
#' @param studyName Character vector of study identifiers, one per
#'   \code{entry}. A study whose entries are already split by block repeats
#'   here, once per entry, with \code{blockId} telling them apart.
#' @param entry A \code{SimpleList} or \code{list} of \code{GRanges}, one per
#'   study (or per block of an already-split study).
#' @param genome Single character string giving the genome build (e.g.,
#'   \code{"hg19"}, \code{"hg38"}). One build per collection.
#' @param ldReference The LD reference: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}, used for every study, or a list with one per
#'   study when they sit on different references. Each distinct reference is
#'   stored once, and the \code{ldReference} column names each row's.
#' @param varY Optional numeric vector of per-study phenotype variances
#'   (\code{NA_real_} entries allowed). Used by the sufficient-statistic
#'   interface; z-score RSS analyses should leave entries as NA.
#' @param nCase,nControl Optional per-study case / control counts. The columns
#'   are attached \strong{only when supplied} (default \code{NULL}), so
#'   quantitative-trait collections keep the original schema. When given, pass
#'   length 1 or length(study) (use \code{NA} for the non-case/control studies
#'   in a mixed collection). For case/control GWAS, downstream consumers (e.g.
#'   \code{\link{colocboostPipeline}}) use the effective sample size \code{4 /
#'   (1/nCase + 1/nControl)} in place of the per-variant \code{N}.
#' @param nSample Optional per-study total sample size (numeric; default
#'   \code{NULL}). Attached only when supplied (length 1 or length(study)). Used
#'   as the study-level fallback for the per-variant \code{N} when a study has
#'   no per-variant \code{N} column and no case/control counts. Named
#'   \code{nSample} to avoid clashing with \code{nSamples()} (the LD-panel
#'   sample size).
#' @param region Optional analysis windows, typically LD blocks: a
#'   \code{GRanges}, or the path of a BED file, read as BED with a 0-based
#'   start (a fourth column, when present, names each window). One set of
#'   windows applies to every study; a list gives one set per study,
#'   \code{NULL} for a study to split by chromosome. Each study's variants are
#'   split into one element per window they fall in, rather than one per
#'   chromosome, and each element's \code{region} is its window, recorded as
#'   given. Its \code{blockId} is the window's name: its \code{names}, else
#'   its \code{blockId} column (a BED file's fourth column), else its
#'   coordinates. cTWAS needs this granularity: its EM estimates parameters
#'   across blocks, so a per-chromosome split is too coarse. Variants outside
#'   every window are dropped with a warning.
#' @param blockId Optional character vector of block names, one per
#'   \code{entry}, for entries that are \strong{already} split. Use it to
#'   carry existing names through a rebuild; without it a rebuild would
#'   re-derive them as seqnames and collapse distinct blocks onto one name. A
#'   name only: it gives no window, which comes from \code{region}. Mutually
#'   exclusive with \code{region}.
#' @param extraCols Optional named list (or \code{DataFrame} /
#'   \code{data.frame}) of additional per-study columns to
#'   attach to the collection's \code{mcols}. Named rather than variadic
#'   because the base class's own constructor, \code{GRangesList(...)}, uses
#'   \code{...} for its \emph{elements}: a bare \code{...} here would read
#'   as adding an entry rather than a metadata column.
#' @param qcOptions A \code{list} recording the QC that ran. Empty
#'   \code{list()} on construction; \code{summaryStatsQc()} fills it.
#'   Fine-mapping / TWAS pipelines reject inputs where
#'   \code{length(qcOptions(x)) == 0}.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{GwasSumStats} object.
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100 * 1:3, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(SNP = paste0("rs", 1:3),
#'   A1 = "A", A2 = "G", Z = rnorm(3), N = 100L)
#' GwasSumStats(studyName = "t1", entry = list(gr), genome = "hg38",
#'   ldReference = panel)
#' @export
GwasSumStats <- function(
    studyName,
    entry,
    genome,
    ldReference = NULL,
    varY = NA_real_,
    nCase = NULL,
    nControl = NULL,
    nSample = NULL,
    qcOptions = list(),
    region = NULL,
    blockId = NULL,
    extraCols = list(),
    metadata = list()
) {
    if (missing(studyName) || missing(entry) || missing(genome)) {
        abort("`studyName`, `entry`, and `genome` are all required.")
    }
    varY <- .gwasValidateArgs(studyName, entry, genome, varY)
    entry <- .ucscEntries(entry)
    cols <- list(
        study = as.character(studyName),
        varY = varY
    ) |>
        .gwasAppendOptional(nCase, nControl, nSample, studyName) |>
        .gwasAppendExtras(extraCols)
    dict <- .ldRefDict(.ldRefPerTuple(ldReference, length(entry)))
    dfArgs <- c(cols, list(ldReference = dict$keys, check.names = FALSE))
    grl <- .gwasElements(entry, dfArgs, region, blockId)
    .withMetadata(
        .sumStatsNewValidated(
            "GwasSumStats",
            grl,
            dict$refs,
            genome,
            qcOptions
        ),
        metadata
    )
}

# The collection's elements with their metadata. The per-study GRanges become
# the ELEMENTS; everything else is per-study metadata and goes in mcols. There
# is no `entry` column. mcols are attached to the GRangesList BEFORE new(),
# because new() validates during initialize() and the validity method needs
# the identity columns to already be there.
# A multi-seqname entry (e.g. a genome-wide GWAS) is split into one element
# per `region` window (or per chromosome when there are none), with its
# metadata row replicated alongside. Splitting is unconditional: a stored
# element always spans exactly one seqname.
# @noRd
.gwasElements <- function(entry, dfArgs, region, blockId) {
    split <- .gwasSplitEntry(entry, region, blockId, length(entry))
    pieces <- .rtlBindElements(split$entry)
    S4Vectors::`mcols<-`(
        pieces,
        value = .gwasPieceMetadata(dfArgs, split, pieces)
    )
}

# Each piece's metadata row: its entry's, replicated across the pieces the
# entry splits into, then the block the piece covers. `blockId` is always
# present, so downstream code (cTWAS in particular) can key regions without
# first asking how the collection was built. So are `region` and
# `regionGiven`, so the block range of a row, and whether the caller gave it,
# need no fallback either.
# @noRd
.gwasPieceMetadata <- function(dfArgs, split, pieces) {
    cbind(
        exec(S4Vectors::DataFrame, !!!dfArgs)[split$fromIdx, , drop = FALSE],
        exec(
            S4Vectors::DataFrame,
            blockId = split$blockId,
            !!!.regionColumns(
                split$given,
                split$hasGiven,
                .rowSpanRanges(pieces)
            )
        )
    )
}

# Build and validate a SumStatsBase subclass from its GRangesList plus the
# three collection-level slots GwasSumStats and QtlSumStats share. Shared
# rather than duplicated so the two constructors cannot drift in how they
# coerce those slots. Used by qtlSumStats.R too. `ldReference` is the named
# list of distinct references (.ldRefDict()), already keyed from the rows.
# @noRd
.sumStatsNewValidated <- function(Class, grl, ldReference, genome, qcOptions) {
    # The build goes into seqinfo, which is where a GRangesList keeps it and
    # where every Bioconductor consumer reads it from. Assigned before new()
    # so validity sees the finished object.
    built <- .rtlNamed(
        GenomeInfoDb::`genome<-`(grl, value = as.character(genome))
    )
    methods::new(
        Class,
        built,
        ldReference = ldReference,
        qcOptions = as.list(qcOptions)
    )
}

# Split across the `region` windows when they are given, else by seqname.
# Both return (entry, fromIdx); this adds each piece's blockId -- its window's
# name, or on the seqname path the seqname it sits on -- and its given window
# (`given`, with `hasGiven` saying which pieces have one).
# @noRd
.gwasSplitEntry <- function(entry, region, blockId, n) {
    if (!is.null(region) && !is.null(blockId)) {
        msg <- glue(
            "pass `region` (its windows name the rows) or `blockId` (supply ",
            "the names), not both."
        )
        abort(msg)
    }
    if (!is.null(region)) {
        return(.gwasSplitByRegion(entry, region, n))
    }
    split <- .rtlSplitBySeqname(entry)
    # Supplied ids are indexed by fromIdx, exactly like the other metadata
    # columns, so an entry that still splits further replicates its id rather
    # than falling out of alignment. This is what lets a rebuild (QC) carry
    # block keys through instead of silently re-deriving them as seqnames.
    # unname(): the seqname splitter names its pieces, and those names would
    # otherwise ride into the mcols column and make it inconsistent with the
    # window path, which produces a bare character vector.
    ids <- if (is.null(blockId)) {
        unname(map_chr(split$entry, .rtlElementSeqname))
    } else {
        .gwasCheckBlockId(blockId, n)[split$fromIdx]
    }
    list_assign(
        split,
        blockId = ids,
        given = .placeholderRanges(length(split$entry)),
        hasGiven = rep(FALSE, length(split$entry))
    )
}

# Each study split across its windows: one set shared by every study, or a
# list with one per study (NULL for a study to split by chromosome), as
# `ldReference` takes one reference or a list.
# @noRd
.gwasSplitByRegion <- function(entry, region, n) {
    if (!.isPerStudyRegion(region)) {
        return(.gwasWindowSplit(entry, .asRegionRanges(region)))
    }
    if (length(region) != n) {
        abort(glue(
            "`region` as a list must have one element per `entry` ",
            "(got {length(region)} vs {n})."
        ))
    }
    .gwasCombineSplits(map2(as.list(entry), as.list(region), .gwasSplitStudy))
}

# @noRd
.isPerStudyRegion <- function(region) {
    methods::is(region, "GRangesList") ||
        (is.list(region) && !is.data.frame(region))
}

# @noRd
.gwasWindowSplit <- function(entry, windows) {
    split <- .rtlSplitByBlocks(entry, windows)
    list_assign(
        split,
        given = .bareRanges(split$blockRange),
        hasGiven = rep(TRUE, length(split$entry)),
        blockRange = zap()
    )
}

# One study against its own windows, or by chromosome when it has none.
# @noRd
.gwasSplitStudy <- function(e, r) {
    if (is.null(r)) {
        return(.gwasSplitEntry(list(e), NULL, NULL, 1L))
    }
    .gwasWindowSplit(list(e), .asRegionRanges(r))
}

# Per-study splits as one: each study's pieces in turn, with `fromIdx`
# pointing back at the study they came from.
# @noRd
.gwasCombineSplits <- function(splits) {
    list(
        entry = list_flatten(map(splits, "entry")),
        fromIdx = rep(seq_along(splits), map_int(splits, .gwasSplitSize)),
        blockId = list_c(map(splits, "blockId"), ptype = character()),
        given = unname(unlist(.rtlBindElements(map(splits, "given")))),
        hasGiven = list_c(map(splits, "hasGiven"), ptype = logical())
    )
}

# @noRd
.gwasSplitSize <- function(split) {
    length(split$entry)
}

# @noRd
.gwasCheckBlockId <- function(blockId, n) {
    if (length(blockId) != n) {
        msg <- glue(
            "`blockId` must have one value per `entry` ",
            "(got {length(blockId)} vs {n})."
        )
        abort(msg)
    }
    as.character(blockId)
}

# Validate genome / entry / length consistency; returns the recycled varY.
# @noRd
.gwasValidateArgs <- function(studyName, entry, genome, varY) {
    if (length(genome) != 1L) {
        msg <- glue(
            "`genome` must be a single character string (one build per ",
            "collection)."
        )
        abort(msg)
    }
    if (!is.list(entry)) {
        abort(
            "`entry` must be a list (or SimpleList) of GRanges, one per study."
        )
    }
    if (length(entry) != length(studyName)) {
        msg <- glue(
            "length(entry) ({length(entry)}) must equal ",
            "length(studyName) ({length(studyName)})."
        )
        abort(msg)
    }
    .recyclePerStudy(varY, "varY", studyName)
}

# Attach the OPTIONAL per-study nCase / nControl / nSample columns (each only
# when supplied; NA for the non-case/control studies in a mixed collection).
# @noRd
.gwasAppendOptional <- function(cols, nCase, nControl, nSample, studyName) {
    c(
        cols,
        compact(list(
            nCase = if (!is.null(nCase)) {
                .recyclePerStudy(nCase, "nCase", studyName)
            },
            nControl = if (!is.null(nControl)) {
                .recyclePerStudy(nControl, "nControl", studyName)
            },
            nSample = if (!is.null(nSample)) {
                .recyclePerStudy(nSample, "nSample", studyName)
            }
        ))
    )
}

# Append any user-supplied extra columns (from `...`).
# @noRd
.gwasAppendExtras <- function(cols, extras) {
    c(cols, .extraColsList(extras))
}


# =============================================================================
# Accessors for the new GwasSumStats collection
# =============================================================================

# Internal: resolve a study selection to a single row index. Errors when
# `study` is missing on a multi-study collection.
# Element indices for one study. Returns a VECTOR, not a scalar: the seqname
# split means one study can own several elements (one per chromosome), and
# sumStats() stitches them back into the single GRanges callers expect.
# @noRd
.gwasSelectStudy <- function(x, studyName) {
    if (nrow(x) == 0L) {
        abort("GwasSumStats has no rows.")
    }
    studies <- as.character(x$study)
    if (missing(studyName) || is.null(studyName)) {
        if (n_distinct(studies) == 1L) {
            return(seq_len(nrow(x)))
        }
        msg <- glue(
            "This GwasSumStats has {n_distinct(studies)} studies. ",
            "Pass `studyName = <name>` to select one. ",
            "Available: {str_flatten(unique(studies), ', ')}"
        )
        abort(msg)
    }
    idx <- which(studies == as.character(studyName))
    if (length(idx) == 0L) {
        msg <- glue(
            "Unknown study: '{studyName}'. ",
            "Available: {str_flatten(unique(studies), ', ')}"
        )
        abort(msg)
    }
    idx
}

#' @title Get a GWAS Study's Summary-Statistic GRanges
#' @description Return the per-variant \code{GRanges} of summary statistics for
#'   one study in a \code{GwasSumStats} collection.
#' @param x A \code{GwasSumStats} object.
#' @param studyName Character (length 1) study identifier. Optional when the
#'   collection has a single row.
#' @param context,trait,annotateSignificance QTL-only selectors, named here so
#'   that supplying one is an error rather than a silently ignored request. A
#'   GWAS has no context or trait axis and carries no significance annotation;
#'   leave them \code{NULL}.
#' @return A \code{GRanges} object.
#' @export
setMethod(
    "sumStats",
    signature(x = "GwasSumStats"),
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL
    ) {
        # A GWAS has no context / trait axis and no significance annotation:
        # those selectors belong to QtlSumStats. They are named here, rather
        # than absorbed by `...`, so asking for one is an error instead of a
        # silently ignored request -- the shared SumStatsBase accessors pass
        # the union of both classes' selectors.
        .gwasRefuseQtlSelectors(context, trait, annotateSignificance)
        .ssStitchElements(x, .gwasSelectStudy(x, studyName))
    }
)

# @noRd
.gwasRefuseQtlSelectors <- function(context, trait, annotateSignificance) {
    given <- c(
        context = !is.null(context),
        trait = !is.null(trait),
        annotateSignificance = !is.null(annotateSignificance)
    )
    if (!any(given)) {
        return(invisible(NULL))
    }
    abort(glue(
        "sumStats(GwasSumStats): {str_flatten(names(given)[given], ', ')} ",
        "{if (sum(given) == 1L) 'is' else 'are'} a QtlSumStats selector. ",
        "A GWAS is selected by `studyName`; narrow it first with ",
        "subsetByOverlaps()."
    ))
}

# z / nSamples / maf / nSnps are provided once by SumStatsBase (AllClasses.R).

#' @rdname varY
#' @export
setMethod("varY", "GwasSumStats", function(x, studyName = NULL) {
    idx <- .gwasSelectStudy(x, studyName)
    # One value per study, held on each of its rows; rows that disagree have
    # no single answer, so that is an error rather than a pick.
    val <- unique(as.numeric(x$varY[idx]))
    if (length(val) > 1L) {
        abort(glue(
            "varY(): the {length(idx)} rows of study ",
            "'{x$study[[idx[[1L]]]]}' disagree on varY ",
            "({str_flatten(val, ', ')})."
        ))
    }
    if (is.na(val)) NULL else val
})


# =============================================================================
# Coercion / converters
# =============================================================================

#' Combine GwasSumStats collections
#'
#' Row-bind two or more \code{\link{GwasSumStats}} collections into one -- e.g.
#' the per-LD-block pieces a block-parallel pipeline writes, back into the
#' single block-keyed collection \code{\link{assembleCtwasInputs}} requires.
#' Per-element metadata (\code{blockId}, \code{varY}, \code{nCase} / ...)
#' is carried through; a collection lacking an optional column is NA-padded.
#'
#' The three collection-level slots are merged rather than taken from the first
#' input: \code{genome} and the \code{\link{summaryStatsQc}} options must agree
#' across the inputs (a mismatch is an error, not a silent first-wins), each
#' row keeps its own QC audit, so \code{\link{qcDiagnostics}} keeps
#' addressing the right element, and every row keeps its own LD reference.
#' References that are views of one source are stored once: genotype panels
#' read from the same files union into one panel over the shared genotype
#' handle (block-parallel pipelines narrow each
#' piece's panel to its own block, so a multi-block collection would otherwise
#' hold one panel per block), and an \code{LdMatrix} or
#' \code{GenotypeMixture} that another input's reference covers gives way to
#' that one. Inputs on unrelated references keep both.
#'
#' @param ... Two or more \code{GwasSumStats} objects, or a single \code{list}
#'   of them.
#' @param ldReference Optional LD reference -- a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{GenotypeMixture} or an
#'   \code{LdMatrix} -- that every row of the combined collection uses
#'   instead of its own. Default \code{NULL}.
#' @return A single combined \code{GwasSumStats}.
#' @seealso \code{\link{combineQtlSumStats}},
#'   \code{\link{combineFineMappingResults}}, \code{\link{combineTwasWeights}}
#' @examples
#' data(gwasSumStatsS4Example)
#' combineGwasSumStats(gwasSumStatsS4Example)
#' @export
combineGwasSumStats <- function(..., ldReference = NULL) {
    parts <- .asCombineList(
        list(...),
        "GwasSumStats",
        "combineGwasSumStats"
    )
    if (length(parts) == 1L) {
        return(parts[[1L]])
    }
    .rbindSumStats(parts, ldReference, "combineGwasSumStats")
}
