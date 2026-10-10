# =============================================================================
# GwasSumStats S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study). Each
# row holds a per-study GRanges of GWAS summary statistics covering a
# single LD block; build a separate collection per block when sweeping
# the genome. Class-level slots ldSketch + genome + qcInfo apply
# uniformly across rows.
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

#' @title GWAS Summary-Statistic Collection
#' @description S4 collection of GWAS summary statistics keyed by the identity
#'   tuple \code{(study)}. Each element is that study's per-variant
#'   \code{GRanges} covering a single LD block, so build one collection per
#'   block when sweeping the genome.
#' @details Required column: \code{study}, unique across rows. The class-level
#'   slots inherited from \code{\linkS4class{SumStatsBase}} --
#'   \code{ldSketch}, \code{genome} and \code{qcInfo} -- apply uniformly to
#'   every row.
#' @seealso \code{\link{GwasSumStats}} for the constructor and
#'   \code{\linkS4class{QtlSumStats}} for the QTL counterpart.
#' @importFrom checkmate makeAssertCollection assertNames assertList
#' @export
setClass(
    "GwasSumStats",
    contains = "SumStatsBase",
    validity = function(object) {
        # The ldSketch slot's class union enforces its type.
        coll <- makeAssertCollection()
        assertNames(
            colnames(mcols(object)) %||% character(0),
            must.include = "study",
            what = "colnames",
            .var.name = "mcols",
            add = coll
        )
        coll$push(.sumStatsCheckGenome(object))
        assertList(object@qcInfo, .var.name = "qcInfo", add = coll)
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
        "GwasSumStats: {nrow(object)} studies, ",
        "genome build {GenomeInfoDb::genome(object)}\n",
        .trim = FALSE
    ))
    ld <- object@ldSketch
    ldSrc <- if (is.null(ld)) {
        "none (LD-free)"
    } else {
        .ldSketchLabel(ld)
    }
    cat(glue("  LD sketch: {ldSrc}\n", .trim = FALSE))
})


#' @title GWAS Summary Statistics Handling
#' @description Constructor, accessors, and converters for \code{GwasSumStats}
#'   (the post-refactor DFrame-subclass collection keyed by \code{study}).
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
#' @description Construct a \code{GwasSumStats} S4 DFrame-subclass collection
#'   from per-study tuple vectors and a list of \code{GRanges} entries (one per
#'   study), plus a single LD sketch handle and a single genome build that apply
#'   to the whole collection.
#'
#' Each \code{GRanges} entry must carry per-variant statistics in its mcols (at
#' minimum \code{SNP}, \code{A1}, \code{A2}, \code{Z}, \code{N}; optionally
#' \code{MAF}, \code{INFO}, \code{BETA}, \code{SE}, \code{P}).
#' @param studyName Character vector of study identifiers (must be unique).
#' @param entry A \code{SimpleList} or \code{list} of \code{GRanges}, one per
#'   study.
#' @param genome Single character string giving the genome build (e.g.,
#'   \code{"hg19"}, \code{"hg38"}). Uniform across the collection because all
#'   entries share the same LD sketch.
#' @param ldSketch A genotype panel (see \code{\link{readGenotypes}})
#'   carrying the LD reference.
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
#' @param ldBlocks Optional LD-block specification: an \code{LdBlocks}, a
#'   \code{GRanges}, a data.frame with \code{chrom}/\code{start}/\code{end}
#'   (plus an optional \code{blockId}), or a path to such a table. When
#'   supplied, each
#'   study's variants are split into one element per block rather than one per
#'   chromosome, and \code{blockId} takes the block's key (its \code{names},
#'   else a \code{blockId} metadata column, else its coordinates). cTWAS needs
#'   this granularity: its EM estimates parameters across blocks, so a
#'   per-chromosome split is too coarse. Variants overlapping no block are
#'   dropped with a warning, because a variant outside every block has no
#'   block-local LD to be fine-mapped against.
#' @param blockId Optional character vector of block keys, one per
#'   \code{entry}, for entries that are \strong{already} split by block. Use
#'   it to carry existing keys through a rebuild; without it a rebuild would
#'   re-derive them as seqnames and collapse distinct blocks onto one key.
#'   Mutually exclusive with \code{ldBlocks}.
#' @param extraCols Optional named list of additional per-study columns to
#'   attach to the collection's \code{mcols}. Named rather than variadic
#'   because the base class's own constructor, \code{GRangesList(...)}, uses
#'   \code{...} for its \emph{elements}: a bare \code{...} here would read
#'   as adding an entry rather than a metadata column.
#' @param qcInfo A \code{list} recording which QC steps ran. Empty \code{list()}
#'   on construction; populated by \code{summaryStatsQc()} with a per-step audit
#'   record. Fine-mapping / TWAS pipelines reject inputs where
#'   \code{length(qcInfo(x)) == 0}.
#' @return A \code{GwasSumStats} object.
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100 * 1:3, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(SNP = paste0("rs", 1:3),
#'   A1 = "A", A2 = "G", Z = rnorm(3), N = 100L)
#' GwasSumStats(studyName = "t1", entry = list(gr), genome = "hg38",
#'   ldSketch = panel)
#' @export
GwasSumStats <- function(
    studyName,
    entry,
    genome,
    ldSketch = NULL,
    varY = NA_real_,
    nCase = NULL,
    nControl = NULL,
    nSample = NULL,
    qcInfo = list(),
    ldBlocks = NULL,
    blockId = NULL,
    extraCols = list()
) {
    if (missing(studyName) || missing(entry) || missing(genome)) {
        abort("`study`, `entry`, and `genome` are all required.")
    }
    varY <- .gwasValidateArgs(studyName, entry, genome, varY)
    cols <- list(
        study = as.character(studyName),
        varY = varY
    ) |>
        .gwasAppendOptional(nCase, nControl, nSample, studyName) |>
        .gwasAppendExtras(extraCols)
    dfArgs <- c(cols, list(check.names = FALSE))
    # The per-study GRanges become the collection's ELEMENTS; everything else
    # is per-study metadata and goes in mcols. There is no `entry` column.
    # mcols are attached to the GRangesList BEFORE new(), because new()
    # validates during initialize() and the validity method needs the identity
    # columns to already be there.
    # A multi-seqname entry (e.g. a genome-wide GWAS) is split into one
    # element per block (or per chromosome when no block manifest is given),
    # with its metadata row replicated alongside. Splitting is unconditional:
    # a stored element always spans exactly one seqname.
    split <- .gwasSplitEntry(entry, ldBlocks, blockId, length(entry))
    # `blockId` is always present, so downstream code (cTWAS in particular) can
    # key regions without first asking how the collection was built.
    md <- cbind(
        exec(S4Vectors::DataFrame, !!!dfArgs)[split$fromIdx, , drop = FALSE],
        S4Vectors::DataFrame(blockId = split$blockId)
    )
    grl <- S4Vectors::`mcols<-`(
        GenomicRanges::GRangesList(split$entry),
        value = md
    )
    .sumStatsNewValidated("GwasSumStats", grl, ldSketch, genome, qcInfo)
}

# Build and validate a SumStatsBase subclass from its GRangesList plus the
# three collection-level slots GwasSumStats and QtlSumStats share. Shared
# rather than duplicated so the two constructors cannot drift in how they
# coerce those slots. Used by qtlSumStats.R too.
# @noRd
.sumStatsNewValidated <- function(Class, grl, ldSketch, genome, qcInfo) {
    # The build goes into seqinfo, which is where a GRangesList keeps it and
    # where every Bioconductor consumer reads it from. Assigned before new()
    # so validity sees the finished object.
    built <- GenomeInfoDb::`genome<-`(grl, value = as.character(genome))
    obj <- methods::new(
        Class,
        built,
        ldSketch = .asLdSketch(ldSketch),
        qcInfo = as.list(qcInfo)
    )
    methods::validObject(obj)
    obj
}

# Split by LD block when a manifest is supplied, else by seqname. Both return
# (entry, fromIdx); this adds the blockId the seqname path does not carry,
# which for that path is just the seqname each piece sits on.
# @noRd
.gwasSplitEntry <- function(entry, ldBlocks, blockId, n) {
    if (!is.null(ldBlocks) && !is.null(blockId)) {
        msg <- glue(
            "pass `ldBlocks` (derive the keys) or `blockId` (supply them), ",
            "not both."
        )
        abort(msg)
    }
    if (!is.null(ldBlocks)) {
        return(.rtlSplitByBlocks(entry, .asLdBlockRanges(ldBlocks)))
    }
    split <- .rtlSplitBySeqname(entry)
    # Supplied ids are indexed by fromIdx, exactly like the other metadata
    # columns, so an entry that still splits further replicates its id rather
    # than falling out of alignment. This is what lets a rebuild (QC) carry
    # block keys through instead of silently re-deriving them as seqnames.
    list_assign(
        split,
        blockId = if (!is.null(blockId)) {
            .gwasCheckBlockId(blockId, n)[split$fromIdx]
        } else {
            # unname(): the seqname splitter names its pieces, and those names
            # would otherwise ride into the mcols column and make it
            # inconsistent with the block path, which produces a bare
            # character vector.
            unname(map_chr(split$entry, .gwasElementSeqname))
        }
    )
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

# The one seqname an already-split element sits on. NA for an empty element,
# which carries no coordinate to name.
# @noRd
.gwasElementSeqname <- function(g) {
    if (length(g) == 0L) {
        return(NA_character_)
    }
    as.character(seqnames(g))[[1L]]
}

# Validate genome / entry / length consistency; returns the recycled varY.
# @noRd
.gwasValidateArgs <- function(studyName, entry, genome, varY) {
    if (length(genome) != 1L) {
        msg <- glue(
            "`genome` must be a single character string (one build per ",
            "collection, because all entries share the LD sketch)."
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
    c(cols, extras)
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
            "Pass `study = <name>` to select one. ",
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
#' @param ranges Optional \code{GRanges} restricting the returned variants to
#'   those it overlaps. \code{NULL} (default) returns the study's full set.
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
        ranges = NULL,
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
        .ssStitchElements(x, .gwasSelectStudy(x, studyName), ranges)
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
        "A GWAS is selected by `study` (and narrowed by `ranges`)."
    ))
}

# z / nSamples / maf / nSnps are provided once by SumStatsBase (AllClasses.R).

#' @title Coerce Summary Statistics to a Data Frame
#' @description Coerce one tuple of a \code{GwasSumStats} or
#'   \code{QtlSumStats} collection to a per-variant \code{data.frame} in the
#'   standardized layout \code{variant_id, chrom, pos, A1, A2, z, beta, se, N,
#'   maf} (optional columns omitted when absent on the entry). This is the
#'   table view of \code{\link{sumStats}}, which returns the same selection as
#'   a \code{GRanges}; coercion rather than a return-type argument is the
#'   Bioconductor idiom for choosing a representation.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param row.names,optional Accepted for compatibility with
#'   \code{as.data.frame} and ignored.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts (\code{QtlSumStats} only).
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits (\code{QtlSumStats} only).
#' @param require Character vector. Columns that must be present (derived if
#'   necessary) in the returned summary-statistics data frame.
#' @param derive Whether to derive missing standard columns from the available
#'   ones; \code{"zFromBetaSe"} recovers \code{z} from \code{beta}/\code{se}.
#' @param keepChrPrefix Logical. If \code{TRUE}, keep the \code{chr} prefix on
#'   chromosome names; otherwise strip it.
#' @param ... Unused.
#' @return A \code{data.frame}.
#' @examples
#' data(qtlSumStatsExample)
#' head(as.data.frame(qtlSumStatsExample))
#' @rdname sumStatsDataFrame
#' @aliases sumStatsDataFrame
#' @export
setMethod(
    "as.data.frame",
    "GwasSumStats",
    function(
        x,
        row.names = NULL,
        optional = FALSE,
        studyName = NULL,
        require = character(0),
        derive = c("none", "zFromBetaSe"),
        keepChrPrefix = TRUE,
        ...
    ) {
        derive <- arg_match(derive)
        gr <- sumStats(x, studyName = studyName)
        .entryToSumstatDf(
            gr,
            require = require,
            derive = derive,
            keepChrPrefix = keepChrPrefix,
            label = glue(
                "GwasSumStats[",
                "{if (is.null(studyName)) '<auto>' else studyName}]"
            )
        )
    }
)


#' @rdname varY
#' @export
setMethod("varY", "GwasSumStats", function(x, studyName = NULL) {
    idx <- .gwasSelectStudy(x, studyName)
    val <- x$varY[[idx]]
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
#' across the inputs (a mismatch is an error, not a silent first-wins), the
#' per-element \code{qcInfo$entryAudit} concatenates in element order so
#' \code{\link{qcDiagnostics}} keeps addressing the right element, and the
#' LD sketches union into one panel over the shared genotype handle. That last
#' one matters: block-parallel pipelines narrow each piece's panel to its own
#' block, so keeping only the first would leave a multi-block collection whose
#' LD reference covers one block.
#'
#' @param ... Two or more \code{GwasSumStats} objects, or a single \code{list}
#'   of them.
#' @param ldSketch Optional genotype panel (see \code{\link{readGenotypes}}) to
#'   attach to the combined collection, overriding the unioned one. Default
#'   \code{NULL}.
#' @return A single combined \code{GwasSumStats}.
#' @seealso \code{\link{combineQtlSumStats}},
#'   \code{\link{combineFineMappingResults}}, \code{\link{combineTwasWeights}}
#' @examples
#' data(gwasSumStatsS4Example)
#' combineGwasSumStats(gwasSumStatsS4Example)
#' @export
combineGwasSumStats <- function(..., ldSketch = NULL) {
    parts <- .asCombineList(
        list(...),
        "GwasSumStats",
        "combineGwasSumStats"
    )
    if (length(parts) == 1L) {
        return(parts[[1L]])
    }
    .rbindSumStats(parts, ldSketch, "combineGwasSumStats")
}
