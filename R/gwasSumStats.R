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

setClass(
    "GwasSumStats",
    contains = "SumStatsBase",
    validity = function(object) {
        # The ldSketch slot's class union enforces its type.
        errors <- character()
        required <- "study"
        missingCols <- setdiff(required, colnames(mcols(object)))
        if (length(missingCols) > 0L) {
            errors <- c(
                errors,
                str_c("missing columns: ", str_flatten(missingCols, ", "))
            )
        }
        if (length(object@genome) != 1L || str_length(object@genome) == 0L) {
            errors <- c(
                errors,
                "'genome' slot must be a single non-empty character string"
            )
        }
        if (!is.list(object@qcInfo)) {
            errors <- c(errors, "'qcInfo' slot must be a list")
        }
        if (length(errors) == 0L) {
            # The elements ARE GRanges by construction now -- the container is
            # a GRangesList -- so the old per-element type and length checks
            # are gone. The one-seqname/one-strand invariant is enforced by
            # RangedTupleList's own validity.
            # Keyed on (study, range), not study alone: a study split across
            # chromosomes contributes one element per seqname, so the study
            # label legitimately repeats.
            if (.ssHasDuplicateKeys(object, "study")) {
                errors <- c(
                    errors,
                    "(study, range) must be unique"
                )
            }
        }
        if (length(errors) == 0L) TRUE else errors
    }
)


#' @rdname show-methods
setMethod("show", "GwasSumStats", function(object) {
    cat(glue(
        "GwasSumStats: {nrow(object)} studies, ",
        "genome build {object@genome}\n",
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
.recyclePerStudy <- function(v, nm, study) {
    if (length(v) == 1L && length(study) > 1L) {
        v <- rep(v, length(study))
    }
    if (length(v) != length(study)) {
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
#' @param study Character vector of study identifiers (must be unique).
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
#'   \code{nSample} to avoid clashing with \code{getNSamples()} (the LD-panel
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
#' @param ... Additional per-study columns to attach to the collection.
#' @param qcInfo A \code{list} recording which QC steps ran. Empty \code{list()}
#'   on construction; populated by \code{summaryStatsQc()} with a per-step audit
#'   record. Fine-mapping / TWAS pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0}.
#' @return A \code{GwasSumStats} object.
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100 * 1:3, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(SNP = paste0("rs", 1:3),
#'   A1 = "A", A2 = "G", Z = rnorm(3), N = 100L)
#' GwasSumStats(study = "t1", entry = list(gr), genome = "hg38",
#'   ldSketch = panel)
#' @export
GwasSumStats <- function(
    study,
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
    ...
) {
    if (missing(study) || missing(entry) || missing(genome)) {
        abort("`study`, `entry`, and `genome` are all required.")
    }
    varY <- .gwasValidateArgs(study, entry, genome, varY)
    cols <- list(
        study = as.character(study),
        varY = varY
    )
    cols <- .gwasAppendOptional(cols, nCase, nControl, nSample, study)
    cols <- .gwasAppendExtras(cols, list(...))
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
    grl <- GenomicRanges::GRangesList(split$entry)
    md <- exec(S4Vectors::DataFrame, !!!dfArgs)
    md <- md[split$fromIdx, , drop = FALSE]
    # `blockId` is always present, so downstream code (cTWAS in particular) can
    # key regions without first asking how the collection was built.
    md$blockId <- split$blockId
    mcols(grl) <- md
    .sumStatsNewValidated("GwasSumStats", grl, ldSketch, genome, qcInfo)
}

# Build and validate a SumStatsBase subclass from its GRangesList plus the
# three collection-level slots GwasSumStats and QtlSumStats share. Shared
# rather than duplicated so the two constructors cannot drift in how they
# coerce those slots. Used by qtlSumStats.R too.
# @noRd
.sumStatsNewValidated <- function(Class, grl, ldSketch, genome, qcInfo) {
    obj <- methods::new(
        Class,
        grl,
        ldSketch = .asLdSketch(ldSketch),
        genome = as.character(genome),
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
    split$blockId <- if (!is.null(blockId)) {
        .gwasCheckBlockId(blockId, n)[split$fromIdx]
    } else {
        # unname(): the seqname splitter names its pieces, and those names
        # would otherwise ride into the mcols column and make it inconsistent
        # with the block path, which produces a bare character vector.
        unname(map_chr(split$entry, .gwasElementSeqname))
    }
    split
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
.gwasValidateArgs <- function(study, entry, genome, varY) {
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
    if (length(entry) != length(study)) {
        msg <- glue(
            "length(entry) ({length(entry)}) must equal ",
            "length(study) ({length(study)})."
        )
        abort(msg)
    }
    .recyclePerStudy(varY, "varY", study)
}

# Attach the OPTIONAL per-study nCase / nControl / nSample columns (each only
# when supplied; NA for the non-case/control studies in a mixed collection).
# @noRd
.gwasAppendOptional <- function(cols, nCase, nControl, nSample, study) {
    if (!is.null(nCase)) {
        cols$nCase <- .recyclePerStudy(nCase, "nCase", study)
    }
    if (!is.null(nControl)) {
        cols$nControl <- .recyclePerStudy(nControl, "nControl", study)
    }
    if (!is.null(nSample)) {
        cols$nSample <- .recyclePerStudy(nSample, "nSample", study)
    }
    cols
}

# Append any user-supplied extra columns (from `...`).
# @noRd
.gwasAppendExtras <- function(cols, extras) {
    for (nm in names(extras)) {
        cols[[nm]] <- extras[[nm]]
    }
    cols
}


# =============================================================================
# Accessors for the new GwasSumStats collection
# =============================================================================

# Internal: resolve a study selection to a single row index. Errors when
# `study` is missing on a multi-study collection.
# Element indices for one study. Returns a VECTOR, not a scalar: the seqname
# split means one study can own several elements (one per chromosome), and
# getSumStats() stitches them back into the single GRanges callers expect.
# @noRd
.gwasSelectStudy <- function(x, study) {
    if (nrow(x) == 0L) {
        abort("GwasSumStats has no rows.")
    }
    studies <- as.character(x$study)
    if (missing(study) || is.null(study)) {
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
    idx <- which(studies == as.character(study))
    if (length(idx) == 0L) {
        msg <- glue(
            "Unknown study: '{study}'. ",
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
#' @param study Character (length 1) study identifier. Optional when the
#'   collection has a single row.
#' @param ranges Optional \code{GRanges} restricting the returned variants to
#'   those it overlaps. \code{NULL} (default) returns the study's full set.
#' @param ... Additional arguments (currently unused).
#' @return A \code{GRanges} object.
#' @export
setMethod(
    "getSumStats",
    signature(x = "GwasSumStats"),
    function(x, study = NULL, ranges = NULL, ...) {
        .ssStitchElements(x, .gwasSelectStudy(x, study), ranges)
    }
)

# getZ / getN / getMaf / nSnps are provided once by SumStatsBase (AllClasses.R).

#' @rdname getSumstatDf
#' @export
setMethod(
    "getSumstatDf",
    "GwasSumStats",
    function(
        x,
        study = NULL,
        require = character(0),
        derive = c("none", "zFromBetaSe"),
        keepChrPrefix = TRUE
    ) {
        derive <- arg_match(derive)
        gr <- getSumStats(x, study = study)
        .entryToSumstatDf(
            gr,
            require = require,
            derive = derive,
            keepChrPrefix = keepChrPrefix,
            label = glue(
                "GwasSumStats[{if (is.null(study)) '<auto>' else study}]"
            )
        )
    }
)


#' @rdname getVarY
#' @export
setMethod("getVarY", "GwasSumStats", function(x, study = NULL) {
    idx <- .gwasSelectStudy(x, study)
    val <- x$varY[[idx]]
    if (is.na(val)) NULL else val
})

# =============================================================================
# Coercion / converters
# =============================================================================

#' @title Convert GwasSumStats to data.frame
#' @description Extracts the per-variant statistics for one study (selected by
#'   \code{study}) into a plain data.frame with columns SNP, CHR, BP, A1, A2, Z,
#'   N (and any optional columns such as MAF, BETA, SE, P).
#' @param x A \code{GwasSumStats} object.
#' @param row.names Ignored (present for S3 generic compatibility).
#' @param optional Ignored.
#' @param study Character (length 1) study identifier. Optional when the
#'   collection has a single row.
#' @param ... Ignored.
#' @return A data.frame.
#' @method as.data.frame GwasSumStats
#' @export
as.data.frame.GwasSumStats <- function(
    x,
    row.names = NULL,
    optional = FALSE,
    study = NULL,
    ...
) {
    gr <- getSumStats(x, study = study)
    mc <- as.data.frame(mcols(gr))
    mc$CHR <- as.character(seqnames(gr))
    mc$BP <- start(gr)
    firstCols <- c("SNP", "CHR", "BP")
    restCols <- setdiff(names(mc), firstCols)
    select(mc, all_of(c(firstCols, restCols)))
}

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
#' \code{\link{getQcDiagnostics}} keeps addressing the right element, and the
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
