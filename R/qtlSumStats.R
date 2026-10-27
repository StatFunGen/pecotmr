# =============================================================================
# QtlSumStats S4 class
# -----------------------------------------------------------------------------
# RangedTupleList collection keyed by the identity tuple (study, context,
# trait). Each row holds a per-tuple GRanges of summary statistics
# (variant_id + per-variant Z/N/MAF mcols). The genome build applies
# uniformly across rows; the ldReference slot holds each distinct LD
# reference once and the ldReference mcols column names the one every row
# uses. The qcOptions slot records the summaryStatsQc() options that ran, and
# the qcAudit column what QC did to each row.
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

#' @title QTL Summary-Statistic Collection
#' @description S4 collection of QTL summary statistics keyed by the identity
#'   tuple \code{(study, context, trait)}. Each element is that tuple's
#'   per-variant \code{GRanges} -- \code{variant_id} plus the per-variant
#'   Z / N / MAF mcols.
#' @details Required columns: \code{study}, \code{context} and \code{trait};
#'   the 3-tuple is unique. The constructor also adds each tuple's analysis
#'   window: \code{region}, a \code{GRanges} with one range per element, and
#'   \code{regionGiven}, \code{TRUE} where the window was given as a
#'   \code{region} and \code{FALSE} where it is the span of the tuple's
#'   variants when the collection was first built, before any QC. A QC pass
#'   carries both unchanged; the span of the variants left is
#'   \code{range(x)}. An optional \code{traitPos} gives the trait's own
#'   coordinates, from which a pipeline's \code{cisWindow} derives a cis
#'   window. The class-level slots inherited from
#'   \code{\linkS4class{SumStatsBase}} -- \code{genome} and
#'   \code{qcOptions} -- apply uniformly to every row, and \code{qcOptions}
#'   records the \code{\link{summaryStatsQc}} options that ran (what QC did
#'   to each row is \code{\link{qcInfo}}); \code{ldReference} holds
#'   each distinct LD reference once, and the \code{ldReference} column names
#'   the one each row uses.
#' @seealso \code{\link{QtlSumStats}} for the constructor and
#'   \code{\linkS4class{GwasSumStats}} for the GWAS counterpart.
#' @export
setClass(
    "QtlSumStats",
    contains = "SumStatsBase",
    validity = function(object) .validateQtlSumStats(object)
)

# ---- QtlSumStats validity helpers ------------------------------------------

# Collect all contract violations (empty vector = valid). The per-entry checks
# run only once the basic slot/column checks pass (they assume those columns).
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames
.validateQtlSumStats <- function(object) {
    coll <- makeAssertCollection()
    # `names(object)` is element names on a RangedTupleList, so the required
    # metadata columns are read from mcols directly.
    assertNames(
        colnames(mcols(object)) %||% character(0),
        must.include = c("study", "context", "trait"),
        what = "colnames",
        .var.name = "mcols",
        add = coll
    )
    coll$push(.qssCheckGenome(object))
    coll$push(.validateTraitPosColumn(object))
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    coll$push(.qssCheckEntries(object))
    coll$getMessages()
}

# The genome build, read from seqinfo (there is no genome slot).
# @noRd
.qssCheckGenome <- function(object) {
    .sumStatsCheckGenome(object)
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
#'   \code{RangedTupleList} collection keyed by \code{(study, context,
#'   trait)}.
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
#' @description Construct a \code{QtlSumStats} S4 collection from per-tuple
#'   vectors and a list of \code{GRanges} entries (one per tuple), plus the LD
#'   reference(s) the rows were harmonized against and a single genome build
#'   that applies to the whole collection. Each \code{GRanges} entry must
#'   carry per-variant statistics in its mcols (\code{SNP}, \code{A1},
#'   \code{A2}, \code{Z}, \code{N}; plus optional \code{MAF}, \code{INFO},
#'   \code{BETA}, \code{SE}, \code{P}).
#' @param studyName Character vector of study identifiers (per tuple).
#' @param context Character vector of context labels (per tuple).
#' @param trait Character vector of trait identifiers (per tuple).
#' @param entry A list / \code{SimpleList} of \code{GRanges}, one per tuple.
#'   Same length as \code{study}, \code{context}, and \code{trait}.
#' @param genome Single character string giving the genome build (e.g.,
#'   \code{"hg19"}, \code{"hg38"}). One build per collection.
#' @param ldReference The LD reference: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}, used for every tuple, or a list with one per
#'   tuple when they sit on different references. Each distinct reference is
#'   stored once, and the \code{ldReference} column names each row's.
#' @param varY Optional numeric vector of per-tuple phenotype variances
#'   (\code{NA_real_} entries allowed).
#' @param nSample Optional per-tuple total sample size (numeric; default
#'   \code{NULL}). Attached only when supplied (length 1 or length(study)). Used
#'   as the study-level fallback for the per-variant \code{N} when a tuple has
#'   no per-variant \code{N} column. Named \code{nSample} to avoid clashing with
#'   \code{nSamples()} (the LD-panel sample size). Unlike GWAS, QTL
#'   collections carry no case/control counts (molecular traits are
#'   quantitative), so only this total-N fallback is exposed.
#' @param extraCols Optional named list (or \code{DataFrame} /
#'   \code{data.frame}) of additional per-tuple columns to
#'   attach to the collection's \code{mcols}. Named rather than variadic
#'   because the base class's own constructor, \code{GRangesList(...)}, uses
#'   \code{...} for its \emph{elements}: a bare \code{...} here would read
#'   as adding an entry rather than a metadata column.
#' @param qcOptions A \code{list} recording the QC that ran. Empty
#'   \code{list()} on construction; \code{summaryStatsQc()} fills it.
#'   Fine-mapping / TWAS pipelines reject inputs where
#'   \code{length(qcOptions(x)) == 0}.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column. A pipeline given a
#'   \code{cisWindow} analyses each tuple over \code{traitPos} \eqn{\pm}
#'   \code{cisWindow}.
#' @param region Optional analysis window of each tuple: a \code{GRanges},
#'   or the path of a BED file, read as BED with a 0-based start. One window
#'   for every tuple, or one per tuple. A tuple's variants outside its window
#'   are dropped with a warning, and its \code{region} is the window,
#'   recorded as given in \code{regionGiven}. Without it, a tuple's
#'   \code{region} is the span of its variants as built, before any QC.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{QtlSumStats} object.
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100 * 1:3, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(SNP = paste0("rs", 1:3),
#'   A1 = "A", A2 = "G", Z = rnorm(3), N = 100L)
#' QtlSumStats(studyName = "s1", context = "brain", trait = "g1", entry =
#'   list(gr), genome = "hg38", ldReference = panel)
#' @export
QtlSumStats <- function(
    studyName,
    context,
    trait,
    entry,
    genome,
    ldReference = NULL,
    varY = NA_real_,
    nSample = NULL,
    qcOptions = list(),
    traitPos = NULL,
    region = NULL,
    extraCols = list(),
    metadata = list()
) {
    if (
        missing(studyName) ||
            missing(context) ||
            missing(trait) ||
            missing(entry) ||
            missing(genome)
    ) {
        msg <- glue(
            "`studyName`, `context`, `trait`, `entry`, and `genome` ",
            "are all required."
        )
        abort(msg)
    }
    n <- length(studyName)
    varY <- .qssValidateArgs(context, trait, entry, genome, varY, n)
    windows <- .qssWindows(region, n)
    entry <- .qssWithinWindows(.ucscEntries(entry), windows) |>
        .qssAddTraitDistances(traitPos, n)
    cols <- .qssBaseCols(studyName, context, trait, varY) |>
        .qssAppendNSample(nSample, n) |>
        .appendTraitPosCol(traitPos, n) |>
        .qssAppendExtras(extraCols)
    dict <- .ldRefDict(.ldRefPerTuple(ldReference, n))
    dfArgs <- c(cols, list(ldReference = dict$keys, check.names = FALSE))
    split <- .rtlSplitBySeqname(entry)
    pieces <- .rtlBindElements(split$entry)
    grl <- S4Vectors::`mcols<-`(
        pieces,
        value = .qssPieceMetadata(dfArgs, split, pieces, windows)
    )
    .withMetadata(
        .sumStatsNewValidated("QtlSumStats", grl, dict$refs, genome, qcOptions),
        metadata
    )
}

# Each piece's metadata row, its tuple's, with the tuple's analysis window:
# the `region` it was given, else the span of its variants as built, before
# any QC. `regionGiven` records which.
#
# The per-tuple GRanges become the collection's ELEMENTS; the tuple keys and
# per-tuple scalars go in mcols. There is no `entry` column. mcols are
# attached to the GRangesList BEFORE new(), because new() validates during
# initialize() and the validity method needs the identity columns to already
# be there. A multi-seqname entry (e.g. a trans tuple) is split into one
# element per chromosome, with its metadata row replicated alongside.
# @noRd
.qssPieceMetadata <- function(dfArgs, split, pieces, windows) {
    n <- length(split$entry)
    given <- if (is.null(windows)) {
        .placeholderRanges(n)
    } else {
        windows[split$fromIdx]
    }
    cbind(
        exec(S4Vectors::DataFrame, !!!dfArgs)[split$fromIdx, , drop = FALSE],
        exec(
            S4Vectors::DataFrame,
            !!!.regionColumns(
                given,
                rep(!is.null(windows), n),
                .rowSpanRanges(pieces)
            )
        )
    )
}

# The window each tuple was given, as bare ranges, one per tuple; NULL when
# none was. A single window stands for every tuple.
# @noRd
.qssWindows <- function(region, n) {
    if (is.null(region)) {
        return(NULL)
    }
    windows <- .asRegionRanges(region)
    if (!is_in(length(windows), c(1L, n))) {
        abort(glue(
            "`region` must give one window, or one per tuple ",
            "(got {length(windows)} for {n} tuples)."
        ))
    }
    .bareRanges(rep(windows, length.out = n))
}

# Each tuple's variants inside its window. Those outside it were not asked
# for, so they are dropped, with a warning saying how many.
# @noRd
.qssWithinWindows <- function(entry, windows) {
    if (is.null(windows)) {
        return(entry)
    }
    kept <- map(
        seq_along(entry),
        .qssWithinWindowAt,
        entry = as.list(entry),
        windows = windows
    )
    dropped <- sum(lengths(entry)) - sum(lengths(kept))
    if (dropped > 0L) {
        warn(glue(
            "{dropped} variant(s) lie outside their tuple's `region` and ",
            "were dropped."
        ))
    }
    kept
}

# @noRd
.qssWithinWindowAt <- function(i, entry, windows) {
    .withinWindow(entry[[i]], windows[i])
}

# Narrow the tuples a pipeline selected -- `rows`, picked by its own
# `contexts` / `traitId` rules -- to their traits' cis windows: traitPos +/-
# cisWindow, within the window each tuple already had. That window becomes
# the tuple's `region`, recorded as given, since the caller asked for it. A
# selected tuple with no trait position has no cis window, which is an error.
# NULL `cisWindow` leaves `x` as it is. The one step behind every pipeline
# that takes a QtlSumStats and a `cisWindow`; `caller` names it in messages.
# @noRd
#' @importFrom checkmate assertNumber
.qssApplyCisWindow <- function(x, cisWindow, caller, rows) {
    if (is.null(cisWindow)) {
        return(x)
    }
    assertNumber(cisWindow, lower = 0, finite = TRUE, .var.name = "cisWindow")
    if (length(rows) == 0L) {
        return(x)
    }
    regions <- .rowRegions(x)
    windows <- .clipWindows(
        .cisWindowsAround(.qssTraitPositions(x, rows, caller), cisWindow),
        regions$region[rows]
    )
    elements <- replace(
        as.list(x),
        rows,
        map(
            seq_along(rows),
            .qssWithinWindowAt,
            entry = as.list(x)[rows],
            windows = windows
        )
    )
    region <- regions$region
    region[rows] <- windows
    .withRegionColumns(
        .rtlRebuild(x, elements, seq_len(nrow(x))),
        list(
            region = region,
            regionGiven = replace(regions$regionGiven, rows, TRUE)
        )
    )
}

# A QtlSumStats run's selection, parsed once (parseContexts() /
# parseTraitIds()), and its data with `cisWindow` laid on the tuples that
# selection keeps -- so the windows cover exactly the tuples the run fits.
# @noRd
.qssSelection <- function(data, contexts, traitId, cisWindow, caller) {
    contexts <- parseContexts(contexts, data)
    traitId <- parseTraitIds(traitId, data)
    rows <- which(.sumStatsRowsSelected(data, contexts, traitId))
    list(
        data = .qssApplyCisWindow(data, cisWindow, caller, rows),
        contexts = contexts,
        traitId = traitId
    )
}

# Each trait position widened by `cisWindow` on both sides, kept on the
# chromosome.
# @noRd
.cisWindowsAround <- function(tp, cisWindow) {
    w <- as.integer(cisWindow)
    GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(tp)),
        IRanges::IRanges(start = pmax(1L, start(tp) - w), end = end(tp) + w)
    )
}

# The trait position of each of `rows`; an error naming the tuples that have
# none, since a cis window is laid around the trait.
# @noRd
.qssTraitPositions <- function(x, rows, caller) {
    tp <- if (is_in("traitPos", .tupleColumnNames(x))) {
        .tupleColumn(x, "traitPos")[rows]
    }
    missing <- if (is.null(tp)) {
        rows
    } else {
        rows[as.character(GenomicRanges::seqnames(tp)) == "chrUn"]
    }
    if (length(missing) > 0L) {
        labels <- str_c(
            as.character(x$study)[missing],
            as.character(x$context)[missing],
            as.character(x$trait)[missing],
            sep = "/"
        )
        abort(c(
            glue(
                "{caller}: `cisWindow` is laid around each tuple's trait ",
                "position, which {length(missing)} tuple(s) do not have: ",
                "{str_flatten(utils::head(labels, 5L), ', ')}."
            ),
            i = str_c(
                "Give one with `traitPos` (or a manifest's traitStart / ",
                "traitEnd columns), or drop `cisWindow` to analyse each ",
                "tuple over its `region`."
            )
        ))
    }
    tp
}

# Each window clipped to its row's `regions` entry. A window that misses it,
# on another chromosome or past its ends, keeps nothing, so it becomes a
# zero-width placeholder at its start.
# @noRd
.clipWindows <- function(windows, regions) {
    same <- as.character(GenomicRanges::seqnames(windows)) ==
        as.character(GenomicRanges::seqnames(regions))
    s <- pmax(start(windows), start(regions))
    e <- pmin(end(windows), end(regions))
    hit <- same & e >= s
    GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(windows)),
        IRanges::IRanges(
            start = if_else(hit, s, start(windows)),
            end = if_else(hit, e, start(windows) - 1L)
        )
    )
}

# Validate genome / entry / length consistency and recycle varY. Returns the
# (possibly recycled) varY.
# @noRd
.qssValidateArgs <- function(context, trait, entry, genome, varY, n) {
    if (length(genome) != 1L) {
        msg <- glue(
            "`genome` must be a single character string (one build per ",
            "collection)."
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
            "`studyName`, `context`, `trait`, and `entry` must all have the ",
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
.qssBaseCols <- function(studyName, context, trait, varY) {
    list(
        study = as.character(studyName),
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
    c(
        cols,
        list(nSample = as.numeric(.qssRecycleTo(nSample, n, "nSample")))
    )
}

# Append any user-supplied extra columns (from `...`).
# @noRd
.qssAppendExtras <- function(cols, extras) {
    c(cols, .extraColsList(extras))
}

# `extraCols` as a named list of columns: a list as it is, a DataFrame or
# data.frame column by column.
# @noRd
.extraColsList <- function(extras) {
    if (methods::is(extras, "DataFrame") || is.data.frame(extras)) {
        return(as.list(extras))
    }
    extras
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
#' @importFrom checkmate assertVector
.qtlSumStatsSelectRow <- function(x, studyName, context, trait) {
    if (nrow(x) == 0L) {
        abort("QtlSumStats has no rows.")
    }
    anyUnset <- missing(studyName) ||
        is.null(studyName) ||
        missing(context) ||
        is.null(context) ||
        missing(trait) ||
        is.null(trait)
    if (anyUnset) {
        if (.qssSingleTuple(x)) {
            return(seq_len(nrow(x)))
        }
        msg <- glue(
            "This QtlSumStats has {nrow(x)} entries. Pass `studyName`, ",
            "`context`, and `trait` to select one."
        )
        abort(msg)
    }
    assertVector(studyName, len = 1L)
    assertVector(context, len = 1L)
    assertVector(trait, len = 1L)
    .qssMatchTuple(x, studyName, context, trait)
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
.qssMatchTuple <- function(x, studyName, context, trait) {
    idx <- .matchTupleRows(
        x,
        list(study = studyName, context = context, trait = trait)
    )
    if (length(idx) == 0L) {
        msg <- glue(
            "No entry for (study='{studyName}', context='{context}', ",
            "trait='{trait}')."
        )
        abort(msg)
    }
    # More than one match is expected, not a violation: uniqueness is keyed on
    # (tuple, range), so a tuple split across chromosomes owns several
    # elements. The caller stitches them.
    idx
}

#' @rdname sumStats
#' @param annotateSignificance Optional correction-method name
#'   (\code{"permutation"} / \code{"bonferroni_original"} /
#'   \code{"bonferroni_filtered"} / \code{"qvalue"}). When set on a QtlSumStats
#'   enriched by \code{\link{qtlAssociationPostprocess}}, a logical
#'   \code{significant} mcol for that method is added to the returned entry (the
#'   significance is derived on the fly, not stored). Take a table of the
#'   whole entry, every mcol included, with \code{as.data.frame()}.
#' @export
setMethod(
    "sumStats",
    signature(x = "QtlSumStats"),
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        annotateSignificance = NULL
    ) {
        idx <- .qtlSumStatsSelectRow(x, studyName, context, trait)
        stitched <- .ssStitchElements(x, idx)
        gr <- if (is.null(annotateSignificance)) {
            stitched
        } else {
            m <- arg_match(
                annotateSignificance,
                c(
                    "permutation",
                    "bonferroni_original",
                    "bonferroni_filtered",
                    "qvalue"
                )
            )
            S4Vectors::`mcols<-`(
                stitched,
                value = `[[<-`(
                    S4Vectors::mcols(stitched, use.names = FALSE),
                    "significant",
                    value = .qapSignificanceMask(x, m)[[idx]]
                )
            )
        }
        gr
    }
)

# z / nSamples / maf / nSnps are provided once by SumStatsBase (AllClasses.R);
# they only delegate to sumStats().

#' @rdname varY
#' @export
setMethod(
    "varY",
    "QtlSumStats",
    function(x, studyName = NULL, context = NULL, trait = NULL) {
        idx <- .qtlSumStatsSelectRow(x, studyName, context, trait)
        val <- x$varY[[idx]]
        if (is.na(val)) NULL else val
    }
)

#' @rdname contexts
#' @export
setMethod("contexts", "QtlSumStats", function(x) {
    unique(as.character(x$context))
})

#' @rdname traitNames
#' @export
setMethod("traitNames", "QtlSumStats", function(x) {
    unique(as.character(x$trait))
})

# Trait position provenance. Optional for a QtlSumStats -- it cannot be inferred
# from summary statistics, so it is only present when the caller supplied it.
# Returns the traitPos GRanges (whole column, or the rows matching `traitId`),
# or a scalar NA when no trait position was supplied.
#' @rdname traitPosition
#' @export
setMethod("traitPosition", "QtlSumStats", function(x, trait = NULL) {
    tp <- .getTraitPosColumn(x)
    if (is.null(trait) || length(tp) == 0L) {
        return(tp)
    }
    idx <- which(is_in(as.character(x$trait), as.character(trait)))
    # The trait position is per-TUPLE, but a tuple split across chromosomes
    # owns several elements and the column is replicated across them, so the
    # matches collapse back to the one position. No match is an empty GRanges.
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
        "genome build {.showGenomeBuild(object)}\n",
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
    ldSrc <- .ldRefsLabel(ldReference(object), "none (LD-free)")
    cat(glue("  LD reference: {ldSrc}\n", .trim = FALSE))
})

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Append tss_distance / tes_distance mcols to entry `i` from traitPos[i] (no-op
# when the entry is empty or the trait position is unset). Existing distance
# mcols are preserved.
# @noRd
.qssAppendTraitDist <- function(i, entry, traitPos) {
    gr <- entry[[i]]
    if (length(gr) == 0L) {
        return(gr)
    }
    # `traitPos[i]` is a length-1 GRanges (an out-of-range `i` errors here
    # rather than yielding an empty one), and a GRanges cannot carry an NA
    # start, so the position is always real.
    tp <- traitPos[i]
    tssPos <- GenomicRanges::start(tp)
    tesPos <- GenomicRanges::end(tp)
    pos <- GenomicRanges::start(gr)
    # No zero-column swap: mcols() on a GRanges is never NULL, and a variant
    # set with no metadata yields a zero-COLUMN DataFrame that already carries
    # one row per variant. Replacing it with DataFrame(row.names = NULL) threw
    # that row count away, so the assignment below died with "n elements in
    # value to replace 0 elements" for any entry without mcols.
    mc <- S4Vectors::mcols(gr)
    added <- compact(list(
        tss_distance = if (is.null(mc[["tss_distance"]])) pos - tssPos,
        tes_distance = if (is.null(mc[["tes_distance"]])) pos - tesPos
    ))
    if (length(added) == 0L) {
        return(gr)
    }
    S4Vectors::`mcols<-`(
        gr,
        value = cbind(mc, exec(S4Vectors::DataFrame, !!!added))
    )
}

#' Combine QtlSumStats collections
#'
#' Row-bind two or more \code{\link{QtlSumStats}} collections into one -- e.g.
#' per-region or per-context association results written by separate pipeline
#' steps. Per-element metadata columns are carried through; a collection
#' lacking an optional column is NA-padded.
#'
#' The three collection-level slots are merged rather than taken from the first
#' input, on the same rules as \code{\link{combineGwasSumStats}}: \code{genome}
#' and the \code{\link{summaryStatsQc}} options must agree, each row keeps
#' its own QC audit, and the rows keep their own LD references, stored the
#' same way.
#'
#' @param ... Two or more \code{QtlSumStats} objects, or a single \code{list}
#'   of them.
#' @param ldReference Optional LD reference -- a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{GenotypeMixture} or an
#'   \code{LdMatrix} -- that every row of the combined collection uses
#'   instead of its own. Default \code{NULL}.
#' @return A single combined \code{QtlSumStats}.
#' @seealso \code{\link{combineGwasSumStats}},
#'   \code{\link{combineFineMappingResults}}, \code{\link{combineTwasWeights}}
#' @examples
#' data(qtlSumStatsExample)
#' combineQtlSumStats(qtlSumStatsExample)
#' @export
combineQtlSumStats <- function(..., ldReference = NULL) {
    parts <- .asCombineList(
        list(...),
        "QtlSumStats",
        "combineQtlSumStats"
    )
    if (length(parts) == 1L) {
        return(parts[[1L]])
    }
    .rbindSumStats(parts, ldReference, "combineQtlSumStats")
}
