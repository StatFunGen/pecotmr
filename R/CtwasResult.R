# =============================================================================
# CtwasResult S4 class
# -----------------------------------------------------------------------------
# RangedTupleList collection of cTWAS runs. A run is keyed by the identity
# tuple (gwasStudy, study, context, method) and, as a genome-wide GWAS fit is
# in a GwasFineMappingResult, split into one element per cTWAS region: the
# element holds that region's gene (and, with keepSnps, SNP) posteriors as
# ranges, and mcols carry the run identity, the region (`blockId`, cTWAS's
# region_id), the region's per-effect susie alpha table (`susieAlpha`), its
# region metadata (`regionInfo`) and the run's jointly-estimated group priors
# (`groupPriors`, repeated on each of the run's rows, as a GWAS fit's
# settings repeat on its blocks).
#
# Unlike the QTL family, `trait` is NOT a key -- cTWAS is multi-gene per run,
# so genes live inside the elements. Optional joint columns (jointStudies,
# jointContexts) tag rows born from a multi-study / multi-context cTWAS run and
# participate in the uniqueness key, exactly as in the TwasWeights /
# FineMappingResultBase families.
# =============================================================================

#' @include AllGenerics.R tupleSelectors.R RangedTupleList.R CtwasResultEntry.R
NULL

#' @title cTWAS Result Collection
#' @description S4 collection of cTWAS runs, each keyed by the identity tuple
#'   \code{(gwasStudy, study, context, method)} and split into one element per
#'   cTWAS region, the way a \code{\linkS4class{GwasFineMappingResult}} holds
#'   a genome-wide fit as one element per LD block. An element is a
#'   \code{GRanges} of that region's gene (and, with \code{keepSnps}, SNP)
#'   posteriors; \code{mcols()} carry the run identity, the region
#'   (\code{blockId}), the region's per-effect susie alpha table
#'   (\code{susieAlpha}), its region metadata (\code{regionInfo}) and the
#'   run's jointly-estimated group priors (\code{groupPriors}), which repeat
#'   on each of the run's rows.
#' @details Unlike the QTL family, \code{trait} is not part of the key: a cTWAS
#'   run is multi-gene, so genes live inside the elements. The optional
#'   \code{jointStudies} / \code{jointContexts} columns tag rows born from a
#'   multi-study or multi-context run and participate in the uniqueness key,
#'   exactly as in the \code{\linkS4class{TwasWeights}} and fine-mapping
#'   families. \code{\link{ctwasPosteriors}} and \code{\link{susieAlpha}}
#'   stack every region back into one \code{GRanges}.
#' @seealso \code{\link{CtwasResult}} for the constructor and
#'   \code{\link{ctwasPipeline}} for the pipeline that builds one.
#' @export
setClass(
    "CtwasResult",
    contains = "RangedTupleList",
    validity = function(object) .validateCtwasResult(object)
)

# The run identity every row carries, and the per-region payload columns.
# @noRd
.CTWAS_ID_COLS <- c("gwasStudy", "study", "context", "method")

# @noRd
.CTWAS_PAYLOAD_COLS <- c("susieAlpha", "regionInfo", "groupPriors")

# ---- CtwasResult validity helpers ------------------------------------------

# Collect all contract violations (TRUE when there are none). The identity and
# uniqueness checks run only once the required columns are present.
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames
.validateCtwasResult <- function(object) {
    cols <- colnames(mcols(object, use.names = FALSE))
    coll <- makeAssertCollection()
    assertNames(
        cols,
        must.include = c(.CTWAS_ID_COLS, "blockId", .CTWAS_PAYLOAD_COLS),
        what = "colnames",
        .var.name = "columns",
        add = coll
    )
    # The checks below read those columns; running them on an object missing
    # them reports the consequence rather than the cause.
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    jointCols <- intersect(c("jointStudies", "jointContexts"), cols)
    problems <- c(
        .ctwasResCheckJointCols(object, jointCols),
        .ctwasResCheckTupleUniqueness(object, jointCols)
    )
    if (length(problems) == 0L) TRUE else problems
}

# @noRd
.ctwasResCheckJointCols <- function(object, jointCols) {
    list_c(compact(map(jointCols, .ctwasResJointColError, object = object)))
}

# @noRd
.ctwasResJointColError <- function(jc, object) {
    value <- .ctwasResColOf(jc, object)
    if (!is.character(value)) {
        return(glue(
            "'{jc}' column must be character (got {class(value)[[1L]]})"
        ))
    }
    NULL
}

# One row per (run, region): the region joins the run identity in the key.
# @noRd
.ctwasResCheckTupleUniqueness <- function(object, jointCols) {
    keyCols <- c(.CTWAS_ID_COLS, jointCols, "blockId")
    keyTbl <- as_tibble(set_names(
        map(keyCols, .ctwasResColOf, object = object),
        keyCols
    ))
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        return(str_c(
            "(gwasStudy, study, context, method[, joint*], region) ",
            "tuple uniqueness violated"
        ))
    }
    NULL
}

# @noRd
.ctwasResColOf <- function(cn, object) {
    mcols(object, use.names = FALSE)[[cn]]
}

# Each joint-provenance column must be one value per run.
# @noRd
.ctwasCheckJointColLength <- function(nm, val, n) {
    if (length(val) == n) {
        return(invisible(NULL))
    }
    msg <- glue("`{nm}` must have the same length as `gwasStudy`.")
    abort(msg)
}

#' @title Create a CtwasResult Collection
#' @description Construct a \code{CtwasResult} collection of cTWAS runs from
#'   per-run identity vectors and a list of \code{\link{CtwasResultEntry}}
#'   payloads (one per run). Each run is split into one element per cTWAS
#'   region -- by the payloads' \code{region_id}, or by chromosome when they
#'   carry none -- so a single-context run is one row per region, and a
#'   multi-context run emits per-context rows sharing the same
#'   \code{jointContexts} set and jointly-estimated group priors. A run whose
#'   payload holds no posteriors is kept as one empty element.
#' @param gwasStudy Character vector of GWAS (disease) study identifiers.
#' @param studyName Character vector of QTL study identifiers.
#' @param context Character vector of context labels (scalar per run).
#' @param method Character vector of weight-method names.
#' @param entry List / \code{SimpleList} of \code{CtwasResultEntry} objects,
#'   one per run; their posteriors and susie alpha tables must be
#'   \code{GRanges}, which the pipeline builds them as.
#' @param jointStudies Optional character vector (length
#'   \code{length(gwasStudy)}) listing the joined studies for a multi-study run,
#'   \code{NA_character_} otherwise. \code{NULL} (default) omits the column.
#' @param jointContexts Optional character vector for multi-context runs. Same
#'   shape as \code{jointStudies}.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{CtwasResult} object.
#' @examples
#' at <- IRanges::IRanges(c(100L, 200L), width = 1L)
#' cre <- CtwasResultEntry(
#'   posteriors = GenomicRanges::GRanges("chr1", at,
#'     id = c("g1", "g2"), region_id = "chr1_1_1000",
#'     susie_pip = c(0.9, 0.1)),
#'   susieAlpha = GenomicRanges::GRanges("chr1", at,
#'     id = c("g1", "g2"), region_id = "chr1_1_1000",
#'     susie_alpha = c(0.9, 0.1)))
#' cr <- CtwasResult(gwasStudy = "gwas1", studyName = "s1",
#'   context = "brain", method = "susie", entry = list(cre))
#' cr
#' @export
CtwasResult <- function(
    gwasStudy,
    studyName,
    context,
    method,
    entry,
    jointStudies = NULL,
    jointContexts = NULL,
    metadata = list()
) {
    n <- length(gwasStudy)
    .ctwasCheckRunLengths(n, studyName, context, method, entry)
    supplied <- compact(list(
        jointStudies = jointStudies,
        jointContexts = jointContexts
    ))
    walk2(names(supplied), supplied, .ctwasCheckJointColLength, n = n)
    entries <- as.list(entry)
    walk(seq_along(entries), .ctwasAssertEntry, entries = entries)
    identity <- exec(
        S4Vectors::DataFrame,
        !!!c(
            list(
                gwasStudy = as.character(gwasStudy),
                study = as.character(studyName),
                context = as.character(context),
                method = as.character(method)
            ),
            map(supplied, as.character),
            list(check.names = FALSE)
        )
    )
    .withMetadata(
        .ctwasNewResult(identity, map(entries, .ctwasSplitRun)),
        metadata
    )
}

# @noRd
.ctwasCheckRunLengths <- function(n, studyName, context, method, entry) {
    if (
        length(studyName) == n &&
            length(context) == n &&
            length(method) == n &&
            length(entry) == n
    ) {
        return(invisible(NULL))
    }
    abort(str_c(
        "`gwasStudy`, `studyName`, `context`, `method`, and `entry` must ",
        "all have the same length."
    ))
}

# A run's payload must be a CtwasResultEntry whose posteriors and susie alpha
# table, when present, are ranges: an element is a GRanges, so a table whose
# ids carry no position has nowhere to go and is refused rather than dropped.
# @noRd
.ctwasAssertEntry <- function(i, entries) {
    e <- entries[[i]]
    if (!methods::is(e, "CtwasResultEntry")) {
        abort(glue(
            "CtwasResult: run {i}'s payload must be a CtwasResultEntry ",
            "(got {class(e)[[1L]]})."
        ))
    }
    parts <- list(posteriors = e@posteriors, susieAlpha = e@susieAlpha)
    tables <- names(parts)[!map_lgl(parts, .ctwasIsRangedOrEmpty)]
    if (length(tables) > 0L) {
        abort(glue(
            "CtwasResult: run {i}'s {str_flatten(tables, ' and ')} carry no ",
            "genomic coordinates, so they cannot be placed in a ranged ",
            "collection. Build them as GRanges."
        ))
    }
    invisible(NULL)
}

# @noRd
.ctwasIsRangedOrEmpty <- function(x) {
    is.null(x) || methods::is(x, "GRanges")
}

# ---- splitting a run into its regions ---------------------------------------

# One run's per-region pieces: the region ids, and for each region its
# posteriors (the element), susie alpha table and region metadata. A run with
# no posteriors at all is one empty region, so it still shows up as a row.
# @noRd
.ctwasSplitRun <- function(e) {
    post <- e@posteriors %||% GenomicRanges::GRanges()
    sa <- e@susieAlpha %||% GenomicRanges::GRanges()
    regions <- unique(c(.ctwasRegionKeys(post), .ctwasRegionKeys(sa)))
    if (length(regions) == 0L) {
        regions <- NA_character_
    }
    list(
        blockId = regions,
        elements = map(regions, .ctwasRegionSlice, gr = post),
        susieAlpha = map(regions, .ctwasRegionSlice, gr = sa),
        regionInfo = map(regions, .ctwasRegionInfoSlice, ri = e@regionInfo),
        groupPriors = rep(list(e@groupPriors), length(regions))
    )
}

# Which region each range belongs to: cTWAS's region_id, or, for a payload
# built without one, the chromosome it sits on.
# @noRd
.ctwasRegionKeys <- function(gr) {
    if (is_in("region_id", colnames(mcols(gr)))) {
        return(as.character(mcols(gr)$region_id))
    }
    as.character(seqnames(gr))
}

# @noRd
.ctwasRegionSlice <- function(region, gr) {
    if (is.na(region)) {
        return(gr[0L])
    }
    gr[.ctwasRegionKeys(gr) == region]
}

# The region's own metadata row, or the run's table as given when it is not
# keyed by region_id.
# @noRd
.ctwasRegionInfoSlice <- function(region, ri) {
    table <- .ctwasRegionInfoTable(ri)
    if (
        is.null(table) || !is_in("region_id", colnames(table)) || is.na(region)
    ) {
        return(ri)
    }
    hit <- as.character(table$region_id) == region
    if (methods::is(ri, "GRanges")) {
        return(ri[hit])
    }
    ri[hit, , drop = FALSE]
}

# Region metadata's columns, whichever form it was built in.
# @noRd
.ctwasRegionInfoTable <- function(ri) {
    if (methods::is(ri, "GRanges")) mcols(ri) else ri
}

# Every run's regions as one collection: the run identity repeated across
# its regions, the region's payload columns beside it.
# @noRd
.ctwasNewResult <- function(identity, runs) {
    reps <- rep(seq_along(runs), map_int(runs, .ctwasRunSize))
    pieces <- c("susieAlpha", "regionInfo", "groupPriors")
    payload <- set_names(map(pieces, .ctwasPayloadColumn, runs = runs), pieces)
    md <- cbind(
        identity[reps, , drop = FALSE],
        exec(
            S4Vectors::DataFrame,
            !!!c(
                list(blockId = list_c(map(runs, "blockId")) %||% character(0)),
                payload,
                list(check.names = FALSE)
            )
        )
    )
    elements <- list_flatten(map(runs, "elements"))
    grl <- S4Vectors::`mcols<-`(.rtlBindElements(unname(elements)), value = md)
    new("CtwasResult", .rtlNamed(grl))
}

# @noRd
.ctwasRunSize <- function(run) {
    length(run$blockId)
}

# @noRd
.ctwasPayloadColumn <- function(piece, runs) {
    S4Vectors::SimpleList(list_flatten(map(runs, piece)))
}

# Before it held one element per region, a CtwasResult was a DFrame of whole
# runs, one CtwasResultEntry per row in an `entry` column; an object saved
# then still carries that table as its `listData` attribute. It is rebuilt
# through the constructor, which splits each run into its regions.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "CtwasResult",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = 'CtwasResult')")
        }
        old <- attr(object, "listData", exact = TRUE)
        if (is.null(old)) {
            return(methods::callNextMethod(object, ..., verbose = verbose))
        }
        CtwasResult(
            gwasStudy = old$gwasStudy,
            studyName = old$study,
            context = old$context,
            method = old$method,
            entry = as.list(old$entry),
            jointStudies = old$jointStudies,
            jointContexts = old$jointContexts
        )
    }
)

#' @rdname methodNames
#' @export
setMethod("methodNames", "CtwasResult", function(x) {
    unique(as.character(x$method))
})

#' @rdname studyName
#' @export
setMethod("studyName", "CtwasResult", function(x) unique(as.character(x$study)))

#' @rdname contexts
#' @export
setMethod("contexts", "CtwasResult", function(x) {
    unique(as.character(x$context))
})

# Stack a per-region ranged payload (the posteriors, or the susie alpha
# tables) across all rows: every row's ranges, its columns kept as mcols and
# prefixed with the row's (gwasStudy, study, context, method) run identity.
# Columns a row lacks are NA-filled.
# @noRd
.ctwasAggregateRows <- function(x, payloads) {
    parts <- compact(
        map(seq_along(payloads), .ctwasAggregateRow, x = x, payloads = payloads)
    )
    if (length(parts) == 0L) {
        return(GenomicRanges::GRanges())
    }
    ranges <- exec(c, !!!map(unname(parts), GenomicRanges::granges))
    out <- S4Vectors::`mcols<-`(
        ranges,
        value = S4Vectors::DataFrame(
            .rbindAligned(map(parts, .ctwasRowTable)),
            check.names = FALSE
        )
    )
    .withGenome(out, .oneBuild(x))
}

#' @rdname ctwasPosteriors
#' @export
setMethod("ctwasPosteriors", "CtwasResult", function(x) {
    .ctwasAggregateRows(x, as.list(x))
})

#' @rdname susieAlpha
#' @export
setMethod("susieAlpha", "CtwasResult", function(x) {
    .ctwasAggregateRows(x, as.list(x$susieAlpha))
})

#' @rdname show-methods
#' @export
setMethod("show", "CtwasResult", function(object) {
    runs <- distinct(as_tibble(
        as.data.frame(mcols(object, use.names = FALSE)[, .CTWAS_ID_COLS])
    ))
    cat(glue(
        "CtwasResult: {nrow(runs)} run(s) over {length(object)} region ",
        "row(s)\n",
        .trim = FALSE
    ))
    if (length(object) > 0L) {
        cat(glue(
            "  {n_distinct(object$gwasStudy)} GWAS stud(y/ies), ",
            "{n_distinct(object$method)} method(s): ",
            "{str_flatten(unique(as.character(object$method)), ', ')}\n",
            .trim = FALSE
        ))
    }
})

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One row's payload with its run identity attached, ranges and all; NULL when
# the row has none.
# @noRd
.ctwasAggregateRow <- function(i, x, payloads) {
    payload <- payloads[[i]]
    if (is.null(payload) || length(payload) == 0L) {
        return(NULL)
    }
    md <- mcols(x, use.names = FALSE)
    identity <- tibble(
        gwasStudy = as.character(md$gwasStudy)[[i]],
        study = as.character(md$study)[[i]],
        context = as.character(md$context)[[i]],
        method = as.character(md$method)[[i]]
    )
    # The row identity is authoritative: a payload column of the same name
    # (cTWAS's own `context`) gives way to it rather than being duplicated.
    own <- .ctwasRowTable(payload)
    S4Vectors::`mcols<-`(
        payload,
        value = S4Vectors::DataFrame(
            bind_cols(
                slice(identity, rep(1L, length(payload))),
                select(own, !any_of(names(identity)))
            ),
            check.names = FALSE
        )
    )
}

# A ranged row part's columns as a table, for stacking across rows.
# @noRd
.ctwasRowTable <- function(part) {
    as_tibble(
        as.data.frame(mcols(part), optional = TRUE),
        .name_repair = "minimal"
    )
}
