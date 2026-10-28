# =============================================================================
# MultiStudyQtlDataset S4 class
# -----------------------------------------------------------------------------
# Multi-study container: a SimpleList of individual-level QtlDataset objects,
# keyed by study, plus an optional QtlSumStats collection for
# summary-statistic-only studies. Used as the entry point for multi-study
# pipelines (fineMappingPipeline, twasWeightsPipeline) that orchestrate joint
# or per-study analyses across multiple cohorts. Being a List, it subsets,
# indexes and maps like one: ms["s1"], ms[["s1"]], names(ms), lapply(ms, f).
# =============================================================================

#' @include QtlDataset.R qtlSumStats.R
NULL

# The summary-statistic-only studies: a QtlSumStats, or NULL when there are
# none. The slot type refuses anything else when the object is built.
setClassUnion("QtlSumStats_OR_NULL", c("QtlSumStats", "NULL"))

#' @title Multi-Study QTL Dataset
#' @description S4 container for a multi-study QTL analysis: a collection of
#'   individual-level \code{QtlDataset} studies, optionally combined with a
#'   \code{QtlSumStats} collection of summary-statistic-only studies. Used as
#'   the input to multi-study fine-mapping and colocboost-style analyses.
#'
#' It is a \code{\link[S4Vectors]{SimpleList}} of the individual-level
#' studies, keyed by study name, so \code{length()}, \code{names()},
#' \code{[[} and \code{[} address them as on any Bioconductor list;
#' \code{[} keeps the summary-statistic studies whole and still requires two
#' studies in total.
#'
#' At least two studies must be present in total, counting \code{qtlDatasets}
#' entries plus the studies in \code{sumStats}.
#'
#' For traits that appear in more than one \code{qtlDatasets} entry, the
#' per-trait genomic positions must agree across the entries (enforced by
#' validity). No cross-checking is performed against \code{sumStats} variant
#' positions, since summary statistics are already computed and cannot be
#' re-aligned at construction time.
#' @slot sumStats An optional \code{QtlSumStats} carrying additional
#'   summary-statistic-only studies. \code{NULL} when absent.
#' @export
setClass(
    "MultiStudyQtlDataset",
    contains = "SimpleList",
    representation(sumStats = "QtlSumStats_OR_NULL"),
    prototype(elementType = "QtlDataset"),
    validity = function(object) .validateMultiStudyQtlDataset(object)
)

# ---- MultiStudyQtlDataset validity helpers ---------------------------------

# @noRd
.validateMultiStudyQtlDataset <- function(object) {
    # Each stage only makes sense once the previous one holds: the study count
    # reads the list the first check just validated, and trait consistency
    # reads across the studies the count check just confirmed. The element
    # type is SimpleList's own check (elementType), which runs first.
    slotErrors <- .msqdCheckDatasets(as.list(object))
    if (length(slotErrors) > 0L) {
        return(slotErrors)
    }
    countErrors <- .msqdCheckStudyCount(object)
    if (length(countErrors) > 0L) {
        return(countErrors)
    }
    traitErrors <- .msqdCheckTraitConsistency(object)
    if (length(traitErrors) == 0L) TRUE else traitErrors
}

# @noRd
#' @importFrom checkmate checkList
.msqdCheckDatasets <- function(qtlDatasets) {
    # names = "unique" subsumes the old named/non-empty/non-NA/unique checks.
    res <- checkList(qtlDatasets, min.len = 1L, names = "unique")
    if (!isTRUE(res)) {
        return(str_c("'qtlDatasets' ", res))
    }
    NULL
}

# SimpleList refuses an element of another class as well, but without naming
# it, so the constructor checks first.
# @noRd
.msqdAssertDatasetTypes <- function(qtlDatasets) {
    bad <- names(qtlDatasets)[!map_lgl(qtlDatasets, .msqdIsQtlDataset)]
    if (length(bad) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "every element of `qtlDatasets` must be a QtlDataset; ",
        "{str_flatten(bad, ', ')} {if (length(bad) == 1L) 'is' else 'are'} ",
        "not."
    ))
}

# @noRd
.msqdIsQtlDataset <- function(d) {
    methods::is(d, "QtlDataset")
}

# At least 2 studies in total (individual-level + summary-statistic).
# @noRd
.msqdCheckStudyCount <- function(object) {
    nQtl <- length(object)
    nSumstats <- if (is.null(object@sumStats)) {
        0L
    } else {
        n_distinct(as.character(object@sumStats$study))
    }
    if ((nQtl + nSumstats) < 2L) {
        return(glue(
            "MultiStudyQtlDataset requires at least 2 studies in total ",
            "(got {nQtl} individual-level + {nSumstats} summary-statistic ",
            "= {nQtl + nSumstats})."
        ))
    }
    NULL
}

# Pairwise trait-position consistency: a trait shared across studies must have
# identical rowRanges.
# @noRd
.msqdCheckTraitConsistency <- function(object) {
    if (length(object) < 2L) {
        return(NULL)
    }
    traitRanges <- map(as.list(object), .msqdTraitRanges)
    pairs <- utils::combn(seq_along(traitRanges), 2L)
    dsNames <- names(object)
    list_c(compact(map(
        seq_len(ncol(pairs)),
        .msqdPairErrors,
        pairs = pairs,
        traitRanges = traitRanges,
        dsNames = dsNames
    )))
}

# Per-dataset trait -> rowRanges map (first occurrence of each trait id).
# @noRd
.msqdTraitRanges <- function(qd) {
    pairs <- .msqdConcat(map(contexts(qd), .msqdContextTraitRanges, qd = qd))
    if (length(pairs) == 0L) {
        return(list())
    }
    # First occurrence of each trait id wins, as the "only set it while still
    # NULL" assignment did.
    tids <- map_chr(pairs, "tid")
    keep <- !duplicated(tids)
    set_names(map(pairs[keep], "range"), tids[keep])
}

# @noRd
.msqdConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# @noRd
.msqdTraitRangePair <- function(i, ids, rr) {
    list(tid = ids[[i]], range = rr[i])
}

# One context's (trait, range) pairs in row order.
# @noRd
.msqdContextTraitRanges <- function(ctx, qd) {
    se <- .qtlContextSe(qd, ctx)
    map(
        seq_along(rownames(se)),
        .msqdTraitRangePair,
        ids = rownames(se),
        rr = SummarizedExperiment::rowRanges(se)
    )
}

# Inconsistency errors for the k-th dataset pair.
# @noRd
.msqdPairErrors <- function(k, pairs, traitRanges, dsNames) {
    a <- traitRanges[[pairs[1L, k]]]
    b <- traitRanges[[pairs[2L, k]]]
    shared <- intersect(names(a), names(b))
    bad <- shared[!map_lgl(shared, .msqdRangesEqual, a = a, b = b)]
    if (length(bad) == 0L) {
        return(NULL)
    }
    map_chr(
        bad,
        .msqdTraitError,
        ni = dsNames[pairs[1L, k]],
        nj = dsNames[pairs[2L, k]]
    )
}

# TRUE when trait `tid`'s ranges match across the two datasets.
# @noRd
.msqdRangesEqual <- function(tid, a, b) {
    identical(
        as.character(GenomicRanges::seqnames(a[[tid]])),
        as.character(GenomicRanges::seqnames(b[[tid]]))
    ) &&
        GenomicRanges::start(a[[tid]]) == GenomicRanges::start(b[[tid]]) &&
        GenomicRanges::end(a[[tid]]) == GenomicRanges::end(b[[tid]])
}

# @noRd
.msqdTraitError <- function(tid, ni, nj) {
    glue(
        "trait '{tid}' has inconsistent rowRanges between studies '{ni}' ",
        "and '{nj}'"
    )
}

# =============================================================================
# MultiStudyQtlDataset constructor and accessors
# =============================================================================

#' @title Create a MultiStudyQtlDataset Object
#' @description Construct a \code{MultiStudyQtlDataset} S4 object from a named
#'   list of \code{QtlDataset} objects (individual-level studies) and an
#'   optional \code{QtlSumStats} of summary-statistic-only studies. The total
#'   study count must be at least two, satisfied by either (a) at least two
#'   \code{qtlDatasets} entries, or (b) at least one \code{qtlDatasets} entry
#'   plus a non-empty \code{sumStats}.
#' @param qtlDatasets A named list of \code{QtlDataset} objects, keyed by study
#'   identifier.
#' @param sumStats An optional \code{QtlSumStats} collection. Default
#'   \code{NULL}.
#' @return A \code{MultiStudyQtlDataset} object.
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr")
#' )
#' rng <- GenomicRanges::GRanges(
#'   "chr22", IRanges::IRanges(14600000L, width = 1000L)
#' )
#' names(rng) <- "ENSG1"
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(expression = matrix(
#'     rnorm(ncol(panel)), 1,
#'     dimnames = list("ENSG1", colnames(panel))
#'   )),
#'   rowRanges = rng
#' )
#' qd1 <- QtlDataset(
#'   studyName = "s1", genotypes = panel, phenotypes = list(brain = se)
#' )
#' qd2 <- QtlDataset(
#'   studyName = "s2", genotypes = panel, phenotypes = list(brain = se)
#' )
#' MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
#' @importFrom checkmate assertList
#' @export
MultiStudyQtlDataset <- function(qtlDatasets, sumStats = NULL) {
    assertList(qtlDatasets, min.len = 1L, names = "unique")
    .msqdAssertDatasetTypes(qtlDatasets)
    new(
        "MultiStudyQtlDataset",
        listData = qtlDatasets,
        sumStats = sumStats
    )
}

# Selecting studies keeps the summary-statistic ones whole, and re-checks the
# result: the stock List subset does not run validity, and a selection below
# two studies in total is not a multi-study dataset.
setMethod("extractROWS", "MultiStudyQtlDataset", function(x, i) {
    ans <- methods::callNextMethod()
    validObject(ans)
    ans
})

# Before it became a SimpleList the studies sat in a `qtlDatasets` slot, which
# an object saved then still carries as an attribute; either way they are
# rebuilt through the constructor, each one updated first.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "MultiStudyQtlDataset",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = 'MultiStudyQtlDataset')")
        }
        studies <- attr(object, "qtlDatasets", exact = TRUE) %||%
            attr(object, "listData", exact = TRUE)
        sumStats <- .msqdSavedSumStats(object)
        MultiStudyQtlDataset(
            qtlDatasets = map(studies, updateObject, ..., verbose = verbose),
            sumStats = if (is.null(sumStats)) {
                NULL
            } else {
                updateObject(sumStats, ..., verbose = verbose)
            }
        )
    }
)

# The summary statistics an object saved under any earlier definition holds,
# NULL when it holds none. slot() rather than attr(): R stores a NULL slot as
# a placeholder symbol, which attr() hands back unconverted.
# @noRd
.msqdSavedSumStats <- function(object) {
    if (is.null(attr(object, "sumStats", exact = TRUE))) {
        return(NULL)
    }
    methods::slot(object, "sumStats")
}

#' @rdname sumStats
#' @export
# No selectors: this returns the embedded QtlSumStats collection whole, and
# one entry is fetched by calling sumStats() on that result. Taking no
# extra formals is what refuses a selector -- R reports the unused argument
# by name.
setMethod("sumStats", "MultiStudyQtlDataset", function(x) x@sumStats)

#' @rdname studyName
#' @export
setMethod("studyName", "MultiStudyQtlDataset", function(x) {
    fromQtl <- names(x)
    fromSs <- if (is.null(x@sumStats)) {
        character(0)
    } else {
        unique(as.character(x@sumStats$study))
    }
    unique(c(fromQtl, fromSs))
})

#' @rdname show-methods
#' @export
setMethod("show", "MultiStudyQtlDataset", function(object) {
    nQtl <- length(object)
    ssEntries <- if (is.null(object@sumStats)) {
        0L
    } else {
        n_distinct(as.character(object@sumStats$study))
    }
    cat(glue(
        "MultiStudyQtlDataset: {nQtl} individual-level + ",
        "{ssEntries} sumstats studies\n",
        .trim = FALSE
    ))
    if (nQtl > 0L) {
        cat(glue(
            "  Individual-level studies: ",
            "{str_flatten(names(object), ', ')}\n",
            .trim = FALSE
        ))
    }
    if (!is.null(object@sumStats)) {
        cat(glue(
            "  Sumstats studies: ",
            "{str_flatten(unique(as.character(object@sumStats$study)), ', ')}",
            "\n",
            .trim = FALSE
        ))
    }
})
