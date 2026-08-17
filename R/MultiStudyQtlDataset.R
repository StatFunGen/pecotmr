# =============================================================================
# MultiStudyQtlDataset S4 class
# -----------------------------------------------------------------------------
# Multi-study container that holds a named list of individual-level
# QtlDataset objects plus an optional QtlSumStats collection for
# summary-statistic-only studies. Used as the entry point for
# multi-study pipelines (fineMappingPipeline, twasWeightsPipeline) that
# orchestrate joint or per-study analyses across multiple cohorts.
# =============================================================================

#' @include QtlDataset.R
NULL

#' @title Multi-Study QTL Dataset
#' @description S4 container for a multi-study QTL analysis: a collection of
#'   individual-level \code{QtlDataset} studies, optionally combined with a
#'   \code{QtlSumStats} collection of summary-statistic-only studies. Used as
#'   the input to multi-study fine-mapping and colocboost-style analyses.
#'
#' At least two studies must be present in total, counting \code{qtlDatasets}
#' entries plus the studies in \code{sumStats}.
#'
#' For traits that appear in more than one \code{qtlDatasets} entry, the
#' per-trait genomic positions must agree across the entries (enforced by
#' validity). No cross-checking is performed against \code{sumStats} variant
#' positions, since summary statistics are already computed and cannot be
#' re-aligned at construction time.
#' @slot qtlDatasets A named list of \code{QtlDataset} objects, keyed by study
#'   identifier.
#' @slot sumStats An optional \code{QtlSumStats} carrying additional
#'   summary-statistic-only studies. \code{NULL} when absent.
#' @export
setClass(
    "MultiStudyQtlDataset",
    representation(
        qtlDatasets = "list",
        sumStats = "ANY"
    ),
    validity = function(object) .validateMultiStudyQtlDataset(object)
)

# ---- MultiStudyQtlDataset validity helpers ---------------------------------

# @noRd
.validateMultiStudyQtlDataset <- function(object) {
    errors <- c(
        .msqdCheckDatasets(object@qtlDatasets),
        .msqdCheckSumStats(object@sumStats)
    )
    if (length(errors) == 0L) {
        errors <- .msqdCheckStudyCount(object)
    }
    if (length(errors) == 0L) {
        errors <- .msqdCheckTraitConsistency(object)
    }
    if (length(errors) == 0L) TRUE else errors
}

# @noRd
.msqdCheckDatasets <- function(qtlDatasets) {
    if (!is.list(qtlDatasets) || length(qtlDatasets) == 0L) {
        return("'qtlDatasets' must be a non-empty named list")
    }
    c(
        .msqdCheckDatasetNames(names(qtlDatasets)),
        .msqdCheckDatasetTypes(qtlDatasets)
    )
}

# @noRd
.msqdCheckDatasetNames <- function(nm) {
    isEmpty <- any(str_length(nm) == 0L, na.rm = TRUE)
    if (is.null(nm) || isEmpty || any(is.na(nm))) {
        return("'qtlDatasets' must be a named list with non-empty names")
    }
    if (n_distinct(nm) < length(nm)) {
        return("names of 'qtlDatasets' must be unique")
    }
    NULL
}

# @noRd
.msqdCheckDatasetTypes <- function(qtlDatasets) {
    if (any(!map_lgl(qtlDatasets, .msqdIsQtlDataset))) {
        return("every element of 'qtlDatasets' must be a QtlDataset")
    }
    NULL
}

# @noRd
.msqdIsQtlDataset <- function(d) {
    methods::is(d, "QtlDataset")
}

# @noRd
.msqdCheckSumStats <- function(sumStats) {
    if (!is.null(sumStats) && !methods::is(sumStats, "QtlSumStats")) {
        return("'sumStats' must be a QtlSumStats object or NULL")
    }
    NULL
}

# At least 2 studies in total (individual-level + summary-statistic).
# @noRd
.msqdCheckStudyCount <- function(object) {
    nQtl <- length(object@qtlDatasets)
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
    if (length(object@qtlDatasets) < 2L) {
        return(NULL)
    }
    traitRanges <- map(object@qtlDatasets, .msqdTraitRanges)
    pairs <- utils::combn(seq_along(traitRanges), 2L)
    dsNames <- names(object@qtlDatasets)
    unlist(compact(map(
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
    out <- list()
    for (ctx in seq_along(qd@phenotypes)) {
        se <- qd@phenotypes[[ctx]]
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- rownames(se)
        for (i in seq_along(ids)) {
            tid <- ids[[i]]
            if (is.null(out[[tid]])) {
                out[[tid]] <- rr[i]
            }
        }
    }
    out
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
    isTRUE(all.equal(
        canonChrom(GenomicRanges::seqnames(a[[tid]])),
        canonChrom(GenomicRanges::seqnames(b[[tid]]))
    )) &&
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
#' gh <- new("GenotypeHandle", path = "toy.gds", format = "gds",
#'   snpInfo = data.frame(SNP = paste0("rs", 1:3), CHR = "1",
#'     BP = c(100L, 200L, 300L), A1 = "A", A2 = "G"),
#'   nSamples = 6L, sampleIds = paste0("s", 1:6), pgenPtr = NULL)
#' rng <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000L, width = 500L))
#' names(rng) <- "ENSG1"
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(expression = matrix(rnorm(6), 1, 6,
#'     dimnames = list("ENSG1", paste0("s", 1:6)))), rowRanges = rng)
#' qd1 <- QtlDataset(study = "s1", genotypes = gh,
#'   phenotypes = list(brain = se),
#'   genotypeCovariates = matrix(0, 6, 0))
#' qd2 <- QtlDataset(study = "s2", genotypes = gh,
#'   phenotypes = list(brain = se),
#'   genotypeCovariates = matrix(0, 6, 0))
#' MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
#' @export
MultiStudyQtlDataset <- function(qtlDatasets, sumStats = NULL) {
    obj <- new(
        "MultiStudyQtlDataset",
        qtlDatasets = qtlDatasets,
        sumStats = sumStats
    )
    validObject(obj)
    obj
}

#' @rdname getQtlDatasets
#' @export
setMethod("getQtlDatasets", "MultiStudyQtlDataset", function(x) x@qtlDatasets)

#' @rdname getSumStats
#' @export
setMethod("getSumStats", "MultiStudyQtlDataset", function(x, ...) {
    if (length(list(...)) > 0L) {
        msg <- glue(
            "getSumStats(MultiStudyQtlDataset) does not accept selection ",
            "arguments; it returns the embedded QtlSumStats collection ",
            "(use getSumStats() on that result to fetch one entry)."
        )
        abort(msg)
    }
    x@sumStats
})

#' @rdname getStudy
#' @export
setMethod("getStudy", "MultiStudyQtlDataset", function(x) {
    fromQtl <- names(x@qtlDatasets)
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
    nQtl <- length(object@qtlDatasets)
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
            "{str_flatten(names(object@qtlDatasets), ', ')}\n",
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
