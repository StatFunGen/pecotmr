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
#' @description S4 container for a multi-study QTL analysis: a collection
#'   of individual-level \code{QtlDataset} studies, optionally combined
#'   with a \code{QtlSumStats} collection of summary-statistic-only
#'   studies. Used as the input to multi-study fine-mapping and
#'   colocboost-style analyses.
#'
#'   At least two studies must be present in total, counting
#'   \code{qtlDatasets} entries plus the studies in \code{sumStats}.
#'
#'   For traits that appear in more than one \code{qtlDatasets} entry,
#'   the per-trait genomic positions must agree across the entries
#'   (enforced by validity). No cross-checking is performed against
#'   \code{sumStats} variant positions, since summary statistics are
#'   already computed and cannot be re-aligned at construction time.
#' @slot qtlDatasets A named list of \code{QtlDataset} objects, keyed by
#'   study identifier.
#' @slot sumStats An optional \code{QtlSumStats} carrying additional
#'   summary-statistic-only studies. \code{NULL} when absent.
#' @export
setClass("MultiStudyQtlDataset",
  representation(
    qtlDatasets = "list",
    sumStats    = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    # qtlDatasets elements
    if (!is.list(object@qtlDatasets) || length(object@qtlDatasets) == 0L) {
      errors <- c(errors, "'qtlDatasets' must be a non-empty named list")
    } else {
      nm <- names(object@qtlDatasets)
      if (is.null(nm) || any(!nzchar(nm)) || any(is.na(nm)))
        errors <- c(errors, "'qtlDatasets' must be a named list with non-empty names")
      else if (anyDuplicated(nm))
        errors <- c(errors, "names of 'qtlDatasets' must be unique")
      bad <- !vapply(object@qtlDatasets,
                    function(d) methods::is(d, "QtlDataset"), logical(1))
      if (any(bad))
        errors <- c(errors,
          "every element of 'qtlDatasets' must be a QtlDataset")
    }
    # sumStats: either NULL or a QtlSumStats
    if (!is.null(object@sumStats) && !methods::is(object@sumStats, "QtlSumStats")) {
      errors <- c(errors,
        "'sumStats' must be a QtlSumStats object or NULL")
    }
    # Total study count >= 2
    nQtl <- length(object@qtlDatasets)
    nSumstats <- if (is.null(object@sumStats)) 0L
                 else length(unique(as.character(object@sumStats$study)))
    if (length(errors) == 0L && (nQtl + nSumstats) < 2L) {
      errors <- c(errors, sprintf(
        paste0("MultiStudyQtlDataset requires at least 2 studies in total ",
               "(got %d individual-level + %d summary-statistic = %d)."),
        nQtl, nSumstats, nQtl + nSumstats))
    }
    # Pairwise trait-position consistency across qtlDatasets
    if (length(errors) == 0L && nQtl >= 2L) {
      traitRanges <- lapply(object@qtlDatasets, function(qd) {
        out <- list()
        for (ctx in seq_along(qd@phenotypes)) {
          se <- qd@phenotypes[[ctx]]
          rr <- SummarizedExperiment::rowRanges(se)
          ids <- rownames(se)
          for (i in seq_along(ids)) {
            tid <- ids[[i]]
            if (is.null(out[[tid]])) out[[tid]] <- rr[i]
          }
        }
        out
      })
      pairs <- utils::combn(seq_along(traitRanges), 2L)
      for (k in seq_len(ncol(pairs))) {
        i <- pairs[1L, k]; j <- pairs[2L, k]
        a <- traitRanges[[i]]; b <- traitRanges[[j]]
        shared <- intersect(names(a), names(b))
        for (tid in shared) {
          if (!isTRUE(all.equal(
                as.character(GenomicRanges::seqnames(a[[tid]])),
                as.character(GenomicRanges::seqnames(b[[tid]])))) ||
              GenomicRanges::start(a[[tid]]) != GenomicRanges::start(b[[tid]]) ||
              GenomicRanges::end(a[[tid]]) != GenomicRanges::end(b[[tid]])) {
            errors <- c(errors, sprintf(
              "trait '%s' has inconsistent rowRanges between studies '%s' and '%s'",
              tid, names(object@qtlDatasets)[i],
              names(object@qtlDatasets)[j]))
          }
        }
      }
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================
# MultiStudyQtlDataset constructor and accessors
# =============================================================================

#' @title Create a MultiStudyQtlDataset Object
#' @description Construct a \code{MultiStudyQtlDataset} S4 object from a
#'   named list of \code{QtlDataset} objects (individual-level studies)
#'   and an optional \code{QtlSumStats} of summary-statistic-only
#'   studies. The total study count must be at least two, satisfied by
#'   either (a) at least two \code{qtlDatasets} entries, or (b) at least
#'   one \code{qtlDatasets} entry plus a non-empty \code{sumStats}.
#' @param qtlDatasets A named list of \code{QtlDataset} objects, keyed
#'   by study identifier.
#' @param sumStats An optional \code{QtlSumStats} collection. Default
#'   \code{NULL}.
#' @return A \code{MultiStudyQtlDataset} object.
#' @export
MultiStudyQtlDataset <- function(qtlDatasets, sumStats = NULL) {
  obj <- new("MultiStudyQtlDataset",
             qtlDatasets = qtlDatasets,
             sumStats    = sumStats)
  validObject(obj)
  obj
}

#' @rdname getQtlDatasets
#' @export
setMethod("getQtlDatasets", "MultiStudyQtlDataset",
          function(x) x@qtlDatasets)

#' @rdname getSumStats
#' @export
setMethod("getSumStats", "MultiStudyQtlDataset",
  function(x, ...) {
    if (length(list(...)) > 0L) {
      stop("getSumStats(MultiStudyQtlDataset) does not accept selection ",
           "arguments; it returns the embedded QtlSumStats collection ",
           "(use getSumStats() on that result to fetch one entry).")
    }
    x@sumStats
  })

#' @rdname getStudy
#' @export
setMethod("getStudy", "MultiStudyQtlDataset", function(x) {
  fromQtl <- names(x@qtlDatasets)
  fromSs  <- if (is.null(x@sumStats)) character(0)
             else unique(as.character(x@sumStats$study))
  unique(c(fromQtl, fromSs))
})


#' @export
setMethod("show", "MultiStudyQtlDataset", function(object) {
  nQtl <- length(object@qtlDatasets)
  ssEntries <- if (is.null(object@sumStats)) 0L
               else length(unique(as.character(object@sumStats$study)))
  cat(sprintf("MultiStudyQtlDataset: %d individual-level + %d sumstats studies\n",
              nQtl, ssEntries))
  if (nQtl > 0L) {
    cat(sprintf("  Individual-level studies: %s\n",
                paste(names(object@qtlDatasets), collapse = ", ")))
  }
  if (!is.null(object@sumStats)) {
    cat(sprintf("  Sumstats studies: %s\n",
                paste(unique(as.character(object@sumStats$study)),
                      collapse = ", ")))
  }
})
