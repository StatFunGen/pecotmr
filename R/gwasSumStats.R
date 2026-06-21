# =============================================================================
# GwasSumStats S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study). Each
# row holds a per-study GRanges of GWAS summary statistics covering a
# single LD block; build a separate collection per block when sweeping
# the genome. Class-level slots ldSketch + genome + qcInfo apply
# uniformly across rows.
# =============================================================================

#' @include SumStatsBase.R tupleSelectors.R
NULL

setClass("GwasSumStats",
  contains = "SumStatsBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(object@genome) != 1L || !nzchar(object@genome))
      errors <- c(errors,
        "'genome' slot must be a single non-empty character string")
    if (!is.list(object@qcInfo))
      errors <- c(errors, "'qcInfo' slot must be a list")
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for GwasSumStats")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "GRanges"), logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a GRanges")
      if (anyDuplicated(as.character(object$study)))
        errors <- c(errors, "`study` must be unique")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)


setMethod("show", "GwasSumStats", function(object) {
  cat(sprintf("GwasSumStats: %d studies, genome build %s\n",
              nrow(object), object@genome))
  cat(sprintf("  LD sketch: %s @ %s\n",
              object@ldSketch@format, object@ldSketch@path))
})


#' @title GWAS Summary Statistics Handling
#' @description Constructor, accessors, and converters for
#'   \code{GwasSumStats} (the post-refactor DFrame-subclass collection
#'   keyed by \code{study}).
#' @name pecotmr-gwas-sumstats
#' @keywords internal
#' @importFrom GenomicRanges GRanges seqnames start
#' @importFrom S4Vectors DataFrame mcols mcols<- SimpleList
#' @importFrom IRanges IRanges
#' @include allGenerics.R
NULL

# =============================================================================
# Constructor
# =============================================================================

#' @title Create a GwasSumStats Collection Object
#' @description Construct a \code{GwasSumStats} S4 DFrame-subclass
#'   collection from per-study tuple vectors and a list of \code{GRanges}
#'   entries (one per study), plus a single LD sketch handle and a
#'   single genome build that apply to the whole collection.
#'
#'   Each \code{GRanges} entry must carry per-variant statistics in its
#'   mcols (at minimum \code{SNP}, \code{A1}, \code{A2}, \code{Z},
#'   \code{N}; optionally \code{MAF}, \code{INFO}, \code{BETA}, \code{SE},
#'   \code{P}).
#' @param study Character vector of study identifiers (must be unique).
#' @param entry A \code{SimpleList} or \code{list} of \code{GRanges},
#'   one per study.
#' @param genome Single character string giving the genome build
#'   (e.g., \code{"hg19"}, \code{"hg38"}). Uniform across the collection
#'   because all entries share the same LD sketch.
#' @param ldSketch A \code{GenotypeHandle} carrying the LD reference.
#' @param varY Optional numeric vector of per-study phenotype variances
#'   (\code{NA_real_} entries allowed). Used by the sufficient-statistic
#'   interface; z-score RSS analyses should leave entries as NA.
#' @param ... Additional per-study columns to attach to the collection.
#' @return A \code{GwasSumStats} object.
#' @export
GwasSumStats <- function(study, entry, genome, ldSketch,
                          varY = NA_real_, qcInfo = list(), ...) {
  if (missing(study) || missing(entry) || missing(genome) || missing(ldSketch)) {
    stop("`study`, `entry`, `genome`, and `ldSketch` are all required.")
  }
  if (length(genome) != 1L) {
    stop("`genome` must be a single character string (one build per ",
         "collection, because all entries share the LD sketch).")
  }
  if (!is.list(entry)) {
    stop("`entry` must be a list (or SimpleList) of GRanges, one per study.")
  }
  if (length(entry) != length(study)) {
    stop("length(entry) (", length(entry),
         ") must equal length(study) (", length(study), ").")
  }
  if (length(varY) == 1L && length(study) > 1L) {
    varY <- rep(varY, length(study))
  }
  if (length(varY) != length(study)) {
    stop("`varY` must have length 1 or length(study).")
  }

  cols <- list(
    study = as.character(study),
    entry = S4Vectors::SimpleList(entry),
    varY  = as.numeric(varY)
  )
  extras <- list(...)
  for (nm in names(extras)) cols[[nm]] <- extras[[nm]]
  df <- do.call(S4Vectors::DataFrame, c(cols, list(check.names = FALSE)))

  obj <- methods::new("GwasSumStats", df,
                     ldSketch = ldSketch,
                     genome   = as.character(genome),
                     qcInfo   = as.list(qcInfo))
  methods::validObject(obj)
  obj
}


# =============================================================================
# Accessors for the new GwasSumStats collection
# =============================================================================

# Internal: resolve a study selection to a single row index. Errors when
# `study` is missing on a multi-study collection.
.gwasSelectStudy <- function(x, study) {
  if (nrow(x) == 0L) stop("GwasSumStats has no rows.")
  if (missing(study) || is.null(study)) {
    if (nrow(x) == 1L) return(1L)
    stop("This GwasSumStats has ", nrow(x),
         " studies. Pass `study = <name>` to select one. ",
         "Available: ", paste(as.character(x$study), collapse = ", "))
  }
  idx <- match(study, as.character(x$study))
  if (is.na(idx)) {
    stop("Unknown study: '", study,
         "'. Available: ",
         paste(as.character(x$study), collapse = ", "))
  }
  idx
}

#' @title Get a GWAS Study's Summary-Statistic GRanges
#' @description Return the per-variant \code{GRanges} of summary
#'   statistics for one study in a \code{GwasSumStats} collection.
#' @param x A \code{GwasSumStats} object.
#' @param study Character (length 1) study identifier. Optional when the
#'   collection has a single row.
#' @return A \code{GRanges} object.
#' @export
setMethod("getSumStats", signature(x = "GwasSumStats"),
  function(x, study = NULL, ...) {
    idx <- .gwasSelectStudy(x, study)
    x$entry[[idx]]
  }
)

#' @rdname getZ
#' @export
setMethod("getZ", "GwasSumStats", function(x, study = NULL) {
  gr <- getSumStats(x, study = study)
  mcols(gr)$Z
})

#' @rdname getN
#' @export
setMethod("getN", "GwasSumStats", function(x, study = NULL) {
  gr <- getSumStats(x, study = study)
  mcols(gr)$N
})

#' @rdname getMaf
#' @export
setMethod("getMaf", "GwasSumStats", function(x, study = NULL, ...) {
  gr <- getSumStats(x, study = study)
  mc <- mcols(gr)
  if ("MAF" %in% colnames(mc)) mc$MAF else NULL
})

#' @rdname getSumstatDf
#' @export
setMethod("getSumstatDf", "GwasSumStats",
  function(x, study = NULL,
           require = character(0),
           derive  = c("none", "zFromBetaSe"),
           keepChrPrefix = TRUE) {
    derive <- match.arg(derive)
    gr <- getSumStats(x, study = study)
    .entryToSumstatDf(gr,
                      require       = require,
                      derive        = derive,
                      keepChrPrefix = keepChrPrefix,
                      label         = sprintf("GwasSumStats[%s]",
                                              if (is.null(study)) "<auto>"
                                              else study))
  })

#' @rdname nSnps
#' @export
setMethod("nSnps", "GwasSumStats", function(x, study = NULL) {
  gr <- getSumStats(x, study = study)
  length(gr)
})

#' @rdname subsetChr
#' @export
setMethod("subsetChr", "GwasSumStats", function(x, chr) {
  chrName <- paste0("chr", sub("^chr", "", as.character(chr)))
  newEntries <- lapply(seq_len(nrow(x)), function(i) {
    gr <- x$entry[[i]]
    idx <- as.character(seqnames(gr)) == chrName
    gr[idx]
  })
  GwasSumStats(
    study    = as.character(x$study),
    entry    = newEntries,
    genome   = x@genome,
    ldSketch = x@ldSketch,
    varY     = as.numeric(x$varY),
    qcInfo   = x@qcInfo)
})

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
#' @description Extracts the per-variant statistics for one study
#'   (selected by \code{study}) into a plain data.frame with columns
#'   SNP, CHR, BP, A1, A2, Z, N (and any optional columns such as MAF,
#'   BETA, SE, P).
#' @param x A \code{GwasSumStats} object.
#' @param row.names Ignored (present for S3 generic compatibility).
#' @param optional Ignored.
#' @param study Character (length 1) study identifier. Optional when the
#'   collection has a single row.
#' @param ... Ignored.
#' @return A data.frame.
#' @method as.data.frame GwasSumStats
#' @export
as.data.frame.GwasSumStats <- function(x, row.names = NULL, optional = FALSE,
                                       study = NULL, ...) {
  gr <- getSumStats(x, study = study)
  mc <- as.data.frame(mcols(gr))
  mc$CHR <- as.character(seqnames(gr))
  mc$BP  <- start(gr)
  firstCols <- c("SNP", "CHR", "BP")
  restCols  <- setdiff(names(mc), firstCols)
  mc[, c(firstCols, restCols), drop = FALSE]
}
