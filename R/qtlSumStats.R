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

setClass("QtlSumStats",
  contains = "SumStatsBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "context", "trait", "entry")
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
          "length(entry) must equal nrow(.) for QtlSumStats")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "GRanges"), logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a GRanges")
      keyDf <- as.data.frame(object[, c("study", "context", "trait")])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, context, trait) tuple uniqueness violated")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title QTL Summary Statistics Handling
#' @description Constructor and accessor methods for \code{QtlSumStats},
#'   the DFrame-subclass collection keyed by
#'   \code{(study, context, trait)}.
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
#' @description Construct a \code{QtlSumStats} S4 DFrame-subclass
#'   collection from per-tuple vectors and a list of \code{GRanges}
#'   entries (one per tuple), plus a single LD sketch handle and a
#'   single genome build that apply to the whole collection. Each
#'   \code{GRanges} entry must carry per-variant statistics in its
#'   mcols (\code{SNP}, \code{A1}, \code{A2}, \code{Z}, \code{N}; plus
#'   optional \code{MAF}, \code{INFO}, \code{BETA}, \code{SE}, \code{P}).
#' @param study Character vector of study identifiers (per tuple).
#' @param context Character vector of context labels (per tuple).
#' @param trait Character vector of trait identifiers (per tuple).
#' @param entry A list / \code{SimpleList} of \code{GRanges}, one per
#'   tuple. Same length as \code{study}, \code{context}, and \code{trait}.
#' @param genome Single character string giving the genome build
#'   (e.g., \code{"hg19"}, \code{"hg38"}). Uniform across the collection
#'   because all entries share the same LD sketch.
#' @param ldSketch A \code{GenotypeHandle} carrying the LD reference.
#' @param varY Optional numeric vector of per-tuple phenotype variances
#'   (\code{NA_real_} entries allowed).
#' @param ... Additional per-tuple columns to attach to the collection.
#' @return A \code{QtlSumStats} object.
#' @export
QtlSumStats <- function(study, context, trait, entry, genome, ldSketch,
                        varY = NA_real_, qcInfo = list(), ...) {
  if (missing(study) || missing(context) || missing(trait) ||
      missing(entry) || missing(genome) || missing(ldSketch)) {
    stop("`study`, `context`, `trait`, `entry`, `genome`, and `ldSketch` ",
         "are all required.")
  }
  if (length(genome) != 1L) {
    stop("`genome` must be a single character string (one build per ",
         "collection, because all entries share the LD sketch).")
  }
  if (!is.list(entry)) {
    stop("`entry` must be a list (or SimpleList) of GRanges, one per tuple.")
  }
  n <- length(study)
  if (length(context) != n || length(trait) != n || length(entry) != n) {
    stop("`study`, `context`, `trait`, and `entry` must all have the ",
         "same length.")
  }
  if (length(varY) == 1L && n > 1L) varY <- rep(varY, n)
  if (length(varY) != n) {
    stop("`varY` must have length 1 or length(study).")
  }

  cols <- list(
    study   = as.character(study),
    context = as.character(context),
    trait   = as.character(trait),
    entry   = S4Vectors::SimpleList(entry),
    varY    = as.numeric(varY)
  )
  extras <- list(...)
  for (nm in names(extras)) cols[[nm]] <- extras[[nm]]
  df <- do.call(S4Vectors::DataFrame, c(cols, list(check.names = FALSE)))

  obj <- methods::new("QtlSumStats", df,
                     ldSketch = ldSketch,
                     genome   = as.character(genome),
                     qcInfo   = as.list(qcInfo))
  methods::validObject(obj)
  obj
}

# =============================================================================
# Accessors
# =============================================================================

# Internal: resolve a (study, context, trait) tuple to a single row index.
.qtlSumStatsSelectRow <- function(x, study, context, trait) {
  if (nrow(x) == 0L) stop("QtlSumStats has no rows.")
  if (missing(study) || is.null(study) ||
      missing(context) || is.null(context) ||
      missing(trait) || is.null(trait)) {
    if (nrow(x) == 1L) return(1L)
    stop("This QtlSumStats has ", nrow(x),
         " entries. Pass `study`, `context`, and `trait` to select one.")
  }
  if (length(study) != 1L || length(context) != 1L || length(trait) != 1L) {
    stop("`study`, `context`, and `trait` must each be length 1.")
  }
  idx <- .matchTupleRows(x, list(study = study, context = context,
                                  trait = trait))
  if (length(idx) == 0L) {
    stop(sprintf(
      "No entry for (study='%s', context='%s', trait='%s').",
      study, context, trait))
  }
  if (length(idx) > 1L) {
    stop(sprintf(
      "Multiple entries match (study='%s', context='%s', trait='%s'); ",
      "tuple uniqueness violation."))
  }
  idx
}

#' @rdname getSumStats
#' @export
setMethod("getSumStats", signature(x = "QtlSumStats"),
  function(x, study = NULL, context = NULL, trait = NULL, ...) {
    idx <- .qtlSumStatsSelectRow(x, study, context, trait)
    x$entry[[idx]]
  }
)

#' @rdname getZ
#' @export
setMethod("getZ", "QtlSumStats",
  function(x, study = NULL, context = NULL, trait = NULL) {
    gr <- getSumStats(x, study = study, context = context, trait = trait)
    mcols(gr)$Z
  }
)

#' @rdname getN
#' @export
setMethod("getN", "QtlSumStats",
  function(x, study = NULL, context = NULL, trait = NULL) {
    gr <- getSumStats(x, study = study, context = context, trait = trait)
    mcols(gr)$N
  }
)

#' @rdname getMaf
#' @export
setMethod("getMaf", "QtlSumStats",
  function(x, study = NULL, context = NULL, trait = NULL, ...) {
    gr <- getSumStats(x, study = study, context = context, trait = trait)
    mc <- mcols(gr)
    if ("MAF" %in% colnames(mc)) mc$MAF else NULL
  }
)

#' @rdname getSumstatDf
#' @export
setMethod("getSumstatDf", "QtlSumStats",
  function(x, study = NULL, context = NULL, trait = NULL,
           require = character(0),
           derive  = c("none", "zFromBetaSe"),
           keepChrPrefix = TRUE) {
    derive <- match.arg(derive)
    gr <- getSumStats(x, study = study, context = context, trait = trait)
    .entryToSumstatDf(gr,
                      require       = require,
                      derive        = derive,
                      keepChrPrefix = keepChrPrefix,
                      label         = sprintf(
                        "QtlSumStats[%s/%s/%s]",
                        if (is.null(study))   "<auto>" else study,
                        if (is.null(context)) "<auto>" else context,
                        if (is.null(trait))   "<auto>" else trait))
  }
)

#' @rdname nSnps
#' @export
setMethod("nSnps", "QtlSumStats",
  function(x, study = NULL, context = NULL, trait = NULL) {
    gr <- getSumStats(x, study = study, context = context, trait = trait)
    length(gr)
  }
)

#' @rdname subsetChr
#' @export
setMethod("subsetChr", "QtlSumStats", function(x, chr) {
  chrName <- paste0("chr", sub("^chr", "", as.character(chr)))
  newEntries <- lapply(seq_len(nrow(x)), function(i) {
    gr <- x$entry[[i]]
    idx <- as.character(seqnames(gr)) == chrName
    gr[idx]
  })
  QtlSumStats(
    study    = as.character(x$study),
    context  = as.character(x$context),
    trait    = as.character(x$trait),
    entry    = newEntries,
    genome   = x@genome,
    ldSketch = x@ldSketch,
    varY     = as.numeric(x$varY),
    qcInfo   = x@qcInfo)
})

#' @rdname getVarY
#' @export
setMethod("getVarY", "QtlSumStats",
  function(x, study = NULL, context = NULL, trait = NULL) {
    idx <- .qtlSumStatsSelectRow(x, study, context, trait)
    val <- x$varY[[idx]]
    if (is.na(val)) NULL else val
  }
)

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlSumStats",
  function(x) unique(as.character(x$context)))

#' @rdname getTraits
#' @export
setMethod("getTraits", "QtlSumStats",
  function(x) unique(as.character(x$trait)))

# =============================================================================
# Show method
# =============================================================================

#' @export
setMethod("show", "QtlSumStats", function(object) {
  cat(sprintf("QtlSumStats: %d entries, genome build %s\n",
              nrow(object), getGenome(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d contexts, %d traits\n",
                length(unique(object$study)),
                length(unique(object$context)),
                length(unique(object$trait))))
  }
  ld <- getLdSketch(object)
  cat(sprintf("  LD sketch: %s @ %s\n", getFormat(ld), getPath(ld)))
})
