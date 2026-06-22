# =============================================================================
# Tuple row matchers
# -----------------------------------------------------------------------------
# Internal helpers shared by the FineMappingResult / TwasWeights /
# SumStats DFrame-subclass collections to resolve a tuple-keyed selection
# to a single row index. Pure R helpers — no S4 dispatch, no exports.
# =============================================================================

# Internal: return integer row indices of `x` where every (column, value)
# pair in `keys` matches as.character(x[[column]]) == value. Shared
# building block for tuple-keyed row selectors and cache lookups
# (.tupleSelectRow, .qtlSumStatsSelectRow, .gwasSelectStudy,
# .fmCacheLookup, .cipFmrHasTuple, etc.). Pure vectorised AND-match with
# character coercion -- no validation, no error reporting.
.matchTupleRows <- function(x, keys) {
  if (length(keys) == 0L) return(seq_len(nrow(x)))
  ok <- rep(TRUE, nrow(x))
  for (k in names(keys)) {
    ok <- ok & as.character(x[[k]]) == keys[[k]]
  }
  which(ok)
}

# Internal: resolve a tuple-keyed selection (study, context, trait,
# method) to a single row index. Used by the QtlFineMappingResult and
# TwasWeights accessors. Returns an error when no row matches; returns
# the single row index when the collection has exactly one row and any
# selector argument was omitted.
.tupleSelectRow <- function(x, study, context, trait, method,
                            cls = "QtlFineMappingResult") {
  if (nrow(x) == 0L) stop(cls, " has no rows.")
  if (missing(study) || is.null(study) ||
      missing(context) || is.null(context) ||
      missing(trait) || is.null(trait) ||
      missing(method) || is.null(method)) {
    if (nrow(x) == 1L) return(1L)
    stop(cls, " has ", nrow(x), " entries. Pass `study`, `context`, ",
         "`trait`, and `method` to select one.")
  }
  if (length(study) != 1L || length(context) != 1L ||
      length(trait) != 1L || length(method) != 1L) {
    stop("`study`, `context`, `trait`, and `method` must each be length 1.")
  }
  idx <- .matchTupleRows(x, list(study = study, context = context,
                                  trait = trait, method = method))
  if (length(idx) == 0L) {
    stop(sprintf(
      "No entry for (study='%s', context='%s', trait='%s', method='%s').",
      study, context, trait, method))
  }
  idx[[1L]]
}

# Internal: resolve a (study, method, region_id) tuple to a single row
# index of a GwasFineMappingResult collection. `region` may be NULL when
# the (study, method) pair maps to a single row; otherwise it disambiguates
# among per-block rows of a genome-wide collection.
.tupleSelectRowGwasFmr <- function(x, study, method, region = NULL) {
  if (nrow(x) == 0L) stop("GwasFineMappingResult has no rows.")
  if (missing(study) || is.null(study) ||
      missing(method) || is.null(method)) {
    if (nrow(x) == 1L) return(1L)
    stop("GwasFineMappingResult has ", nrow(x), " entries. Pass `study` ",
         "and `method` to select one.")
  }
  if (length(study) != 1L || length(method) != 1L)
    stop("`study` and `method` must each be length 1.")
  if (!is.null(region) && length(region) != 1L)
    stop("`region` must be length 1 when supplied.")
  keys <- list(study = study, method = method)
  if (!is.null(region)) keys$region_id <- region
  idx <- .matchTupleRows(x, keys)
  if (length(idx) == 0L) {
    stop(sprintf(
      "No entry for (study='%s', method='%s'%s).",
      study, method,
      if (is.null(region)) "" else sprintf(", region='%s'", region)))
  }
  if (length(idx) > 1L) {
    stop(sprintf(
      "GwasFineMappingResult has %d rows matching (study='%s', method='%s'); ",
      length(idx), study, method),
      "pass `region` to disambiguate (available: ",
      paste(shQuote(as.character(x$region_id[idx])), collapse = ", "),
      ").")
  }
  idx[[1L]]
}
