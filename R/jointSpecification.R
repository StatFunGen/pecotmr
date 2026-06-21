# =============================================================================
# Joint-specification grammar and ragged input-argument parsing for the
# fineMapping / twasWeights pipelines. Pure validation + normalization;
# no pipeline dispatch or fits live here.
# =============================================================================

# -----------------------------------------------------------------------------
# Internal scope helpers — what (study, context, trait, dataForm) tuples does
# the input cover? `dataForm` is "individual" for QtlDataset-located studies
# and "sumstats" for QtlSumStats-located studies.
# -----------------------------------------------------------------------------

# Return character vector of all studies present in `data`.
# @noRd
.spListStudies <- function(data) {
  if (is(data, "QtlDataset"))           return(data@study)
  if (is(data, "QtlSumStats"))          return(unique(as.character(data$study)))
  if (is(data, "MultiStudyQtlDataset")) {
    indStudies <- names(getQtlDatasets(data))
    ss <- getSumStats(data)
    ssStudies <- if (is.null(ss)) character(0)
                 else unique(as.character(ss$study))
    return(unique(c(indStudies, ssStudies)))
  }
  stop(".spListStudies: unsupported class: ", class(data)[[1L]])
}

# Return "individual" or "sumstats" for a single study in `data`. Errors if
# the study is not present.
# @noRd
.spStudyDataForm <- function(data, study) {
  if (is(data, "QtlDataset")) {
    if (!identical(study, data@study))
      stop(".spStudyDataForm: study '", study,
           "' not in QtlDataset (study='", data@study, "')")
    return("individual")
  }
  if (is(data, "QtlSumStats")) {
    if (!(study %in% unique(as.character(data$study))))
      stop(".spStudyDataForm: study '", study, "' not in QtlSumStats")
    return("sumstats")
  }
  if (is(data, "MultiStudyQtlDataset")) {
    if (study %in% names(getQtlDatasets(data))) return("individual")
    ss <- getSumStats(data)
    if (!is.null(ss) && study %in% unique(as.character(ss$study)))
      return("sumstats")
    stop(".spStudyDataForm: study '", study, "' not in MultiStudyQtlDataset")
  }
  stop(".spStudyDataForm: unsupported class: ", class(data)[[1L]])
}

# Return character vector of contexts in `data` (across all studies when
# `study = NULL`, or for one study otherwise).
# @noRd
.spListContexts <- function(data, study = NULL) {
  if (is(data, "QtlDataset")) {
    if (!is.null(study) && !identical(study, data@study))
      return(character(0))
    return(names(data@phenotypes))
  }
  if (is(data, "QtlSumStats")) {
    if (is.null(study))
      return(unique(as.character(data$context)))
    return(unique(as.character(
      data$context[as.character(data$study) == study])))
  }
  if (is(data, "MultiStudyQtlDataset")) {
    indDatasets <- getQtlDatasets(data)
    ss <- getSumStats(data)
    if (is.null(study)) {
      out <- character(0)
      for (qd in indDatasets) out <- c(out, names(qd@phenotypes))
      if (!is.null(ss)) out <- c(out, unique(as.character(ss$context)))
      return(unique(out))
    }
    if (study %in% names(indDatasets))
      return(names(indDatasets[[study]]@phenotypes))
    if (!is.null(ss) && study %in% unique(as.character(ss$study)))
      return(unique(as.character(
        ss$context[as.character(ss$study) == study])))
    return(character(0))
  }
  stop(".spListContexts: unsupported class: ", class(data)[[1L]])
}

# Return character vector of traits in `data` (filtered by study and/or
# context when supplied).
# @noRd
.spListTraits <- function(data, study = NULL, context = NULL) {
  if (is(data, "QtlDataset")) {
    if (!is.null(study) && !identical(study, data@study))
      return(character(0))
    if (is.null(context))
      return(unique(unlist(lapply(data@phenotypes, rownames),
                           use.names = FALSE)))
    se <- data@phenotypes[[context]]
    if (is.null(se)) return(character(0))
    return(rownames(se))
  }
  if (is(data, "QtlSumStats")) {
    keep <- rep(TRUE, nrow(data))
    if (!is.null(study))   keep <- keep & as.character(data$study)   == study
    if (!is.null(context)) keep <- keep & as.character(data$context) == context
    return(unique(as.character(data$trait[keep])))
  }
  if (is(data, "MultiStudyQtlDataset")) {
    indDatasets <- getQtlDatasets(data)
    ss <- getSumStats(data)
    if (!is.null(study) && study %in% names(indDatasets))
      return(.spListTraits(indDatasets[[study]], context = context))
    if (!is.null(ss) &&
        (is.null(study) || study %in% unique(as.character(ss$study))))
      return(.spListTraits(ss, study = study, context = context))
    if (is.null(study)) {
      out <- character(0)
      for (qd in indDatasets) out <- c(out, .spListTraits(qd))
      if (!is.null(ss)) out <- c(out, .spListTraits(ss))
      return(unique(out))
    }
    return(character(0))
  }
  stop(".spListTraits: unsupported class: ", class(data)[[1L]])
}


# -----------------------------------------------------------------------------
# parseJointSpecification — normalize the user-supplied joint spec into a
# canonical list of `list(axes = <character>, scope = <named list or NULL>)`
# entries. Validates axes ⊂ {study, context, trait}, no per-spec duplicates,
# scope keys and values present in `data`.
# -----------------------------------------------------------------------------

.spValidJointAxes <- c("study", "context", "trait")

# @noRd
parseJointSpecification <- function(jointSpecification, data) {
  if (is.null(jointSpecification)) return(list())

  # Auto-wrap a top-level character vector as a single spec
  if (is.character(jointSpecification)) {
    jointSpecification <- list(jointSpecification)
  }
  if (!is.list(jointSpecification)) {
    stop("`jointSpecification` must be NULL, a character vector of axes, ",
         "or a list of joint specs.")
  }

  lapply(seq_along(jointSpecification), function(i) {
    spec <- jointSpecification[[i]]
    label <- sprintf("jointSpecification[[%d]]", i)
    if (is.character(spec)) {
      axes <- spec
      scope <- NULL
    } else if (is.list(spec)) {
      if (!"axes" %in% names(spec))
        stop(label, ": missing `axes` element")
      axes <- spec$axes
      scope <- spec$scope
      extras <- setdiff(names(spec), c("axes", "scope"))
      if (length(extras) > 0L)
        stop(label, ": unknown element(s): ",
             paste(extras, collapse = ", "))
    } else {
      stop(label, ": each spec must be a character vector or a named list ",
           "with `axes` (and optional `scope`)")
    }

    if (!is.character(axes) || length(axes) == 0L)
      stop(label, ": `axes` must be a non-empty character vector")
    badAxes <- setdiff(axes, .spValidJointAxes)
    if (length(badAxes) > 0L)
      stop(label, ": unknown axes: ",
           paste(badAxes, collapse = ", "),
           ". Valid axes: ", paste(.spValidJointAxes, collapse = ", "))
    if (anyDuplicated(axes))
      stop(label, ": duplicate axes in `axes`")

    if (!is.null(scope)) {
      if (!is.list(scope) || is.null(names(scope)) ||
          any(!nzchar(names(scope))))
        stop(label, ": `scope` must be a named list keyed by ",
             "study / context / trait")
      badKeys <- setdiff(names(scope), .spValidJointAxes)
      if (length(badKeys) > 0L)
        stop(label, ": unknown scope key(s): ",
             paste(badKeys, collapse = ", "))
      for (k in names(scope)) {
        v <- scope[[k]]
        if (!is.character(v) || length(v) == 0L)
          stop(label, ": scope$", k, " must be a non-empty character vector")
        available <- switch(k,
          study   = .spListStudies(data),
          context = .spListContexts(data),
          trait   = .spListTraits(data))
        missing <- setdiff(v, available)
        if (length(missing) > 0L)
          stop(label, ": scope$", k, " contains values not in data: ",
               paste(missing, collapse = ", "))
      }
    }
    list(axes = axes, scope = scope)
  })
}


# -----------------------------------------------------------------------------
# parseContexts — normalize the user-supplied `contexts` argument to a named
# list keyed by every study in `data`, with each entry the character vector
# of selected contexts. NULL input is preserved as NULL ("all contexts").
# -----------------------------------------------------------------------------

# @noRd
parseContexts <- function(contexts, data) {
  if (is.null(contexts)) return(NULL)
  studies <- .spListStudies(data)

  # Vector form: applied uniformly to every study; filter to each study's
  # availability, warn on missing.
  isPlainCharVec <- is.character(contexts) &&
    (is.null(names(contexts)) || all(names(contexts) == ""))
  if (isPlainCharVec) {
    if (length(contexts) == 0L)
      stop("`contexts` must be NULL or a non-empty character vector ",
           "(or named list).")
    out <- list()
    for (s in studies) {
      avail <- .spListContexts(data, s)
      missing <- setdiff(contexts, avail)
      if (length(missing) > 0L)
        warning(sprintf(
          "parseContexts: study '%s' is missing requested context(s): %s",
          s, paste(missing, collapse = ", ")))
      out[[s]] <- intersect(contexts, avail)
    }
    return(out)
  }

  # Named-list form: explicit per-study selection. Studies not in the list
  # default to all available contexts.
  if (is.list(contexts)) {
    if (is.null(names(contexts)) || any(!nzchar(names(contexts))))
      stop("`contexts` must be NULL, a character vector, or a named list ",
           "keyed by study.")
    badStudies <- setdiff(names(contexts), studies)
    if (length(badStudies) > 0L)
      stop("`contexts` references unknown studies: ",
           paste(badStudies, collapse = ", "))
    out <- list()
    for (s in studies) {
      avail <- .spListContexts(data, s)
      if (s %in% names(contexts)) {
        requested <- as.character(contexts[[s]])
        if (length(requested) == 0L)
          stop(sprintf("contexts[['%s']] must be a non-empty character vector",
                       s))
        missing <- setdiff(requested, avail)
        if (length(missing) > 0L)
          stop(sprintf(
            "contexts[['%s']] contains unknown contexts: %s",
            s, paste(missing, collapse = ", ")))
        out[[s]] <- requested
      } else {
        out[[s]] <- avail
      }
    }
    return(out)
  }

  stop("`contexts` must be NULL, a character vector, or a named list ",
       "keyed by study.")
}


# -----------------------------------------------------------------------------
# parseTraitIds — normalize the user-supplied `traitId` argument. Accepts a
# character vector (applied uniformly), a study-keyed list, or a doubly-
# nested study→context list. Returns NULL when input is NULL (= use all
# available traits). Validates IDs against `.spListTraits` lookups.
# -----------------------------------------------------------------------------

# @noRd
parseTraitIds <- function(traitId, data) {
  if (is.null(traitId)) return(NULL)
  studies <- .spListStudies(data)

  isPlainCharVec <- is.character(traitId) &&
    (is.null(names(traitId)) || all(names(traitId) == ""))
  if (isPlainCharVec) {
    if (length(traitId) == 0L)
      stop("`traitId` must be NULL or a non-empty character vector ",
           "(or named list).")
    return(as.character(traitId))
  }

  if (!is.list(traitId)) {
    stop("`traitId` must be NULL, a character vector, or a named list ",
         "keyed by study (optionally nested by context).")
  }
  if (is.null(names(traitId)) || any(!nzchar(names(traitId))))
    stop("`traitId` (list form) must be named by study.")
  badStudies <- setdiff(names(traitId), studies)
  if (length(badStudies) > 0L)
    stop("`traitId` references unknown studies: ",
         paste(badStudies, collapse = ", "))

  out <- list()
  for (s in names(traitId)) {
    val <- traitId[[s]]
    if (is.character(val)) {
      if (length(val) == 0L)
        stop(sprintf("traitId[['%s']] must be a non-empty character vector", s))
      availTraits <- .spListTraits(data, study = s)
      missing <- setdiff(val, availTraits)
      if (length(missing) > 0L)
        stop(sprintf("traitId[['%s']] contains unknown traits: %s",
                     s, paste(missing, collapse = ", ")))
      out[[s]] <- as.character(val)
    } else if (is.list(val)) {
      if (is.null(names(val)) || any(!nzchar(names(val))))
        stop(sprintf("traitId[['%s']] (list form) must be named by context", s))
      badContexts <- setdiff(names(val), .spListContexts(data, s))
      if (length(badContexts) > 0L)
        stop(sprintf("traitId[['%s']] references unknown contexts: %s",
                     s, paste(badContexts, collapse = ", ")))
      sub <- list()
      for (cx in names(val)) {
        v2 <- val[[cx]]
        if (!is.character(v2) || length(v2) == 0L)
          stop(sprintf("traitId[['%s']][['%s']] must be a non-empty character vector",
                       s, cx))
        availTraits <- .spListTraits(data, study = s, context = cx)
        missing <- setdiff(v2, availTraits)
        if (length(missing) > 0L)
          stop(sprintf("traitId[['%s']][['%s']] contains unknown traits: %s",
                       s, cx, paste(missing, collapse = ", ")))
        sub[[cx]] <- as.character(v2)
      }
      out[[s]] <- sub
    } else {
      stop(sprintf("traitId[['%s']] must be a character vector or a named list keyed by context",
                   s))
    }
  }
  out
}


# -----------------------------------------------------------------------------
# parseMethods — normalize and validate the `methods` argument with optional
# `sumStatsMethods` / `qtlDatasetMethods` overrides. Validates:
#   * mutual exclusivity (methods XOR split-by-data-form)
#   * nested list structure (vector OR named list at each level; never both)
#   * method names against the capability table
#   * multi-axis methods may NOT appear at per-context or per-trait levels
#   * mr.mash and mvsusie pipeline scope
#
# `caps` is the capability table for the pipeline (see `.fineMappingMethodCapabilities`
#   and `.twasMethodCapabilities`).
# `multivariateMethods` is the subset of tokens whose `multivariate = TRUE`;
#   used for per-context / per-trait placement rejection.
# `rejectedAtUser` is a character vector of tokens forbidden as user-requested
#   methods on this pipeline (e.g. "mrmash" in fineMapping, "mvsusie" in twas).
#
# Returns a list with components:
#   methods            (NULL if not given)
#   sumStatsMethods    (NULL if not given)
#   qtlDatasetMethods  (NULL if not given)
#   shape              "primary" if `methods` was given, "split" otherwise
# -----------------------------------------------------------------------------

# Walk a (possibly nested) methods spec and return the depth at which leaf
# vectors live: 1 = top-level vector, 2 = per-study, 3 = per-(study,context),
# 4 = per-(study,context,trait). Returns a tibble-like list of (path, vec).
# `levelNames` is c("study", "context", "trait"); the leaf level is where
# the vector lives.
# @noRd
.spWalkMethods <- function(spec, label = "methods", depth = 0L,
                           maxDepth = 3L, path = character(0)) {
  if (is.character(spec)) {
    return(list(list(depth = depth, path = path, methods = unique(spec))))
  }
  if (!is.list(spec))
    stop(label, ": every node must be a character vector or a named list ",
         "(got class '", class(spec)[[1L]], "')")
  if (depth >= maxDepth)
    stop(label, ": cannot nest below the trait level (depth ",
         maxDepth, " is the deepest a vector may appear at).")
  if (is.null(names(spec)) || any(!nzchar(names(spec))))
    stop(label, ": named-list nodes must have non-empty names at depth ",
         depth + 1L)
  if (length(spec) == 0L)
    stop(label, ": empty named list at depth ", depth + 1L)
  out <- list()
  for (nm in names(spec)) {
    out <- c(out, .spWalkMethods(spec[[nm]], label = label,
                                  depth = depth + 1L,
                                  maxDepth = maxDepth,
                                  path = c(path, nm)))
  }
  out
}

# @noRd
parseMethods <- function(methods,
                         sumStatsMethods   = NULL,
                         qtlDatasetMethods = NULL,
                         data,
                         caps,
                         multivariateMethods,
                         rejectedAtUser = character(0)) {
  primaryGiven <- !is.null(methods)
  splitGiven   <- !is.null(sumStatsMethods) || !is.null(qtlDatasetMethods)

  if (primaryGiven && splitGiven)
    stop("Use either `methods` or (`sumStatsMethods` + `qtlDatasetMethods`), ",
         "not both.")
  if (!primaryGiven && !splitGiven)
    stop("Specify `methods`, or both `sumStatsMethods` and ",
         "`qtlDatasetMethods`.")
  if (splitGiven) {
    if (is.null(sumStatsMethods) || is.null(qtlDatasetMethods))
      stop("`sumStatsMethods` and `qtlDatasetMethods` must be given together.")
    if (!is.character(sumStatsMethods) || length(sumStatsMethods) == 0L)
      stop("`sumStatsMethods` must be a non-empty character vector.")
    if (!is.character(qtlDatasetMethods) || length(qtlDatasetMethods) == 0L)
      stop("`qtlDatasetMethods` must be a non-empty character vector.")
  }

  validateLeafVec <- function(vec, label) {
    if (!is.character(vec) || length(vec) == 0L)
      stop(label, ": method vector must be a non-empty character vector")
    bad <- setdiff(vec, names(caps))
    if (length(bad) > 0L)
      stop(label, ": unknown method token(s): ",
           paste(bad, collapse = ", "),
           ". Known tokens: ", paste(names(caps), collapse = ", "))
    rejected <- intersect(vec, rejectedAtUser)
    if (length(rejected) > 0L)
      stop(label, ": method(s) cannot be user-requested on this pipeline: ",
           paste(rejected, collapse = ", "))
    invisible(NULL)
  }

  if (splitGiven) {
    validateLeafVec(sumStatsMethods,   "sumStatsMethods")
    validateLeafVec(qtlDatasetMethods, "qtlDatasetMethods")
  } else {
    walked <- .spWalkMethods(methods, label = "methods", maxDepth = 3L)
    studyNames <- .spListStudies(data)
    for (leaf in walked) {
      lab <- sprintf("methods[[%s]]",
                     paste0("'", leaf$path, "'", collapse = "$"))
      if (length(leaf$path) == 0L) lab <- "methods"
      validateLeafVec(leaf$methods, lab)
      # Multi-axis methods may not appear at per-context or per-trait levels.
      if (leaf$depth >= 2L) {
        bad <- intersect(leaf$methods, multivariateMethods)
        if (length(bad) > 0L)
          stop(lab,
               ": multi-axis method(s) ",
               paste(bad, collapse = ", "),
               " cannot be assigned at the ",
               c("per-study", "per-context", "per-trait")[[leaf$depth]],
               " level (multi-axis methods operate across axes).")
      }
      # Study-keyed nodes must reference valid studies
      if (leaf$depth >= 1L) {
        s <- leaf$path[[1L]]
        if (!(s %in% studyNames))
          stop(lab, ": unknown study '", s, "'")
      }
      # Context-keyed nodes must reference valid contexts for the study
      if (leaf$depth >= 2L) {
        s <- leaf$path[[1L]]; cx <- leaf$path[[2L]]
        avail <- .spListContexts(data, s)
        if (!(cx %in% avail))
          stop(lab, ": unknown context '", cx, "' for study '", s, "'")
      }
      # Trait-keyed nodes must reference valid traits for the (study, context)
      if (leaf$depth >= 3L) {
        s <- leaf$path[[1L]]; cx <- leaf$path[[2L]]; tr <- leaf$path[[3L]]
        avail <- .spListTraits(data, study = s, context = cx)
        if (!(tr %in% avail))
          stop(lab, ": unknown trait '", tr, "' for (study '", s,
               "', context '", cx, "')")
      }
    }
  }

  list(
    methods           = methods,
    sumStatsMethods   = sumStatsMethods,
    qtlDatasetMethods = qtlDatasetMethods,
    shape             = if (primaryGiven) "primary" else "split")
}


# -----------------------------------------------------------------------------
# validateMethodsVsJointSpec — cross-validation. A per-axis method assignment
# at or below the axis being jointed in any spec contradicts user intent.
# E.g. axes = "context" + per-context methods = contradiction. Joint flags
# operate on axes that haven't been pinned to per-axis methods.
#
# The rule: for each jointSpec, every axis in `axes` must NOT be a level at
# which the methods list nests. Concretely:
#   - "study"   in axes -> methods must not be a named list keyed by study
#                          (i.e. methods must be a top-level vector OR the
#                          split form). Per-study methods would mean different
#                          methods per study, incompatible with cross-study
#                          joints.
#   - "context" in axes -> no per-(study, context) nesting at any study.
#   - "trait"   in axes -> no per-(study, context, trait) nesting at any
#                          (study, context).
# -----------------------------------------------------------------------------

# @noRd
validateMethodsVsJointSpec <- function(methodsParsed, jointSpecParsed) {
  # Split-form methods are flat per-data-form vectors — nothing to check.
  if (methodsParsed$shape == "split") return(invisible(NULL))
  if (length(jointSpecParsed) == 0L) return(invisible(NULL))
  methods <- methodsParsed$methods
  if (is.character(methods)) return(invisible(NULL))  # top-level vector OK

  walked <- .spWalkMethods(methods, label = "methods", maxDepth = 3L)
  # depth observed at leaves; max depth in the spec reflects nesting level.
  maxDepth <- max(vapply(walked, function(L) L$depth, integer(1)))

  for (i in seq_along(jointSpecParsed)) {
    axes <- jointSpecParsed[[i]]$axes
    lab  <- sprintf("jointSpecification[[%d]]", i)
    if ("study" %in% axes && maxDepth >= 1L)
      stop(lab, ": `axes` includes 'study' but `methods` nests per-study; ",
           "remove per-study method assignment when joining over studies.")
    if ("context" %in% axes && maxDepth >= 2L)
      stop(lab, ": `axes` includes 'context' but `methods` nests per-context; ",
           "remove per-context method assignment when joining over contexts.")
    if ("trait" %in% axes && maxDepth >= 3L)
      stop(lab, ": `axes` includes 'trait' but `methods` nests per-trait; ",
           "remove per-trait method assignment when joining over traits.")
  }
  invisible(NULL)
}
