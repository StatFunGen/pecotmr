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

# =============================================================================
# Joint-specification dispatchers (merged from former R/jointDispatchers.R)
# =============================================================================


# =============================================================================
# Shared helpers
# =============================================================================

# Resolve which studies / contexts / traits participate in `spec` given
# `data`. Filters data scope through the spec's `scope` and any explicit
# pipeline-level `contexts` / `traitIds` arguments. Returns a list with
# `studies` (character), `contexts` (named list keyed by study), `traits`
# (named list keyed by study).
# @noRd
.fmResolveSpecScope <- function(spec, data, contexts = NULL,
                                traitIds = NULL) {
  scope <- spec$scope
  studies <- .spListStudies(data)
  if (!is.null(scope$study))
    studies <- intersect(studies, scope$study)

  contextsOut <- list()
  traitsOut <- list()
  for (s in studies) {
    ctxAvail <- .spListContexts(data, s)
    if (!is.null(scope$context))
      ctxAvail <- intersect(ctxAvail, scope$context)
    if (!is.null(contexts)) {
      if (is.list(contexts) && s %in% names(contexts))
        ctxAvail <- intersect(ctxAvail, contexts[[s]])
      else if (is.character(contexts))
        ctxAvail <- intersect(ctxAvail, contexts)
    }
    contextsOut[[s]] <- ctxAvail

    trAvail <- .spListTraits(data, study = s)
    if (!is.null(scope$trait))
      trAvail <- intersect(trAvail, scope$trait)
    if (!is.null(traitIds)) {
      if (is.character(traitIds))
        trAvail <- intersect(trAvail, traitIds)
      else if (is.list(traitIds) && s %in% names(traitIds)) {
        tv <- traitIds[[s]]
        if (is.character(tv)) trAvail <- intersect(trAvail, tv)
      }
    }
    traitsOut[[s]] <- trAvail
  }
  list(studies = studies, contexts = contextsOut, traits = traitsOut)
}


# Build a (variants × tupleRows) Z matrix from a QtlSumStats subset,
# requiring all rows to share an identical SNP order (the post-
# summaryStatsQc contract). Returns list(Z, nVec, variantIds).
# `errorLabel` is woven into the SNP-order error to identify the caller.
# @noRd
.buildJointSumstatZMatrix <- function(data, tupleRows, colLabels, errorLabel) {
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  firstDf <- getSumstatDf(data,
                           study   = studyCol[[tupleRows[[1L]]]],
                           context = contextCol[[tupleRows[[1L]]]],
                           trait   = traitCol[[tupleRows[[1L]]]],
                           require = c("SNP", "Z", "N"))
  variantIds <- firstDf$variant_id
  Z <- matrix(NA_real_, nrow = length(variantIds), ncol = length(tupleRows),
              dimnames = list(variantIds, colLabels))
  nVec <- numeric(length(tupleRows))
  for (kk in seq_along(tupleRows)) {
    i <- tupleRows[[kk]]
    d <- getSumstatDf(data,
                       study   = studyCol[[i]],
                       context = contextCol[[i]],
                       trait   = traitCol[[i]],
                       require = c("SNP", "Z", "N"))
    if (!identical(d$variant_id, variantIds))
      stop(sprintf("%s: every entry in a joint group must share an identical SNP order after summaryStatsQc().",
                   errorLabel))
    Z[, kk] <- d$z
    nVec[kk] <- stats::median(d$N, na.rm = TRUE)
  }
  list(Z = Z, nVec = nVec, variantIds = variantIds)
}


# Build a multi-context Y matrix for a single (study, trait) from an
# individual-level QtlDataset. Returns list(X, Y, perTraitContexts) or
# NULL when fewer than 2 contexts carry `tid` or the sample / complete-Y
# subset is too small to fit.
# @noRd
.buildIndividualCrossContextXY <- function(data, tid, scopedContexts,
                                           cisWindow, verbose, label,
                                           region = NULL) {
  perTraitContexts <- character(0)
  for (cx in scopedContexts) {
    se <- getPhenotypes(data, contexts = cx)
    if (tid %in% rownames(se))
      perTraitContexts <- c(perTraitContexts, cx)
  }
  if (length(perTraitContexts) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: trait '%s' present in %d scoped context(s); skipping.",
        label, tid, length(perTraitContexts)))
    return(NULL)
  }
  X <- if (is.null(region)) {
    .fmResidGeno(data, contexts = perTraitContexts, traitId = tid,
                 cisWindow = cisWindow)
  } else {
    .fmResidGeno(data, contexts = perTraitContexts, region = region)
  }
  Yres <- .fmResidPheno(
    data, contexts = perTraitContexts, traitId = tid)
  commonSamples <- Reduce(intersect,
    c(list(rownames(X)), lapply(Yres, rownames)))
  if (length(commonSamples) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: trait '%s' has too few shared samples across contexts; skipping.",
        label, tid))
    return(NULL)
  }
  X <- X[commonSamples, , drop = FALSE]
  Y <- do.call(cbind, lapply(perTraitContexts, function(cx) {
    ym <- Yres[[cx]][commonSamples, , drop = FALSE]
    colnames(ym) <- cx
    ym
  }))
  keep <- stats::complete.cases(Y)
  if (sum(keep) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: trait '%s' has too few complete-Y subjects; skipping.",
        label, tid))
    return(NULL)
  }
  list(X = X[keep, , drop = FALSE], Y = Y[keep, , drop = FALSE],
       perTraitContexts = perTraitContexts)
}


# Build a multi-trait Y matrix for a single (study, context) from an
# individual-level QtlDataset. Returns list(X, Y, traitsHere, se) or NULL
# when fewer than 2 traits live in the context or the sample / complete-Y
# subset is too small.
# @noRd
.buildIndividualCrossTraitXY <- function(data, cx, scopedTraits,
                                         cisWindow, verbose, label, study,
                                         region = NULL) {
  se <- getPhenotypes(data, contexts = cx)
  traitsHere <- intersect(scopedTraits, rownames(se))
  if (length(traitsHere) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: context '%s' (study '%s') has %d scoped trait(s); skipping.",
        label, cx, study, length(traitsHere)))
    return(NULL)
  }
  X <- if (is.null(region)) {
    .fmResidGeno(data, contexts = cx, traitId = traitsHere,
                 cisWindow = cisWindow)
  } else {
    .fmResidGeno(data, contexts = cx, region = region)
  }
  Y <- .fmResidPheno(
    data, contexts = cx, traitId = traitsHere)
  common <- intersect(rownames(X), rownames(Y))
  if (length(common) < 2L) return(NULL)
  X <- X[common, , drop = FALSE]; Y <- Y[common, , drop = FALSE]
  keep <- stats::complete.cases(Y)
  if (sum(keep) < 2L) return(NULL)
  list(X = X[keep, , drop = FALSE], Y = Y[keep, , drop = FALSE],
       traitsHere = traitsHere, se = se)
}


# Build a composed-axes (context, trait) X/Y for individual-level
# QtlDataset. Returns list(X, Y, tuples) or NULL.
# @noRd
.buildComposedIndividualXY <- function(data, scope, study, cisWindow,
                                       verbose, label, region = NULL) {
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  tuples <- list()
  for (cx in scopedContexts) {
    se <- getPhenotypes(data, contexts = cx)
    for (tid in intersect(scopedTraits, rownames(se))) {
      tuples[[length(tuples) + 1L]] <- list(context = cx, trait = tid)
    }
  }
  if (length(tuples) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: study '%s' has %d (context, trait) tuple(s) in scope; skipping.",
        label, study, length(tuples)))
    return(NULL)
  }
  allContexts <- unique(vapply(tuples, function(t) t$context, character(1L)))
  allTraits   <- unique(vapply(tuples, function(t) t$trait,   character(1L)))
  X <- if (is.null(region)) {
    .fmResidGeno(data, contexts = allContexts, traitId = allTraits,
                 cisWindow = cisWindow)
  } else {
    .fmResidGeno(data, contexts = allContexts, region = region)
  }
  YresList <- .fmResidPheno(
    data, contexts = allContexts, traitId = allTraits)
  if (length(allContexts) == 1L) YresList <- setNames(list(YresList), allContexts)
  commonSamples <- Reduce(intersect,
    c(list(rownames(X)), lapply(YresList, rownames)))
  if (length(commonSamples) < 2L) return(NULL)
  X <- X[commonSamples, , drop = FALSE]
  yCols <- list(); colLabels <- character(0)
  for (t in tuples) {
    ym <- YresList[[t$context]]
    if (!(t$trait %in% colnames(ym))) next
    col <- ym[commonSamples, t$trait, drop = FALSE]
    colnames(col) <- paste(t$context, t$trait, sep = ":")
    yCols[[length(yCols) + 1L]] <- col
    colLabels <- c(colLabels, paste(t$context, t$trait, sep = ":"))
  }
  if (length(yCols) < 2L) return(NULL)
  Y <- do.call(cbind, yCols)
  keep <- stats::complete.cases(Y)
  if (sum(keep) < 2L) return(NULL)
  list(X = X[keep, , drop = FALSE], Y = Y[keep, , drop = FALSE],
       tuples = tuples)
}


# Enumerate composed-axes row groups for a QtlSumStats input. Returns the
# list of (rowIdx) per group along with the per-axis identity columns
# needed to label the output row. Groups containing fewer than 2 rows
# are returned unfiltered; the caller decides whether to skip.
# @noRd
.enumerateComposedSumstatGroups <- function(spec, data, scope) {
  axes <- spec$axes
  complement <- setdiff(c("study", "context", "trait"), axes)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  inScope <- vapply(seq_len(nrow(data)), function(i) {
    s <- studyCol[i]; cx <- contextCol[i]; tr <- traitCol[i]
    (s %in% scope$studies) &&
    (cx %in% scope$contexts[[s]]) &&
    (tr %in% scope$traits[[s]])
  }, logical(1L))
  rowIdx <- which(inScope)
  if (length(rowIdx) == 0L) return(NULL)
  groupKey <- if (length(complement) == 0L) {
    rep("__all__", length(rowIdx))
  } else {
    do.call(paste, c(lapply(complement, function(a)
      switch(a, study = studyCol[rowIdx],
                context = contextCol[rowIdx],
                trait = traitCol[rowIdx])),
      sep = "||"))
  }
  groups <- split(rowIdx, groupKey)
  list(groups = groups, axes = axes,
       studyCol = studyCol, contextCol = contextCol, traitCol = traitCol)
}


# =============================================================================
# Fine-mapping dispatchers
# =============================================================================

# Cross-context joint dispatcher for QtlDataset. For each trait in scope
# with >= 2 contexts in scope, fits mvsusieR::mvsusie on the multi-column
# Y matrix and emits ONE result row with context = "joint" and
# jointContexts = "ctx1;ctx2;...".
# @noRd
.fmDispatchCrossContextQtlDataset <- function(spec, data, methods,
                                               contexts, traitIds,
                                               cisWindow,
                                               coverage, secondaryCoverage,
                                               signalCutoff, minAbsCorr,
                                               verbose,
                                               methodArgs = list(),
                                               region = NULL,
                                               twasWeights = NULL,
                                               dataDrivenPriorWeightsCutoff = 1e-10) {
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  if (length(scopedContexts) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext: study '%s' has %d context(s) in scope; skipping cross-context fits.",
        study, length(scopedContexts)))
    return(NULL)
  }

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (tid in scopedTraits) {
    xy <- .buildIndividualCrossContextXY(
      data, tid, scopedContexts, cisWindow, verbose,
      label = "jointCrossContext", region = region)
    if (is.null(xy)) next

    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext: fitting mvsusie for (study='%s', trait='%s') across contexts (%s) ...",
        study, tid, paste(xy$perTraitContexts, collapse = ", ")))
    # Reweighted prior from a cross-context mr.mash joint twas run, keyed on
    # (study, trait=tid, context="joint"); conditions are the contexts.
    mvPrior <- .buildMvsusieReweightedPrior(
      .fmLookupMrmashFit(twasWeights, study, tid, context = "joint"),
      colnames(xy$Y), dataDrivenPriorWeightsCutoff)
    mvBaseArgs <- list(
      X = xy$X, Y = xy$Y,
      prior_variance = mvPrior$priorVariance,
      coverage = coverage)
    if (!is.null(mvPrior$residualVariance))
      mvBaseArgs$residual_variance <- mvPrior$residualVariance
    fit <- do.call(fitMvsusie,
                   .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                    methodArgs[["mvsusie"]]))
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    entry <- .fmPostprocessOne(
      fit = fit, method = "mvsusie",
      dataX = xy$X, dataY = NULL,
      coverage = coverage,
      secondaryCoverage = secondaryCoverage,
      signalCutoff = signalCutoff,
      minAbsCorr = minAbsCorr,
      csInput = "X")
    rowStudy   <- c(rowStudy,   study)
    rowContext <- c(rowContext, "joint")
    rowTrait   <- c(rowTrait,   tid)
    rowMethod  <- c(rowMethod,  "mvsusie")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointContexts <- c(rowJointContexts,
                          paste(xy$perTraitContexts, collapse = ";"))
  }

  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = NULL)
}


# Cross-trait joint dispatcher for QtlDataset. Per (study, context), fits
# mvsusieR::mvsusie or fsusieR::susiF (when in `methods`) jointly across
# the scoped traits within that context. Emits ONE result row per
# (study, context, method) with trait = "joint" and jointTraits populated.
# @noRd
.fmDispatchCrossTraitQtlDataset <- function(spec, data, methods,
                                             contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose,
                                             methodArgs = list(),
                                             region = NULL,
                                             twasWeights = NULL,
                                             dataDrivenPriorWeightsCutoff = 1e-10) {
  jointMethods <- intersect(methods, c("mvsusie", "fsusie"))
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (cx in scopedContexts) {
    xy <- .buildIndividualCrossTraitXY(
      data, cx, scopedTraits, cisWindow, verbose,
      label = "jointCrossTrait", study = study, region = region)
    if (is.null(xy)) next

    for (mm in jointMethods) {
      if (verbose >= 1)
        message(sprintf(
          "jointCrossTrait: fitting %s for (study='%s', context='%s') across traits (%s) ...",
          mm, study, cx, paste(xy$traitsHere, collapse = ", ")))
      if (mm == "mvsusie") {
        # Reweighted prior from a cross-trait mr.mash joint twas run, keyed on
        # (study, trait="joint", context=cx); conditions are the traits.
        mvPrior <- .buildMvsusieReweightedPrior(
          .fmLookupMrmashFit(twasWeights, study, "joint", context = cx),
          colnames(xy$Y), dataDrivenPriorWeightsCutoff)
        mvBaseArgs <- list(
          X = xy$X, Y = xy$Y,
          prior_variance = mvPrior$priorVariance,
          coverage = coverage)
        if (!is.null(mvPrior$residualVariance))
          mvBaseArgs$residual_variance <- mvPrior$residualVariance
        fit <- do.call(fitMvsusie,
                       .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                        methodArgs[["mvsusie"]]))
        fit <- .setFinemappingFitClass(fit, "mvsusie")
        entry <- .fmPostprocessOne(
          fit = fit, method = "mvsusie",
          dataX = xy$X, dataY = NULL,
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          csInput = "X")
      } else {
        rr <- SummarizedExperiment::rowRanges(xy$se)
        ord <- match(colnames(xy$Y), rownames(xy$se))
        rr <- rr[ord]
        pos <- (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2
        fit <- do.call(fitFsusie,
                       .fmMergeUserArgs(list(X = xy$X, Y = xy$Y, pos = pos),
                                        "fsusie", methodArgs[["fsusie"]]))
        fit <- .setFinemappingFitClass(fit, "fsusie")
        entry <- .fmPostprocessOne(
          fit = fit, method = "fsusie",
          dataX = xy$X, dataY = NULL,
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          csInput = "fsusie")
      }
      rowStudy   <- c(rowStudy,   study)
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   "joint")
      rowMethod  <- c(rowMethod,  mm)
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointTraits <- c(rowJointTraits,
                          paste(xy$traitsHere, collapse = ";"))
    }
  }

  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = NULL)
}


# Cross-context joint dispatcher for QtlSumStats input. Groups the
# selected sumstats rows by (study, trait); each group with >= 2 contexts
# in scope produces one mvsusie_rss fit and one result row with context =
# "joint" and jointContexts populated.
# @noRd
.fmDispatchCrossContextQtlSumStats <- function(spec, data, methods,
                                                contexts, traitIds,
                                                coverage, secondaryCoverage,
                                                signalCutoff, minAbsCorr,
                                                verbose,
                                                methodArgs = list(),
                                                twasWeights = NULL,
                                                dataDrivenPriorWeightsCutoff = 1e-10) {
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    if (length(scopedContexts) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (QtlSumStats): study '%s' has %d context(s) in scope; skipping.",
          s, length(scopedContexts)))
      next
    }
    for (tid in scopedTraits) {
      tupleRows <- which(studyCol == s & traitCol == tid &
                         contextCol %in% scopedContexts)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossContext (QtlSumStats): (study='%s', trait='%s') has %d scoped context(s); skipping.",
            s, tid, length(tupleRows)))
        next
      }
      ctxNames <- contextCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, ctxNames,
        errorLabel = "jointCrossContext (QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (QtlSumStats): fitting mvsusie_rss for (study='%s', trait='%s', %d contexts) ...",
          s, tid, length(ctxNames)))
      # Reweighted prior from a cross-context mr.mash joint twas run, keyed on
      # (study, trait=tid, context="joint"); conditions are the contexts.
      mvPrior <- .buildMvsusieReweightedPrior(
        .fmLookupMrmashFit(twasWeights, s, tid, context = "joint"),
        colnames(jz$Z), dataDrivenPriorWeightsCutoff)
      mvBaseArgs <- list(
        Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
        prior_variance = mvPrior$priorVariance,
        coverage = coverage)
      if (!is.null(mvPrior$residualVariance))
        mvBaseArgs$residual_variance <- mvPrior$residualVariance
      fit <- do.call(fitMvsusieRss,
                     .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                      methodArgs[["mvsusie"]]))
      fit <- .setFinemappingFitClass(fit, "mvsusie")
      entry <- .fmPostprocessOne(
        fit = fit, method = "mvsusie",
        dataX = ldMat, dataY = NULL,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        csInput = "Xcorr")
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, "joint")
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mvsusie")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointContexts <- c(rowJointContexts,
                            paste(ctxNames, collapse = ";"))
    }
  }

  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = ldSketch)
}


# Cross-trait joint dispatcher for QtlSumStats: groups by (study, context),
# requires >= 2 scoped traits per group. mvsusie_rss only -- no RSS fsusie.
# @noRd
.fmDispatchCrossTraitQtlSumStats <- function(spec, data, methods,
                                              contexts, traitIds,
                                              coverage, secondaryCoverage,
                                              signalCutoff, minAbsCorr,
                                              verbose,
                                              methodArgs = list(),
                                              twasWeights = NULL,
                                              dataDrivenPriorWeightsCutoff = 1e-10) {
  if ("fsusie" %in% methods)
    stop("jointCrossTrait (QtlSumStats): fsusie has no RSS variant; ",
         "fsusie cannot participate in sumstats-based joint fits.")
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    for (cx in scopedContexts) {
      tupleRows <- which(studyCol == s & contextCol == cx &
                         traitCol %in% scopedTraits)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossTrait (QtlSumStats): (study='%s', context='%s') has %d scoped trait(s); skipping.",
            s, cx, length(tupleRows)))
        next
      }
      trNames <- traitCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, trNames,
        errorLabel = "jointCrossTrait (QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossTrait (QtlSumStats): fitting mvsusie_rss for (study='%s', context='%s', %d traits) ...",
          s, cx, length(trNames)))
      # Reweighted prior from a cross-trait mr.mash joint twas run, keyed on
      # (study, trait="joint", context=cx); conditions are the traits.
      mvPrior <- .buildMvsusieReweightedPrior(
        .fmLookupMrmashFit(twasWeights, s, "joint", context = cx),
        colnames(jz$Z), dataDrivenPriorWeightsCutoff)
      mvBaseArgs <- list(
        Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
        prior_variance = mvPrior$priorVariance,
        coverage = coverage)
      if (!is.null(mvPrior$residualVariance))
        mvBaseArgs$residual_variance <- mvPrior$residualVariance
      fit <- do.call(fitMvsusieRss,
                     .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                      methodArgs[["mvsusie"]]))
      fit <- .setFinemappingFitClass(fit, "mvsusie")
      entry <- .fmPostprocessOne(
        fit = fit, method = "mvsusie",
        dataX = ldMat, dataY = NULL,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        csInput = "Xcorr")
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   "joint")
      rowMethod  <- c(rowMethod,  "mvsusie")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointTraits <- c(rowJointTraits,
                          paste(trNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = ldSketch)
}


# Cross-study joint dispatcher for QtlSumStats: groups by (context, trait),
# requires >= 2 scoped studies per group. Sumstats-only by definition;
# individual-level studies are excluded with a message at the caller.
# mvsusie_rss only.
# @noRd
.fmDispatchCrossStudyQtlSumStats <- function(spec, data, methods,
                                              contexts, traitIds,
                                              coverage, secondaryCoverage,
                                              signalCutoff, minAbsCorr,
                                              verbose,
                                              methodArgs = list(),
                                              twasWeights = NULL,
                                              dataDrivenPriorWeightsCutoff = 1e-10) {
  if ("fsusie" %in% methods)
    stop("jointCrossStudy: fsusie cannot participate (no RSS variant).")
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  allCtxs <- unique(unlist(scope$contexts, use.names = FALSE))
  allTrs  <- unique(unlist(scope$traits,   use.names = FALSE))

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointStudies <- character(0)

  for (cx in allCtxs) {
    for (tid in allTrs) {
      tupleRows <- which(contextCol == cx & traitCol == tid &
                         studyCol %in% scope$studies)
      keep <- logical(length(tupleRows))
      for (k in seq_along(tupleRows)) {
        s <- studyCol[tupleRows[k]]
        keep[k] <- (cx %in% scope$contexts[[s]]) &&
                   (tid %in% scope$traits[[s]])
      }
      tupleRows <- tupleRows[keep]
      if (length(tupleRows) < 2L) {
        if (length(tupleRows) > 0L && verbose >= 1)
          message(sprintf(
            "jointCrossStudy: (context='%s', trait='%s') has %d study(ies) in scope; skipping.",
            cx, tid, length(tupleRows)))
        next
      }
      stNames <- studyCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, stNames,
        errorLabel = "jointCrossStudy")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossStudy: fitting mvsusie_rss for (context='%s', trait='%s', %d studies) ...",
          cx, tid, length(stNames)))
      mvBaseArgs <- list(
        Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
        # TODO(mvsusie-prior): cross-study lookup key undecided; canonical prior for now
        prior_variance = mvsusieR::create_mixture_prior(R = ncol(jz$Z)),
        coverage = coverage)
      fit <- do.call(fitMvsusieRss,
                     .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                      methodArgs[["mvsusie"]]))
      fit <- .setFinemappingFitClass(fit, "mvsusie")
      entry <- .fmPostprocessOne(
        fit = fit, method = "mvsusie",
        dataX = ldMat, dataY = NULL,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        csInput = "Xcorr")
      rowStudy   <- c(rowStudy,   "joint")
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mvsusie")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointStudies <- c(rowJointStudies,
                           paste(stNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study        = rowStudy,
    context      = rowContext,
    trait        = rowTrait,
    method       = rowMethod,
    entry        = rowEntries,
    jointStudies = rowJointStudies,
    ldSketch     = ldSketch)
}


# Composed multi-axis joint dispatcher for QtlDataset. Only axes =
# c("context", "trait") is meaningful for a single-study individual-
# level input. Iterates per (study) (just one), enumerates the
# (context, trait) tuples in scope where the trait exists in the
# context, and fits one mvsusie joint over those tuples.
# @noRd
.fmDispatchComposedQtlDataset <- function(spec, data, methods,
                                           contexts, traitIds, cisWindow,
                                           coverage, secondaryCoverage,
                                           signalCutoff, minAbsCorr,
                                           verbose,
                                           methodArgs = list(),
                                           region = NULL,
                                           twasWeights = NULL,
                                           dataDrivenPriorWeightsCutoff = 1e-10) {
  axes <- spec$axes
  if ("study" %in% axes)
    stop("composed jointSpecification (QtlDataset): axes including 'study' require sumstats input.")
  if (!setequal(axes, c("context", "trait")))
    stop(sprintf(
      "composed jointSpecification (QtlDataset): unsupported axes (%s) for individual-level input.",
      paste(axes, collapse = ", ")))
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  xy <- .buildComposedIndividualXY(data, scope, study, cisWindow,
                                    verbose,
                                    label = "composed joint (QtlDataset)",
                                    region = region)
  if (is.null(xy)) return(NULL)

  if (verbose >= 1)
    message(sprintf(
      "composed joint (QtlDataset): fitting mvsusie for study='%s' over %d (context, trait) columns ...",
      study, ncol(xy$Y)))
  # Reweighted prior from a composed mr.mash joint twas run, keyed on
  # (study, trait="joint", context="joint"); conditions are the (context,trait)
  # columns of the composed design.
  mvPrior <- .buildMvsusieReweightedPrior(
    .fmLookupMrmashFit(twasWeights, study, "joint", context = "joint"),
    colnames(xy$Y), dataDrivenPriorWeightsCutoff)
  mvBaseArgs <- list(
    X = xy$X, Y = xy$Y,
    prior_variance = mvPrior$priorVariance,
    coverage = coverage)
  if (!is.null(mvPrior$residualVariance))
    mvBaseArgs$residual_variance <- mvPrior$residualVariance
  fit <- do.call(fitMvsusie,
                 .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                  methodArgs[["mvsusie"]]))
  fit <- .setFinemappingFitClass(fit, "mvsusie")
  entry <- .fmPostprocessOne(
    fit = fit, method = "mvsusie",
    dataX = xy$X, dataY = NULL,
    coverage = coverage,
    secondaryCoverage = secondaryCoverage,
    signalCutoff = signalCutoff,
    minAbsCorr = minAbsCorr,
    csInput = "X")
  QtlFineMappingResult(
    study         = study,
    context       = "joint",
    trait         = "joint",
    method        = "mvsusie",
    entry         = list(entry),
    jointContexts = paste(vapply(xy$tuples, function(t) t$context,
                                  character(1)), collapse = ";"),
    jointTraits   = paste(vapply(xy$tuples, function(t) t$trait,
                                  character(1)), collapse = ";"),
    ldSketch      = NULL)
}


# Composed multi-axis joint dispatcher for QtlSumStats. Handles any
# `axes` subset of {study, context, trait} of size >= 2 by iterating the
# complement-axis Cartesian product and emitting one joint fit per
# iteration unit.
# @noRd
.fmDispatchComposedQtlSumStats <- function(spec, data, methods,
                                            contexts, traitIds,
                                            coverage, secondaryCoverage,
                                            signalCutoff, minAbsCorr,
                                            verbose,
                                            methodArgs = list(),
                                            twasWeights = NULL,
                                            dataDrivenPriorWeightsCutoff = 1e-10) {
  if ("fsusie" %in% methods)
    stop("composed jointSpecification (QtlSumStats): fsusie has no RSS variant.")
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  groupInfo <- .enumerateComposedSumstatGroups(spec, data, scope)
  if (is.null(groupInfo)) return(NULL)
  axes <- groupInfo$axes
  studyCol <- groupInfo$studyCol
  contextCol <- groupInfo$contextCol
  traitCol <- groupInfo$traitCol

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list()
  rowJointStudies  <- character(0)
  rowJointContexts <- character(0)
  rowJointTraits   <- character(0)

  for (gIdx in groupInfo$groups) {
    if (length(gIdx) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "composed joint (QtlSumStats): group has %d row(s); skipping.",
          length(gIdx)))
      next
    }
    colLabels <- vapply(gIdx, function(i)
      paste(studyCol[i], contextCol[i], traitCol[i], sep = ":"),
      character(1L))
    jz <- .buildJointSumstatZMatrix(
      data, gIdx, colLabels,
      errorLabel = "composed joint (QtlSumStats)")
    ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
    if (verbose >= 1)
      message(sprintf(
        "composed joint (QtlSumStats): fitting mvsusie_rss for axes=(%s), %d columns ...",
        paste(axes, collapse = ", "), length(gIdx)))
    # Reweighted prior from a composed mr.mash joint twas run, keyed on
    # (study, trait="joint", context="joint"); conditions are the joint columns.
    mvPrior <- .buildMvsusieReweightedPrior(
      .fmLookupMrmashFit(twasWeights, studyCol[gIdx[[1L]]], "joint",
                         context = "joint"),
      colnames(jz$Z), dataDrivenPriorWeightsCutoff)
    mvBaseArgs <- list(
      Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
      prior_variance = mvPrior$priorVariance,
      coverage = coverage)
    if (!is.null(mvPrior$residualVariance))
      mvBaseArgs$residual_variance <- mvPrior$residualVariance
    fit <- do.call(fitMvsusieRss,
                   .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                    methodArgs[["mvsusie"]]))
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    entry <- .fmPostprocessOne(
      fit = fit, method = "mvsusie",
      dataX = ldMat, dataY = NULL,
      coverage = coverage,
      secondaryCoverage = secondaryCoverage,
      signalCutoff = signalCutoff,
      minAbsCorr = minAbsCorr,
      csInput = "Xcorr")

    repStudy   <- if ("study"   %in% axes) "joint" else studyCol[gIdx[[1L]]]
    repContext <- if ("context" %in% axes) "joint" else contextCol[gIdx[[1L]]]
    repTrait   <- if ("trait"   %in% axes) "joint" else traitCol[gIdx[[1L]]]
    rowStudy   <- c(rowStudy,   repStudy)
    rowContext <- c(rowContext, repContext)
    rowTrait   <- c(rowTrait,   repTrait)
    rowMethod  <- c(rowMethod,  "mvsusie")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointStudies <- c(rowJointStudies,
      if ("study" %in% axes) paste(studyCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointContexts <- c(rowJointContexts,
      if ("context" %in% axes) paste(contextCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointTraits <- c(rowJointTraits,
      if ("trait" %in% axes) paste(traitCol[gIdx], collapse = ";")
      else NA_character_)
  }

  if (length(rowStudy) == 0L) return(NULL)
  jsArg <- if (all(is.na(rowJointStudies))) NULL else rowJointStudies
  jcArg <- if (all(is.na(rowJointContexts))) NULL else rowJointContexts
  jtArg <- if (all(is.na(rowJointTraits))) NULL else rowJointTraits
  QtlFineMappingResult(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointStudies  = jsArg,
    jointContexts = jcArg,
    jointTraits   = jtArg,
    ldSketch      = ldSketch)
}


# Top-level joint dispatcher for fineMappingPipeline(QtlDataset).
# @noRd
# Merge per-region QtlFineMappingResult collections (same keys across regions)
# into one by merging each (study, context, trait, method) row's entry via
# .fmMergeEntries (per-region susieFit list + renumbered credible sets).
# @noRd
.fmMergeResultsByKey <- function(results) {
  base <- results[[1L]]
  n <- nrow(base)
  if (n == 0L) return(base)
  keyOf <- function(r) paste(as.character(r$study), as.character(r$context),
                             as.character(r$trait), as.character(r$method),
                             sep = "\r")
  baseKeys <- keyOf(base)
  mergedEntries <- lapply(seq_len(n), function(i) {
    perRegion <- lapply(results, function(r) {
      hit <- which(keyOf(r) == baseKeys[[i]])
      if (length(hit)) r$entry[[hit[[1L]]]] else NULL
    })
    .fmMergeEntries(Filter(Negate(is.null), perRegion))
  })
  QtlFineMappingResult(
    study = as.character(base$study), context = as.character(base$context),
    trait = as.character(base$trait), method = as.character(base$method),
    entry = mergedEntries,
    jointStudies  = if ("jointStudies"  %in% names(base)) base$jointStudies  else NULL,
    jointContexts = if ("jointContexts" %in% names(base)) base$jointContexts else NULL,
    jointTraits   = if ("jointTraits"   %in% names(base)) base$jointTraits   else NULL,
    ldSketch = NULL)
}

.fmDispatchJointSpecsQtlDataset <- function(parsedJointSpec, data,
                                             methods, contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose,
                                             methodArgs = list(),
                                             xRegions = list(NULL),
                                             twasWeights = NULL,
                                             dataDrivenPriorWeightsCutoff = 1e-10) {
  # Run the joint dispatch once per region block, then merge per
  # (study, context, trait, method) across regions. A single block (cis or
  # jointRegions=TRUE concatenated) returns its result directly.
  perRegion <- lapply(xRegions, function(rg) {
    .fmDispatchJointSpecsQtlDatasetOneRegion(
      parsedJointSpec, data, methods, contexts, traitIds, cisWindow,
      coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
      methodArgs = methodArgs, region = rg,
      twasWeights = twasWeights,
      dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff)
  })
  perRegion <- Filter(Negate(is.null), perRegion)
  if (length(perRegion) == 0L) return(NULL)
  if (length(perRegion) == 1L) return(perRegion[[1L]])
  .fmMergeResultsByKey(perRegion)
}

.fmDispatchJointSpecsQtlDatasetOneRegion <- function(parsedJointSpec, data,
                                             methods, contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose,
                                             methodArgs = list(),
                                             region = NULL,
                                             twasWeights = NULL,
                                             dataDrivenPriorWeightsCutoff = 1e-10) {
  # Bundle the data-driven mvSuSiE prior pass-through args once; every leaf
  # dispatcher accepts the same pair.
  priorArgs <- list(twasWeights = twasWeights,
                    dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff)
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- do.call(.fmDispatchComposedQtlDataset, c(list(
        spec, data, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs, region = region), priorArgs))
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindFineMappingResult(out, res, ldSketch = NULL)
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = do.call(.fmDispatchCrossContextQtlDataset, c(list(
        spec, data, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs, region = region), priorArgs)),
      trait = do.call(.fmDispatchCrossTraitQtlDataset, c(list(
        spec, data, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs, region = region), priorArgs)),
      study = stop(
        "fineMappingPipeline(QtlDataset): jointSpecification with axes = 'study' requires sumstats input. ",
        "QtlDataset represents a single individual-level study; cross-study joints operate on the sumstats slot of MultiStudyQtlDataset or on QtlSumStats directly."),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindFineMappingResult(out, res, ldSketch = NULL)
  }
  out
}


# Top-level joint dispatcher for fineMappingPipeline(QtlSumStats).
# @noRd
.fmDispatchJointSpecsQtlSumStats <- function(parsedJointSpec, data,
                                              methods, contexts, traitIds,
                                              coverage, secondaryCoverage,
                                              signalCutoff, minAbsCorr,
                                              verbose,
                                              methodArgs = list(),
                                              twasWeights = NULL,
                                              dataDrivenPriorWeightsCutoff = 1e-10) {
  # Bundle the data-driven mvSuSiE prior pass-through args once; every leaf
  # dispatcher accepts the same pair.
  priorArgs <- list(twasWeights = twasWeights,
                    dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff)
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- do.call(.fmDispatchComposedQtlSumStats, c(list(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs), priorArgs))
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindFineMappingResult(out, res, ldSketch = getLdSketch(data))
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = do.call(.fmDispatchCrossContextQtlSumStats, c(list(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs), priorArgs)),
      trait = do.call(.fmDispatchCrossTraitQtlSumStats, c(list(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs), priorArgs)),
      study = do.call(.fmDispatchCrossStudyQtlSumStats, c(list(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs), priorArgs)),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindFineMappingResult(out, res,
                                          ldSketch = getLdSketch(data))
  }
  out
}


# Top-level joint dispatcher for fineMappingPipeline(MultiStudyQtlDataset).
# Routes per-component AND per-axis: a spec with `axes = "study"` only
# touches the sumStats slot; `axes = "context"` and `axes = "trait"` run
# on every component.
# @noRd
.fmDispatchJointSpecsMultiStudy <- function(parsedJointSpec, data,
                                             methods, contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose,
                                             methodArgs = list(),
                                             xRegions = list(NULL),
                                             twasWeights = NULL,
                                             dataDrivenPriorWeightsCutoff = 1e-10) {
  out <- NULL
  embeddedLd <- NULL
  qtlDatasets <- getQtlDatasets(data)
  sumStats <- getSumStats(data)

  studyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) "study" %in% s$axes, logical(1L))]
  nonStudyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) !("study" %in% s$axes), logical(1L))]

  if (length(studyAxisSpecs) > 0L && length(qtlDatasets) > 0L && verbose >= 1) {
    message(sprintf(
      "jointCrossStudy: excluding individual-level studies (%s) from cross-study fits (no LD sketch available); sumstats studies participate.",
      paste(names(qtlDatasets), collapse = ", ")))
  }

  if (length(nonStudyAxisSpecs) > 0L) {
    for (qdName in names(qtlDatasets)) {
      qd <- qtlDatasets[[qdName]]
      qdRes <- .fmDispatchJointSpecsQtlDataset(
        nonStudyAxisSpecs, qd, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs, xRegions = xRegions,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff)
      if (!is.null(qdRes))
        out <- if (is.null(out)) qdRes
               else .rbindFineMappingResult(out, qdRes, ldSketch = NULL)
    }
  }

  if (!is.null(sumStats)) {
    ssRes <- .fmDispatchJointSpecsQtlSumStats(
      parsedJointSpec, sumStats, methods, contexts, traitIds,
      coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
      methodArgs = methodArgs,
      twasWeights = twasWeights,
      dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff)
    if (!is.null(ssRes)) {
      embeddedLd <- getLdSketch(ssRes)
      out <- if (is.null(out)) ssRes
             else .rbindFineMappingResult(out, ssRes,
                                          ldSketch = embeddedLd)
    }
  } else if (length(studyAxisSpecs) > 0L && verbose >= 1) {
    message("jointCrossStudy: no sumStats slot present on this MultiStudyQtlDataset; cross-study specs produce no result.")
  }
  out
}


# =============================================================================
# TWAS-weights dispatchers
# =============================================================================

# Cross-context joint dispatcher for QtlDataset (twas). Mr.mash across
# scoped contexts per (study, trait).
# @noRd
.twasDispatchCrossContextQtlDataset <- function(spec, data, methods,
                                                 contexts, traitIds,
                                                 cisWindow, dataType,
                                                 verbose, region = NULL,
                                                 retainFit = TRUE,
                                                 retainFitDetail = "slim") {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  if (length(scopedContexts) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext (twas QtlDataset): study '%s' has %d context(s) in scope; skipping.",
        study, length(scopedContexts)))
    return(NULL)
  }

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (tid in scopedTraits) {
    xy <- .buildIndividualCrossContextXY(
      data, tid, scopedContexts, cisWindow, verbose,
      label = "jointCrossContext (twas QtlDataset)", region = region)
    if (is.null(xy)) next

    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext (twas QtlDataset): fitting mr.mash for (study='%s', trait='%s') across contexts (%s) ...",
        study, tid, paste(xy$perTraitContexts, collapse = ", ")))
    weights <- mrmashWeights(X = xy$X, Y = xy$Y,
                             retainFit = retainFit, fitDetail = retainFitDetail)
    if (is.null(rownames(weights))) rownames(weights) <- colnames(xy$X)
    entry <- TwasWeightsEntry(
      variantIds   = rownames(weights),
      weights      = weights,
      fits         = attr(weights, "fit"),
      standardized = FALSE,
      dataType     = dataType)
    rowStudy   <- c(rowStudy,   study)
    rowContext <- c(rowContext, "joint")
    rowTrait   <- c(rowTrait,   tid)
    rowMethod  <- c(rowMethod,  "mrmash")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointContexts <- c(rowJointContexts,
                          paste(xy$perTraitContexts, collapse = ";"))
  }

  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = NULL)
}


# Cross-trait joint dispatcher for QtlDataset (twas). Mr.mash per
# (study, context) across scoped traits.
# @noRd
.twasDispatchCrossTraitQtlDataset <- function(spec, data, methods,
                                               contexts, traitIds,
                                               cisWindow, dataType,
                                               verbose, region = NULL,
                                               retainFit = TRUE,
                                               retainFitDetail = "slim") {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (cx in scopedContexts) {
    xy <- .buildIndividualCrossTraitXY(
      data, cx, scopedTraits, cisWindow, verbose,
      label = "jointCrossTrait (twas)", study = study, region = region)
    if (is.null(xy)) next

    if (verbose >= 1)
      message(sprintf(
        "jointCrossTrait (twas): fitting mr.mash for (study='%s', context='%s') across traits (%s) ...",
        study, cx, paste(xy$traitsHere, collapse = ", ")))
    weights <- mrmashWeights(X = xy$X, Y = xy$Y,
                             retainFit = retainFit, fitDetail = retainFitDetail)
    if (is.null(rownames(weights))) rownames(weights) <- colnames(xy$X)
    entry <- TwasWeightsEntry(
      variantIds   = rownames(weights),
      weights      = weights,
      fits         = attr(weights, "fit"),
      standardized = FALSE,
      dataType     = dataType)
    rowStudy   <- c(rowStudy,   study)
    rowContext <- c(rowContext, cx)
    rowTrait   <- c(rowTrait,   "joint")
    rowMethod  <- c(rowMethod,  "mrmash")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointTraits <- c(rowJointTraits,
                        paste(xy$traitsHere, collapse = ";"))
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = NULL)
}


# Cross-context joint dispatcher for QtlSumStats (twas). Mr.mash.rss per
# (study, trait).
# @noRd
.twasDispatchCrossContextQtlSumStats <- function(spec, data, methods,
                                                  contexts, traitIds,
                                                  dataType, verbose,
                                                  retainFit = TRUE,
                                                  retainFitDetail = "slim") {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    if (length(scopedContexts) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (twas QtlSumStats): study '%s' has %d context(s) in scope; skipping.",
          s, length(scopedContexts)))
      next
    }
    for (tid in scopedTraits) {
      tupleRows <- which(studyCol == s & traitCol == tid &
                         contextCol %in% scopedContexts)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossContext (twas QtlSumStats): (study='%s', trait='%s') has %d scoped context(s); skipping.",
            s, tid, length(tupleRows)))
        next
      }
      ctxNames <- contextCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, ctxNames,
        errorLabel = "jointCrossContext (twas QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      stat <- list(z = jz$Z, N = jz$nVec)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (twas QtlSumStats): fitting mr.mash.rss for (study='%s', trait='%s', %d contexts) ...",
          s, tid, length(ctxNames)))
      weights <- mrmashRssWeights(stat = stat, LD = ldMat,
                                  retainFit = retainFit,
                                  fitDetail = retainFitDetail)
      if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
      entry <- TwasWeightsEntry(
        variantIds   = rownames(weights),
        weights      = weights,
        fits         = attr(weights, "fit"),
        standardized = TRUE,
        dataType     = dataType)
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, "joint")
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mrmash")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointContexts <- c(rowJointContexts,
                            paste(ctxNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = ldSketch)
}


# Cross-trait joint dispatcher for QtlSumStats (twas). Mr.mash.rss per
# (study, context).
# @noRd
.twasDispatchCrossTraitQtlSumStats <- function(spec, data, methods,
                                                contexts, traitIds,
                                                dataType, verbose,
                                                retainFit = TRUE,
                                                retainFitDetail = "slim") {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    for (cx in scopedContexts) {
      tupleRows <- which(studyCol == s & contextCol == cx &
                         traitCol %in% scopedTraits)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossTrait (twas QtlSumStats): (study='%s', context='%s') has %d scoped trait(s); skipping.",
            s, cx, length(tupleRows)))
        next
      }
      trNames <- traitCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, trNames,
        errorLabel = "jointCrossTrait (twas QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      stat <- list(z = jz$Z, N = jz$nVec)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossTrait (twas QtlSumStats): fitting mr.mash.rss for (study='%s', context='%s', %d traits) ...",
          s, cx, length(trNames)))
      weights <- mrmashRssWeights(stat = stat, LD = ldMat,
                                  retainFit = retainFit,
                                  fitDetail = retainFitDetail)
      if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
      entry <- TwasWeightsEntry(
        variantIds   = rownames(weights),
        weights      = weights,
        fits         = attr(weights, "fit"),
        standardized = TRUE,
        dataType     = dataType)
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   "joint")
      rowMethod  <- c(rowMethod,  "mrmash")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointTraits <- c(rowJointTraits,
                          paste(trNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = ldSketch)
}


# Cross-study joint dispatcher for QtlSumStats (twas). Mr.mash.rss per
# (context, trait).
# @noRd
.twasDispatchCrossStudyQtlSumStats <- function(spec, data, methods,
                                                contexts, traitIds,
                                                dataType, verbose,
                                                retainFit = TRUE,
                                                retainFitDetail = "slim") {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  allCtxs <- unique(unlist(scope$contexts, use.names = FALSE))
  allTrs  <- unique(unlist(scope$traits,   use.names = FALSE))

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointStudies <- character(0)

  for (cx in allCtxs) {
    for (tid in allTrs) {
      tupleRows <- which(contextCol == cx & traitCol == tid &
                         studyCol %in% scope$studies)
      keep <- logical(length(tupleRows))
      for (k in seq_along(tupleRows)) {
        s <- studyCol[tupleRows[k]]
        keep[k] <- (cx %in% scope$contexts[[s]]) &&
                   (tid %in% scope$traits[[s]])
      }
      tupleRows <- tupleRows[keep]
      if (length(tupleRows) < 2L) {
        if (length(tupleRows) > 0L && verbose >= 1)
          message(sprintf(
            "jointCrossStudy (twas): (context='%s', trait='%s') has %d study(ies) in scope; skipping.",
            cx, tid, length(tupleRows)))
        next
      }
      stNames <- studyCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, stNames,
        errorLabel = "jointCrossStudy (twas)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      stat <- list(z = jz$Z, N = jz$nVec)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossStudy (twas): fitting mr.mash.rss for (context='%s', trait='%s', %d studies) ...",
          cx, tid, length(stNames)))
      weights <- mrmashRssWeights(stat = stat, LD = ldMat,
                                  retainFit = retainFit,
                                  fitDetail = retainFitDetail)
      if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
      entry <- TwasWeightsEntry(
        variantIds   = rownames(weights),
        weights      = weights,
        fits         = attr(weights, "fit"),
        standardized = TRUE,
        dataType     = dataType)
      rowStudy   <- c(rowStudy,   "joint")
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mrmash")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointStudies <- c(rowJointStudies,
                           paste(stNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study        = rowStudy,
    context      = rowContext,
    trait        = rowTrait,
    method       = rowMethod,
    entry        = rowEntries,
    jointStudies = rowJointStudies,
    ldSketch     = ldSketch)
}


# Composed multi-axis joint dispatcher for QtlSumStats (twas).
# @noRd
.twasDispatchComposedQtlSumStats <- function(spec, data, methods,
                                              contexts, traitIds,
                                              dataType, verbose,
                                              retainFit = TRUE,
                                              retainFitDetail = "slim") {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  groupInfo <- .enumerateComposedSumstatGroups(spec, data, scope)
  if (is.null(groupInfo)) return(NULL)
  axes <- groupInfo$axes
  studyCol <- groupInfo$studyCol
  contextCol <- groupInfo$contextCol
  traitCol <- groupInfo$traitCol

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list()
  rowJointStudies  <- character(0)
  rowJointContexts <- character(0)
  rowJointTraits   <- character(0)

  for (gIdx in groupInfo$groups) {
    if (length(gIdx) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "composed joint (twas QtlSumStats): group has %d row(s); skipping.",
          length(gIdx)))
      next
    }
    colLabels <- vapply(gIdx, function(i)
      paste(studyCol[i], contextCol[i], traitCol[i], sep = ":"),
      character(1L))
    jz <- .buildJointSumstatZMatrix(
      data, gIdx, colLabels,
      errorLabel = "composed joint (twas QtlSumStats)")
    ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
    stat <- list(z = jz$Z, N = jz$nVec)
    if (verbose >= 1)
      message(sprintf(
        "composed joint (twas QtlSumStats): fitting mr.mash.rss for axes=(%s), %d columns ...",
        paste(axes, collapse = ", "), length(gIdx)))
    weights <- mrmashRssWeights(stat = stat, LD = ldMat,
                                retainFit = retainFit,
                                fitDetail = retainFitDetail)
    if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
    entry <- TwasWeightsEntry(
      variantIds   = rownames(weights),
      weights      = weights,
      fits         = attr(weights, "fit"),
      standardized = TRUE,
      dataType     = dataType)

    repStudy   <- if ("study"   %in% axes) "joint" else studyCol[gIdx[[1L]]]
    repContext <- if ("context" %in% axes) "joint" else contextCol[gIdx[[1L]]]
    repTrait   <- if ("trait"   %in% axes) "joint" else traitCol[gIdx[[1L]]]
    rowStudy   <- c(rowStudy,   repStudy)
    rowContext <- c(rowContext, repContext)
    rowTrait   <- c(rowTrait,   repTrait)
    rowMethod  <- c(rowMethod,  "mrmash")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointStudies <- c(rowJointStudies,
      if ("study" %in% axes) paste(studyCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointContexts <- c(rowJointContexts,
      if ("context" %in% axes) paste(contextCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointTraits <- c(rowJointTraits,
      if ("trait" %in% axes) paste(traitCol[gIdx], collapse = ";")
      else NA_character_)
  }
  if (length(rowStudy) == 0L) return(NULL)
  jsArg <- if (all(is.na(rowJointStudies))) NULL else rowJointStudies
  jcArg <- if (all(is.na(rowJointContexts))) NULL else rowJointContexts
  jtArg <- if (all(is.na(rowJointTraits))) NULL else rowJointTraits
  TwasWeights(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointStudies  = jsArg,
    jointContexts = jcArg,
    jointTraits   = jtArg,
    ldSketch      = ldSketch)
}


# Composed multi-axis joint dispatcher for QtlDataset (twas). axes =
# c("context", "trait") only.
# @noRd
.twasDispatchComposedQtlDataset <- function(spec, data, methods,
                                             contexts, traitIds, cisWindow,
                                             dataType, verbose,
                                             region = NULL,
                                             retainFit = TRUE,
                                             retainFitDetail = "slim") {
  axes <- spec$axes
  if ("study" %in% axes)
    stop("composed jointSpecification (twas QtlDataset): axes including 'study' require sumstats input.")
  if (!setequal(axes, c("context", "trait")))
    stop(sprintf("composed jointSpecification (twas QtlDataset): unsupported axes (%s) for individual-level input.",
                 paste(axes, collapse = ", ")))
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  xy <- .buildComposedIndividualXY(data, scope, study, cisWindow,
                                    verbose,
                                    label = "composed joint (twas QtlDataset)",
                                    region = region)
  if (is.null(xy)) return(NULL)

  if (verbose >= 1)
    message(sprintf(
      "composed joint (twas QtlDataset): fitting mr.mash for study='%s' over %d (context, trait) columns ...",
      study, ncol(xy$Y)))
  weights <- mrmashWeights(X = xy$X, Y = xy$Y,
                           retainFit = retainFit, fitDetail = retainFitDetail)
  if (is.null(rownames(weights))) rownames(weights) <- colnames(xy$X)
  entry <- TwasWeightsEntry(
    variantIds   = rownames(weights),
    weights      = weights,
    fits         = attr(weights, "fit"),
    standardized = FALSE,
    dataType     = dataType)
  TwasWeights(
    study         = study,
    context       = "joint",
    trait         = "joint",
    method        = "mrmash",
    entry         = list(entry),
    jointContexts = paste(vapply(xy$tuples, function(t) t$context,
                                  character(1)), collapse = ";"),
    jointTraits   = paste(vapply(xy$tuples, function(t) t$trait,
                                  character(1)), collapse = ";"),
    ldSketch      = NULL)
}


# Top-level joint dispatcher for twasWeightsPipeline(QtlDataset).
# @noRd
# Merge per-region TwasWeights collections (same keys across regions) into one
# by concatenating each (study, context, trait, method) row's entry via
# .twasMergeRegionEntries (stacked weights + flat per-region cvPerformance).
# @noRd
.twasMergeResultsByKey <- function(results, regionLabels) {
  base <- results[[1L]]
  n <- length(base$method)
  if (n == 0L) return(base)
  keyOf <- function(r) paste(as.character(r$study), as.character(r$context),
                             as.character(r$trait), as.character(r$method),
                             sep = "\r")
  baseKeys <- keyOf(base)
  mergedEntries <- lapply(seq_len(n), function(i) {
    perRegion <- lapply(results, function(r) {
      hit <- which(keyOf(r) == baseKeys[[i]])
      if (length(hit)) r$entry[[hit[[1L]]]] else NULL
    })
    keep <- !vapply(perRegion, is.null, logical(1))
    .twasMergeRegionEntries(perRegion[keep], regionLabels[keep])
  })
  TwasWeights(
    study = as.character(base$study), context = as.character(base$context),
    trait = as.character(base$trait), method = as.character(base$method),
    entry = mergedEntries)
}

.twasDispatchJointSpecsQtlDataset <- function(parsedJointSpec, data,
                                               methods, contexts, traitIds,
                                               cisWindow, dataType,
                                               verbose, xRegions = list(NULL),
                                               retainFit = TRUE,
                                               retainFitDetail = "slim") {
  # Run the joint dispatch once per region block, then merge per
  # (study, context, trait, method) across regions. A single block (cis or
  # jointRegions=TRUE concatenated) returns its result directly.
  perRegion <- lapply(xRegions, function(rg) {
    .twasDispatchJointSpecsQtlDatasetOneRegion(
      parsedJointSpec, data, methods, contexts, traitIds, cisWindow, dataType,
      verbose, region = rg,
      retainFit = retainFit, retainFitDetail = retainFitDetail)
  })
  labs <- vapply(xRegions, .twasRegionLabel, character(1))
  keep <- !vapply(perRegion, is.null, logical(1))
  perRegion <- perRegion[keep]; labs <- labs[keep]
  if (length(perRegion) == 0L) return(NULL)
  if (length(perRegion) == 1L) return(perRegion[[1L]])
  .twasMergeResultsByKey(perRegion, labs)
}

.twasDispatchJointSpecsQtlDatasetOneRegion <- function(parsedJointSpec, data,
                                               methods, contexts, traitIds,
                                               cisWindow, dataType,
                                               verbose, region = NULL,
                                               retainFit = TRUE,
                                               retainFitDetail = "slim") {
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- .twasDispatchComposedQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow, dataType, verbose,
        region = region,
        retainFit = retainFit, retainFitDetail = retainFitDetail)
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindTwasWeights(out, res, ldSketch = NULL)
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = .twasDispatchCrossContextQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow, dataType, verbose,
        region = region,
        retainFit = retainFit, retainFitDetail = retainFitDetail),
      trait = .twasDispatchCrossTraitQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow, dataType, verbose,
        region = region,
        retainFit = retainFit, retainFitDetail = retainFitDetail),
      study = stop(
        "twasWeightsPipeline(QtlDataset): jointSpecification with axes = 'study' requires sumstats input."),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindTwasWeights(out, res, ldSketch = NULL)
  }
  out
}


# Top-level joint dispatcher for twasWeightsPipeline(QtlSumStats).
# @noRd
.twasDispatchJointSpecsQtlSumStats <- function(parsedJointSpec, data,
                                                methods, contexts, traitIds,
                                                dataType, verbose,
                                                retainFit = TRUE,
                                                retainFitDetail = "slim") {
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- .twasDispatchComposedQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose,
        retainFit = retainFit, retainFitDetail = retainFitDetail)
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindTwasWeights(out, res, ldSketch = getLdSketch(data))
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = .twasDispatchCrossContextQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose,
        retainFit = retainFit, retainFitDetail = retainFitDetail),
      trait = .twasDispatchCrossTraitQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose,
        retainFit = retainFit, retainFitDetail = retainFitDetail),
      study = .twasDispatchCrossStudyQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose,
        retainFit = retainFit, retainFitDetail = retainFitDetail),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindTwasWeights(out, res, ldSketch = getLdSketch(data))
  }
  out
}


# Top-level joint dispatcher for twasWeightsPipeline(MultiStudyQtlDataset).
# @noRd
.twasDispatchJointSpecsMultiStudy <- function(parsedJointSpec, data,
                                               methods, contexts, traitIds,
                                               cisWindow, dataType, verbose,
                                               xRegions = list(NULL),
                                               retainFit = TRUE,
                                               retainFitDetail = "slim") {
  out <- NULL
  embeddedLd <- NULL
  qtlDatasets <- getQtlDatasets(data)
  sumStats <- getSumStats(data)

  studyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) "study" %in% s$axes, logical(1L))]
  nonStudyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) !("study" %in% s$axes), logical(1L))]

  if (length(studyAxisSpecs) > 0L && length(qtlDatasets) > 0L && verbose >= 1) {
    message(sprintf(
      "jointCrossStudy (twas): excluding individual-level studies (%s) from cross-study fits; sumstats studies participate.",
      paste(names(qtlDatasets), collapse = ", ")))
  }

  if (length(nonStudyAxisSpecs) > 0L) {
    for (qdName in names(qtlDatasets)) {
      qd <- qtlDatasets[[qdName]]
      qdRes <- .twasDispatchJointSpecsQtlDataset(
        nonStudyAxisSpecs, qd, methods, contexts, traitIds, cisWindow,
        dataType, verbose, xRegions = xRegions,
        retainFit = retainFit, retainFitDetail = retainFitDetail)
      if (!is.null(qdRes))
        out <- if (is.null(out)) qdRes
               else .rbindTwasWeights(out, qdRes, ldSketch = NULL)
    }
  }

  if (!is.null(sumStats)) {
    ssRes <- .twasDispatchJointSpecsQtlSumStats(
      parsedJointSpec, sumStats, methods, contexts, traitIds, dataType,
      verbose,
      retainFit = retainFit, retainFitDetail = retainFitDetail)
    if (!is.null(ssRes)) {
      embeddedLd <- getLdSketch(ssRes)
      out <- if (is.null(out)) ssRes
             else .rbindTwasWeights(out, ssRes, ldSketch = embeddedLd)
    }
  } else if (length(studyAxisSpecs) > 0L && verbose >= 1) {
    message("jointCrossStudy (twas): no sumStats slot present on this MultiStudyQtlDataset; cross-study specs produce no result.")
  }
  out
}
