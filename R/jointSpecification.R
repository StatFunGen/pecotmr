# =============================================================================
# Joint-specification grammar and ragged input-argument parsing for the
# fineMapping / twasWeights pipelines. Pure validation + normalization;
# no pipeline dispatch or fits live here.
# =============================================================================

# -----------------------------------------------------------------------------
# Internal scope helpers -- what (study, context, trait, dataForm) tuples does
# the input cover? `dataForm` is "individual" for QtlDataset-located studies
# and "sumstats" for QtlSumStats-located studies.
# -----------------------------------------------------------------------------

# Return character vector of all studies present in `data`.
# @noRd
.spListStudies <- function(data) {
    if (is(data, "QtlDataset")) {
        return(getStudy(data))
    }
    if (is(data, "QtlSumStats")) {
        return(unique(as.character(data$study)))
    }
    if (is(data, "MultiStudyQtlDataset")) {
        indStudies <- names(getQtlDatasets(data))
        ss <- getSumStats(data)
        ssStudies <- if (is.null(ss)) {
            character(0)
        } else {
            unique(as.character(ss$study))
        }
        return(unique(c(indStudies, ssStudies)))
    }
    cls <- class(data)[[1L]]
    msg <- glue(".spListStudies: unsupported class: {cls}")
    abort(msg)
}

# Return "individual" or "sumstats" for a single study in `data`. Errors if
# the study is not present.
# @noRd
.spStudyDataForm <- function(data, study) {
    if (is(data, "QtlDataset")) {
        if (!identical(study, getStudy(data))) {
            dataStudy <- getStudy(data)
            msg <- glue(
                ".spStudyDataForm: study '{study}' not in QtlDataset ",
                "(study='{dataStudy}')"
            )
            abort(msg)
        }
        return("individual")
    }
    if (is(data, "QtlSumStats")) {
        if (!is_in(study, unique(as.character(data$study)))) {
            msg <- glue(".spStudyDataForm: study '{study}' not in QtlSumStats")
            abort(msg)
        }
        return("sumstats")
    }
    if (is(data, "MultiStudyQtlDataset")) {
        if (is_in(study, names(getQtlDatasets(data)))) {
            return("individual")
        }
        ss <- getSumStats(data)
        if (!is.null(ss) && is_in(study, unique(as.character(ss$study)))) {
            return("sumstats")
        }
        msg <- glue(
            ".spStudyDataForm: study '{study}' not in MultiStudyQtlDataset"
        )
        abort(msg)
    }
    cls <- class(data)[[1L]]
    msg <- glue(".spStudyDataForm: unsupported class: {cls}")
    abort(msg)
}

# Return character vector of contexts in `data` (across all studies when
# `study = NULL`, or for one study otherwise).
# @noRd
# Concatenate a list of pieces into one, the way `c()` did, but safe on the
# empty case (`list_c()` has no zero-length identity to return).
# @noRd
.spConcat <- function(pieces, empty = list()) {
    if (length(pieces) == 0L) {
        return(empty)
    }
    list_c(pieces)
}

.spListContexts <- function(data, study = NULL) {
    if (is(data, "QtlDataset")) {
        if (!is.null(study) && !identical(study, getStudy(data))) {
            return(character(0))
        }
        return(getContexts(data))
    }
    if (is(data, "QtlSumStats")) {
        if (is.null(study)) {
            return(unique(as.character(data$context)))
        }
        return(unique(as.character(
            data$context[as.character(data$study) == study]
        )))
    }
    if (is(data, "MultiStudyQtlDataset")) {
        indDatasets <- getQtlDatasets(data)
        ss <- getSumStats(data)
        if (is.null(study)) {
            fromInd <- .spConcat(
                map(indDatasets, getContexts),
                empty = character(0)
            )
            fromSs <- if (is.null(ss)) {
                character(0)
            } else {
                unique(as.character(ss$context))
            }
            return(unique(c(fromInd, fromSs)))
        }
        if (is_in(study, names(indDatasets))) {
            return(getContexts(indDatasets[[study]]))
        }
        if (!is.null(ss) && is_in(study, unique(as.character(ss$study)))) {
            return(unique(as.character(
                ss$context[as.character(ss$study) == study]
            )))
        }
        return(character(0))
    }
    cls <- class(data)[[1L]]
    msg <- glue(".spListContexts: unsupported class: {cls}")
    abort(msg)
}

# Return character vector of traits in `data` (filtered by study and/or
# context when supplied).
# @noRd
# The trait ids present in one context of a QtlDataset.
# @noRd
.spTraitsInContext <- function(context, data) {
    rownames(getPhenotypes(data, context))
}

# Traits available in a single individual-level QtlDataset (optionally scoped).
.spListTraitsQtlDataset <- function(data, study, context) {
    if (!is.null(study) && !identical(study, getStudy(data))) {
        return(character(0))
    }
    if (is.null(context)) {
        return(unique(unname(list_c(
            map(getContexts(data), .spTraitsInContext, data = data)
        ))))
    }
    # Checked here rather than left to the accessor: `.spListTraits` answers
    # "which traits are in this scope", and an absent context is an empty
    # answer, not an error. The raw slot read this replaced returned NULL for
    # an unknown context; getPhenotypes() rejects one, so the empty case has to
    # be stated instead of falling out of a NULL.
    if (!is_in(context, getContexts(data))) {
        return(character(0))
    }
    rownames(getPhenotypes(data, context))
}

# Traits across a MultiStudyQtlDataset (individual studies + sumstats).
.spListTraitsMultiStudy <- function(data, study, context) {
    indDatasets <- getQtlDatasets(data)
    ss <- getSumStats(data)
    # No study filter: aggregate across every component (individual + sumstats).
    # Must precede per-study branches so a present sumStats slot does not shadow
    # the individual-level studies' traits when study = NULL.
    if (is.null(study)) {
        fromInd <- .spConcat(
            map(indDatasets, .spListTraits, context = context),
            empty = character(0)
        )
        fromSs <- if (is.null(ss)) {
            character(0)
        } else {
            .spListTraits(ss, context = context)
        }
        return(unique(c(fromInd, fromSs)))
    }
    if (is_in(study, names(indDatasets))) {
        return(.spListTraits(indDatasets[[study]], context = context))
    }
    if (!is.null(ss) && is_in(study, unique(as.character(ss$study)))) {
        return(.spListTraits(ss, study = study, context = context))
    }
    character(0)
}

.spListTraits <- function(data, study = NULL, context = NULL) {
    if (is(data, "QtlDataset")) {
        return(.spListTraitsQtlDataset(data, study, context))
    }
    if (is(data, "QtlSumStats")) {
        byStudy <- if (is.null(study)) {
            rep(TRUE, nrow(data))
        } else {
            as.character(data$study) == study
        }
        byContext <- if (is.null(context)) {
            TRUE
        } else {
            as.character(data$context) == context
        }
        keep <- byStudy & byContext
        return(unique(as.character(data$trait[keep])))
    }
    if (is(data, "MultiStudyQtlDataset")) {
        return(.spListTraitsMultiStudy(data, study, context))
    }
    cls <- class(data)[[1L]]
    msg <- glue(".spListTraits: unsupported class: {cls}")
    abort(msg)
}


# -----------------------------------------------------------------------------
# parseJointSpecification -- normalize the user-supplied joint spec into a
# canonical list of `list(axes = <character>, scope = <named list or NULL>)`
# entries. Validates axes subset of {study, context, trait}, no per-spec
# duplicates,
# scope keys and values present in `data`.
# -----------------------------------------------------------------------------

.spValidJointAxes <- c("study", "context", "trait")

# @noRd
# --- parseJointSpecification helpers ----------------------------------------

# Extract (axes, scope) from a spec (character vector or named list).
.parseJointSpecExtract <- function(spec, label) {
    if (is.character(spec)) {
        return(list(axes = spec, scope = NULL))
    }
    if (is.list(spec)) {
        if (!is_in("axes", names(spec))) {
            msg <- glue("{label}: missing `axes` element")
            abort(msg)
        }
        extras <- setdiff(names(spec), c("axes", "scope"))
        if (length(extras) > 0L) {
            extraStr <- str_flatten(extras, ", ")
            msg <- glue("{label}: unknown element(s): {extraStr}")
            abort(msg)
        }
        return(list(axes = spec$axes, scope = spec$scope))
    }
    msg <- glue(
        "{label}: each spec must be a character vector or a named list ",
        "with `axes` (and optional `scope`)"
    )
    abort(msg)
}

# Axes must be a non-empty, duplicate-free vector drawn from the valid axes.
.parseJointSpecCheckAxes <- function(axes, label) {
    if (!is.character(axes) || length(axes) == 0L) {
        msg <- glue("{label}: `axes` must be a non-empty character vector")
        abort(msg)
    }
    badAxes <- setdiff(axes, .spValidJointAxes)
    if (length(badAxes) > 0L) {
        badStr <- str_flatten(badAxes, ", ")
        validStr <- str_flatten(.spValidJointAxes, ", ")
        msg <- glue(
            "{label}: unknown axes: {badStr}. Valid axes: {validStr}"
        )
        abort(msg)
    }
    if (n_distinct(axes) < length(axes)) {
        msg <- glue("{label}: duplicate axes in `axes`")
        abort(msg)
    }
}

# One scope entry must be a non-empty vector of values present in the data.
.parseJointSpecCheckScopeKey <- function(k, v, label, data) {
    if (!is.character(v) || length(v) == 0L) {
        msg <- glue(
            "{label}: scope${k} must be a non-empty character vector"
        )
        abort(msg)
    }
    available <- switch(
        k,
        study = .spListStudies(data),
        context = .spListContexts(data),
        trait = .spListTraits(data)
    )
    missing <- setdiff(v, available)
    if (length(missing) > 0L) {
        missingStr <- str_flatten(missing, ", ")
        msg <- glue(
            "{label}: scope${k} contains values not in data: {missingStr}"
        )
        abort(msg)
    }
}

# Scope (optional) must be a named list keyed by study/context/trait.
.parseJointSpecCheckScope <- function(scope, label, data) {
    if (is.null(scope)) {
        return(invisible(NULL))
    }
    if (
        !is.list(scope) ||
            is.null(names(scope)) ||
            any(str_length(names(scope)) == 0L, na.rm = TRUE)
    ) {
        msg <- glue(
            "{label}: `scope` must be a named list keyed by ",
            "study / context / trait"
        )
        abort(msg)
    }
    badKeys <- setdiff(names(scope), .spValidJointAxes)
    if (length(badKeys) > 0L) {
        keyStr <- str_flatten(badKeys, ", ")
        msg <- glue("{label}: unknown scope key(s): {keyStr}")
        abort(msg)
    }
    for (k in names(scope)) {
        .parseJointSpecCheckScopeKey(k, scope[[k]], label, data)
    }
}

# Validate one joint specification entry.
.parseOneJointSpec <- function(spec, i, data) {
    label <- glue("jointSpecification[[{i}]]")
    extracted <- .parseJointSpecExtract(spec, label)
    axes <- extracted$axes
    scope <- extracted$scope
    .parseJointSpecCheckAxes(axes, label)
    .parseJointSpecCheckScope(scope, label, data)
    list(axes = axes, scope = scope)
}

parseJointSpecification <- function(jointSpecification, data) {
    if (is.null(jointSpecification)) {
        return(list())
    }
    # Auto-wrap a top-level character vector as a single spec.
    if (is.character(jointSpecification)) {
        jointSpecification <- list(jointSpecification)
    }
    if (!is.list(jointSpecification)) {
        msg <- glue(
            "`jointSpecification` must be NULL, a character vector of axes, ",
            "or a list of joint specs."
        )
        abort(msg)
    }
    map(
        seq_along(jointSpecification),
        .parseJointSpecAt,
        jointSpecification = jointSpecification,
        data = data
    )
}


# -----------------------------------------------------------------------------
# parseContexts -- normalize the user-supplied `contexts` argument to a named
# list keyed by every study in `data`, with each entry the character vector
# of selected contexts. NULL input is preserved as NULL ("all contexts").
# -----------------------------------------------------------------------------

# @noRd
# --- parseContexts helpers --------------------------------------------------

# Vector form: apply uniformly to every study, filtering to availability.
.parseContextsVec <- function(contexts, studies, data) {
    if (length(contexts) == 0L) {
        msg <- glue(
            "`contexts` must be NULL or a non-empty character vector ",
            "(or named list)."
        )
        abort(msg)
    }
    set_names(
        map(studies, .parseContextsAvail, data = data, contexts = contexts),
        studies
    )
}

# One study's requested contexts, warning about any it does not have.
# @noRd
.parseContextsAvail <- function(s, data, contexts) {
    avail <- .spListContexts(data, s)
    missing <- setdiff(contexts, avail)
    if (length(missing) > 0L) {
        missingStr <- str_flatten(missing, ", ")
        msg <- glue(
            "parseContexts: study '{s}' is missing requested ",
            "context(s): {missingStr}"
        )
        warn(msg)
    }
    intersect(contexts, avail)
}

# One study's contexts from the named-list form: its own entry when given,
# otherwise everything it has.
# @noRd
.parseContextsListAt <- function(s, data, contexts) {
    avail <- .spListContexts(data, s)
    if (is_in(s, names(contexts))) {
        .parseContextsStudy(as.character(contexts[[s]]), s, avail)
    } else {
        avail
    }
}

# Validate one study's explicitly-requested contexts against availability.
.parseContextsStudy <- function(requested, s, avail) {
    if (length(requested) == 0L) {
        msg <- glue("contexts[['{s}']] must be a non-empty character vector")
        abort(msg)
    }
    missing <- setdiff(requested, avail)
    if (length(missing) > 0L) {
        missingStr <- str_flatten(missing, ", ")
        msg <- glue("contexts[['{s}']] contains unknown contexts: {missingStr}")
        abort(msg)
    }
    requested
}

# Named-list form: explicit per-study selection; unlisted studies get all.
.parseContextsList <- function(contexts, studies, data) {
    ctxNm <- names(contexts)
    if (is.null(ctxNm) || any(str_length(ctxNm) == 0L, na.rm = TRUE)) {
        msg <- glue(
            "`contexts` must be NULL, a character vector, or a named list ",
            "keyed by study."
        )
        abort(msg)
    }
    badStudies <- setdiff(names(contexts), studies)
    if (length(badStudies) > 0L) {
        badStr <- str_flatten(badStudies, ", ")
        msg <- glue("`contexts` references unknown studies: {badStr}")
        abort(msg)
    }
    set_names(
        map(studies, .parseContextsListAt, data = data, contexts = contexts),
        studies
    )
}

parseContexts <- function(contexts, data) {
    if (is.null(contexts)) {
        return(NULL)
    }
    studies <- .spListStudies(data)
    isPlainCharVec <- is.character(contexts) &&
        (is.null(names(contexts)) || all(names(contexts) == ""))
    if (isPlainCharVec) {
        return(.parseContextsVec(contexts, studies, data))
    }
    if (is.list(contexts)) {
        return(.parseContextsList(contexts, studies, data))
    }
    msg <- glue(
        "`contexts` must be NULL, a character vector, or a named list ",
        "keyed by study."
    )
    abort(msg)
}


# -----------------------------------------------------------------------------
# parseTraitIds -- normalize the user-supplied `traitId` argument. Accepts a
# character vector (applied uniformly), a study-keyed list, or a doubly-
# nested study->context list. Returns NULL when input is NULL (= use all
# available traits). Validates IDs against `.spListTraits` lookups.
# -----------------------------------------------------------------------------

# @noRd
# --- parseTraitIds helpers --------------------------------------------------

# Validate a per-(study, context) trait vector against the data.
.parseTraitIdContext <- function(v2, s, cx, data) {
    if (!is.character(v2) || length(v2) == 0L) {
        msg <- glue(
            "traitId[['{s}']][['{cx}']] must be a non-empty ",
            "character vector"
        )
        abort(msg)
    }
    missing <- setdiff(v2, .spListTraits(data, study = s, context = cx))
    if (length(missing) > 0L) {
        missingStr <- str_flatten(missing, ", ")
        msg <- glue(
            "traitId[['{s}']][['{cx}']] contains unknown traits: ",
            "{missingStr}"
        )
        abort(msg)
    }
    as.character(v2)
}

# Validate a per-study character vector of traits.
.parseTraitIdStudyChar <- function(val, s, data) {
    if (length(val) == 0L) {
        msg <- glue("traitId[['{s}']] must be a non-empty character vector")
        abort(msg)
    }
    missing <- setdiff(val, .spListTraits(data, study = s))
    if (length(missing) > 0L) {
        missingStr <- str_flatten(missing, ", ")
        msg <- glue("traitId[['{s}']] contains unknown traits: {missingStr}")
        abort(msg)
    }
    as.character(val)
}

# Validate a per-study context-keyed list of trait vectors.
.parseTraitIdStudyList <- function(val, s, data) {
    valNm <- names(val)
    if (is.null(valNm) || any(str_length(valNm) == 0L, na.rm = TRUE)) {
        msg <- glue("traitId[['{s}']] (list form) must be named by context")
        abort(msg)
    }
    badContexts <- setdiff(names(val), .spListContexts(data, s))
    if (length(badContexts) > 0L) {
        badStr <- str_flatten(badContexts, ", ")
        msg <- glue("traitId[['{s}']] references unknown contexts: {badStr}")
        abort(msg)
    }
    set_names(
        map(names(val), .parseTraitIdContextAt, val = val, s = s, data = data),
        names(val)
    )
}

# One context's trait spec within a study.
# @noRd
.parseTraitIdContextAt <- function(cx, val, s, data) {
    .parseTraitIdContext(val[[cx]], s, cx, data)
}

# Dispatch one study's trait spec (character vector or context-keyed list).
.parseTraitIdStudy <- function(val, s, data) {
    if (is.character(val)) {
        return(.parseTraitIdStudyChar(val, s, data))
    }
    if (is.list(val)) {
        return(.parseTraitIdStudyList(val, s, data))
    }
    msg <- glue(
        "traitId[['{s}']] must be a character vector or a named ",
        "list keyed by context"
    )
    abort(msg)
}

parseTraitIds <- function(traitId, data) {
    if (is.null(traitId)) {
        return(NULL)
    }
    studies <- .spListStudies(data)
    isPlainCharVec <- is.character(traitId) &&
        (is.null(names(traitId)) || all(names(traitId) == ""))
    if (isPlainCharVec) {
        if (length(traitId) == 0L) {
            msg <- glue(
                "`traitId` must be NULL or a non-empty character vector ",
                "(or named list)."
            )
            abort(msg)
        }
        return(as.character(traitId))
    }
    if (!is.list(traitId)) {
        msg <- glue(
            "`traitId` must be NULL, a character vector, or a named list ",
            "keyed by study (optionally nested by context)."
        )
        abort(msg)
    }
    trNm <- names(traitId)
    if (is.null(trNm) || any(str_length(trNm) == 0L, na.rm = TRUE)) {
        abort("`traitId` (list form) must be named by study.")
    }
    badStudies <- setdiff(names(traitId), studies)
    if (length(badStudies) > 0L) {
        badStr <- str_flatten(badStudies, ", ")
        msg <- glue("`traitId` references unknown studies: {badStr}")
        abort(msg)
    }
    set_names(
        map(
            names(traitId),
            .parseTraitIdStudyAt,
            traitId = traitId,
            data = data
        ),
        names(traitId)
    )
}


# One study's parsed trait spec.
# @noRd
.parseTraitIdStudyAt <- function(s, traitId, data) {
    .parseTraitIdStudy(traitId[[s]], s, data)
}

# -----------------------------------------------------------------------------
# parseMethods -- normalize and validate the `methods` argument with optional
# `sumStatsMethods` / `qtlDatasetMethods` overrides. Validates:
#   * mutual exclusivity (methods XOR split-by-data-form)
#   * nested list structure (vector OR named list at each level; never both)
#   * method names against the capability table
#   * multi-axis methods may NOT appear at per-context or per-trait levels
#   * mr.mash and mvsusie pipeline scope
#
# `caps` is the capability table for the pipeline (see
# `.fineMappingMethodCapabilities` and `.twasMethodCapabilities`).
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
# Validate a named-list method node before recursing into it.
.spWalkValidate <- function(spec, label, depth, maxDepth) {
    if (!is.list(spec)) {
        cls <- class(spec)[[1L]]
        msg <- glue(
            "{label}: every node must be a character vector or a named ",
            "list (got class '{cls}')"
        )
        abort(msg)
    }
    if (depth >= maxDepth) {
        msg <- glue(
            "{label}: cannot nest below the trait level (depth ",
            "{maxDepth} is the deepest a vector may appear at)."
        )
        abort(msg)
    }
    specNm <- names(spec)
    if (is.null(specNm) || any(str_length(specNm) == 0L, na.rm = TRUE)) {
        d <- depth + 1L
        msg <- glue(
            "{label}: named-list nodes must have non-empty names at ",
            "depth {d}"
        )
        abort(msg)
    }
    if (length(spec) == 0L) {
        d <- depth + 1L
        msg <- glue("{label}: empty named list at depth {d}")
        abort(msg)
    }
}

.spWalkMethods <- function(
    spec,
    label = "methods",
    depth = 0L,
    maxDepth = 3L,
    path = character(0)
) {
    if (is.character(spec)) {
        return(list(list(depth = depth, path = path, methods = unique(spec))))
    }
    .spWalkValidate(spec, label, depth, maxDepth)
    out <- .spConcat(map(
        names(spec),
        .spWalkMethodsAt,
        spec = spec,
        label = label,
        depth = depth,
        maxDepth = maxDepth,
        path = path
    ))
    out
}

# Walk one named branch of a methods spec.
# @noRd
.spWalkMethodsAt <- function(nm, spec, label, depth, maxDepth, path) {
    .spWalkMethods(
        spec[[nm]],
        label = label,
        depth = depth + 1L,
        maxDepth = maxDepth,
        path = c(path, nm)
    )
}

# Validate one leaf method vector: non-empty character, all tokens known (in
# `caps`), and none in `rejectedAtUser`.
# @noRd
.jointValidateLeafVec <- function(vec, label, caps, rejectedAtUser) {
    if (!is.character(vec) || length(vec) == 0L) {
        msg <- glue(
            "{label}: method vector must be a non-empty character vector"
        )
        abort(msg)
    }
    bad <- setdiff(vec, names(caps))
    if (length(bad) > 0L) {
        badStr <- str_flatten(bad, ", ")
        knownStr <- str_flatten(names(caps), ", ")
        msg <- glue(
            "{label}: unknown method token(s): {badStr}. ",
            "Known tokens: {knownStr}"
        )
        abort(msg)
    }
    rejected <- intersect(vec, rejectedAtUser)
    if (length(rejected) > 0L) {
        rejectedStr <- str_flatten(rejected, ", ")
        msg <- glue(
            "{label}: method(s) cannot be user-requested on this ",
            "pipeline: {rejectedStr}"
        )
        abort(msg)
    }
    invisible(NULL)
}

# @noRd
# --- parseMethods helpers ---------------------------------------------------

# Validate mutual exclusivity of primary vs split method specs.
#' @importFrom checkmate assertCharacter
.parseMethodsValidateArgs <- function(
    primaryGiven,
    splitGiven,
    sumStatsMethods,
    qtlDatasetMethods
) {
    if (primaryGiven && splitGiven) {
        msg <- glue(
            "Use either `methods` or (`sumStatsMethods` + ",
            "`qtlDatasetMethods`), not both."
        )
        abort(msg)
    }
    if (!primaryGiven && !splitGiven) {
        msg <- glue(
            "Specify `methods`, or both `sumStatsMethods` and ",
            "`qtlDatasetMethods`."
        )
        abort(msg)
    }
    if (splitGiven) {
        if (is.null(sumStatsMethods) || is.null(qtlDatasetMethods)) {
            msg <- glue(
                "`sumStatsMethods` and `qtlDatasetMethods` must be given ",
                "together."
            )
            abort(msg)
        }
        assertCharacter(sumStatsMethods, min.len = 1L)
        assertCharacter(qtlDatasetMethods, min.len = 1L)
    }
}

# Multi-axis methods may not appear below the per-study level.
.parseMethodsCheckMultiAxis <- function(leaf, lab, multivariateMethods) {
    if (leaf$depth < 2L) {
        return(invisible(NULL))
    }
    bad <- intersect(leaf$methods, multivariateMethods)
    if (length(bad) > 0L) {
        badStr <- str_flatten(bad, ", ")
        levelName <- c("per-study", "per-context", "per-trait")[[leaf$depth]]
        msg <- glue(
            "{lab}: multi-axis method(s) {badStr} cannot be assigned at ",
            "the {levelName} level (multi-axis methods operate across ",
            "axes)."
        )
        abort(msg)
    }
}

# Study/context/trait keys of a method leaf must reference valid entities.
.parseMethodsCheckPath <- function(leaf, lab, data, studyNames) {
    # The path is (study, context, trait); each level is only checked once the
    # leaf is deep enough to name it -- a shallower leaf has no such entry to
    # read, so the depth guard has to come before the extraction.
    if (leaf$depth < 1L) {
        return(invisible(NULL))
    }
    s <- leaf$path[[1L]]
    if (!is_in(s, studyNames)) {
        msg <- glue("{lab}: unknown study '{s}'")
        abort(msg)
    }
    if (leaf$depth < 2L) {
        return(invisible(NULL))
    }
    cx <- leaf$path[[2L]]
    if (!is_in(cx, .spListContexts(data, s))) {
        msg <- glue("{lab}: unknown context '{cx}' for study '{s}'")
        abort(msg)
    }
    if (leaf$depth < 3L) {
        return(invisible(NULL))
    }
    tr <- leaf$path[[3L]]
    if (!is_in(tr, .spListTraits(data, study = s, context = cx))) {
        msg <- glue(
            "{lab}: unknown trait '{tr}' for (study '{s}', ",
            "context '{cx}')"
        )
        abort(msg)
    }
}

# Validate one method leaf (leaf vector + multi-axis + path references).
.parseMethodsCheckLeaf <- function(
    leaf,
    data,
    caps,
    multivariateMethods,
    rejectedAtUser,
    studyNames
) {
    lab <- if (length(leaf$path) == 0L) {
        "methods"
    } else {
        pathStr <- str_flatten(str_c("'", leaf$path, "'"), "$")
        glue("methods[[{pathStr}]]")
    }
    .jointValidateLeafVec(leaf$methods, lab, caps, rejectedAtUser)
    .parseMethodsCheckMultiAxis(leaf, lab, multivariateMethods)
    .parseMethodsCheckPath(leaf, lab, data, studyNames)
}

# Walk a primary `methods` tree and validate every leaf.
.parseMethodsWalked <- function(
    methods,
    data,
    caps,
    multivariateMethods,
    rejectedAtUser
) {
    walked <- .spWalkMethods(methods, label = "methods", maxDepth = 3L)
    studyNames <- .spListStudies(data)
    for (leaf in walked) {
        .parseMethodsCheckLeaf(
            leaf,
            data,
            caps,
            multivariateMethods,
            rejectedAtUser,
            studyNames
        )
    }
}

parseMethods <- function(
    methods,
    sumStatsMethods = NULL,
    qtlDatasetMethods = NULL,
    data,
    caps,
    multivariateMethods,
    rejectedAtUser = character(0)
) {
    primaryGiven <- !is.null(methods)
    splitGiven <- !is.null(sumStatsMethods) || !is.null(qtlDatasetMethods)
    .parseMethodsValidateArgs(
        primaryGiven,
        splitGiven,
        sumStatsMethods,
        qtlDatasetMethods
    )
    if (splitGiven) {
        .jointValidateLeafVec(
            sumStatsMethods,
            "sumStatsMethods",
            caps,
            rejectedAtUser
        )
        .jointValidateLeafVec(
            qtlDatasetMethods,
            "qtlDatasetMethods",
            caps,
            rejectedAtUser
        )
    } else {
        .parseMethodsWalked(
            methods,
            data,
            caps,
            multivariateMethods,
            rejectedAtUser
        )
    }
    list(
        methods = methods,
        sumStatsMethods = sumStatsMethods,
        qtlDatasetMethods = qtlDatasetMethods,
        shape = if (primaryGiven) "primary" else "split"
    )
}


# -----------------------------------------------------------------------------
# validateMethodsVsJointSpec -- cross-validation. A per-axis method assignment
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
    # Split-form methods are flat per-data-form vectors -- nothing to check.
    if (methodsParsed$shape == "split") {
        return(invisible(NULL))
    }
    if (length(jointSpecParsed) == 0L) {
        return(invisible(NULL))
    }
    methods <- methodsParsed$methods
    if (is.character(methods)) {
        return(invisible(NULL))
    } # top-level vector OK

    walked <- .spWalkMethods(methods, label = "methods", maxDepth = 3L)
    # depth observed at leaves; max depth in the spec reflects nesting level.
    maxDepth <- max(map_int(walked, "depth"))

    for (i in seq_along(jointSpecParsed)) {
        axes <- jointSpecParsed[[i]]$axes
        lab <- glue("jointSpecification[[{i}]]")
        if (is_in("study", axes) && maxDepth >= 1L) {
            msg <- glue(
                "{lab}: `axes` includes 'study' but `methods` nests ",
                "per-study; remove per-study method assignment when ",
                "joining over studies."
            )
            abort(msg)
        }
        if (is_in("context", axes) && maxDepth >= 2L) {
            msg <- glue(
                "{lab}: `axes` includes 'context' but `methods` nests ",
                "per-context; remove per-context method assignment when ",
                "joining over contexts."
            )
            abort(msg)
        }
        if (is_in("trait", axes) && maxDepth >= 3L) {
            msg <- glue(
                "{lab}: `axes` includes 'trait' but `methods` nests ",
                "per-trait; remove per-trait method assignment when ",
                "joining over traits."
            )
            abort(msg)
        }
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
.fmResolveSpecScope <- function(spec, data, contexts = NULL, traitIds = NULL) {
    scope <- spec$scope
    allStudies <- .spListStudies(data)
    studies <- if (is.null(scope$study)) {
        allStudies
    } else {
        intersect(allStudies, scope$study)
    }

    contextsOut <- set_names(
        map(
            studies,
            .fmScopeAxis,
            avail = .spListContexts,
            scopeLimit = scope$context,
            userSpec = contexts,
            requireCharacter = FALSE,
            data = data
        ),
        studies
    )
    traitsOut <- set_names(
        map(
            studies,
            .fmScopeAxis,
            avail = .fmScopeTraits,
            scopeLimit = scope$trait,
            userSpec = traitIds,
            requireCharacter = TRUE,
            data = data
        ),
        studies
    )
    list(studies = studies, contexts = contextsOut, traits = traitsOut)
}

# `.spListTraits` takes the study by name, so give it the same (data, s)
# shape the context lister has.
# @noRd
.fmScopeTraits <- function(data, s) {
    .spListTraits(data, study = s)
}

# What a pipeline-level `contexts` / `traitIds` argument narrows study `s` to,
# or NULL when it says nothing about it. `requireCharacter` mirrors the trait
# axis, which ignores a non-character per-study entry.
# @noRd
.fmScopeUserLimit <- function(spec, s, requireCharacter) {
    if (is.null(spec)) {
        return(NULL)
    }
    if (is.character(spec)) {
        return(spec)
    }
    if (!is.list(spec) || !is_in(s, names(spec))) {
        return(NULL)
    }
    entry <- spec[[s]]
    if (!requireCharacter || is.character(entry)) entry else NULL
}

# One axis for one study: everything it has, narrowed by the spec's scope and
# then by the caller's argument. Each narrowing is an intersect, so the chain
# is a fold rather than a variable rewritten in place.
# @noRd
.fmScopeAxis <- function(
    s,
    avail,
    scopeLimit,
    userSpec,
    requireCharacter,
    data
) {
    limits <- compact(list(
        scopeLimit,
        .fmScopeUserLimit(userSpec, s, requireCharacter)
    ))
    reduce(limits, intersect, .init = avail(data, s))
}


# Build a (variants x tupleRows) Z matrix from a QtlSumStats subset,
# requiring all rows to share an identical SNP order (the post-
# summaryStatsQc contract). Returns list(Z, nVec, variantIds).
# `errorLabel` is woven into the SNP-order error to identify the caller.
# @noRd
.buildJointSumstatZMatrix <- function(
    data,
    tupleRows,
    colLabels,
    errorLabel,
    ldSketch = NULL,
    cutoffs = NULL
) {
    cols <- list(
        study = as.character(data$study),
        context = as.character(data$context),
        trait = as.character(data$trait)
    )
    allDf <- getSumStatsDf(
        data,
        study = cols$study[[tupleRows[[1L]]]],
        context = cols$context[[tupleRows[[1L]]]],
        trait = cols$trait[[tupleRows[[1L]]]],
        require = c("SNP", "Z", "N")
    )
    # Every entry shares one SNP order (asserted below), so filtering the first
    # entry's ids filters the group: Z is built against this vector and each
    # entry is checked against it.
    firstDf <- allDf[
        .panelKeepMask(allDf$variant_id, ldSketch, cutoffs, errorLabel),
        ,
        drop = FALSE
    ]
    variantIds <- firstDf$variant_id
    filled <- .jointFillSumstatZ(
        data,
        tupleRows,
        cols,
        variantIds,
        colLabels,
        errorLabel
    )
    list(Z = filled$Z, nVec = filled$nVec, variantIds = variantIds)
}

# One tuple's z column and median N, checked against the shared SNP order.
# @noRd
.jointSumstatEntryAt <- function(
    kk,
    data,
    tupleRows,
    cols,
    variantIds,
    errorLabel
) {
    i <- tupleRows[[kk]]
    d <- getSumStatsDf(
        data,
        study = cols$study[[i]],
        context = cols$context[[i]],
        trait = cols$trait[[i]],
        require = c("SNP", "Z", "N")
    )
    # Narrowed to the same set before the order check, or that check fires
    # on a length mismatch the panel filter itself created.
    kept <- d[is_in(d$variant_id, variantIds), , drop = FALSE]
    .jointCheckSnpOrder(kept$variant_id, variantIds, errorLabel)
    list(z = kept$z, n = stats::median(kept$N, na.rm = TRUE))
}

# Read each tuple's z / N into one matrix. The fitters index z by column
# position, so every entry must present the same SNPs in the same order --
# which is what lets the columns simply be laid side by side here.
# @noRd
.jointFillSumstatZ <- function(
    data,
    tupleRows,
    cols,
    variantIds,
    colLabels,
    errorLabel
) {
    entries <- map(
        seq_along(tupleRows),
        .jointSumstatEntryAt,
        data = data,
        tupleRows = tupleRows,
        cols = cols,
        variantIds = variantIds,
        errorLabel = errorLabel
    )
    Z <- matrix(
        unname(list_c(map(entries, "z"))),
        nrow = length(variantIds),
        ncol = length(tupleRows),
        dimnames = list(variantIds, colLabels)
    )
    list(Z = Z, nVec = map_dbl(entries, "n"))
}

# @noRd
.jointCheckSnpOrder <- function(got, want, errorLabel) {
    if (identical(got, want)) {
        return(invisible(NULL))
    }
    abort(glue(
        "{errorLabel}: every entry in a joint group must share an ",
        "identical SNP order after summaryStatsQc()."
    ))
}


# Build a multi-context Y matrix for a single (study, trait) from an
# individual-level QtlDataset. Returns list(X, Y, perTraitContexts) or
# NULL when fewer than 2 contexts carry `tid` or the sample / complete-Y
# subset is too small to fit.
# @noRd
# --- .buildIndividualCrossContextXy helpers ---------------------------------

# Scoped contexts in which a trait is present; NULL if fewer than 2.
# Does context `cx` carry trait `tid`?
# @noRd
.ccContextHasTrait <- function(cx, data, tid) {
    is_in(tid, rownames(getPhenotypes(data, contexts = cx)))
}

.crossContextPerTrait <- function(data, tid, scopedContexts, verbose, label) {
    perTraitContexts <- keep(
        scopedContexts,
        .ccContextHasTrait,
        data = data,
        tid = tid
    )
    if (length(perTraitContexts) < 2L) {
        if (verbose >= 1) {
            nCtx <- length(perTraitContexts)
            msg <- glue(
                "{label}: trait '{tid}' present in {nCtx} scoped ",
                "context(s); skipping."
            )
            inform(msg)
        }
        return(NULL)
    }
    perTraitContexts
}

# Cross-context response matrix (one column per context) on the shared samples.
.crossContextY <- function(Yres, perTraitContexts, commonSamples) {
    yCols <- map(
        perTraitContexts,
        .crossContextYCol,
        Yres = Yres,
        commonSamples = commonSamples
    )
    exec(cbind, !!!yCols)
}

# Intersect samples, build the response matrix, and drop incomplete rows.
.crossContextAssemble <- function(
    X,
    Yres,
    perTraitContexts,
    verbose,
    label,
    tid
) {
    commonSamples <- reduce(
        c(list(rownames(X)), map(Yres, rownames)),
        intersect
    )
    if (length(commonSamples) < 2L) {
        if (verbose >= 1) {
            inform(glue(
                "{label}: trait '{tid}' has too few shared samples across ",
                "contexts; skipping."
            ))
        }
        return(NULL)
    }
    shared <- X[commonSamples, , drop = FALSE]
    Y <- .crossContextY(Yres, perTraitContexts, commonSamples)
    keep <- stats::complete.cases(Y)
    if (sum(keep) < 2L) {
        if (verbose >= 1) {
            inform(glue(
                "{label}: trait '{tid}' has too few complete-Y subjects; ",
                "skipping."
            ))
        }
        return(NULL)
    }
    list(
        X = shared[keep, , drop = FALSE],
        Y = Y[keep, , drop = FALSE],
        perTraitContexts = perTraitContexts
    )
}

.buildIndividualCrossContextXy <- function(
    data,
    tid,
    scopedContexts,
    cisWindow,
    verbose,
    label,
    region = NULL,
    residualizationArgs
) {
    perTraitContexts <- .crossContextPerTrait(
        data,
        tid,
        scopedContexts,
        verbose,
        label
    )
    if (is.null(perTraitContexts)) {
        return(NULL)
    }
    X <- .buildResidGeno(
        data,
        perTraitContexts,
        tid,
        cisWindow,
        region,
        residualizationArgs
    )
    Yres <- .fmResidPheno(
        data,
        contexts = perTraitContexts,
        traitId = tid,
        residualizationArgs = residualizationArgs
    )
    .crossContextAssemble(X, Yres, perTraitContexts, verbose, label, tid)
}


# Subset `traits` to those whose phenotype coordinates overlap `region` (the
# genes at a locus). region = NULL -> all `traits` unchanged (gene/cisWindow
# mode does not region-filter). Mirrors fineMappingPipeline's univariate region
# trait selection (ids[overlapsAny(rowRanges(se), region)]) so the joint-engine
# region path joins the same gene set.
# @noRd
.fmTraitsInRegion <- function(se, traits, region) {
    if (is.null(region) || length(traits) == 0L) {
        return(traits)
    }
    rr <- SummarizedExperiment::rowRanges(se)
    traits[IRanges::overlapsAny(rr[traits], region)]
}

# Build a multi-trait Y matrix for a single (study, context) from an
# individual-level QtlDataset. Returns list(X, Y, traitsHere, se) or NULL
# when fewer than 2 traits live in the context or the sample / complete-Y
# subset is too small.
# @noRd
# Warn + signal skip when a context has fewer than 2 scoped traits.
.crossTraitTooFew <- function(traitsHere, cx, study, verbose, label) {
    if (length(traitsHere) >= 2L) {
        return(FALSE)
    }
    if (verbose >= 1) {
        nTraits <- length(traitsHere)
        msg <- glue(
            "{label}: context '{cx}' (study '{study}') has {nTraits} ",
            "scoped trait(s); skipping."
        )
        inform(msg)
    }
    TRUE
}

.buildIndividualCrossTraitXy <- function(
    data,
    cx,
    scopedTraits,
    cisWindow,
    verbose,
    label,
    study,
    region = NULL,
    residualizationArgs
) {
    se <- getPhenotypes(data, contexts = cx)
    # scopedTraits is already region-restricted upstream (.runJointSpecs) when
    # region mode is used without an explicit traitId.
    traitsHere <- intersect(scopedTraits, rownames(se))
    if (.crossTraitTooFew(traitsHere, cx, study, verbose, label)) {
        return(NULL)
    }
    allX <- .buildResidGeno(
        data,
        cx,
        traitsHere,
        cisWindow,
        region,
        residualizationArgs
    )
    allY <- .fmResidPheno(
        data,
        contexts = cx,
        traitId = traitsHere,
        residualizationArgs = residualizationArgs
    )
    common <- intersect(rownames(allX), rownames(allY))
    if (length(common) < 2L) {
        return(NULL)
    }
    X <- allX[common, , drop = FALSE]
    Y <- allY[common, , drop = FALSE]
    keep <- stats::complete.cases(Y)
    if (sum(keep) < 2L) {
        return(NULL)
    }
    list(
        X = X[keep, , drop = FALSE],
        Y = Y[keep, , drop = FALSE],
        traitsHere = traitsHere,
        se = se
    )
}


# Build a composed-axes (context, trait) X/Y for individual-level
# QtlDataset. Returns list(X, Y, tuples) or NULL.
# @noRd
# Residualized genotype for the cross-* / composed builders (cis-window or an
# explicit region). Shared by the individual multi-axis Xy builders.
.buildResidGeno <- function(
    data,
    contexts,
    traitId,
    cisWindow,
    region,
    residualizationArgs
) {
    if (is.null(region)) {
        .fmResidGeno(
            data,
            contexts = contexts,
            traitId = traitId,
            cisWindow = cisWindow,
            residualizationArgs = residualizationArgs
        )
    } else {
        .fmResidGeno(
            data,
            contexts = contexts,
            region = region,
            residualizationArgs = residualizationArgs
        )
    }
}

# --- .buildComposedIndividualXy helpers -------------------------------------

# Enumerate the in-scope (context, trait) tuples for a study; NULL if < 2.
# @noRd
.composedTuple <- function(tid, cx) {
    list(context = cx, trait = tid)
}

# The (context, trait) tuples context `cx` contributes.
# @noRd
.composedTuplesForContext <- function(cx, data, scopedTraits) {
    se <- getPhenotypes(data, contexts = cx)
    # scopedTraits is region-restricted upstream (.runJointSpecs) if needed.
    map(intersect(scopedTraits, rownames(se)), .composedTuple, cx = cx)
}

.composedTuples <- function(data, scope, study, verbose, label) {
    scopedContexts <- scope$contexts[[study]]
    scopedTraits <- scope$traits[[study]]
    tuples <- .spConcat(map(
        scopedContexts,
        .composedTuplesForContext,
        data = data,
        scopedTraits = scopedTraits
    ))
    if (length(tuples) < 2L) {
        if (verbose >= 1) {
            nTuples <- length(tuples)
            msg <- glue(
                "{label}: study '{study}' has {nTuples} (context, trait) ",
                "tuple(s) in scope; skipping."
            )
            inform(msg)
        }
        return(NULL)
    }
    tuples
}

# Assemble the composed response matrix (one column per tuple); NULL if < 2.
# One tuple's phenotype column, labelled "context:trait", or NULL when that
# context does not carry the trait. `colnames<-` applied as a function returns
# a relabelled copy rather than renaming a binding in place.
# @noRd
.composedYCol <- function(t, YresList, commonSamples) {
    ym <- YresList[[t$context]]
    if (!is_in(t$trait, colnames(ym))) {
        return(NULL)
    }
    `colnames<-`(
        ym[commonSamples, t$trait, drop = FALSE],
        str_c(t$context, t$trait, sep = ":")
    )
}

.composedYCols <- function(YresList, tuples, commonSamples) {
    yCols <- compact(map(
        tuples,
        .composedYCol,
        YresList = YresList,
        commonSamples = commonSamples
    ))
    if (length(yCols) < 2L) {
        return(NULL)
    }
    exec(cbind, !!!yCols)
}

# The residualized genotypes and the per-context residualized phenotypes
# for one composed block, both over the same context/trait span.
# @noRd
.composedResidualized <- function(
    data,
    allContexts,
    allTraits,
    cisWindow,
    region,
    residualizationArgs
) {
    # Genotypes first, as the inline version did: reading them is what
    # fails on a missing panel file, and that error must stay the one the
    # caller sees rather than a downstream lm.fit() complaint.
    allX <- .buildResidGeno(
        data,
        allContexts,
        allTraits,
        cisWindow,
        region,
        residualizationArgs
    )
    resid <- .fmResidPheno(
        data,
        contexts = allContexts,
        traitId = allTraits,
        residualizationArgs = residualizationArgs
    )
    list(
        X = allX,
        # A single context returns the bare matrix rather than a named list.
        YresList = if (length(allContexts) == 1L) {
            set_names(list(resid), allContexts)
        } else {
            resid
        }
    )
}

.buildComposedIndividualXy <- function(
    data,
    scope,
    study,
    cisWindow,
    verbose,
    label,
    region = NULL,
    residualizationArgs
) {
    tuples <- .composedTuples(data, scope, study, verbose, label)
    if (is.null(tuples)) {
        return(NULL)
    }
    allContexts <- unique(map_chr(tuples, "context"))
    allTraits <- unique(map_chr(tuples, "trait"))
    sides <- .composedResidualized(
        data,
        allContexts,
        allTraits,
        cisWindow,
        region,
        residualizationArgs
    )
    allX <- sides$X
    YresList <- sides$YresList
    commonSamples <- reduce(
        c(list(rownames(allX)), map(YresList, rownames)),
        intersect
    )
    if (length(commonSamples) < 2L) {
        return(NULL)
    }
    X <- allX[commonSamples, , drop = FALSE]
    Y <- .composedYCols(YresList, tuples, commonSamples)
    if (is.null(Y)) {
        return(NULL)
    }
    keep <- stats::complete.cases(Y)
    if (sum(keep) < 2L) {
        return(NULL)
    }
    list(
        X = X[keep, , drop = FALSE],
        Y = Y[keep, , drop = FALSE],
        tuples = tuples
    )
}


# Enumerate composed-axes row groups for a QtlSumStats input. Returns the
# list of (rowIdx) per group along with the per-axis identity columns
# needed to label the output row. Groups containing fewer than 2 rows
# are returned unfiltered; the caller decides whether to skip.
# @noRd
.enumerateComposedSumstatGroups <- function(spec, data, scope) {
    axes <- spec$axes
    complement <- setdiff(c("study", "context", "trait"), axes)
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    inScope <- map_lgl(
        seq_len(nrow(data)),
        .composedRowInScope,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol,
        scope = scope
    )
    rowIdx <- which(inScope)
    if (length(rowIdx) == 0L) {
        return(NULL)
    }
    groupKey <- if (length(complement) == 0L) {
        rep("__all__", length(rowIdx))
    } else {
        pasteArgs <- c(
            map(
                complement,
                .composedAxisCol,
                studyCol = studyCol,
                contextCol = contextCol,
                traitCol = traitCol,
                rowIdx = rowIdx
            ),
            list(sep = "||")
        )
        exec(paste, !!!pasteArgs)
    }
    groups <- split(rowIdx, groupKey)
    list(
        groups = groups,
        axes = axes,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol
    )
}


# =============================================================================
# Fine-mapping dispatchers
# =============================================================================

# Identity-tuple key (study/context/trait/method joined by "\r"), used to align
# per-region result entries when merging. Shared by the fm/twas mergers.
# @noRd
.mergeResultKeyOf <- function(r) {
    str_c(
        as.character(r$study),
        as.character(r$context),
        as.character(r$trait),
        as.character(r$method),
        sep = "\r"
    )
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
    if (n == 0L) {
        return(base)
    }
    baseKeys <- .mergeResultKeyOf(base)
    mergedEntries <- map(
        seq_len(n),
        .fmMergedEntryAt,
        results = results,
        baseKeys = baseKeys
    )
    qfmrArgs <- c(
        list(
            study = as.character(base$study),
            context = as.character(base$context),
            trait = as.character(base$trait),
            method = as.character(base$method),
            entry = mergedEntries
        ),
        .jointCols(base),
        list(ldSketch = NULL)
    )
    exec(QtlFineMappingResult, !!!qfmrArgs)
}

# One passthrough column of a joint result row as a character vector (the joint-
# key columns), or NULL when absent.
# @noRd
.jointStrCol <- function(nm, df) {
    if (is_in(nm, names(df))) as.character(df[[nm]]) else NULL
}

# One passthrough column carried through UNCOERCED (GRanges provenance columns
# region / traitPos), or NULL when absent.
# @noRd
.jointRawCol <- function(nm, df) if (is_in(nm, names(df))) df[[nm]] else NULL

# The optional passthrough columns of a per-tuple result row, as a named list
# (NULL for any absent column): the three joint-key columns (jointStudies /
# jointContexts / jointTraits) plus the GRanges provenance columns (region /
# traitPos). Spliced into the QtlFineMappingResult / TwasWeights constructors so
# by-key / cross-study rebuilds preserve them.
.jointCols <- function(df) {
    # traitPos is a GRanges provenance column: carry it through uncoerced so
    # by-key / cross-study rebuilds keep the trait position instead of silently
    # dropping it. There is no `region` counterpart -- the fine-mapping window
    # retired (section 4.4) because the element's own span is the region and a
    # stored window had no correct update rule under subsetRegion().
    list(
        jointStudies = .jointStrCol("jointStudies", df),
        jointContexts = .jointStrCol("jointContexts", df),
        jointTraits = .jointStrCol("jointTraits", df),
        traitPos = .jointRawCol("traitPos", df)
    )
}

# Shared tail of the MultiStudyQtlDataset fineMapping / twasWeights pipeline
# methods: recurse into each embedded QtlDataset then the embedded QtlSumStats,
# row-bind the per-study results, build the per-tuple result object, and combine
# it with any joint-specification result. `perStudyFn` / `sumStatsFn` are the
# (method-specific) single-study recursions; `rbindFn` / `resultCtor` /
# `pipelineName` supply the per-pipeline pieces. Each method keeps its own
# (divergent) method-gating / joint-dispatch preamble and passes the computed
# `jointResult` in.
# @noRd
# --- .multiStudyPipelineDriver helpers --------------------------------------

# Accumulate per-study + sumstats results (tracking the embedded LD sketch).
# rbind two results with no LD sketch -- the per-study combine.
# @noRd
.msRbindNoLd <- function(acc, x, rbindFn) {
    rbindFn(acc, x, ldSketch = NULL)
}

# Fold per-study results together, skipping studies that produced nothing.
# NULL when nothing did, which is what the callers' `out` started as.
# @noRd
.msFoldStudyResults <- function(results, rbindFn) {
    kept <- compact(results)
    if (length(kept) == 0L) {
        return(NULL)
    }
    reduce(kept, .msRbindNoLd, rbindFn = rbindFn)
}

.msDriverAccumulate <- function(
    qtlDatasets,
    sumStats,
    perStudyFn,
    sumStatsFn,
    cfg,
    rbindFn
) {
    fromStudies <- .msFoldStudyResults(
        map(qtlDatasets, perStudyFn, cfg),
        rbindFn
    )
    ssRes <- if (is.null(sumStats)) NULL else sumStatsFn(sumStats, cfg)
    if (is.null(ssRes)) {
        return(list(out = fromStudies, embeddedLd = NULL))
    }
    embeddedLd <- getLdSketch(ssRes)
    out <- if (is.null(fromStudies)) {
        ssRes
    } else {
        rbindFn(fromStudies, ssRes, ldSketch = embeddedLd)
    }
    list(out = out, embeddedLd = embeddedLd)
}

# Reconstruct the per-tuple result object from the accumulated rows.
.msDriverPerTuple <- function(out, resultCtor, embeddedLd) {
    # ldSketch: NULL if all studies were individual-level; the embedded
    # sumStats's ldSketch otherwise.
    ctorArgs <- c(
        list(
            study = as.character(out$study),
            context = as.character(out$context),
            trait = as.character(out$trait),
            method = as.character(out$method),
            entry = .collectionEntries(out)
        ),
        .jointCols(out),
        list(ldSketch = embeddedLd)
    )
    exec(resultCtor, !!!ctorArgs)
}

.multiStudyPipelineDriver <- function(
    data,
    jointResult,
    perStudyFn,
    sumStatsFn,
    cfg,
    rbindFn,
    resultCtor,
    pipelineName,
    noun = "a result"
) {
    acc <- .msDriverAccumulate(
        getQtlDatasets(data),
        getSumStats(data),
        perStudyFn,
        sumStatsFn,
        cfg,
        rbindFn
    )
    out <- acc$out
    embeddedLd <- acc$embeddedLd
    perTupleResult <- if (!is.null(out)) {
        .msDriverPerTuple(out, resultCtor, embeddedLd)
    } else {
        NULL
    }
    if (is.null(jointResult)) {
        if (is.null(perTupleResult)) {
            msg <- glue(
                "{pipelineName}(MultiStudyQtlDataset): no entries produced ",
                "{noun}."
            )
            abort(msg)
        }
        return(perTupleResult)
    }
    if (is.null(perTupleResult)) {
        return(jointResult)
    }
    rbindFn(perTupleResult, jointResult, ldSketch = embeddedLd)
}

# Synthesize a jointSpecification for the AUTO-DETECTION path (no explicit
# jointSpecification supplied): route mvsusie / fsusie over the data's natural
# multi-axis shape through the SAME engine as an explicit jointSpecification,
# matching the historical mvJobs / runMultivariate detection. >= 2 traits ->
# cross-trait (covers multi-trait single-context AND both-multi, since the
# cross-trait enumerator iterates contexts -> per-context multi-trait fits);
# else >= 2 contexts (single trait) -> cross-context. Single context & single
# trait -> no joint (the caller's multivariate guard already rejects mvsusie /
# fsusie there). Returns a list of parsed specs (full scope) or list().
# @noRd
.fmSynthesizeJointSpec <- function(nCtx, nTraits) {
    if (nTraits >= 2L) {
        list(list(axes = "trait", scope = NULL))
    } else if (nCtx >= 2L) {
        list(list(axes = "context", scope = NULL))
    } else {
        list()
    }
}


.fmDispatchJointSpecsQtlDataset <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    credibleSetArgs = CredibleSetParam(includeAllCs = FALSE),
    verbose,
    methodArgs = list(),
    xRegions = list(NULL),
    mrmashPrior = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    crossValidationArgs = CrossValidationParam(),
    residualizationArgs = ResidualizationParam(),
    pipCutoffToSkip = 0,
    fineMappingResult = NULL,
    fitRetention = c("slim", "none", "full"),
    seed = NULL
) {
    fitRetention <- arg_match(fitRetention)
    # Run the joint dispatch once per region block, then merge per
    # (study, context, trait, method) across regions. A single block (cis or
    # jointRegions=TRUE concatenated) returns its result directly.
    # xRegions is deliberately absent: each region is supplied per call.
    perRegion <- map(
        xRegions,
        .fmDispatchJointSpecRegion,
        parsedJointSpec = parsedJointSpec,
        data = data,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        credibleSetArgs = credibleSetArgs,
        verbose = verbose,
        methodArgs = methodArgs,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        crossValidationArgs = crossValidationArgs,
        residualizationArgs = residualizationArgs,
        pipCutoffToSkip = pipCutoffToSkip,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        seed = seed
    ) |>
        compact()
    .fmCombineRegionResults(perRegion)
}

# Merge the per-region joint results by (study, context, trait, method). A
# single block -- cis, or jointRegions = TRUE concatenated -- is returned
# directly rather than merged with itself.
# @noRd
.fmCombineRegionResults <- function(perRegion) {
    if (length(perRegion) == 0L) {
        return(NULL)
    }
    if (length(perRegion) == 1L) {
        return(perRegion[[1L]])
    }
    .fmMergeResultsByKey(perRegion)
}

# FmJointPipeline for individual-level fine-mapping, built from the call params.
.fmJointPipeline <- function(
    credibleSetArgs = CredibleSetParam(includeAllCs = FALSE),
    dataDrivenPriorWeightsCutoff,
    crossValidationArgs,
    residualizationArgs,
    verbose,
    fitRetention = c("slim", "none", "full"),
    seed
) {
    fitRetention <- arg_match(fitRetention)
    new(
        "FmJointPipeline",
        config = list(
            credibleSetArgs = credibleSetArgs,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            crossValidationArgs = crossValidationArgs,
            residualizationArgs = residualizationArgs,
            verbose = verbose,
            fitRetention = fitRetention,
            seed = seed,
            ldSketch = NULL
        )
    )
}

# The joint pipeline record, once the spec has been checked against the
# individual-level form -- a study axis has no meaning there, and building
# a pipeline for a spec that names one would defer the error past the fit.
# @noRd
.fmJointPipelineChecked <- function(
    parsedJointSpec,
    credibleSetArgs,
    dataDrivenPriorWeightsCutoff,
    crossValidationArgs,
    residualizationArgs,
    verbose,
    fitRetention,
    seed
) {
    .jointRejectStudyOnIndividual(parsedJointSpec)
    .fmJointPipeline(
        credibleSetArgs = credibleSetArgs,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        crossValidationArgs = crossValidationArgs,
        residualizationArgs = residualizationArgs,
        verbose = verbose,
        fitRetention = fitRetention,
        seed = seed
    )
}

.fmDispatchJointSpecsQtlDatasetOneRegion <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    verbose,
    methodArgs = list(),
    region = NULL,
    mrmashPrior = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    crossValidationArgs = CrossValidationParam(),
    residualizationArgs = ResidualizationParam(),
    pipCutoffToSkip = 0,
    fineMappingResult = NULL,
    seed = NULL,
    credibleSetArgs,
    fitRetention
) {
    # Engine routing (jointEngine.R); one region block (the caller loops
    # regions).
    pipeline <- .fmJointPipelineChecked(
        parsedJointSpec,
        credibleSetArgs = credibleSetArgs,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        crossValidationArgs = crossValidationArgs,
        residualizationArgs = residualizationArgs,
        verbose = verbose,
        fitRetention = fitRetention,
        seed = seed
    )
    .runJointSpecs(
        parsedJointSpec,
        data,
        dataForm = "individual",
        pipeline = pipeline,
        jointMethods = intersect(methods, c("mvsusie", "fsusie")),
        contexts = contexts,
        traitIds = traitIds,
        args = list(
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            methodArgs = methodArgs,
            cisWindow = cisWindow,
            region = region,
            verbose = verbose,
            pipCutoffToSkip = pipCutoffToSkip,
            cache = fineMappingResult
        )
    )
}


# Top-level joint dispatcher for fineMappingPipeline(QtlSumStats).
# @noRd
# FmJointPipeline for summary-statistics fine-mapping (RSS: no sample folds;
# LD sketch drawn from the data).
.fmSumStatsPipeline <- function(
    data,
    credibleSetArgs,
    dataDrivenPriorWeightsCutoff,
    verbose,
    fitRetention
) {
    new(
        "FmJointPipeline",
        config = list(
            credibleSetArgs = credibleSetArgs,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            verbose = verbose,
            fitRetention = fitRetention,
            crossValidationArgs = .cvResolve(CrossValidationParam()),
            ldSketch = getLdSketch(data)
        )
    )
}

.fmDispatchJointSpecsQtlSumStats <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    verbose,
    methodArgs = list(),
    mrmashPrior = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    fineMappingResult = NULL,
    panelFilterArgs = PanelFilterParam(),
    credibleSetArgs,
    fitRetention
) {
    # Engine routing (jointEngine.R): the dispatch table + .runJointCell replace
    # the per-axis switch + the cross-context/trait/study/composed leaf
    # dispatchers.
    pipeline <- .fmSumStatsPipeline(
        data = data,
        credibleSetArgs = credibleSetArgs,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        verbose = verbose,
        fitRetention = fitRetention
    )
    .runJointSpecs(
        parsedJointSpec,
        data,
        dataForm = "sumstats",
        pipeline = pipeline,
        jointMethods = intersect(methods, c("mvsusie", "fsusie")),
        contexts = contexts,
        traitIds = traitIds,
        args = list(
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            methodArgs = methodArgs,
            verbose = verbose,
            cache = fineMappingResult,
            cutoffs = .panelCutoffs(panelFilterArgs)
        )
    )
}


# Top-level joint dispatcher for fineMappingPipeline(MultiStudyQtlDataset).
# Routes per-component AND per-axis: a spec with `axes = "study"` only
# touches the sumStats slot; `axes = "context"` and `axes = "trait"` run
# on every component.
# @noRd
# --- .fmDispatchJointSpecsMultiStudy helpers --------------------------------

# Partition joint specs into those with a `study` axis and the rest.
.fmSplitStudyAxisSpecs <- function(parsedJointSpec) {
    hasStudy <- map_lgl(parsedJointSpec, .jsHasStudyAxis)
    list(
        study = parsedJointSpec[hasStudy],
        nonStudy = parsedJointSpec[!hasStudy]
    )
}

# Note that individual-level studies are excluded from cross-study fits.
.fmMultiStudyWarnExcluded <- function(studyAxisSpecs, qtlDatasets, verbose) {
    if (
        length(studyAxisSpecs) > 0L && length(qtlDatasets) > 0L && verbose >= 1
    ) {
        qdNames <- str_flatten(names(qtlDatasets), ", ")
        msg <- glue(
            "jointCrossStudy: excluding individual-level studies ",
            "({qdNames}) from cross-study fits (no LD sketch available); ",
            "sumstats studies participate."
        )
        inform(msg)
    }
}

# Fine-map the non-study-axis specs on each individual-level QtlDataset.
.fmMultiStudyQtlLoop <- function(
    nonStudyAxisSpecs,
    qtlDatasets,
    methods,
    contexts,
    traitIds,
    cisWindow,
    credibleSetArgs,
    verbose,
    methodArgs,
    xRegions,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff
) {
    if (length(nonStudyAxisSpecs) == 0L) {
        return(NULL)
    }
    .msFoldStudyResults(
        map(
            qtlDatasets,
            .fmDispatchForStudy,
            nonStudyAxisSpecs = nonStudyAxisSpecs,
            methods = methods,
            contexts = contexts,
            traitIds = traitIds,
            cisWindow = cisWindow,
            credibleSetArgs = credibleSetArgs,
            verbose = verbose,
            methodArgs = methodArgs,
            xRegions = xRegions,
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
        ),
        .rbindFineMappingResult
    )
}

# `map()` hands the dataset first; the dispatcher wants the specs first.
# @noRd
.fmDispatchForStudy <- function(
    qd,
    nonStudyAxisSpecs,
    methods,
    contexts,
    traitIds,
    cisWindow,
    verbose,
    methodArgs,
    xRegions,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff,
    credibleSetArgs
) {
    .fmDispatchJointSpecsQtlDataset(
        nonStudyAxisSpecs,
        qd,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        credibleSetArgs = credibleSetArgs,
        verbose = verbose,
        methodArgs = methodArgs,
        xRegions = xRegions,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
}

# Fine-map all specs on the sumstats collection; rbind onto `out`.
.fmMultiStudySumStats <- function(
    parsedJointSpec,
    sumStats,
    studyAxisSpecs,
    out,
    methods,
    contexts,
    traitIds,
    verbose,
    methodArgs,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff,
    credibleSetArgs
) {
    if (is.null(sumStats)) {
        if (length(studyAxisSpecs) > 0L && verbose >= 1) {
            msg <- glue(
                "jointCrossStudy: no sumStats slot present on this ",
                "MultiStudyQtlDataset; cross-study specs produce no result."
            )
            inform(msg)
        }
        return(out)
    }
    ssRes <- .fmDispatchJointSpecsQtlSumStats(
        parsedJointSpec,
        sumStats,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        credibleSetArgs = credibleSetArgs,
        verbose = verbose,
        methodArgs = methodArgs,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
    if (is.null(ssRes)) {
        return(out)
    }
    embeddedLd <- getLdSketch(ssRes)
    if (is.null(out)) {
        ssRes
    } else {
        .rbindFineMappingResult(out, ssRes, ldSketch = embeddedLd)
    }
}

.fmDispatchJointSpecsMultiStudy <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    credibleSetArgs,
    verbose,
    methodArgs = list(),
    xRegions = list(NULL),
    mrmashPrior = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10
) {
    qtlDatasets <- getQtlDatasets(data)
    sumStats <- getSumStats(data)
    specs <- .fmSplitStudyAxisSpecs(parsedJointSpec)
    .fmMultiStudyWarnExcluded(specs$study, qtlDatasets, verbose)
    out <- .fmMultiStudyQtlLoop(
        specs$nonStudy,
        qtlDatasets,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        credibleSetArgs = credibleSetArgs,
        verbose = verbose,
        methodArgs = methodArgs,
        xRegions = xRegions,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
    .fmMultiStudySumStats(
        parsedJointSpec,
        sumStats,
        specs$study,
        out,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        verbose = verbose,
        methodArgs = methodArgs,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
}


# =============================================================================
# TWAS-weights dispatchers
# =============================================================================

# Top-level joint dispatcher for twasWeightsPipeline(QtlDataset).
# @noRd
# Merge per-region TwasWeights collections (same keys across regions) into one
# by concatenating each (study, context, trait, method) row's entry via
# .twasMergeRegionEntries (stacked weights + flat per-region cvResult).
# @noRd
.twasMergeResultsByKey <- function(results, regionLabels) {
    base <- results[[1L]]
    n <- length(base$method)
    if (n == 0L) {
        return(base)
    }
    baseKeys <- .mergeResultKeyOf(base)
    mergedEntries <- map(
        seq_len(n),
        .twasMergedEntryAt,
        results = results,
        baseKeys = baseKeys,
        regionLabels = regionLabels
    )
    # Passthrough columns (joint keys + region + traitPos) are per-row
    # properties of `base`, which aligns row-for-row with mergedEntries; splice
    # them so a multi-region merge preserves provenance instead of dropping it.
    twArgs <- c(
        list(
            study = as.character(base$study),
            context = as.character(base$context),
            trait = as.character(base$trait),
            method = as.character(base$method),
            entry = mergedEntries
        ),
        .jointCols(base)
    )
    exec(TwasWeights, !!!twArgs)
}

.twasDispatchJointSpecsQtlDataset <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    dataType,
    verbose,
    xRegions = list(NULL),
    fitRetention = c("slim", "none", "full"),
    seed = NULL
) {
    fitRetention <- arg_match(fitRetention)
    # Run the joint dispatch once per region block, then merge per
    # (study, context, trait, method) across regions. A single block (cis or
    # jointRegions=TRUE concatenated) returns its result directly.
    allRegions <- map(
        xRegions,
        .twasDispatchJointSpecRegion,
        parsedJointSpec = parsedJointSpec,
        data = data,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        dataType = dataType,
        verbose = verbose,
        fitRetention = fitRetention,
        seed = seed
    )
    allLabs <- map_chr(xRegions, .twasRegionLabel)
    keep <- !map_lgl(allRegions, is.null)
    perRegion <- allRegions[keep]
    labs <- allLabs[keep]
    if (length(perRegion) == 0L) {
        return(NULL)
    }
    if (length(perRegion) == 1L) {
        return(perRegion[[1L]])
    }
    .twasMergeResultsByKey(perRegion, labs)
}

.twasDispatchJointSpecsQtlDatasetOneRegion <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    dataType,
    verbose,
    region = NULL,
    fitRetention = c("slim", "none", "full"),
    seed = NULL
) {
    fitRetention <- arg_match(fitRetention)
    # Engine routing (jointEngine.R); one region block (the caller loops
    # regions).
    .jointRejectStudyOnIndividual(parsedJointSpec)
    pipeline <- new(
        "TwasJointPipeline",
        config = list(
            dataType = dataType,
            crossValidationArgs = .cvResolve(CrossValidationParam()),
            fitFullData = TRUE,
            standardized = FALSE,
            seed = seed,
            ldSketch = NULL
        )
    )
    .runJointSpecs(
        parsedJointSpec,
        data,
        dataForm = "individual",
        pipeline = pipeline,
        jointMethods = intersect(methods, "mrmash"),
        contexts = contexts,
        traitIds = traitIds,
        args = list(
            methodArgs = list(),
            cisWindow = cisWindow,
            region = region,
            verbose = verbose
        )
    )
}


# Top-level joint dispatcher for twasWeightsPipeline(QtlSumStats).
# @noRd
.twasDispatchJointSpecsQtlSumStats <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    dataType,
    verbose,
    fitRetention = c("slim", "none", "full"),
    panelFilterArgs = PanelFilterParam()
) {
    fitRetention <- arg_match(fitRetention)
    # Engine routing (jointEngine.R).
    pipeline <- new(
        "TwasJointPipeline",
        config = list(
            dataType = dataType,
            crossValidationArgs = .cvResolve(CrossValidationParam()),
            fitFullData = TRUE,
            standardized = TRUE,
            ldSketch = getLdSketch(data)
        )
    )
    .runJointSpecs(
        parsedJointSpec,
        data,
        dataForm = "sumstats",
        pipeline = pipeline,
        jointMethods = intersect(methods, "mrmash"),
        contexts = contexts,
        traitIds = traitIds,
        args = list(
            methodArgs = list(),
            verbose = verbose,
            cutoffs = .panelCutoffs(panelFilterArgs)
        )
    )
}


# Top-level joint dispatcher for twasWeightsPipeline(MultiStudyQtlDataset).
# @noRd
# --- .twasDispatchJointSpecsMultiStudy helpers ------------------------------

# Note that individual-level studies are excluded from cross-study TWAS fits.
.twasMultiStudyWarnExcluded <- function(studyAxisSpecs, qtlDatasets, verbose) {
    if (
        length(studyAxisSpecs) > 0L && length(qtlDatasets) > 0L && verbose >= 1
    ) {
        qdNames <- str_flatten(names(qtlDatasets), ", ")
        msg <- glue(
            "jointCrossStudy (twas): excluding individual-level ",
            "studies ({qdNames}) from cross-study fits; sumstats studies ",
            "participate."
        )
        inform(msg)
    }
}

# Learn weights for the non-study-axis specs on each individual-level dataset.
.twasMultiStudyQtlLoop <- function(
    nonStudyAxisSpecs,
    qtlDatasets,
    methods,
    contexts,
    traitIds,
    cisWindow,
    dataType,
    verbose,
    xRegions,
    fitRetention,
    seed
) {
    if (length(nonStudyAxisSpecs) == 0L) {
        return(NULL)
    }
    .msFoldStudyResults(
        map(
            qtlDatasets,
            .twasDispatchForStudy,
            nonStudyAxisSpecs = nonStudyAxisSpecs,
            methods = methods,
            contexts = contexts,
            traitIds = traitIds,
            cisWindow = cisWindow,
            dataType = dataType,
            verbose = verbose,
            xRegions = xRegions,
            fitRetention = fitRetention,
            seed = seed
        ),
        .rbindTwasWeights
    )
}

# `map()` hands the dataset first; the dispatcher wants the specs first.
# @noRd
.twasDispatchForStudy <- function(
    qd,
    nonStudyAxisSpecs,
    methods,
    contexts,
    traitIds,
    cisWindow,
    dataType,
    verbose,
    xRegions,
    fitRetention,
    seed
) {
    .twasDispatchJointSpecsQtlDataset(
        nonStudyAxisSpecs,
        qd,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        dataType = dataType,
        verbose = verbose,
        xRegions = xRegions,
        fitRetention = fitRetention,
        seed = seed
    )
}

# Learn weights for all specs on the sumstats collection; rbind onto `out`.
.twasMultiStudySumStats <- function(
    parsedJointSpec,
    sumStats,
    studyAxisSpecs,
    out,
    methods,
    contexts,
    traitIds,
    dataType,
    verbose,
    fitRetention
) {
    if (is.null(sumStats)) {
        if (length(studyAxisSpecs) > 0L && verbose >= 1) {
            msg <- glue(
                "jointCrossStudy (twas): no sumStats slot present on this ",
                "MultiStudyQtlDataset; cross-study specs produce no result."
            )
            inform(msg)
        }
        return(out)
    }
    ssRes <- .twasDispatchJointSpecsQtlSumStats(
        parsedJointSpec,
        sumStats,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        dataType = dataType,
        verbose = verbose,
        fitRetention = fitRetention
    )
    if (is.null(ssRes)) {
        return(out)
    }
    embeddedLd <- getLdSketch(ssRes)
    if (is.null(out)) {
        ssRes
    } else {
        .rbindTwasWeights(out, ssRes, ldSketch = embeddedLd)
    }
}

.twasDispatchJointSpecsMultiStudy <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    dataType,
    verbose,
    xRegions = list(NULL),
    fitRetention = c("slim", "none", "full"),
    seed = NULL
) {
    fitRetention <- arg_match(fitRetention)
    qtlDatasets <- getQtlDatasets(data)
    sumStats <- getSumStats(data)
    specs <- .fmSplitStudyAxisSpecs(parsedJointSpec)
    .twasMultiStudyWarnExcluded(specs$study, qtlDatasets, verbose)
    out <- .twasMultiStudyQtlLoop(
        specs$nonStudy,
        qtlDatasets,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        dataType = dataType,
        verbose = verbose,
        xRegions = xRegions,
        fitRetention = fitRetention,
        seed = seed
    )
    .twasMultiStudySumStats(
        parsedJointSpec,
        sumStats,
        specs$study,
        out,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        dataType = dataType,
        verbose = verbose,
        fitRetention = fitRetention
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Parse the joint spec at position `i` (element + its index carry to the
# parser).
# @noRd
.parseJointSpecAt <- function(i, jointSpecification, data) {
    .parseOneJointSpec(jointSpecification[[i]], i, data)
}

# One context's response column on the shared samples, named by the context.
# @noRd
.crossContextYCol <- function(cx, Yres, commonSamples) {
    `colnames<-`(Yres[[cx]][commonSamples, , drop = FALSE], cx)
}

# TRUE when sumstats row `i`'s (study, context, trait) is entirely in scope.
# @noRd
.composedRowInScope <- function(i, studyCol, contextCol, traitCol, scope) {
    s <- studyCol[i]
    cx <- contextCol[i]
    tr <- traitCol[i]
    is_in(s, scope$studies) &&
        is_in(cx, scope$contexts[[s]]) &&
        is_in(tr, scope$traits[[s]])
}

# The scoped-row identity column for complement axis `a` (study/context/trait).
# @noRd
.composedAxisCol <- function(a, studyCol, contextCol, traitCol, rowIdx) {
    switch(
        a,
        study = studyCol[rowIdx],
        context = contextCol[rowIdx],
        trait = traitCol[rowIdx]
    )
}

# The entry in result `r` matching merge key `key`, or NULL when absent. Shared
# by the FM and TWAS cross-region merges.
# @noRd
.mergeEntryForKey <- function(r, key) {
    hit <- which(.mergeResultKeyOf(r) == key)
    if (!length(hit)) {
        return(NULL)
    }
    .rowParts(r, hit[[1L]])
}

# The merged FM entry for base row `i`: gather that key's entry from every
# region
# and fold them together.
# @noRd
.fmMergedEntryAt <- function(i, results, baseKeys) {
    perRegion <- map(results, .mergeEntryForKey, key = baseKeys[[i]])
    .fmMergeEntries(compact(perRegion))
}

# The merged TWAS entry for base row `i`: gather + concatenate that key's entry
# across regions, keeping region labels aligned to the surviving entries.
# @noRd
.twasMergedEntryAt <- function(i, results, baseKeys, regionLabels) {
    perRegion <- map(results, .mergeEntryForKey, key = baseKeys[[i]])
    keep <- !map_lgl(perRegion, is.null)
    .twasMergeRegionEntries(perRegion[keep], regionLabels[keep])
}

# TRUE when a parsed joint spec includes the `study` axis.
# @noRd
.jsHasStudyAxis <- function(s) {
    is_in("study", s$axes)
}

# One region's FM joint-spec dispatch; `args` bundles the shared call arguments.
# @noRd
.fmDispatchJointSpecRegion <- function(
    rg,
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    credibleSetArgs,
    verbose,
    methodArgs,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff,
    crossValidationArgs,
    residualizationArgs,
    pipCutoffToSkip,
    fineMappingResult,
    fitRetention,
    seed
) {
    .fmDispatchJointSpecsQtlDatasetOneRegion(
        parsedJointSpec = parsedJointSpec,
        data = data,
        methods = methods,
        contexts = contexts,
        traitIds = traitIds,
        cisWindow = cisWindow,
        credibleSetArgs = credibleSetArgs,
        verbose = verbose,
        methodArgs = methodArgs,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        crossValidationArgs = crossValidationArgs,
        residualizationArgs = residualizationArgs,
        pipCutoffToSkip = pipCutoffToSkip,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        seed = seed,
        region = rg
    )
}

# One region's TWAS joint-spec dispatch over a QtlDataset.
# @noRd
.twasDispatchJointSpecRegion <- function(
    rg,
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    dataType,
    verbose,
    fitRetention,
    seed = NULL
) {
    .twasDispatchJointSpecsQtlDatasetOneRegion(
        parsedJointSpec,
        data,
        methods,
        contexts,
        traitIds,
        cisWindow,
        dataType,
        verbose,
        region = rg,
        fitRetention = fitRetention,
        seed = seed
    )
}
