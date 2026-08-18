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
        return(data@study)
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
        if (!identical(study, data@study)) {
            dataStudy <- data@study
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
.spListContexts <- function(data, study = NULL) {
    if (is(data, "QtlDataset")) {
        if (!is.null(study) && !identical(study, data@study)) {
            return(character(0))
        }
        return(names(data@phenotypes))
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
            out <- character(0)
            for (qd in indDatasets) {
                out <- c(out, names(qd@phenotypes))
            }
            if (!is.null(ss)) {
                out <- c(out, unique(as.character(ss$context)))
            }
            return(unique(out))
        }
        if (is_in(study, names(indDatasets))) {
            return(names(indDatasets[[study]]@phenotypes))
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
# Traits available in a single individual-level QtlDataset (optionally scoped).
.spListTraitsQtlDataset <- function(data, study, context) {
    if (!is.null(study) && !identical(study, data@study)) {
        return(character(0))
    }
    if (is.null(context)) {
        return(unique(unlist(
            map(data@phenotypes, rownames),
            use.names = FALSE
        )))
    }
    se <- data@phenotypes[[context]]
    if (is.null(se)) {
        return(character(0))
    }
    rownames(se)
}

# Traits across a MultiStudyQtlDataset (individual studies + sumstats).
.spListTraitsMultiStudy <- function(data, study, context) {
    indDatasets <- getQtlDatasets(data)
    ss <- getSumStats(data)
    # No study filter: aggregate across every component (individual + sumstats).
    # Must precede per-study branches so a present sumStats slot does not shadow
    # the individual-level studies' traits when study = NULL.
    if (is.null(study)) {
        out <- character(0)
        for (qd in indDatasets) {
            out <- c(out, .spListTraits(qd, context = context))
        }
        if (!is.null(ss)) {
            out <- c(out, .spListTraits(ss, context = context))
        }
        return(unique(out))
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
        keep <- rep(TRUE, nrow(data))
        if (!is.null(study)) {
            keep <- keep & as.character(data$study) == study
        }
        if (!is.null(context)) {
            keep <- keep & as.character(data$context) == context
        }
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
    out <- list()
    for (s in studies) {
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
        out[[s]] <- intersect(contexts, avail)
    }
    out
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
    out <- list()
    for (s in studies) {
        avail <- .spListContexts(data, s)
        out[[s]] <- if (is_in(s, names(contexts))) {
            .parseContextsStudy(as.character(contexts[[s]]), s, avail)
        } else {
            avail
        }
    }
    out
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
    sub <- list()
    for (cx in names(val)) {
        sub[[cx]] <- .parseTraitIdContext(val[[cx]], s, cx, data)
    }
    sub
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
    out <- list()
    for (s in names(traitId)) {
        out[[s]] <- .parseTraitIdStudy(traitId[[s]], s, data)
    }
    out
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
    out <- list()
    for (nm in names(spec)) {
        out <- c(
            out,
            .spWalkMethods(
                spec[[nm]],
                label = label,
                depth = depth + 1L,
                maxDepth = maxDepth,
                path = c(path, nm)
            )
        )
    }
    out
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
        if (!is.character(sumStatsMethods) || length(sumStatsMethods) == 0L) {
            abort("`sumStatsMethods` must be a non-empty character vector.")
        }
        if (
            !is.character(qtlDatasetMethods) || length(qtlDatasetMethods) == 0L
        ) {
            abort("`qtlDatasetMethods` must be a non-empty character vector.")
        }
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
    if (leaf$depth >= 1L) {
        s <- leaf$path[[1L]]
        if (!is_in(s, studyNames)) {
            msg <- glue("{lab}: unknown study '{s}'")
            abort(msg)
        }
    }
    if (leaf$depth >= 2L) {
        s <- leaf$path[[1L]]
        cx <- leaf$path[[2L]]
        if (!is_in(cx, .spListContexts(data, s))) {
            msg <- glue("{lab}: unknown context '{cx}' for study '{s}'")
            abort(msg)
        }
    }
    if (leaf$depth >= 3L) {
        s <- leaf$path[[1L]]
        cx <- leaf$path[[2L]]
        tr <- leaf$path[[3L]]
        if (!is_in(tr, .spListTraits(data, study = s, context = cx))) {
            msg <- glue(
                "{lab}: unknown trait '{tr}' for (study '{s}', ",
                "context '{cx}')"
            )
            abort(msg)
        }
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
    studies <- .spListStudies(data)
    if (!is.null(scope$study)) {
        studies <- intersect(studies, scope$study)
    }

    contextsOut <- list()
    traitsOut <- list()
    for (s in studies) {
        ctxAvail <- .spListContexts(data, s)
        if (!is.null(scope$context)) {
            ctxAvail <- intersect(ctxAvail, scope$context)
        }
        if (!is.null(contexts)) {
            if (is.list(contexts) && is_in(s, names(contexts))) {
                ctxAvail <- intersect(ctxAvail, contexts[[s]])
            } else if (is.character(contexts)) {
                ctxAvail <- intersect(ctxAvail, contexts)
            }
        }
        contextsOut[[s]] <- ctxAvail

        trAvail <- .spListTraits(data, study = s)
        if (!is.null(scope$trait)) {
            trAvail <- intersect(trAvail, scope$trait)
        }
        if (!is.null(traitIds)) {
            if (is.character(traitIds)) {
                trAvail <- intersect(trAvail, traitIds)
            } else if (is.list(traitIds) && is_in(s, names(traitIds))) {
                tv <- traitIds[[s]]
                if (is.character(tv)) trAvail <- intersect(trAvail, tv)
            }
        }
        traitsOut[[s]] <- trAvail
    }
    list(studies = studies, contexts = contextsOut, traits = traitsOut)
}


# Build a (variants x tupleRows) Z matrix from a QtlSumStats subset,
# requiring all rows to share an identical SNP order (the post-
# summaryStatsQc contract). Returns list(Z, nVec, variantIds).
# `errorLabel` is woven into the SNP-order error to identify the caller.
# @noRd
.buildJointSumstatZMatrix <- function(data, tupleRows, colLabels, errorLabel) {
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    firstDf <- getSumstatDf(
        data,
        study = studyCol[[tupleRows[[1L]]]],
        context = contextCol[[tupleRows[[1L]]]],
        trait = traitCol[[tupleRows[[1L]]]],
        require = c("SNP", "Z", "N")
    )
    variantIds <- firstDf$variant_id
    Z <- matrix(
        NA_real_,
        nrow = length(variantIds),
        ncol = length(tupleRows),
        dimnames = list(variantIds, colLabels)
    )
    nVec <- numeric(length(tupleRows))
    for (kk in seq_along(tupleRows)) {
        i <- tupleRows[[kk]]
        d <- getSumstatDf(
            data,
            study = studyCol[[i]],
            context = contextCol[[i]],
            trait = traitCol[[i]],
            require = c("SNP", "Z", "N")
        )
        if (!identical(d$variant_id, variantIds)) {
            msg <- glue(
                "{errorLabel}: every entry in a joint group must share an ",
                "identical SNP order after summaryStatsQc()."
            )
            abort(msg)
        }
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
# --- .buildIndividualCrossContextXy helpers ---------------------------------

# Scoped contexts in which a trait is present; NULL if fewer than 2.
.crossContextPerTrait <- function(data, tid, scopedContexts, verbose, label) {
    perTraitContexts <- character(0)
    for (cx in scopedContexts) {
        se <- getPhenotypes(data, contexts = cx)
        if (is_in(tid, rownames(se))) {
            perTraitContexts <- c(perTraitContexts, cx)
        }
    }
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
            msg <- glue(
                "{label}: trait '{tid}' has too few shared samples across ",
                "contexts; skipping."
            )
            inform(msg)
        }
        return(NULL)
    }
    X <- X[commonSamples, , drop = FALSE]
    Y <- .crossContextY(Yres, perTraitContexts, commonSamples)
    keep <- stats::complete.cases(Y)
    if (sum(keep) < 2L) {
        if (verbose >= 1) {
            msg <- glue(
                "{label}: trait '{tid}' has too few complete-Y subjects; ",
                "skipping."
            )
            inform(msg)
        }
        return(NULL)
    }
    list(
        X = X[keep, , drop = FALSE],
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
    region = NULL
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
    X <- .buildResidGeno(data, perTraitContexts, tid, cisWindow, region)
    Yres <- .fmResidPheno(data, contexts = perTraitContexts, traitId = tid)
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
    region = NULL
) {
    se <- getPhenotypes(data, contexts = cx)
    # scopedTraits is already region-restricted upstream (.runJointSpecs) when
    # region mode is used without an explicit traitId.
    traitsHere <- intersect(scopedTraits, rownames(se))
    if (.crossTraitTooFew(traitsHere, cx, study, verbose, label)) {
        return(NULL)
    }
    X <- .buildResidGeno(data, cx, traitsHere, cisWindow, region)
    Y <- .fmResidPheno(data, contexts = cx, traitId = traitsHere)
    common <- intersect(rownames(X), rownames(Y))
    if (length(common) < 2L) {
        return(NULL)
    }
    X <- X[common, , drop = FALSE]
    Y <- Y[common, , drop = FALSE]
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
.buildResidGeno <- function(data, contexts, traitId, cisWindow, region) {
    if (is.null(region)) {
        .fmResidGeno(
            data,
            contexts = contexts,
            traitId = traitId,
            cisWindow = cisWindow
        )
    } else {
        .fmResidGeno(data, contexts = contexts, region = region)
    }
}

# --- .buildComposedIndividualXy helpers -------------------------------------

# Enumerate the in-scope (context, trait) tuples for a study; NULL if < 2.
.composedTuples <- function(data, scope, study, verbose, label) {
    scopedContexts <- scope$contexts[[study]]
    scopedTraits <- scope$traits[[study]]
    tuples <- list()
    for (cx in scopedContexts) {
        se <- getPhenotypes(data, contexts = cx)
        # scopedTraits is region-restricted upstream (.runJointSpecs) if needed.
        for (tid in intersect(scopedTraits, rownames(se))) {
            tuples[[length(tuples) + 1L]] <- list(context = cx, trait = tid)
        }
    }
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
.composedYCols <- function(YresList, tuples, commonSamples) {
    yCols <- list()
    for (t in tuples) {
        ym <- YresList[[t$context]]
        if (!is_in(t$trait, colnames(ym))) {
            next
        }
        col <- ym[commonSamples, t$trait, drop = FALSE]
        colnames(col) <- str_c(t$context, t$trait, sep = ":")
        yCols[[length(yCols) + 1L]] <- col
    }
    if (length(yCols) < 2L) {
        return(NULL)
    }
    exec(cbind, !!!yCols)
}

.buildComposedIndividualXy <- function(
    data,
    scope,
    study,
    cisWindow,
    verbose,
    label,
    region = NULL
) {
    tuples <- .composedTuples(data, scope, study, verbose, label)
    if (is.null(tuples)) {
        return(NULL)
    }
    allContexts <- unique(map_chr(tuples, "context"))
    allTraits <- unique(map_chr(tuples, "trait"))
    X <- .buildResidGeno(data, allContexts, allTraits, cisWindow, region)
    YresList <- .fmResidPheno(data, contexts = allContexts, traitId = allTraits)
    if (length(allContexts) == 1L) {
        YresList <- set_names(list(YresList), allContexts)
    }
    commonSamples <- reduce(
        c(list(rownames(X)), map(YresList, rownames)),
        intersect
    )
    if (length(commonSamples) < 2L) {
        return(NULL)
    }
    X <- X[commonSamples, , drop = FALSE]
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
    # region / traitPos are GRanges provenance columns: carry them through
    # uncoerced so by-key / cross-study rebuilds keep the fine-mapping window
    # and trait position instead of silently dropping them.
    list(
        jointStudies = .jointStrCol("jointStudies", df),
        jointContexts = .jointStrCol("jointContexts", df),
        jointTraits = .jointStrCol("jointTraits", df),
        region = .jointRawCol("region", df),
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
.msDriverAccumulate <- function(
    qtlDatasets,
    sumStats,
    perStudyFn,
    sumStatsFn,
    cfg,
    rbindFn
) {
    out <- NULL
    embeddedLd <- NULL
    for (qdName in names(qtlDatasets)) {
        res <- perStudyFn(qtlDatasets[[qdName]], cfg)
        if (!is.null(res)) {
            out <- if (is.null(out)) res else rbindFn(out, res, ldSketch = NULL)
        }
    }
    if (!is.null(sumStats)) {
        ssRes <- sumStatsFn(sumStats, cfg)
        if (!is.null(ssRes)) {
            embeddedLd <- getLdSketch(ssRes)
            out <- if (is.null(out)) {
                ssRes
            } else {
                rbindFn(out, ssRes, ldSketch = embeddedLd)
            }
        }
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
            entry = as.list(out$entry)
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
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    methodArgs = list(),
    xRegions = list(NULL),
    twasWeights = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    cvFolds = 0,
    cvThreads = 1,
    samplePartition = NULL,
    pipCutoffToSkip = 0,
    fineMappingResult = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE,
    seed = NULL
) {
    # Run the joint dispatch once per region block, then merge per
    # (study, context, trait, method) across regions. A single block (cis or
    # jointRegions=TRUE concatenated) returns its result directly.
    args <- as.list(environment())
    args$xRegions <- NULL
    perRegion <- map(xRegions, .fmDispatchJointSpecRegion, args = args)
    perRegion <- compact(perRegion)
    if (length(perRegion) == 0L) {
        return(NULL)
    }
    if (length(perRegion) == 1L) {
        return(perRegion[[1L]])
    }
    .fmMergeResultsByKey(perRegion)
}

# FmJointPipeline for individual-level fine-mapping, built from the call params.
.fmJointPipeline <- function(args) {
    new(
        "FmJointPipeline",
        config = c(
            args[c(
                "coverage",
                "secondaryCoverage",
                "signalCutoff",
                "minAbsCorr",
                "dataDrivenPriorWeightsCutoff",
                "cvFolds",
                "cvThreads",
                "samplePartition",
                "verbose",
                "fullFit",
                "fullFitAlphaOnly",
                "includeAllCs",
                "seed"
            )],
            list(ldSketch = NULL)
        )
    )
}

.fmDispatchJointSpecsQtlDatasetOneRegion <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    cisWindow,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    methodArgs = list(),
    region = NULL,
    twasWeights = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    cvFolds = 0,
    cvThreads = 1,
    samplePartition = NULL,
    pipCutoffToSkip = 0,
    fineMappingResult = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE,
    seed = NULL
) {
    # Engine routing (jointEngine.R); one region block (the caller loops
    # regions).
    .jointRejectStudyOnIndividual(parsedJointSpec)
    pipeline <- .fmJointPipeline(as.list(environment()))
    .runJointSpecs(
        parsedJointSpec,
        data,
        dataForm = "individual",
        pipeline = pipeline,
        jointMethods = intersect(methods, c("mvsusie", "fsusie")),
        contexts = contexts,
        traitIds = traitIds,
        args = list(
            twasWeights = twasWeights,
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
.fmSumStatsPipeline <- function(args) {
    new(
        "FmJointPipeline",
        config = c(
            args[c(
                "coverage",
                "secondaryCoverage",
                "signalCutoff",
                "minAbsCorr",
                "dataDrivenPriorWeightsCutoff",
                "verbose",
                "fullFit",
                "fullFitAlphaOnly",
                "includeAllCs"
            )],
            list(cvFolds = 0L, ldSketch = getLdSketch(args$data))
        )
    )
}

.fmDispatchJointSpecsQtlSumStats <- function(
    parsedJointSpec,
    data,
    methods,
    contexts,
    traitIds,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    methodArgs = list(),
    twasWeights = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    fineMappingResult = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE
) {
    # Engine routing (jointEngine.R): the dispatch table + .runJointCell replace
    # the per-axis switch + the cross-context/trait/study/composed leaf
    # dispatchers.
    pipeline <- .fmSumStatsPipeline(as.list(environment()))
    .runJointSpecs(
        parsedJointSpec,
        data,
        dataForm = "sumstats",
        pipeline = pipeline,
        jointMethods = intersect(methods, c("mvsusie", "fsusie")),
        contexts = contexts,
        traitIds = traitIds,
        args = list(
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            methodArgs = methodArgs,
            verbose = verbose,
            cache = fineMappingResult
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
.fmMultiStudyQtlLoop <- function(nonStudyAxisSpecs, qtlDatasets, args) {
    out <- NULL
    if (length(nonStudyAxisSpecs) == 0L) {
        return(out)
    }
    fwd <- args[c(
        "methods",
        "contexts",
        "traitIds",
        "cisWindow",
        "coverage",
        "secondaryCoverage",
        "signalCutoff",
        "minAbsCorr",
        "verbose",
        "methodArgs",
        "xRegions",
        "twasWeights",
        "dataDrivenPriorWeightsCutoff"
    )]
    for (qdName in names(qtlDatasets)) {
        qdArgs <- c(
            list(nonStudyAxisSpecs, qtlDatasets[[qdName]]),
            fwd
        )
        qdRes <- exec(.fmDispatchJointSpecsQtlDataset, !!!qdArgs)
        if (!is.null(qdRes)) {
            out <- if (is.null(out)) {
                qdRes
            } else {
                .rbindFineMappingResult(out, qdRes, ldSketch = NULL)
            }
        }
    }
    out
}

# Fine-map all specs on the sumstats collection; rbind onto `out`.
.fmMultiStudySumStats <- function(
    parsedJointSpec,
    sumStats,
    studyAxisSpecs,
    out,
    args,
    verbose
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
    fwd <- args[c(
        "methods",
        "contexts",
        "traitIds",
        "coverage",
        "secondaryCoverage",
        "signalCutoff",
        "minAbsCorr",
        "verbose",
        "methodArgs",
        "twasWeights",
        "dataDrivenPriorWeightsCutoff"
    )]
    ssArgs <- c(
        list(parsedJointSpec, sumStats),
        fwd
    )
    ssRes <- exec(.fmDispatchJointSpecsQtlSumStats, !!!ssArgs)
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
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    methodArgs = list(),
    xRegions = list(NULL),
    twasWeights = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10
) {
    args <- as.list(environment())
    qtlDatasets <- getQtlDatasets(data)
    sumStats <- getSumStats(data)
    specs <- .fmSplitStudyAxisSpecs(parsedJointSpec)
    .fmMultiStudyWarnExcluded(specs$study, qtlDatasets, verbose)
    out <- .fmMultiStudyQtlLoop(specs$nonStudy, qtlDatasets, args)
    .fmMultiStudySumStats(
        parsedJointSpec,
        sumStats,
        specs$study,
        out,
        args,
        verbose
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
    retainFit = TRUE,
    retainFitDetail = "slim",
    seed = NULL
) {
    # Run the joint dispatch once per region block, then merge per
    # (study, context, trait, method) across regions. A single block (cis or
    # jointRegions=TRUE concatenated) returns its result directly.
    perRegion <- map(
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
        retainFit = retainFit,
        retainFitDetail = retainFitDetail,
        seed = seed
    )
    labs <- map_chr(xRegions, .twasRegionLabel)
    keep <- !map_lgl(perRegion, is.null)
    perRegion <- perRegion[keep]
    labs <- labs[keep]
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
    retainFit = TRUE,
    retainFitDetail = "slim",
    seed = NULL
) {
    # Engine routing (jointEngine.R); one region block (the caller loops
    # regions).
    .jointRejectStudyOnIndividual(parsedJointSpec)
    pipeline <- new(
        "TwasJointPipeline",
        config = list(
            retainFitDetail = retainFitDetail,
            dataType = dataType,
            cvFolds = 0L,
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
    retainFit = TRUE,
    retainFitDetail = "slim"
) {
    # Engine routing (jointEngine.R).
    pipeline <- new(
        "TwasJointPipeline",
        config = list(
            retainFitDetail = retainFitDetail,
            dataType = dataType,
            cvFolds = 0L,
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
        args = list(methodArgs = list(), verbose = verbose)
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
.twasMultiStudyQtlLoop <- function(nonStudyAxisSpecs, qtlDatasets, args) {
    out <- NULL
    if (length(nonStudyAxisSpecs) == 0L) {
        return(out)
    }
    fwd <- args[c(
        "methods",
        "contexts",
        "traitIds",
        "cisWindow",
        "dataType",
        "verbose",
        "xRegions",
        "retainFit",
        "retainFitDetail",
        "seed"
    )]
    for (qdName in names(qtlDatasets)) {
        qdArgs <- c(
            list(nonStudyAxisSpecs, qtlDatasets[[qdName]]),
            fwd
        )
        qdRes <- exec(.twasDispatchJointSpecsQtlDataset, !!!qdArgs)
        if (!is.null(qdRes)) {
            out <- if (is.null(out)) {
                qdRes
            } else {
                .rbindTwasWeights(out, qdRes, ldSketch = NULL)
            }
        }
    }
    out
}

# Learn weights for all specs on the sumstats collection; rbind onto `out`.
.twasMultiStudySumStats <- function(
    parsedJointSpec,
    sumStats,
    studyAxisSpecs,
    out,
    args,
    verbose
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
    fwd <- args[c(
        "methods",
        "contexts",
        "traitIds",
        "dataType",
        "verbose",
        "retainFit",
        "retainFitDetail"
    )]
    ssArgs <- c(
        list(parsedJointSpec, sumStats),
        fwd
    )
    ssRes <- exec(.twasDispatchJointSpecsQtlSumStats, !!!ssArgs)
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
    retainFit = TRUE,
    retainFitDetail = "slim",
    seed = NULL
) {
    args <- as.list(environment())
    qtlDatasets <- getQtlDatasets(data)
    sumStats <- getSumStats(data)
    specs <- .fmSplitStudyAxisSpecs(parsedJointSpec)
    .twasMultiStudyWarnExcluded(specs$study, qtlDatasets, verbose)
    out <- .twasMultiStudyQtlLoop(specs$nonStudy, qtlDatasets, args)
    .twasMultiStudySumStats(
        parsedJointSpec,
        sumStats,
        specs$study,
        out,
        args,
        verbose
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
    ym <- Yres[[cx]][commonSamples, , drop = FALSE]
    colnames(ym) <- cx
    ym
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
    if (length(hit)) r$entry[[hit[[1L]]]] else NULL
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
.fmDispatchJointSpecRegion <- function(rg, args) {
    regionArgs <- c(args, list(region = rg))
    exec(.fmDispatchJointSpecsQtlDatasetOneRegion, !!!regionArgs)
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
    retainFit,
    retainFitDetail,
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
        retainFit = retainFit,
        retainFitDetail = retainFitDetail,
        seed = seed
    )
}
