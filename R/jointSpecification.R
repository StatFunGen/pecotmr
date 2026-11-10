# =============================================================================
# Joint-specification grammar and ragged input-argument parsing for the
# fineMapping / twasWeights pipelines, plus the joint-specification
# dispatchers and multi-region / multi-study merges they share.
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
        return(studyName(data))
    }
    if (is(data, "QtlSumStats")) {
        return(unique(as.character(data$study)))
    }
    if (is(data, "MultiStudyQtlDataset")) {
        indStudies <- names(data)
        ss <- sumStats(data)
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
.spStudyDataForm <- function(data, studyName) {
    if (is(data, "QtlDataset")) {
        if (!identical(studyName, studyName(data))) {
            dataStudy <- studyName(data)
            msg <- glue(
                ".spStudyDataForm: study '{studyName}' not in QtlDataset ",
                "(study='{dataStudy}')"
            )
            abort(msg)
        }
        return("individual")
    }
    if (is(data, "QtlSumStats")) {
        if (!is_in(studyName, unique(as.character(data$study)))) {
            msg <- glue(
                ".spStudyDataForm: study '{studyName}' not in QtlSumStats"
            )
            abort(msg)
        }
        return("sumstats")
    }
    if (is(data, "MultiStudyQtlDataset")) {
        if (is_in(studyName, names(data))) {
            return("individual")
        }
        ss <- sumStats(data)
        if (!is.null(ss) && is_in(studyName, unique(as.character(ss$study)))) {
            return("sumstats")
        }
        msg <- glue(
            ".spStudyDataForm: study '{studyName}' not in MultiStudyQtlDataset"
        )
        abort(msg)
    }
    cls <- class(data)[[1L]]
    msg <- glue(".spStudyDataForm: unsupported class: {cls}")
    abort(msg)
}

# Concatenate a list of pieces into one, the way `c()` did, but safe on the
# empty case (`list_c()` has no zero-length identity to return).
# @noRd
.spConcat <- function(pieces, empty = list()) {
    if (length(pieces) == 0L) {
        return(empty)
    }
    list_c(pieces)
}

# Return character vector of contexts in `data` (across all studies when
# `study = NULL`, or for one study otherwise).
# @noRd
.spListContexts <- function(data, studyName = NULL) {
    if (is(data, "QtlDataset")) {
        if (!is.null(studyName) && !identical(studyName, studyName(data))) {
            return(character(0))
        }
        return(contexts(data))
    }
    if (is(data, "QtlSumStats")) {
        if (is.null(studyName)) {
            return(unique(as.character(data$context)))
        }
        return(unique(as.character(
            data$context[as.character(data$study) == studyName]
        )))
    }
    if (is(data, "MultiStudyQtlDataset")) {
        indDatasets <- as.list(data)
        ss <- sumStats(data)
        if (is.null(studyName)) {
            fromInd <- .spConcat(
                map(indDatasets, contexts),
                empty = character(0)
            )
            fromSs <- if (is.null(ss)) {
                character(0)
            } else {
                unique(as.character(ss$context))
            }
            return(unique(c(fromInd, fromSs)))
        }
        if (is_in(studyName, names(indDatasets))) {
            return(contexts(indDatasets[[studyName]]))
        }
        if (!is.null(ss) && is_in(studyName, unique(as.character(ss$study)))) {
            return(unique(as.character(
                ss$context[as.character(ss$study) == studyName]
            )))
        }
        return(character(0))
    }
    cls <- class(data)[[1L]]
    msg <- glue(".spListContexts: unsupported class: {cls}")
    abort(msg)
}

# The trait ids present in one context of a QtlDataset.
# @noRd
.spTraitsInContext <- function(context, data) {
    rownames(.qtlContextSe(data, context))
}

# Traits available in a single individual-level QtlDataset (optionally scoped).
.spListTraitsQtlDataset <- function(data, studyName, context) {
    if (!is.null(studyName) && !identical(studyName, studyName(data))) {
        return(character(0))
    }
    if (is.null(context)) {
        return(unique(unname(list_c(
            map(contexts(data), .spTraitsInContext, data = data)
        ))))
    }
    # Checked here rather than left to the accessor: `.spListTraits` answers
    # "which traits are in this scope", and an absent context is an empty
    # answer, not an error. The raw slot read this replaced returned NULL for
    # an unknown context; .qtlContextSe() rejects one, so the empty case
    # has to
    # be stated instead of falling out of a NULL.
    if (!is_in(context, contexts(data))) {
        return(character(0))
    }
    rownames(.qtlContextSe(data, context))
}

# Traits across a MultiStudyQtlDataset (individual studies + sumstats).
.spListTraitsMultiStudy <- function(data, studyName, context) {
    indDatasets <- as.list(data)
    ss <- sumStats(data)
    # No study filter: aggregate across every component (individual + sumstats).
    # Must precede per-study branches so a present sumStats slot does not shadow
    # the individual-level studies' traits when study = NULL.
    if (is.null(studyName)) {
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
    if (is_in(studyName, names(indDatasets))) {
        return(.spListTraits(indDatasets[[studyName]], context = context))
    }
    if (!is.null(ss) && is_in(studyName, unique(as.character(ss$study)))) {
        return(.spListTraits(ss, studyName = studyName, context = context))
    }
    character(0)
}

# Return character vector of traits in `data` (filtered by study and/or
# context when supplied).
# @noRd
.spListTraits <- function(data, studyName = NULL, context = NULL) {
    if (is(data, "QtlDataset")) {
        return(.spListTraitsQtlDataset(data, studyName, context))
    }
    if (is(data, "QtlSumStats")) {
        byStudy <- if (is.null(studyName)) {
            rep(TRUE, nrow(data))
        } else {
            as.character(data$study) == studyName
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
        return(.spListTraitsMultiStudy(data, studyName, context))
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

# @noRd
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

# --- parseContexts helpers --------------------------------------------------

# Vector form: apply uniformly to every study, filtering to availability. A
# context no study has is a mistake rather than a gap -- with one study it
# is the only way to name a missing context -- so that is an error; one that
# only some studies lack is dropped for those, with a warning.
.parseContextsVec <- function(contexts, studies, data) {
    if (length(contexts) == 0L) {
        msg <- glue(
            "`contexts` must be NULL or a non-empty character vector ",
            "(or named list)."
        )
        abort(msg)
    }
    unknown <- setdiff(contexts, .spListContexts(data))
    if (length(unknown) > 0L) {
        unknownStr <- str_flatten(unknown, ", ")
        msg <- glue("`contexts` names context(s) no study has: {unknownStr}")
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

# @noRd
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
    missing <- setdiff(v2, .spListTraits(data, studyName = s, context = cx))
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
    missing <- setdiff(val, .spListTraits(data, studyName = s))
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

# @noRd
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
# Reading a parsed selection. The pipelines normalize `contexts` / `traitId`
# once, at entry, with parseContexts() / parseTraitIds(); these answer "what
# does it select here?" for one study, context or row. A study or context the
# selection does not mention keeps everything it has.
# -----------------------------------------------------------------------------

# The contexts selected for study `s` (possibly none), or NULL for all.
# @noRd
.selectedContexts <- function(contexts, s) {
    if (is.null(contexts)) NULL else contexts[[s]]
}

# The traits selected for study `s` in context `cx`, or NULL for all.
# @noRd
.selectedTraits <- function(traitId, s, cx) {
    if (!is.list(traitId)) {
        return(traitId)
    }
    entry <- traitId[[s]]
    if (is.list(entry)) entry[[cx]] else entry
}

# TRUE when trait `tid` is selected for study `s` in context `cx`.
# @noRd
.traitSelectedAt <- function(tid, traitId, s, cx) {
    sel <- .selectedTraits(traitId, s, cx)
    is.null(sel) || is_in(tid, sel)
}

# The traits study `s` takes part with across `studyContexts`, the contexts in
# play for it -- what a fit spanning those contexts draws on: the union of
# their selections, or NULL when any of them keeps all of its traits.
# @noRd
.selectedStudyTraits <- function(traitId, s, studyContexts) {
    entry <- if (is.list(traitId)) traitId[[s]] else traitId
    if (!is.list(entry)) {
        return(entry)
    }
    if (!all(is_in(studyContexts, names(entry)))) {
        return(NULL)
    }
    unique(as.character(list_c(unname(entry[studyContexts]))))
}

# TRUE when some study's traits are chosen context by context. Only then can
# a (context, trait) pair inside a study's selected traits still be left out.
# @noRd
.hasNestedTraitSelection <- function(traitId) {
    is.list(traitId) && any(map_lgl(traitId, is.list))
}

# Which rows of a summary-statistics collection the parsed `contexts` /
# `traitId` selection keeps.
# @noRd
.sumStatsRowsSelected <- function(data, contexts, traitId) {
    map_lgl(
        seq_len(nrow(data)),
        .sumStatsRowSelectedAt,
        studyCol = as.character(data$study),
        contextCol = as.character(data$context),
        traitCol = as.character(data$trait),
        contexts = contexts,
        traitId = traitId
    )
}

# Whether the selection keeps row `i`.
# @noRd
.sumStatsRowSelectedAt <- function(
    i,
    studyCol,
    contextCol,
    traitCol,
    contexts,
    traitId
) {
    s <- studyCol[[i]]
    cx <- contextCol[[i]]
    ctxSel <- .selectedContexts(contexts, s)
    (is.null(ctxSel) || is_in(cx, ctxSel)) &&
        .traitSelectedAt(traitCol[[i]], traitId, s, cx)
}

# Study `s`'s share of a parsed `traitId`, in a form the single-study entry
# point parses back to the same selection: unchanged when it is one vector
# for every study, else that study's entry, keyed by the study when it is
# chosen context by context. NULL when the study keeps all of its traits.
# @noRd
.studyTraitIdArg <- function(traitId, s) {
    if (!is.list(traitId)) {
        return(traitId)
    }
    entry <- traitId[[s]]
    if (is.list(entry)) set_names(list(entry), s) else entry
}

# The traits a run selects in context `cx` of study `qd`: its share of the
# parsed `traitId` (all of them when the selection leaves the context out),
# else -- when no trait was named at all -- those `region` overlaps, else all.
# @noRd
.contextSelectedTraits <- function(cx, qd, traitId, region) {
    se <- .qtlContextSe(qd, context = cx)
    ids <- rownames(se)
    if (!is.null(traitId)) {
        intersect(ids, .selectedTraits(traitId, studyName(qd), cx) %||% ids)
    } else if (!is.null(region)) {
        rr <- SummarizedExperiment::rowRanges(se)
        ids[IRanges::overlapsAny(rr, region)]
    } else {
        ids
    }
}

# Whether context `cx` of study `qd` has a trait selected.
# @noRd
.contextHasSelectedTraits <- function(cx, qd, traitId, region) {
    length(.contextSelectedTraits(cx, qd, traitId, region)) > 0L
}

# Whether study `qd` has a trait selected in a context it runs on. In a
# multi-study run a study without one sits out, as one with none of its
# contexts selected does, rather than failing the run.
# @noRd
.studyHasSelectedTraits <- function(qd, contexts, traitId, region) {
    useCtx <- .selectedContexts(contexts, studyName(qd)) %||% contexts(qd)
    any(map_lgl(
        useCtx,
        .contextHasSelectedTraits,
        qd = qd,
        traitId = traitId,
        region = region
    ))
}

# A multi-study run in which no component has a selected trait is refused:
# each one without sits out, so the run would otherwise come back empty. An
# embedded summary-statistics row counts when the selection keeps it. `cfg`
# is the run's record, carrying the parsed selection and its `region`.
# @noRd
.msAssertTraitsSelected <- function(data, cfg, pipeline) {
    inStudies <- map_lgl(
        as.list(data),
        .studyHasSelectedTraits,
        contexts = cfg$contexts,
        traitId = cfg$traitId,
        region = cfg$region
    )
    ss <- sumStats(data)
    inSumStats <- !is.null(ss) &&
        any(.sumStatsRowsSelected(ss, cfg$contexts, cfg$traitId))
    if (!any(inStudies) && !inSumStats) {
        abort(glue("{pipeline}(MultiStudyQtlDataset): no traits selected."))
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
            data = data
        ),
        studies
    )
    traitsOut <- set_names(
        map(
            studies,
            .fmScopeStudyTraits,
            scopeLimit = scope$trait,
            traitIds = traitIds,
            contextsOut = contextsOut,
            data = data
        ),
        studies
    )
    # A joint fit spans a study's contexts, so its traits are the union of
    # the per-context selections; when the selection is per context, the
    # cell runner then drops each (context, trait) pair it leaves out.
    list(
        studies = studies,
        contexts = contextsOut,
        traits = traitsOut,
        traitSelection = if (.hasNestedTraitSelection(traitIds)) traitIds
    )
}

# `.spListTraits` takes the study by name, so give it the same (data, s)
# shape the context lister has.
# @noRd
.fmScopeTraits <- function(data, s) {
    .spListTraits(data, studyName = s)
}

# What a pipeline-level `contexts` argument narrows study `s` to, or NULL when
# it says nothing about it: one vector for every study, or parseContexts()'s
# per-study list.
# @noRd
.fmScopeUserLimit <- function(spec, s) {
    if (is.list(spec)) .selectedContexts(spec, s) else spec
}

# One axis for one study: everything it has, narrowed by the spec's scope and
# then by the caller's argument. Each narrowing is an intersect, so the chain
# is a fold rather than a variable rewritten in place.
# @noRd
.fmScopeAxis <- function(s, avail, scopeLimit, userSpec, data) {
    limits <- compact(list(scopeLimit, .fmScopeUserLimit(userSpec, s)))
    reduce(limits, intersect, .init = avail(data, s))
}

# Study `s`'s traits: everything it has, narrowed by the spec's scope and by
# the caller's `traitId` over the contexts in play for it.
# @noRd
.fmScopeStudyTraits <- function(s, scopeLimit, traitIds, contextsOut, data) {
    limits <- compact(list(
        scopeLimit,
        .selectedStudyTraits(traitIds, s, contextsOut[[s]])
    ))
    reduce(limits, intersect, .init = .fmScopeTraits(data, s))
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
    ldReference = NULL,
    cutoffs = NULL
) {
    cols <- list(
        study = as.character(data$study),
        context = as.character(data$context),
        trait = as.character(data$trait)
    )
    allDf <- .sumStatsDf(
        data,
        studyName = cols$study[[tupleRows[[1L]]]],
        context = cols$context[[tupleRows[[1L]]]],
        trait = cols$trait[[tupleRows[[1L]]]],
        require = c("SNP", "Z", "N")
    )
    # Every entry shares one SNP order (asserted below), so filtering the first
    # entry's ids filters the group: Z is built against this vector and each
    # entry is checked against it.
    firstDf <- allDf[
        .panelKeepMask(allDf$variant_id, ldReference, cutoffs, errorLabel),
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

# The z-scores of studies each on their own LD reference
# (.enumCrossStudyGroup()), over the union of their variants: each study's
# variants kept by its own reference's panel filters, the union matched by
# allele (.xsVariantUnion()), and each study's z signed into the union's
# allele frame, NA at a variant it lacks. `read` holds the data, the panel
# cutoffs and the label for messages.
# @noRd
.buildJointSumstatZUnion <- function(tupleRows, cols, refs, read) {
    entries <- map2(tupleRows, refs, .jointUnionEntry, cols = cols, read = read)
    ids <- .xsVariantUnion(map(entries, "variantId"))
    Z <- exec(cbind, !!!map(entries, .jointUnionColumn, ids = ids))
    list(
        Z = `dimnames<-`(Z, list(ids, cols$study[tupleRows])),
        nVec = map_dbl(entries, "n"),
        variantIds = ids
    )
}

# One study's variants, z-scores and median N, its variants kept by its own
# reference's panel filters.
# @noRd
.jointUnionEntry <- function(row, ref, cols, read) {
    d <- .sumStatsDf(
        read$data,
        studyName = cols$study[[row]],
        context = cols$context[[row]],
        trait = cols$trait[[row]],
        require = c("SNP", "Z", "N")
    )
    kept <- d[
        .panelKeepMask(d$variant_id, ref, read$cutoffs, read$label),
        ,
        drop = FALSE
    ]
    list(
        variantId = kept$variant_id,
        z = kept$z,
        n = stats::median(kept$N, na.rm = TRUE)
    )
}

# A study's z-scores over the union `ids`, signed into its allele frame and
# NA at a variant the study lacks.
# @noRd
.jointUnionColumn <- function(entry, ids) {
    m <- matchVariants(entry$variantId, ids, removeStrandAmbiguous = FALSE)
    at <- match(seq_along(ids), m$idxB)
    entry$z[m$idxA[at]] * m$sign[at]
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
    d <- .sumStatsDf(
        data,
        studyName = cols$study[[i]],
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


# --- .crossContextPlan helpers ----------------------------------------------

# Does context `cx` carry trait `tid`?
# @noRd
.ccContextHasTrait <- function(cx, data, tid) {
    is_in(tid, rownames(.qtlContextSe(data, context = cx)))
}

# Scoped contexts in which a trait is present; NULL if fewer than 2.
# @noRd
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

# Skip a trait with too few of something to fit: NULL, and a message saying so.
# `report` holds the caller's verbose, label and tid.
# @noRd
.crossContextSkip <- function(tooFew, report) {
    if (report$verbose >= 1) {
        inform(glue(
            "{report$label}: trait '{report$tid}' has too few {tooFew}; ",
            "skipping."
        ))
    }
    NULL
}

# The samples a cross-context fit is on: every context's Y holds them and none
# is missing a value. NULL, with a message, when fewer than two remain.
# @noRd
.crossContextFitSamples <- function(Yres, perTraitContexts, report) {
    shared <- reduce(map(Yres, rownames), intersect)
    if (length(shared) < 2L) {
        return(.crossContextSkip("shared samples across contexts", report))
    }
    Y <- .crossContextY(Yres, perTraitContexts, shared)
    keep <- stats::complete.cases(Y)
    if (sum(keep) < 2L) {
        return(.crossContextSkip("complete-Y subjects", report))
    }
    shared[keep]
}

# X and the response matrix on the fit's samples that X holds.
# @noRd
.crossContextAssemble <- function(
    X,
    Yres,
    perTraitContexts,
    fitSamples,
    report
) {
    common <- intersect(rownames(X), fitSamples)
    if (length(common) < 2L) {
        return(.crossContextSkip("shared samples across contexts", report))
    }
    list(
        X = X[common, , drop = FALSE],
        Y = .crossContextY(Yres, perTraitContexts, common),
        perTraitContexts = perTraitContexts
    )
}

# Build a multi-context Y matrix for a single (study, trait) from an
# individual-level QtlDataset, before any genotype is read. Returns the trait's
# contexts, their residualized phenotypes and the samples the fit is on, or NULL
# when fewer than 2 contexts carry `tid` or the sample / complete-Y subset is
# too small to fit.
# @noRd
.crossContextPlan <- function(
    data,
    tid,
    scopedContexts,
    verbose,
    label,
    residualizationParam,
    phenotypeFilter = NULL
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
    Yres <- .fmResidPheno(
        data,
        contexts = perTraitContexts,
        traitId = tid,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam
    )
    report <- list(verbose = verbose, label = label, tid = tid)
    fitSamples <- .crossContextFitSamples(Yres, perTraitContexts, report)
    if (is.null(fitSamples)) {
        return(NULL)
    }
    list(
        perTraitContexts = perTraitContexts,
        Yres = Yres,
        fitSamples = fitSamples,
        report = report
    )
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

# Warn + signal skip when a context has fewer than 2 scoped traits.
# @noRd
.crossTraitTooFew <- function(traitsHere, cx, studyName, verbose, label) {
    if (length(traitsHere) >= 2L) {
        return(FALSE)
    }
    if (verbose >= 1) {
        nTraits <- length(traitsHere)
        msg <- glue(
            "{label}: context '{cx}' (study '{studyName}') has {nTraits} ",
            "scoped trait(s); skipping."
        )
        inform(msg)
    }
    TRUE
}

# Build a multi-trait Y matrix for a single (study, context) from an
# individual-level QtlDataset. Returns list(X, Y, traitsHere, se) or NULL
# when fewer than 2 traits live in the context or the sample / complete-Y
# subset is too small. Before any genotype is read.
# @noRd
.crossTraitPlan <- function(
    data,
    cx,
    scopedTraits,
    studyName,
    verbose,
    label,
    residualizationParam,
    phenotypeFilter = NULL
) {
    se <- .qtlContextSe(data, context = cx)
    # scopedTraits is already region-restricted upstream (.runJointSpecs) when
    # region mode is used without an explicit traitId.
    traitsHere <- intersect(scopedTraits, rownames(se))
    if (.crossTraitTooFew(traitsHere, cx, studyName, verbose, label)) {
        return(NULL)
    }
    allY <- .fmResidPheno(
        data,
        contexts = cx,
        traitId = traitsHere,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam
    )
    # The samples observed for every trait, which are the ones the fit is on.
    fitSamples <- rownames(allY)[stats::complete.cases(allY)]
    if (length(fitSamples) < 2L) {
        return(NULL)
    }
    list(
        se = se,
        traitsHere = traitsHere,
        allY = allY,
        fitSamples = fitSamples
    )
}

# The cross-trait X / Y on the samples they share, keeping only the samples
# observed for every trait; NULL when fewer than two remain at either step.
# @noRd
.crossTraitAssemble <- function(allX, allY, traitsHere, se) {
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


# --- .composedPlan helpers ---------------------------------------------------

# Enumerate the in-scope (context, trait) tuples for a study; NULL if < 2.
# @noRd
.composedTuple <- function(tid, cx) {
    list(context = cx, trait = tid)
}

# The (context, trait) tuples context `cx` contributes.
# @noRd
.composedTuplesForContext <- function(cx, data, scopedTraits) {
    se <- .qtlContextSe(data, context = cx)
    # scopedTraits is region-restricted upstream (.runJointSpecs) if needed.
    map(intersect(scopedTraits, rownames(se)), .composedTuple, cx = cx)
}

.composedTuples <- function(data, scope, studyName, verbose, label) {
    scopedContexts <- scope$contexts[[studyName]]
    scopedTraits <- scope$traits[[studyName]]
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
                "{label}: study '{studyName}' has {nTuples} (context, trait) ",
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

# The samples a composed fit is on: every context's Y holds them and no tuple's
# column misses them. Empty when the tuples cannot be fit together.
# @noRd
.composedFitSamples <- function(YresList, tuples) {
    shared <- reduce(map(YresList, rownames), intersect)
    Y <- if (length(shared) >= 2L) {
        .composedYCols(YresList, tuples, shared)
    }
    if (is.null(Y)) {
        return(character(0))
    }
    shared[stats::complete.cases(Y)]
}

# A composed-axes (context, trait) fit on an individual-level QtlDataset, before
# any genotype is read: its tuples, the per-context residualized phenotypes and
# the samples it is on, or NULL when too few samples could be fit to read any
# genotypes for them.
# @noRd
.composedPlan <- function(
    data,
    scope,
    studyName,
    verbose,
    label,
    residualizationParam,
    phenotypeFilter = NULL
) {
    tuples <- .composedTuples(data, scope, studyName, verbose, label)
    if (is.null(tuples)) {
        return(NULL)
    }
    allContexts <- unique(map_chr(tuples, "context"))
    allTraits <- unique(map_chr(tuples, "trait"))
    resid <- .fmResidPheno(
        data,
        contexts = allContexts,
        traitId = allTraits,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam
    )
    # A single context returns the bare matrix rather than a named list.
    YresList <- if (length(allContexts) == 1L) {
        set_names(list(resid), allContexts)
    } else {
        resid
    }
    fitSamples <- .composedFitSamples(YresList, tuples)
    if (length(fitSamples) < 2L) {
        return(NULL)
    }
    list(
        tuples = tuples,
        allContexts = allContexts,
        allTraits = allTraits,
        YresList = YresList,
        fitSamples = fitSamples
    )
}

# X and the tuples' Y columns on the samples both sides share, without the
# samples any Y column misses; NULL when fewer than two samples remain.
# @noRd
.composedAlignXy <- function(sides, tuples) {
    commonSamples <- reduce(
        c(list(rownames(sides$X)), map(sides$YresList, rownames)),
        intersect
    )
    if (length(commonSamples) < 2L) {
        return(NULL)
    }
    X <- sides$X[commonSamples, , drop = FALSE]
    Y <- .composedYCols(sides$YresList, tuples, commonSamples)
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

# Identity key (study/context/trait/method and the joint key, joined by "\r"),
# used to align per-region result entries when merging: two joint groups
# sharing a tuple stay apart. Shared by the fm/twas mergers.
# @noRd
.mergeResultKeyOf <- function(r) {
    joint <- map(.JOINT_COLUMNS, .fmrIdentityColumn, x = r)
    exec(
        str_c,
        as.character(r$study),
        as.character(r$context),
        as.character(r$trait),
        as.character(r$method),
        !!!map(joint, coalesce, ""),
        sep = "\r"
    )
}

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
        # Merged entries sit on the rows they were keyed to, references and
        # all.
        list(ldReference = .rowLdRefs(base))
    )
    .fmWithMergedJointIndex(
        exec(QtlFineMappingResult, !!!qfmrArgs),
        results,
        baseKeys
    )
}

# `merged` with each row's position in its joint fit's condition axis: the
# one its regions agree on, else NA (a region's screen can drop different
# conditions), as .twasMergedJointIndex() keeps it. No column when no region
# records one.
# @noRd
.fmWithMergedJointIndex <- function(merged, results, baseKeys) {
    if (!any(map_lgl(results, .fmrHasJointIndex))) {
        return(merged)
    }
    `$<-`(
        merged,
        "jointIndex",
        map_int(baseKeys, .fmMergedJointIndexFor, results = results)
    )
}

# @noRd
.fmrHasJointIndex <- function(r) {
    is_in("jointIndex", colnames(mcols(r)))
}

# @noRd
.fmMergedJointIndexFor <- function(key, results) {
    v <- unique(discard(
        map_int(results, .fmJointIndexForKey, key = key),
        is.na
    ))
    if (length(v) == 1L) v else NA_integer_
}

# The position recorded by the row of result `r` matching merge key `key`; NA
# when it has none.
# @noRd
.fmJointIndexForKey <- function(r, key) {
    hit <- which(.mergeResultKeyOf(r) == key)
    if (length(hit) == 0L) {
        return(NA_integer_)
    }
    as.integer(.fmrJointIndexAt(r, hit[[1L]]))
}

# One passthrough column of a joint result row as a character vector (the joint-
# key columns), or NULL when absent.
# @noRd
.jointStrCol <- function(nm, rows) {
    if (is_in(nm, colnames(mcols(rows)))) as.character(mcols(rows)[[nm]])
}

# One passthrough column carried through UNCOERCED (GRanges provenance columns
# region / traitPos), or NULL when absent.
# @noRd
.jointRawCol <- function(nm, rows) {
    if (is_in(nm, colnames(mcols(rows)))) mcols(rows)[[nm]]
}

# The optional passthrough columns of a per-tuple result row, as a named list
# (NULL for any absent column): the three joint-key columns (jointStudies /
# jointContexts / jointTraits) plus the GRanges provenance columns (region /
# traitPos). Spliced into the QtlFineMappingResult / TwasWeights constructors so
# by-key / cross-study rebuilds preserve them.
.jointCols <- function(rows) {
    # traitPos is a GRanges provenance column: carry it through uncoerced so
    # by-key / cross-study rebuilds keep the trait position instead of silently
    # dropping it. There is no `region` counterpart -- the fine-mapping window
    # retired (section 4.4) because the element's own span is the region and a
    # stored window had no correct update rule when variants are narrowed.
    list(
        jointStudies = .jointStrCol("jointStudies", rows),
        jointContexts = .jointStrCol("jointContexts", rows),
        jointTraits = .jointStrCol("jointTraits", rows),
        traitPos = .jointRawCol("traitPos", rows)
    )
}

# --- .multiStudyPipelineDriver helpers --------------------------------------

# rbind two results with no LD sketch -- the per-study combine.
# @noRd
.msRbindNoLd <- function(acc, x, rbindFn) {
    rbindFn(acc, x, ldReference = NULL)
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

# Accumulate per-study + sumstats results; every row keeps its own LD
# reference (none for an individual-level study's rows).
# @noRd
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
        return(fromStudies)
    }
    if (is.null(fromStudies)) {
        return(ssRes)
    }
    rbindFn(fromStudies, ssRes)
}

# Shared tail of the MultiStudyQtlDataset fineMapping / twasWeights pipeline
# methods: recurse into each embedded QtlDataset then the embedded QtlSumStats,
# row-bind the per-study results, and combine them with any
# joint-specification result. `perStudyFn` / `sumStatsFn` are the
# (method-specific) single-study recursions; `rbindFn` / `pipelineName` supply
# the per-pipeline pieces. Each method keeps its own
# (divergent) method-gating / joint-dispatch preamble and passes the computed
# `jointResult` in.
# @noRd
.multiStudyPipelineDriver <- function(
    data,
    jointResult,
    perStudyFn,
    sumStatsFn,
    cfg,
    rbindFn,
    pipelineName,
    noun = "a result"
) {
    perTupleResult <- .msDriverAccumulate(
        as.list(data),
        sumStats(data),
        perStudyFn,
        sumStatsFn,
        cfg,
        rbindFn
    )
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
    rbindFn(perTupleResult, jointResult)
}

# The jointSpecification of a run on individual-level data that names none:
# both pipelines group every joint method (mvsusie, fsusie, mr.mash) by it,
# so a TWAS run reads the very fits a fine-mapping run made, through the SAME
# engine as an explicit jointSpecification. >= 2 traits -> cross-trait (covers
# multi-trait single-context AND both-multi, since the cross-trait enumerator
# iterates contexts -> per-context multi-trait fits); else >= 2 contexts
# (single trait) -> cross-context. Single context & single trait -> no joint
# (the callers' multivariate guards already reject the joint methods there).
# Returns a list of parsed specs (full scope) or list().
# @noRd
.defaultJointSpec <- function(nCtx, nTraits) {
    if (nTraits >= 2L) {
        list(list(axes = "trait", scope = NULL))
    } else if (nCtx >= 2L) {
        list(list(axes = "context", scope = NULL))
    } else {
        list()
    }
}

# The same for summary statistics: across contexts.
# @noRd
.defaultSumStatsJointSpec <- function() {
    list(list(axes = "context", scope = NULL))
}
# ---- fine-mapping joint dispatch --------------------------------------------
#
# Every fine-mapping joint dispatch hands the engine a pipeline config and an
# `args` list built from ONE record: the pipeline's own joint record
# (.fmIndividualJointCfg() / .fmQssConfigs()), or, for a component of a
# MultiStudyQtlDataset, that record with the component's data in it.
# Re-listing the settings at each dispatch layer is how they went missing --
# the per-study dispatch dropped cross-validation, residualization, the
# screen, fitRetention and the seed, and passed the rest positionally into
# the wrong slots.

# Run joint specs on a QtlDataset, once per region block, merging per (study,
# context, trait, method) across the blocks. A study axis has no meaning on
# individual-level data, and rejecting it up front keeps the error ahead of
# the fit.
# @noRd
.fmDispatchJointSpecsQtlDataset <- function(
    parsedJointSpec,
    tokens,
    methodArgs,
    cfg
) {
    .qtlRunWork(
        .fmJointSpecsWork(parsedJointSpec, tokens, methodArgs, cfg),
        cfg$data,
        cfg$genotypeFilter,
        cfg$residualizationParam,
        withAf = TRUE
    )
}

# .fmDispatchJointSpecsQtlDataset()'s fits, planned, for a run with others to
# fit alongside.
# @noRd
.fmJointSpecsWork <- function(parsedJointSpec, tokens, methodArgs, cfg) {
    .jointRejectStudyOnIndividual(parsedJointSpec)
    .fmJointRegionsWork(parsedJointSpec, tokens, methodArgs, cfg)
}

# Individual-level joint specs' fits, planned, one region block at a time.
# @noRd
.fmJointRegionsWork <- function(parsedJointSpec, tokens, methodArgs, cfg) {
    .qtlWorkOf(
        map(
            cfg$xRegions,
            .fmJointSpecRegionWork,
            parsedJointSpec = parsedJointSpec,
            tokens = tokens,
            methodArgs = methodArgs,
            cfg = cfg
        ),
        .fmJointRegionsResult
    )
}

# @noRd
.fmJointRegionsResult <- function(perRegion) {
    .fmCombineRegionResults(compact(perRegion))
}

# One region block of .fmJointSpecsWork().
# @noRd
.fmJointSpecRegionWork <- function(
    region,
    parsedJointSpec,
    tokens,
    methodArgs,
    cfg
) {
    .jointSpecsWork(
        parsedJointSpec,
        cfg$data,
        dataForm = "individual",
        pipeline = .fmJointPipeline(cfg, cfg$crossValidationParam),
        jointMethods = intersect(tokens, c("mvsusie", "fsusie")),
        contexts = cfg$contexts,
        traitIds = cfg$traitId,
        args = .fmJointIndividualArgs(cfg, methodArgs, region)
    )
}

# The FmJointPipeline a joint run fits under: the settings the fitters read
# off the pipeline. Summary statistics pass the no-CV record (there are no
# samples to hold out); individual-level data pass the caller's resolved CV
# record. The LD a group is fitted against is the group's own
# (.jointRowsLdRef()), not a setting of the run.
# @noRd
.fmJointPipeline <- function(cfg, crossValidationParam) {
    new(
        "FmJointPipeline",
        config = list(
            credibleSetParam = cfg$credibleSetParam,
            dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff,
            crossValidationParam = crossValidationParam,
            verbose = cfg$verbose,
            fitRetention = cfg$fitRetention,
            seed = cfg$seed
        )
    )
}

# The engine `args` of an individual-level joint run over one region block.
# Residualization rides here rather than on the pipeline config: the group
# enumerators read it while extracting X and Y.
# @noRd
.fmJointIndividualArgs <- function(cfg, methodArgs, region) {
    list(
        mrmashPrior = cfg$mrmashPrior,
        dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff,
        methodArgs = methodArgs,
        cisWindow = cfg$cisWindow,
        region = region,
        verbose = cfg$verbose,
        pipCutoffToSkip = cfg$screen,
        cache = cfg$fineMappingResult,
        residualizationParam = cfg$residualizationParam,
        genotypeFilter = cfg$genotypeFilter,
        phenotypeFilter = cfg$phenotypeFilter
    )
}

# Run joint specs on a QtlSumStats collection. `cfg` is the summary-statistics
# joint record (.fmQssConfigs()$joint), or a MultiStudy record with the
# collection as its data.
# @noRd
.fmDispatchJointSpecsQtlSumStats <- function(
    parsedJointSpec,
    tokens,
    methodArgs,
    cfg
) {
    .runJointSpecs(
        parsedJointSpec,
        cfg$data,
        dataForm = "sumstats",
        pipeline = .fmJointPipeline(cfg, .cvResolve(CrossValidationParam())),
        jointMethods = intersect(tokens, c("mvsusie", "fsusie")),
        contexts = cfg$contexts,
        traitIds = cfg$traitId,
        args = list(
            mrmashPrior = cfg$mrmashPrior,
            dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff,
            methodArgs = methodArgs,
            verbose = cfg$verbose,
            cache = cfg$fineMappingResult,
            cutoffs = .panelCutoffs(
                cfg$panelFilterParam %||% PanelFilterParam()
            )
        )
    )
}

# Joint dispatch for a MultiStudyQtlDataset. Per component AND per axis: a
# spec with a `study` axis alone joins the individual-level studies with one
# another (.fmDispatchStudySpecsIndividual()), while the other specs run on
# every individual-level component -- each under the same settings, with the
# component's own data (its genotype-filter overrides applied, as the
# per-study recursion applies them). Every spec also runs on the summary
# statistics, a composed one including the study only there. `cfg` is the
# individual-level joint record plus `genotypeFilterParam` and
# `panelFilterParam`, with the collection as data. `components = FALSE` leaves
# the components' own fits to their own runs (.fmPerStudy()), which have other
# fits to share their reads with; a cross-study fit belongs to no one run.
# @noRd
.fmDispatchJointSpecsMultiStudy <- function(
    parsedJointSpec,
    tokens,
    methodArgs,
    cfg,
    components = TRUE
) {
    qtlDatasets <- as.list(cfg$data)
    specs <- .fmSplitStudyAxisSpecs(parsedJointSpec)
    .fmMultiStudyWarnExcluded(specs$study, qtlDatasets, cfg$verbose)
    perStudy <- if (length(specs$nonStudy) == 0L || !components) {
        NULL
    } else {
        map(
            qtlDatasets,
            .fmStudyJointResult,
            specs = specs$nonStudy,
            tokens = tokens,
            methodArgs = methodArgs,
            cfg = cfg
        )
    }
    acrossStudies <- .fmDispatchStudySpecsIndividual(
        keep(specs$study, .jsIsStudyAxisOnly),
        tokens,
        methodArgs,
        cfg
    )
    .fmMultiStudySumStats(
        parsedJointSpec,
        specs$study,
        .msFoldStudyResults(
            c(perStudy, list(acrossStudies)),
            .rbindFineMappingResult
        ),
        tokens,
        methodArgs,
        cfg
    )
}

# The study-axis specs on a MultiStudyQtlDataset's individual-level studies:
# one mvsusie fit per (context, trait) two or more of them carry, their samples
# stacked with each study's phenotype missing outside its own rows
# (.enumCrossStudyIndividual()), once per region block. Each study is read with
# its genotype-filter overrides applied, as its own fits read it. NULL when no
# spec applies or fewer than two individual-level studies are present. fsusie
# has no cross-study fit: its functional domain is a study's traits.
# @noRd
.fmDispatchStudySpecsIndividual <- function(
    studySpecs,
    tokens,
    methodArgs,
    cfg
) {
    jointTokens <- intersect(tokens, "mvsusie")
    if (length(studySpecs) == 0L || length(cfg$data) < 2L) {
        return(NULL)
    }
    if (is_in("fsusie", tokens) && (cfg$verbose %||% 1) >= 1) {
        inform(str_c(
            "jointCrossStudy: fsusie has no cross-study fit, so the ",
            "individual-level studies are fitted across studies by mvsusie ",
            "alone."
        ))
    }
    if (length(jointTokens) == 0L) {
        return(NULL)
    }
    studies <- MultiStudyQtlDataset(
        qtlDatasets = map(
            as.list(cfg$data),
            .qtlRestrictSamples,
            keepSamples = cfg$genotypeFilter$keepSamples
        )
    )
    .qtlRunWork(
        .fmJointRegionsWork(
            studySpecs,
            jointTokens,
            methodArgs,
            list_assign(cfg, data = studies)
        ),
        studies,
        cfg$genotypeFilter,
        cfg$residualizationParam,
        withAf = TRUE
    )
}

# One individual-level component's joint fits: the shared record with the
# component as its data.
# @noRd
.fmStudyJointResult <- function(qd, specs, tokens, methodArgs, cfg) {
    studyCfg <- list_assign(
        cfg,
        data = .qtlRestrictSamples(qd, cfg$genotypeFilter$keepSamples)
    )
    .fmDispatchJointSpecsQtlDataset(specs, tokens, methodArgs, studyCfg)
}

# A MultiStudyQtlDataset's summary statistics as its joint fits read them:
# with `cfg$cisWindow` laid on the tuples the parsed selection keeps, as the
# per-study recursion has its QtlSumStats run lay it. The joint fits read the
# collection directly rather than recursing, so they lay it here.
# @noRd
.msJointSumStats <- function(ss, cfg, caller) {
    rows <- which(.sumStatsRowsSelected(ss, cfg$contexts, cfg$traitId))
    .qssApplyCisWindow(ss, cfg$cisWindow, caller, rows)
}

# Fine-map every spec on the summary-statistics component and rbind onto the
# individual-level results `out`. Without one, a composed spec including the
# study axis has nothing to fit: the individual-level studies join one another
# on the study axis alone.
# @noRd
.fmMultiStudySumStats <- function(
    parsedJointSpec,
    studyAxisSpecs,
    out,
    tokens,
    methodArgs,
    cfg
) {
    sumStats <- sumStats(cfg$data)
    if (is.null(sumStats)) {
        composed <- discard(studyAxisSpecs, .jsIsStudyAxisOnly)
        if (length(composed) > 0L && cfg$verbose >= 1) {
            msg <- glue(
                "jointCrossStudy: no sumStats slot present on this ",
                "MultiStudyQtlDataset; composed specs including the study ",
                "axis produce no result."
            )
            inform(msg)
        }
        return(out)
    }
    ssRes <- .fmDispatchJointSpecsQtlSumStats(
        parsedJointSpec,
        tokens,
        methodArgs,
        list_assign(
            cfg,
            data = .msJointSumStats(sumStats, cfg, "fineMappingPipeline")
        )
    )
    if (is.null(ssRes)) {
        return(out)
    }
    if (is.null(out)) {
        return(ssRes)
    }
    .rbindFineMappingResult(out, ssRes)
}


# Merge the per-region joint results by (study, context, trait, method) and
# joint key. A single block -- cis, or jointRegions = TRUE concatenated -- is
# returned directly rather than merged with itself.
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


# --- .fmDispatchJointSpecsMultiStudy helpers --------------------------------

# The joint specs a component's own run fits, each scope narrowed to what the
# component holds, as a scope only names what some component does. A spec that
# leaves it nothing does not apply to it.
# @noRd
.msStudySpecs <- function(specs, qd) {
    compact(map(specs, .msStudySpec, qd = qd))
}

# @noRd
.msStudySpec <- function(spec, qd) {
    if (is.null(spec$scope)) {
        return(spec)
    }
    scope <- imap(spec$scope, .msScopeValuesIn, qd = qd)
    if (any(lengths(scope) == 0L)) {
        return(NULL)
    }
    list_assign(spec, scope = scope)
}

# @noRd
.msScopeValuesIn <- function(values, key, qd) {
    intersect(
        values,
        switch(
            key,
            study = .spListStudies(qd),
            context = .spListContexts(qd),
            trait = .spListTraits(qd)
        )
    )
}

# Partition joint specs into those with a `study` axis and the rest.
.fmSplitStudyAxisSpecs <- function(parsedJointSpec) {
    hasStudy <- map_lgl(parsedJointSpec, .jsHasStudyAxis)
    list(
        study = parsedJointSpec[hasStudy],
        nonStudy = parsedJointSpec[!hasStudy]
    )
}

# Note that individual-level studies are excluded from composed specs
# including the study axis: they join one another on the study axis alone.
.fmMultiStudyWarnExcluded <- function(studyAxisSpecs, qtlDatasets, verbose) {
    composed <- discard(studyAxisSpecs, .jsIsStudyAxisOnly)
    if (length(composed) > 0L && length(qtlDatasets) > 0L && verbose >= 1) {
        qdNames <- str_flatten(names(qtlDatasets), ", ")
        msg <- glue(
            "jointCrossStudy: excluding individual-level studies ",
            "({qdNames}) from composed specs including the study axis; ",
            "sumstats studies participate. Individual-level studies join one ",
            "another on the study axis alone."
        )
        inform(msg)
    }
}


# =============================================================================
# TWAS-weights dispatchers
# =============================================================================

# Merge per-region TwasWeights collections into one by concatenating each
# entry across the regions that fitted it (.twasMergeRegionEntries: stacked
# weights + flat per-region cvResult). Entries are matched on (study,
# context, trait, joint*, method) over EVERY region, so one that a region
# lacks is kept and two joint groups stay apart. A merged entry keeps its
# trait position and LD reference, and its window spans the windows it was
# fitted over.
# @noRd
.twasMergeResultsByKey <- function(results, regionLabels) {
    inputs <- map(results, .twasAllCellInputs)
    stacked <- list(
        keys = bind_rows(map(inputs, "keys")),
        method = list_c(map(inputs, "method"), ptype = character()),
        entry = list_flatten(map(inputs, "entry")),
        traitPos = .twasJoinRanges(map(inputs, "traitPos"), inputs),
        region = .twasJoinRanges(map(inputs, "region"), inputs),
        ldRefs = list_flatten(map(inputs, "ldRefs")),
        jointIndex = .twasJoinJointIndex(inputs),
        block = rep(seq_along(inputs), map_int(inputs, .twasInputCells))
    )
    groups <- unname(split(
        seq_along(stacked$method),
        as.integer(vctrs::vec_group_id(
            bind_cols(stacked$keys, .method = stacked$method)
        ))
    ))
    first <- map_int(groups, .twasFirstCell)
    merged <- list(
        keys = stacked$keys[first, ],
        method = stacked$method[first],
        entry = map(
            groups,
            .twasMergedGroupEntry,
            stacked = stacked,
            regionLabels = regionLabels
        ),
        traitPos = stacked$traitPos[first],
        region = .twasMergedWindows(stacked$region, groups),
        ldRefs = stacked$ldRefs[first],
        jointIndex = .twasMergedJointIndex(stacked$jointIndex, groups)
    )
    .withMetadata(.twasNewValidated(merged), S4Vectors::metadata(results[[1L]]))
}

# The merged entry of one group of stacked cells (one entry, several regions).
# @noRd
.twasMergedGroupEntry <- function(group, stacked, regionLabels) {
    .twasMergeRegionEntries(
        stacked$entry[group],
        regionLabels[stacked$block[group]]
    )
}

# Each group's joint-fit position: the one its regions agree on, else NA (a
# region's screen can drop different conditions). NULL when no cell has one.
# @noRd
.twasMergedJointIndex <- function(jointIndex, groups) {
    if (is.null(jointIndex)) {
        return(NULL)
    }
    map_int(groups, .twasGroupJointIndex, jointIndex = jointIndex)
}

# @noRd
.twasGroupJointIndex <- function(group, jointIndex) {
    v <- unique(jointIndex[group])
    if (length(v) == 1L) v else NA_integer_
}

# Each group's window: the span of its regions' windows (NULL for none).
# @noRd
.twasMergedWindows <- function(region, groups) {
    if (is.null(region)) {
        return(NULL)
    }
    windowSpans <- map(groups, .twasWindowSpan, region = region)
    suppressWarnings(exec(c, !!!windowSpans))
}

# @noRd
.twasWindowSpan <- function(group, region) {
    gr <- region[group]
    known <- gr[!.twasIsPlaceholder(gr)]
    if (length(known) == 0L) {
        return(gr[1L])
    }
    span <- range(known, ignore.strand = TRUE)
    if (length(span) != 1L) {
        abort(glue(
            "TwasWeights: an entry was fitted over windows on ",
            "{length(span)} chromosomes; a row's window is one range."
        ))
    }
    span
}

# ---- TWAS joint dispatch ----------------------------------------------------
#
# An explicit jointSpecification fits under the same pipeline marker and the
# same engine `args` as the per-tuple multivariate grid (.twasJointSetup() /
# .twasRunJointSpecGrid()), so cross-validation, residualization,
# fitRetention, the per-method options and the mash prior reach an explicit
# joint fit exactly as they reach the rest. The explicit path used to build a
# bare pipeline record of its own -- default cross-validation, no
# fitRetention -- and to hand the engine no residualization and no options.

# Explicit joint specs (mr.mash) on a QtlDataset.
# @noRd
.twasDispatchJointSpecsQtlDataset <- function(
    parsedJointSpec,
    tokens,
    marker,
    ctx,
    contexts,
    traitIds
) {
    .qtlRunWork(
        .twasJointSpecsWork(
            parsedJointSpec,
            tokens,
            marker,
            ctx,
            contexts,
            traitIds
        ),
        ctx$data,
        ctx$genotypeFilter,
        ctx$residualizationParam
    )
}

# .twasDispatchJointSpecsQtlDataset()'s fits, planned, for a run with others to
# fit alongside: the joint methods' (.twasJointTokens()), and those of the
# methods also fitted jointly on a group of two or more conditions
# (`multitraitTokens`, .twasMultitraitTokens()).
# @noRd
.twasJointSpecsWork <- function(
    parsedJointSpec,
    tokens,
    marker,
    ctx,
    contexts,
    traitIds
) {
    .jointRejectStudyOnIndividual(parsedJointSpec)
    .twasJointSpecGridWork(
        parsedJointSpec,
        union(.twasJointTokens(tokens), .jpConfig(marker)$multitraitTokens),
        contexts,
        traitIds,
        marker,
        ctx
    )
}

# Joint specs on a QtlSumStats collection: mr.mash RSS fitted on each group,
# mvsusie read off its fine-mapping fit of the group. `cfg` is the
# summary-statistics joint record, or a MultiStudy record with the collection
# as its data; `methodArgs` are the summary-statistics path's per-method
# options.
# @noRd
.twasDispatchJointSpecsQtlSumStats <- function(
    parsedJointSpec,
    tokens,
    methodArgs,
    cfg
) {
    pipeline <- new(
        "TwasJointPipeline",
        config = list(
            crossValidationParam = .cvResolve(CrossValidationParam()),
            fitRetention = cfg$fitRetention,
            standardized = TRUE,
            verbose = cfg$verbose,
            seed = cfg$seed,
            studyCoupling = cfg$studyCoupling
        )
    )
    args <- list(
        methodArgs = methodArgs,
        verbose = cfg$verbose,
        fineMappingResult = cfg$fineMappingResult,
        cache = cfg$twasWeights,
        cutoffs = .panelCutoffs(cfg$panelFilterParam %||% PanelFilterParam())
    )
    joint <- .runJointSpecs(
        parsedJointSpec,
        cfg$data,
        dataForm = "sumstats",
        pipeline = pipeline,
        jointMethods = .twasJointTokens(tokens),
        contexts = cfg$contexts,
        traitIds = cfg$traitId,
        args = args
    )
    .msFoldStudyResults(
        list(joint, .twasCoupledSumStats(parsedJointSpec, tokens, cfg, args)),
        .rbindTwasWeights
    )
}

# The methods fitted across the summary-statistics studies
# (.jointTwasCoupledRss()) -- the penalized ones under the penalty on the
# difference between their effects, PRS-CS as JointPRS and SDPR as SDPRX fit
# ancestries -- and the ensemble stacking the studies' own fits
# (EnsembleParam(acrossStudies = TRUE), .twasSumstatEnsembleLayer()), on the
# specs with the study axis alone; their per-study fits run with the
# per-tuple ones. NULL when none applies.
# @noRd
.twasCoupledSumStats <- function(parsedJointSpec, tokens, cfg, args) {
    studySpecs <- .twasCoupledSpecs(parsedJointSpec, tokens, cfg)
    if (length(studySpecs) == 0L) {
        return(NULL)
    }
    pipeline <- new(
        "TwasJointPipeline",
        config = list(
            crossValidationParam = .cvResolve(CrossValidationParam()),
            fitRetention = cfg$fitRetention,
            standardized = TRUE,
            verbose = cfg$verbose,
            seed = cfg$seed,
            studyCoupling = cfg$studyCoupling,
            ensembleParam = list_assign(
                .ensembleResolve(cfg$ensemble),
                stackTokens = .twasStackTokens(tokens, cfg$ensemble, "sumstat")
            ),
            h2 = cfg$h2Spec
        )
    )
    .runJointSpecs(
        studySpecs,
        cfg$data,
        dataForm = "sumstats",
        pipeline = pipeline,
        jointMethods = .twasAcrossStudyTokens(tokens, cfg),
        contexts = cfg$contexts,
        traitIds = cfg$traitId,
        args = args
    )
}

# The specs a run fits methods across studies on: those with the study axis
# alone, when it couples a method across studies or stacks the studies' own
# fits.
# @noRd
.twasCoupledSpecs <- function(parsedJointSpec, tokens, cfg) {
    if (length(.twasAcrossStudyTokens(tokens, cfg)) == 0L) {
        return(list())
    }
    keep(.fmSplitStudyAxisSpecs(parsedJointSpec)$study, .jsIsStudyAxisOnly)
}

# The summary-statistics methods a run fits across studies: those it couples
# (`studyCoupling` given) and those its ensemble stacks across them.
# @noRd
.twasAcrossStudyTokens <- function(tokens, cfg) {
    coupled <- if (!is.null(cfg$studyCoupling)) {
        intersect(tokens, .TWAS_COUPLED_SUMSTATS)
    }
    union(coupled, .twasStackTokens(tokens, cfg$ensemble, "sumstat"))
}

# Joint dispatch for a MultiStudyQtlDataset, mirroring
# .fmDispatchJointSpecsMultiStudy(). Non-study specs run on every
# individual-level component, each fitting under the setup its own QtlDataset
# run would use (.twasJointSetup() on the component, its genotype-filter
# overrides applied); a spec with a `study` axis alone reads mvsusie off the
# fits across the individual-level studies
# (.twasDispatchStudySpecsIndividual()); every spec runs on the
# summary-statistics component. `cfg` is .twasJointCfg()'s record plus the
# collection, `genotypeFilterParam` and `panelFilterParam`; `methods` is the
# translated methods record.
# @noRd
.twasDispatchJointSpecsMultiStudy <- function(
    parsedJointSpec,
    tokens,
    methods,
    cfg,
    components = TRUE
) {
    qtlDatasets <- as.list(cfg$data)
    specs <- .fmSplitStudyAxisSpecs(parsedJointSpec)
    .twasMultiStudyWarnExcluded(specs$study, qtlDatasets, cfg$verbose)
    perStudy <- if (length(specs$nonStudy) == 0L || !components) {
        NULL
    } else {
        map(
            qtlDatasets,
            .twasStudyJointResult,
            specs = specs$nonStudy,
            tokens = tokens,
            methods = methods,
            cfg = cfg
        )
    }
    studyOnly <- keep(specs$study, .jsIsStudyAxisOnly)
    acrossStudies <- .msFoldStudyResults(
        list(
            .twasDispatchStudySpecsIndividual(studyOnly, tokens, cfg),
            .twasDispatchStudySpecsCoupled(studyOnly, tokens, methods, cfg)
        ),
        .rbindTwasWeights
    )
    .twasMultiStudySumStats(
        parsedJointSpec,
        specs$study,
        .msFoldStudyResults(
            c(perStudy, list(acrossStudies)),
            .rbindTwasWeights
        ),
        tokens,
        methods,
        cfg
    )
}

# The study-axis specs on a MultiStudyQtlDataset's individual-level studies:
# mvsusie's weights over each group fineMappingPipeline() fits across them,
# read off those fits (.twasAcrossStudyWeights()), each study taken with its
# genotype-filter overrides applied, as fine-mapping took it. NULL when no
# spec applies, fewer than two individual-level studies are present, or
# mvsusie is not asked for: mr.mash fits no individual-level group across
# studies, and fsusie has none.
# @noRd
.twasDispatchStudySpecsIndividual <- function(studySpecs, tokens, cfg) {
    if (
        length(studySpecs) == 0L ||
            length(cfg$data) < 2L ||
            !is_in("mvsusie", tokens)
    ) {
        return(NULL)
    }
    studies <- MultiStudyQtlDataset(
        qtlDatasets = map(
            as.list(cfg$data),
            .qtlRestrictSamples,
            keepSamples = cfg$genotypeFilter$keepSamples
        )
    )
    .twasAcrossStudyWeights(
        .twasStudyAxisGroups(studySpecs, studies, cfg),
        cfg
    )
}

# The study-axis specs on a MultiStudyQtlDataset's individual-level studies,
# for the univariate methods fitted across them under the penalty on the
# difference between their effects (lasso, enet, ridge, l0learn;
# .lassoCoupledWeights()), for the methods fitted across them jointly
# (bayesC and bayesN on BGLR; .twasMultitraitTokens()), and for the ensemble
# across studies (EnsembleParam(acrossStudies = TRUE)), which stacks each
# study's own fits of the methods it names (.twasStackTokens()): one group
# per (context, trait) two or more studies carry, their samples stacked as
# fineMappingPipeline() stacks them (.enumCrossStudyIndividual()), per region
# block, each study read with its genotype-filter overrides applied. The
# methods' per-study fits run as well, in the per-component recursion. NULL
# when no spec applies, fewer than two studies are present, or none of the
# coupling (`studyCoupling` and a coupled method), a joint method and the
# ensemble across studies is asked for.
# @noRd
.twasDispatchStudySpecsCoupled <- function(studySpecs, tokens, methods, cfg) {
    coupled <- if (!is.null(cfg$studyCoupling)) {
        intersect(tokens, names(.TWAS_COUPLED_IMPL))
    }
    stacked <- .twasStackTokens(tokens, cfg$ensembleParam, "individual")
    fitted <- union(
        union(coupled, stacked),
        .twasMultitraitTokens(tokens, cfg$backend)
    )
    if (
        length(studySpecs) == 0L ||
            length(cfg$data) < 2L ||
            length(fitted) == 0L
    ) {
        return(NULL)
    }
    studies <- MultiStudyQtlDataset(
        qtlDatasets = map(
            as.list(cfg$data),
            .qtlRestrictSamples,
            keepSamples = cfg$genotypeFilter$keepSamples
        )
    )
    setup <- .twasJointSetup(
        studies,
        .twasNormWithBackends(methods, cfg$backend),
        list_assign(cfg, stackTokens = stacked)
    )
    .qtlRunWork(
        .twasJointSpecGridWork(
            studySpecs,
            fitted,
            cfg$contexts,
            cfg$traitId,
            setup$marker,
            setup$ctx
        ),
        studies,
        cfg$genotypeFilter,
        cfg$residualizationParam,
        withAf = TRUE
    )
}

# The methods a run stacks across studies (EnsembleParam(acrossStudies =
# TRUE)): its univariate ones fitted from their own data on `path`
# ("individual" or "sumstat"), every fine-mapping method excepted, whose
# weights come from fits a group across studies does not hold. None unless
# the ensemble is on and across studies.
# @noRd
.twasStackTokens <- function(tokens, ensembleParam, path) {
    if (
        !isTRUE(ensembleParam$enabled) ||
            !isTRUE(ensembleParam$acrossStudies)
    ) {
        return(character(0))
    }
    keep(tokens, .twasStacksAcross, path = path)
}

# @noRd
.twasStacksAcross <- function(tk, path) {
    impl <- .twasMethodCapabilities[[tk]][[str_c(path, "Impl")]]
    !is.null(impl) &&
        is.null(.twasFineMappingMethodAdapters[[tk]]) &&
        !.twasIsMultivariateToken(tk)
}

# One individual-level component's joint fits.
# @noRd
.twasStudyJointResult <- function(qd, specs, tokens, methods, cfg) {
    setup <- .twasJointSetup(
        .qtlRestrictSamples(qd, cfg$genotypeFilter$keepSamples),
        .twasNormWithBackends(methods, cfg$backend),
        cfg
    )
    .twasDispatchJointSpecsQtlDataset(
        specs,
        tokens,
        setup$marker,
        setup$ctx,
        cfg$contexts,
        cfg$traitId
    )
}

# Fit every spec on the summary-statistics component and rbind onto the
# individual-level results `out`.
# @noRd
.twasMultiStudySumStats <- function(
    parsedJointSpec,
    studyAxisSpecs,
    out,
    tokens,
    methods,
    cfg
) {
    sumStats <- sumStats(cfg$data)
    if (is.null(sumStats)) {
        .twasMultiStudyNoSumStats(studyAxisSpecs, tokens, cfg$verbose)
        return(out)
    }
    ssData <- .msJointSumStats(sumStats, cfg, "twasWeightsPipeline")
    ssRes <- .twasDispatchJointSpecsQtlSumStats(
        parsedJointSpec,
        tokens,
        .twasQssSetMrmashPrior(
            .twasSumStatsMethodTokens(methods)$methodArgs,
            cfg$mashPrior,
            tokens
        ),
        list_assign(
            cfg,
            data = ssData,
            ensemble = .ensembleResolve(cfg$ensembleParam),
            h2Spec = .twasQssH2Spec(
                cfg$h2,
                list(
                    data = ssData,
                    contexts = cfg$contexts,
                    traitId = cfg$traitId
                )
            )
        )
    )
    if (is.null(ssRes)) {
        return(out)
    }
    if (is.null(out)) {
        return(ssRes)
    }
    .rbindTwasWeights(out, ssRes)
}


# --- .twasDispatchJointSpecsMultiStudy helpers ------------------------------

# Note that individual-level studies are excluded from cross-study TWAS fits.
.twasMultiStudyWarnExcluded <- function(studyAxisSpecs, qtlDatasets, verbose) {
    composed <- discard(studyAxisSpecs, .jsIsStudyAxisOnly)
    if (length(composed) > 0L && length(qtlDatasets) > 0L && verbose >= 1) {
        qdNames <- str_flatten(names(qtlDatasets), ", ")
        msg <- glue(
            "jointCrossStudy (twas): excluding individual-level studies ",
            "({qdNames}) from composed specs including the study axis; ",
            "sumstats studies participate. Individual-level studies join one ",
            "another on the study axis alone."
        )
        inform(msg)
    }
}

# Without summary statistics, a composed spec including the study axis has
# nothing to fit, and neither has mr.mash on the study axis: it fits no
# group of individual-level studies, which join one another only on the
# study axis alone (mvsusie read off fine-mapping fits, the coupled and
# joint methods, the ensemble across studies).
# @noRd
.twasMultiStudyNoSumStats <- function(studyAxisSpecs, tokens, verbose) {
    composed <- discard(studyAxisSpecs, .jsIsStudyAxisOnly)
    mrmashAcross <- length(composed) < length(studyAxisSpecs) &&
        is_in("mrmash", tokens)
    if ((length(composed) == 0L && !mrmashAcross) || verbose < 1) {
        return(invisible(NULL))
    }
    inform(glue(
        "jointCrossStudy (twas): no sumStats slot present on this ",
        "MultiStudyQtlDataset; composed specs including the study axis, and ",
        "mr.mash across studies, produce no result."
    ))
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
    .collectionEntry(r, hit[[1L]])
}

# The merged FM entry for base row `i`: gather that key's entry from every
# region
# and fold them together.
# @noRd
.fmMergedEntryAt <- function(i, results, baseKeys) {
    perRegion <- map(results, .mergeEntryForKey, key = baseKeys[[i]])
    .fmMergeEntries(compact(perRegion))
}

# TRUE when a parsed joint spec includes the `study` axis.
# @noRd
.jsHasStudyAxis <- function(s) {
    is_in("study", s$axes)
}

# TRUE when a parsed joint spec's one axis is the study.
# @noRd
.jsIsStudyAxisOnly <- function(s) {
    identical(s$axes, "study")
}
