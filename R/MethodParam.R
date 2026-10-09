# =============================================================================
# MethodParam S4 virtual base
# -----------------------------------------------------------------------------
# The other half of the argument-management split.
#
# A *Options() record is an argument bag bound for an external engine. Nothing
# in pecotmr defines what those arguments mean, so the only honest check is
# that each NAME exists in the engine's live formals; the values travel as
# given. That is MethodOptions.
#
# A *Param() record is the opposite case: a setting pecotmr itself reads. No
# upstream signature defines it, so there is nothing to validate names
# against -- but equally, nothing stops pecotmr from declaring each field's
# type and checking it, which is what Bioconductor's own *Param classes do.
# ScanBamParam, PileupParam and ScanVcfParam are each a plain setClass with
# typed slots and a validity method, inheriting from nothing.
#
# So: Options check names against someone else's signature, Params check
# types against their own.
#
# Why a shared virtual base rather than 26 unrelated classes: these records
# are read list-style throughout the pipelines -- `args$coverage`,
# `names(args)`, `as.list(args)` spliced into an engine. Declaring the slots
# typed and keeping the list view as methods on one base buys the typing
# without rewriting every read site.
# =============================================================================

#' @include AllGenerics.R
#' @include MethodOptions.R
#' @importFrom methods slot slotNames new is validObject
#' @importClassesFrom S4Vectors character_OR_NULL
NULL

# --- nullable slot types ---------------------------------------------------
#
# A typed slot cannot hold NULL, and "unset" is a real state for most of
# these fields -- `medianAbsCorr = NULL` means "do not apply one", not "zero".
# A union per base type is the standard S4 way to say "this type, or unset",
# and it keeps the type declared rather than retreating to "ANY".
#
# `character_OR_NULL` is not defined here: S4Vectors exports it already,
# with the same definition, so it is imported rather than shadowed.

setClassUnion("numeric_OR_NULL", c("numeric", "NULL"))
setClassUnion("logical_OR_NULL", c("logical", "NULL"))
setClassUnion("list_OR_NULL", c("list", "NULL"))

# The one cross-family slot: SusieRssParam() carries a SusieRssControlOptions()
# record, which is a MethodOptions rather than a Param -- `control` really is
# someone else's argument list, nested inside a pecotmr setting.
setClassUnion("MethodOptions_OR_NULL", c("MethodOptions", "NULL"))

#' @title Typed Settings For One Pipeline Stage
#' @description The virtual base of pecotmr's \code{*Param()} records --
#'   \code{\link{CredibleSetParam}}, \code{\link{SusieRssParam}},
#'   \code{\link{PanelFilterParam}} and the rest. Each one holds settings
#'   \strong{pecotmr itself reads}, declared as typed slots and checked when
#'   the record is built.
#'
#'   This is the counterpart to \code{\link{MethodOptions}}, which carries
#'   arguments bound for an \emph{external} engine. The distinction is
#'   deliberate and visible in the name: a \code{*Param()} is a pecotmr
#'   setting whose meaning and type pecotmr defines, so pecotmr validates it;
#'   a \code{*Options()} is somebody else's argument list, so pecotmr checks
#'   only that each name exists in that engine's live formals and otherwise
#'   does not attempt to mirror its semantics.
#' @section List interface:
#'   A \code{MethodParam} reads like a list: \code{$} and \code{[[} fetch a
#'   field, and \code{names()}, \code{length()} and \code{as.list()} describe
#'   it.
#'
#'   Those last three report only the fields that are \strong{set}. A slot
#'   holding \code{NULL} is absent from all three, while \code{$} and
#'   \code{[[} still return \code{NULL} for it. This is not cosmetic:
#'   \code{as.list()} on one of these records is spliced into an engine call,
#'   and an engine that distinguishes an explicit \code{NULL} argument from an
#'   omitted one would otherwise see a different call.
#'
#'   Unlike a list, \code{$} on a name that is not a field is an error rather
#'   than \code{NULL}. Silently reading a misspelled setting as "unset" is the
#'   failure these classes exist to prevent.
#' @aliases MethodParam
#' @name MethodParam-class
#' @docType class
#' @export
setClass("MethodParam", representation("VIRTUAL"))

# The fields a record actually carries, in declaration order: every slot
# except those holding NULL. The single definition behind names(), length()
# and as.list(), so the three can never disagree about what is set.
# @noRd
.paramSetFields <- function(x) {
    nms <- slotNames(x)
    keep(nms, .paramFieldIsSet, x = x)
}

# @noRd
.paramFieldIsSet <- function(nm, x) {
    !is.null(slot(x, nm))
}

# Reject a field name the class does not declare, naming the ones it does.
# @noRd
.paramAssertField <- function(x, name) {
    nms <- slotNames(x)
    if (is_in(name, nms)) {
        return(invisible(NULL))
    }
    abort(glue(
        "`{name}` is not a setting of {class(x)[[1L]]}. ",
        "It has: {str_flatten(nms, ', ')}."
    ))
}

#' @rdname MethodParam-class
#' @param x A \code{MethodParam} object.
#' @param name The name of a setting.
#' @return \code{$} returns that setting's value, \code{NULL} when it is
#'   unset.
#' @export
setMethod("$", "MethodParam", function(x, name) {
    .paramAssertField(x, name)
    slot(x, name)
})

#' @rdname MethodParam-class
#' @param i A setting's name, or its position among the settings that are
#'   set.
#' @param j Unused, for consistency with the generic.
#' @param ... Unused, for consistency with the generic.
#' @return \code{[[} returns that setting's value.
#' @export
setMethod("[[", "MethodParam", function(x, i, j, ...) {
    if (is.character(i)) {
        .paramAssertField(x, i)
        return(slot(x, i))
    }
    # A positional index counts over the fields that are set, matching what
    # names() and as.list() report.
    set <- .paramSetFields(x)
    if (
        !is.numeric(i) ||
            length(i) != 1L ||
            is.na(i) ||
            i < 1L ||
            i > length(set)
    ) {
        abort(glue(
            "subscript out of bounds: {class(x)[[1L]]} has ",
            "{length(set)} setting(s) set."
        ))
    }
    slot(x, set[[as.integer(i)]])
})

#' @rdname MethodParam-class
#' @return \code{names} returns the settings that are set.
#' @export
setMethod("names", "MethodParam", function(x) {
    .paramSetFields(x)
})

#' @rdname MethodParam-class
#' @return \code{length} returns how many settings are set.
#' @export
setMethod("length", "MethodParam", function(x) {
    length(.paramSetFields(x))
})

#' @rdname MethodParam-class
#' @return \code{as.list} returns the settings that are set, as a named list.
#' @export
setMethod("as.list", "MethodParam", function(x, ...) {
    set <- .paramSetFields(x)
    set_names(map(set, slot, object = x), set)
})

#' @rdname MethodParam-class
#' @param object A \code{MethodParam} object.
#' @return \code{show} is called for its side effect and returns
#'   \code{invisible(NULL)}.
#' @export
setMethod("show", "MethodParam", function(object) {
    cat(glue("<{class(object)[[1L]]}>"), "\n", sep = "")
    set <- .paramSetFields(object)
    unset <- setdiff(slotNames(object), set)
    if (length(set) == 0L) {
        cat("  (nothing set)\n")
    }
    for (nm in set) {
        cat(sprintf("  %-24s %s\n", nm, .paramSlotText(slot(object, nm))))
    }
    if (length(unset) > 0L) {
        cat("  unset: ", str_flatten(unset, ", "), "\n", sep = "")
    }
    invisible(NULL)
})

# =============================================================================
# Method-selection Params
# -----------------------------------------------------------------------------
# Which methods a pipeline runs, and each one's engine arguments.
#
# Three slots rather than one because a method running on both data paths
# reaches a DIFFERENT engine on each -- lasso is glmnet on individual data
# and lassosum on summary statistics -- so one bag of arguments per method
# cannot serve a MultiStudyQtlDataset carrying both. The split is declared
# by the user rather than guessed:
#
#   methods             entries whose input path need not be stated: a
#                       single-type run, or a method with nothing
#                       path-specific to configure (the fit-derived ones
#                       carry no arguments at all).
#   qtlDatasetMethods   entries bound for the individual-level path.
#   qtlSumStatsMethods  entries bound for the summary-statistics path.
#
# A method appears in exactly ONE of the three. That is per-method, not
# per-slot: a dual-type run tuning `lasso` differently per path while
# selecting `susie` with defaults has something to say in all three places,
# and forcing `susie` into both path slots would imply two configurations
# where there is only one fit.
# =============================================================================

#' @title Which Methods A Pipeline Runs, And How
#' @description The virtual base of
#'   \code{\link{FineMappingMethodsParam}} and
#'   \code{\link{TwasWeightsMethodsParam}}: which methods a pipeline runs,
#'   and each one's engine arguments.
#'
#'   Three slots rather than one because a method running on both data paths
#'   reaches a different engine on each --- lasso is glmnet on
#'   individual-level data and lassosum on summary statistics --- so one bag
#'   of arguments per method cannot serve a \code{MultiStudyQtlDataset}
#'   carrying both. \code{methods} holds entries whose input path need not
#'   be stated; \code{qtlDatasetMethods} and \code{qtlSumStatsMethods} hold
#'   the ones that must state it. Naming a method selects it.
#' @slot methods Named list of per-method options whose input path need not
#'   be stated.
#' @slot qtlDatasetMethods Named list of per-method options for the
#'   individual-level path.
#' @slot qtlSumStatsMethods Named list of per-method options for the
#'   summary-statistics path.
#' @seealso \code{\link{FineMappingMethodsParam}},
#'   \code{\link{TwasWeightsMethodsParam}}
#' @aliases MethodsSelectionParam-class
#' @name MethodsSelectionParam-class
#' @docType class
#' @exportClass MethodsSelectionParam
setClass(
    "MethodsSelectionParam",
    contains = c("MethodParam", "VIRTUAL"),
    slots = c(
        methods = "list_OR_NULL",
        qtlDatasetMethods = "list_OR_NULL",
        qtlSumStatsMethods = "list_OR_NULL"
    )
)

# The input path a slot is bound for; NULL for `methods`, which is the slot
# that does not state one.
# @noRd
.methodsSlotInput <- function(slotName) {
    switch(
        slotName,
        qtlDatasetMethods = "QtlDataset",
        qtlSumStatsMethods = "QtlSumStats",
        NULL
    )
}

# Which methods each family offers on a given path. Keyed on the class so
# the two concrete subclasses differ in nothing but their registry.
# @noRd
.methodsParamAvailable <- function(cls, inputKind) {
    if (identical(cls, "FineMappingMethodsParam")) {
        return(.fmMethodsFor(inputKind))
    }
    .twasMethodsFor(inputKind)
}

# Every method the family knows, on any path -- the vocabulary `methods`
# itself is checked against, since that slot does not name a path.
# @noRd
.methodsParamAnyMethod <- function(cls) {
    unique(unlist(map(
        c("QtlDataset", "QtlSumStats", "GwasSumStats"),
        .methodsParamAvailable,
        cls = cls
    )))
}

.methodsSlotNames <- c("methods", "qtlDatasetMethods", "qtlSumStatsMethods")

setValidity("MethodsSelectionParam", function(object) {
    cls <- class(object)[[1L]]
    problems <- c(
        unlist(map(
            .methodsSlotNames,
            .methodsSlotProblems,
            object = object,
            cls = cls
        )),
        .methodsParamOverlapProblem(object),
        .methodsParamSharedInputProblem(object)
    )
    if (length(problems) == 0L) TRUE else as.character(problems)
})

# Every method a Param names, across all three slots. The answer to
# token-level questions that do not care which path a method came from --
# "is mrmash requested", "are these methods available at all".
# @noRd
.methodsParamTokens <- function(param) {
    unique(unlist(map(
        .methodsSlotNames,
        function(s) names(slot(param, s)) %||% character(0)
    )))
}

# The same Param with some methods removed from every slot -- the joint
# phase handles mrmash itself, so the per-method loop must not also run it.
# Returns a Param rather than a flat list so the path structure survives
# into the per-component recursion.
# @noRd
.methodsParamDrop <- function(param, drop) {
    args <- set_names(
        map(.methodsSlotNames, .methodsSlotWithout, param = param, drop = drop),
        .methodsSlotNames
    )
    exec(match.fun(class(param)[[1L]]), !!!args)
}

# @noRd
.methodsSlotWithout <- function(slotName, param, drop) {
    v <- slot(param, slotName)
    if (is.null(v)) {
        return(NULL)
    }
    kept <- v[!is_in(names(v), drop)]
    if (length(kept) == 0L) NULL else kept
}

# The methods to run on one input class, and each one's arguments.
#
# `methods` applies to every input class -- it is the slot that states no
# path -- so it is always included; the path slot for this class is added to
# it. The other path slot is simply not this run's business.
# @noRd
.methodsParamResolve <- function(param, inputKind) {
    entries <- c(
        slot(param, "methods") %||% list(),
        slot(param, .methodsPathSlotFor(inputKind)) %||% list()
    )
    list(
        tokens = names(entries),
        methodArgs = map(entries, as.list)
    )
}

# Which path slot serves an input class. GwasSumStats is a
# summary-statistics class, so it reads the same slot as QtlSumStats.
# @noRd
.methodsPathSlotFor <- function(inputKind) {
    if (identical(inputKind, "QtlDataset")) {
        "qtlDatasetMethods"
    } else {
        "qtlSumStatsMethods"
    }
}

# The same, for a MultiStudyQtlDataset -- which may carry individual-level
# studies AND summary statistics, so one input class does not describe it.
#
# Two extra list shapes, named for the data they target and not mixable:
#
#   selection  list("susie", qtlDataset = "fsusie", qtlSumStats = "ser")
#              bare strings run on whatever the dataset holds; a path key
#              takes a character vector of methods for that half only.
#   override   list(susie = list(L = 20),
#                   qtlDataset = list(susie = list(...)),
#                   qtlSumStats = list(susie = list(...)))
#              a method named at the top carries arguments SHARED by both
#              halves, distributed into both path slots so each is checked
#              against the engine that half actually runs; a path key takes
#              a list of per-method argument lists for that half only.
#
# Mixing them is refused rather than guessed at: `list("susie", lasso =
# list(...))` could mean either and nothing says which.
# @noRd
.methodsParamForMulti <- function(value, cls, label, hasSumStats) {
    param <- if (is(value, cls)) {
        value
    } else if (is(value, "MethodsSelectionParam")) {
        abort(glue(
            "{label}: `methods` is a {class(value)[[1L]]}, which configures ",
            "a different pipeline. Use {cls}()."
        ))
    } else if (is.character(value)) {
        exec(match.fun(cls), methods = value)
    } else if (is.list(value)) {
        .methodsMultiFromList(value, cls, label)
    } else {
        abort(glue(
            "{label}: `methods` must be a character vector, a named list of ",
            "per-method options, or a {cls}() record; got ",
            "{class(value)[[1L]]}"
        ))
    }
    .methodsAssertSumStats(param, hasSumStats, label)
    param
}

# @noRd
.methodsMultiFromList <- function(value, cls, label) {
    nms <- names(value) %||% rep("", length(value))
    pathAt <- is_in(nms, c("qtlDataset", "qtlSumStats"))
    rest <- value[!pathAt]
    restNames <- nms[!pathAt]
    .methodsAssertOneShape(rest, restNames, label)
    shared <- if (length(rest) > 0L && all(nzchar(restNames))) {
        rest
    } else {
        list()
    }
    selected <- if (length(rest) > 0L && !all(nzchar(restNames))) {
        unlist(rest, use.names = FALSE)
    } else {
        character(0)
    }
    exec(
        match.fun(cls),
        methods = if (length(selected) > 0L) selected else NULL,
        qtlDatasetMethods = .methodsMultiSlot(
            value[nms == "qtlDataset"],
            shared,
            label,
            "qtlDataset"
        ),
        qtlSumStatsMethods = .methodsMultiSlot(
            value[nms == "qtlSumStats"],
            shared,
            label,
            "qtlSumStats"
        )
    )
}

# A path key's own entries, with the shared ones merged in beneath them: a
# per-path value wins over the shared value of the same name.
# @noRd
.methodsMultiSlot <- function(entry, shared, label, key) {
    own <- if (length(entry) == 0L) {
        list()
    } else {
        .methodsMultiPathEntries(entry[[1L]], label, key)
    }
    if (length(shared) == 0L && length(own) == 0L) {
        return(NULL)
    }
    merged <- shared
    for (nm in names(own)) {
        merged[[nm]] <- if (is.list(own[[nm]]) && is.list(shared[[nm]])) {
            list_modify(shared[[nm]], !!!own[[nm]])
        } else {
            own[[nm]]
        }
    }
    merged
}

# A path key takes either a character vector of methods or a named list of
# per-method arguments; a character vector becomes empty entries, which is
# what "run it with its defaults" looks like everywhere else.
# @noRd
.methodsMultiPathEntries <- function(v, label, key) {
    if (is.character(v)) {
        return(set_names(rep(list(list()), length(v)), v))
    }
    if (is.list(v)) {
        return(v)
    }
    abort(glue(
        "{label}: `{key}` must be a character vector of methods or a named ",
        "list of per-method options; got {class(v)[[1L]]}"
    ))
}

# The two list shapes are not mixable.
# @noRd
.methodsAssertOneShape <- function(rest, restNames, label) {
    if (length(rest) == 0L) {
        return(invisible(NULL))
    }
    named <- nzchar(restNames)
    if (all(named) || !any(named)) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: `methods` mixes two shapes -- bare method names, which ",
        "select, and named entries, which configure. Use one: a character ",
        "vector to select, or a named list to configure."
    ))
}

# Summary-statistics options are an instruction that cannot be carried out
# on a dataset holding none.
# @noRd
.methodsAssertSumStats <- function(param, hasSumStats, label) {
    if (hasSumStats) {
        return(invisible(NULL))
    }
    entries <- slot(param, "qtlSumStatsMethods")
    if (length(entries) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: `qtlSumStats` names ",
        "{str_flatten(names(entries), ', ')}, but this ",
        "MultiStudyQtlDataset carries no summary statistics. Remove them, ",
        "or supply a dataset whose sumStats slot is set."
    ))
}

# Whatever the caller gave a pipeline as `methods`, as a Param of this
# family.
#
# A pipeline knows its input class by dispatch, which is what makes a plain
# named list routable here but not in the constructor: the overrides go
# straight into the slot for this run's path, so they are checked against
# the engine that will actually receive them. That is why
# `methods = list(susie = list(L = 20))` stays as short at a single-type
# pipeline as it ever was.
# @noRd
.methodsParamFor <- function(value, inputKind, cls, label) {
    if (is(value, cls)) {
        return(value)
    }
    if (is(value, "MethodsSelectionParam")) {
        abort(glue(
            "{label}: `methods` is a {class(value)[[1L]]}, which configures ",
            "a different pipeline. Use {cls}()."
        ))
    }
    ctor <- match.fun(cls)
    if (is.character(value)) {
        return(ctor(methods = value))
    }
    args <- list()
    args[[.methodsPathSlotFor(inputKind)]] <- value
    exec(ctor, !!!args)
}

# One slot's value as it is stored: a named list of records.
#
# A bare character vector selects methods with their defaults, which is the
# same thing `methods = c("susie", "lasso")` means at a pipeline. A plain
# list entry is spliced into that method's own constructor, so it gets the
# same checking as writing the constructor out -- but only in a path slot,
# where the path is known. Under `methods` the path is not known until the
# pipeline dispatches, so an entry there must either be an already-built
# record (which carries its own input type) or list(); the pipelines route
# plain lists themselves, where the input class is settled.
# @noRd
.methodsNormalizeSlot <- function(value, slotName, cls) {
    if (is.null(value) || length(value) == 0L) {
        return(NULL)
    }
    if (is.character(value)) {
        return(set_names(rep(list(list()), length(value)), value))
    }
    if (!is.list(value)) {
        abort(glue(
            "{cls}: `{slotName}` must be a named list of per-method ",
            "options, or a character vector of method names; got ",
            "{class(value)[[1L]]}"
        ))
    }
    if (
        length(value) > 0L &&
            (is.null(names(value)) || any(!nzchar(names(value))))
    ) {
        abort(glue(
            "{cls}: every `{slotName}` entry must be named for its method"
        ))
    }
    # Detected with a map but reported OUTSIDE one: aborting from inside any
    # purrr call wraps the message in "In index: 1. With name: x." framing,
    # which is the caller's mistake dressed as an internal one. Same reason
    # .nestedAssertEngines() checks up front.
    imap(value, .methodsNormalizeEntry, slotName = slotName, cls = cls)
}

# An entry is stored as the caller wrote it. A plain list is NOT coerced
# into the engine's Options constructor: a per-method bundle is legitimately
# mixed -- the pecotmr wrapper's own formals (`dataDrivenPriorMatrices`,
# `fitRetention`) alongside the engine's arguments -- and .splitMethodArgs
# separates them at call time. Coercing would reject the wrapper half.
#
# Its NAMES are still checked here when the slot names a path, against the
# whole chain that path runs; see .methodsEntryNameProblem.
# @noRd
.methodsNormalizeEntry <- function(value, key, slotName, cls) {
    value
}

# Every name one method accepts on one input path: the wrapper chain's
# formals plus the engine's arguments. NULL when that cannot be determined,
# which is the existing convention for "do not reject anything here".
# @noRd
.methodsPathAccepted <- function(cls, token, inputKind) {
    if (identical(cls, "FineMappingMethodsParam")) {
        callee <- .fmMethodCalleeFor(token, inputKind)
        if (is.null(callee)) {
            return(NULL)
        }
        accepted <- .engineAcceptedNames(callee)
        if (is.null(accepted)) {
            return(NULL)
        }
        return(unique(c(accepted, .fmSeededArgNames())))
    }
    .twasChainAccepted(token, inputKind)
}

# One slot's own problems: shape, vocabulary, element class, and -- for the
# two path slots -- that each record is actually usable on that path.
#
# The first two stages stop on the first problem because the stages are
# ordered: a slot whose names are malformed has nothing to compare against
# the method vocabulary, and a name that is not a method of this pipeline
# cannot be asked whether it is usable on a given input path.
# @noRd
.methodsSlotProblems <- function(slotName, object, cls) {
    entries <- slot(object, slotName)
    if (is.null(entries) || length(entries) == 0L) {
        return(character(0))
    }
    nms <- names(entries)
    shape <- .methodsSlotShapeProblem(slotName, nms)
    if (length(shape) > 0L) {
        return(shape)
    }
    inputKind <- .methodsSlotInput(slotName)
    vocabulary <- .methodsSlotVocabularyProblem(
        slotName,
        nms,
        cls,
        inputKind
    )
    if (length(vocabulary) > 0L) {
        return(vocabulary)
    }
    .methodsSlotEntryProblems(nms, entries, slotName, cls, inputKind)
}

# Every entry named, exactly once.
# @noRd
.methodsSlotShapeProblem <- function(slotName, nms) {
    if (is.null(nms) || any(!nzchar(nms))) {
        return(glue("`{slotName}`: every entry must be named for its method"))
    }
    dup <- unique(nms[duplicated(nms)])
    if (length(dup) == 0L) {
        return(character(0))
    }
    glue(
        "`{slotName}`: {str_flatten(dup, ', ')} given more than once; ",
        "one set of options per method"
    )
}

# Every name a method this pipeline runs -- on `inputKind` for a path slot,
# anywhere for the shared one.
# @noRd
.methodsSlotVocabularyProblem <- function(slotName, nms, cls, inputKind) {
    allowed <- if (is.null(inputKind)) {
        .methodsParamAnyMethod(cls)
    } else {
        .methodsParamAvailable(cls, inputKind)
    }
    unknown <- setdiff(nms, allowed)
    if (length(unknown) == 0L) {
        return(character(0))
    }
    glue(
        "`{slotName}`: {str_flatten(unknown, ', ')} ",
        "{if (length(unknown) == 1L) 'is not a method' else ",
        "'are not methods'} this pipeline runs",
        "{if (is.null(inputKind)) '' else glue(' on {inputKind}')}. ",
        "Available: {str_flatten(allowed, ', ')}."
    )
}

# Each entry's own problems. Unlike the two stages above these are collected
# rather than short-circuited: one bad entry says nothing about the others.
# @noRd
.methodsSlotEntryProblems <- function(
    nms,
    entries,
    slotName,
    cls,
    inputKind
) {
    c(
        unlist(map(
            nms,
            .methodsEntryClassProblem,
            entries = entries,
            slotName = slotName
        )),
        unlist(map(
            nms,
            .methodsEntryOptionlessProblem,
            entries = entries,
            slotName = slotName,
            cls = cls
        )),
        if (is.null(inputKind)) {
            NULL
        } else {
            .methodsSlotPathProblems(nms, entries, slotName, inputKind, cls)
        }
    )
}

# The checks that only mean anything on a path slot: the record has to be
# usable on `inputKind`, name arguments that method accepts, and declare an
# engine that path provides.
# @noRd
.methodsSlotPathProblems <- function(
    nms,
    entries,
    slotName,
    inputKind,
    cls
) {
    c(
        unlist(map(
            nms,
            .methodsEntryPathProblem,
            entries = entries,
            slotName = slotName,
            inputKind = inputKind
        )),
        unlist(map(
            nms,
            .methodsEntryNameProblem,
            entries = entries,
            slotName = slotName,
            inputKind = inputKind,
            cls = cls
        )),
        unlist(map(
            nms,
            .methodsEntryEngineProblem,
            entries = entries,
            slotName = slotName,
            inputKind = inputKind,
            cls = cls
        ))
    )
}

# An entry is a constructor-built record or a plain list of that method's
# arguments; list() means "run it with its defaults". A plain list is
# allowed because the bundle is legitimately mixed -- the wrapper's own
# formals and the engine's arguments together -- so only its NAMES can be
# judged, which .methodsEntryNameProblem does wherever the path is known.
# @noRd
.methodsEntryClassProblem <- function(nm, entries, slotName) {
    v <- entries[[nm]]
    if (.isMethodOptions(v) || is.list(v)) {
        return(character(0))
    }
    glue(
        "`{slotName}${nm}` must be that method's Options record or a list ",
        "of its arguments; got {class(v)[[1L]]}"
    )
}

# A method with nothing to configure, given something to configure.
#
# The fit-derived methods extract weights from a fit a prior
# fineMappingPipeline() run produced, so they have no arguments of their own
# on either path -- which is why naming one with list() SELECTS it and is
# the only thing it can say. A non-empty entry is a misunderstanding worth
# naming, and .twasTokenNoArgsReason already says where the setting lives.
#
# Path-independent, so this applies to `methods` as well: the answer does
# not depend on which input class runs.
# @noRd
.methodsEntryOptionlessProblem <- function(nm, entries, slotName, cls) {
    v <- entries[[nm]]
    if (.isMethodOptions(v) || !is.list(v) || length(v) == 0L) {
        return(character(0))
    }
    reason <- .methodsNoOptionsReason(cls, nm)
    if (is.null(reason)) {
        return(character(0))
    }
    glue("`{slotName}${nm}`: {reason}")
}

# @noRd
.methodsNoOptionsReason <- function(cls, token) {
    if (identical(cls, "FineMappingMethodsParam")) {
        return(NULL)
    }
    .twasTokenNoArgsReason(token)
}

# The engine one method's options must have been built for on one path.
#
# Not the token: `SusieRssOptions()` records engine "susieRss" while the
# token is "susie". Derived by taking that token's OWN constructors and
# picking the one whose callee matches this path's -- which is also what
# tells `susie` from `susieInf`, since both run `susieR::susie` on
# individual data and only the per-token pair distinguishes them.
# @noRd
.methodsExpectedEngine <- function(cls, token, inputKind) {
    if (identical(cls, "FineMappingMethodsParam")) {
        return(.fmExpectedEngine(token, inputKind))
    }
    ctor <- .twasMethodCtorFor(token, inputKind)
    if (is.null(ctor)) {
        return(NULL)
    }
    .methodsRecordEngine(ctor)
}

# @noRd
.fmExpectedEngine <- function(token, inputKind) {
    ctors <- .fmTokenCtorsByPath()[[token]]
    want <- .fmMethodCalleeFor(token, inputKind)
    if (is.null(ctors) || is.null(want)) {
        return(NULL)
    }
    hit <- keep(ctors, .methodsCtorHasCallee, callee = want)
    if (length(hit) == 0L) {
        return(NULL)
    }
    .methodsRecordEngine(hit[[1L]])
}

# @noRd
.methodsCtorHasCallee <- function(ctor, callee) {
    rec <- tryCatch(ctor(), error = function(cnd) NULL)
    !is.null(rec) && identical(metadata(rec)$callee, callee)
}

# @noRd
.methodsRecordEngine <- function(ctor) {
    rec <- tryCatch(ctor(), error = function(cnd) NULL)
    if (is.null(rec)) {
        return(NULL)
    }
    metadata(rec)$engine
}

# An Options record keyed under one method must have been built by THAT
# method's constructor. The input-type check cannot see this: `lasso` and
# `bayesA` both run on individual data, so `GlmnetOptions()` and
# `QggOptions()` are indistinguishable by path and only the engine tells
# them apart.
# @noRd
.methodsEntryEngineProblem <- function(nm, entries, slotName, inputKind, cls) {
    v <- entries[[nm]]
    if (!.isMethodOptions(v)) {
        return(character(0))
    }
    want <- .methodsExpectedEngine(cls, nm, inputKind)
    got <- metadata(v)$engine
    if (is.null(want) || is.null(got) || identical(got, want)) {
        return(character(0))
    }
    glue(
        "`{slotName}${nm}` was built with the constructor for '{got}', but ",
        "{nm} on {inputKind} uses the '{want}' constructor"
    )
}

# A plain list's names, against everything that method accepts on this path
# -- the wrapper chain plus the engine. Only in a path slot: `methods`
# names no path, so the pipeline checks it once the input class is settled.
# @noRd
.methodsEntryNameProblem <- function(nm, entries, slotName, inputKind, cls) {
    v <- entries[[nm]]
    if (.isMethodOptions(v) || !is.list(v) || length(v) == 0L) {
        return(character(0))
    }
    accepted <- .methodsPathAccepted(cls, nm, inputKind)
    if (is.null(accepted)) {
        return(character(0))
    }
    unknown <- setdiff(names(v), accepted)
    if (length(unknown) == 0L) {
        return(character(0))
    }
    glue(
        "`{slotName}${nm}`: unknown argument(s) ",
        "{str_flatten(unknown, ', ')} for {nm} on {inputKind}."
    )
}

# A record carries the input classes it is for (see .engineInputType), so a
# summary-statistics record under `qtlDatasetMethods` is caught here rather
# than failing inside the engine.
# @noRd
.methodsEntryPathProblem <- function(nm, entries, slotName, inputKind) {
    it <- .optionsInputType(entries[[nm]])
    if (is.null(it) || is_in(inputKind, it)) {
        return(character(0))
    }
    glue(
        "`{slotName}${nm}`: these options are for ",
        "{str_flatten(sort(it), ' / ')}, but `{slotName}` is forwarded to ",
        "{inputKind}"
    )
}

# A method named in `methods` must not also be named in a path slot: the
# first says its path need not be stated, the second states one, and
# nothing decides which wins.
#
# Appearing in BOTH path slots is not a conflict -- it is the whole point of
# having two, and how a dual-type run gives one method different arguments
# per path.
# @noRd
.methodsParamOverlapProblem <- function(object) {
    pathless <- names(slot(object, "methods")) %||% character(0)
    paths <- unique(c(
        names(slot(object, "qtlDatasetMethods")) %||% character(0),
        names(slot(object, "qtlSumStatsMethods")) %||% character(0)
    ))
    clash <- intersect(pathless, paths)
    if (length(clash) == 0L) {
        return(character(0))
    }
    glue(
        "{str_flatten(clash, ', ')} ",
        "{if (length(clash) == 1L) 'is' else 'are'} named in `methods` and ",
        "also under an input path. Put each method in one or the other: ",
        "`methods` for options that need no path, the path slots for ones ",
        "that do."
    )
}

# @noRd
.methodsInputTypeLabel <- function(key, types) {
    sprintf("%s -> %s", key, str_flatten(sort(types[[key]]), "/"))
}

# `methods` states no path, so everything in it has to be usable on one
# single input class -- otherwise the slot describes a run that cannot
# happen.
# @noRd
.methodsParamSharedInputProblem <- function(object) {
    entries <- slot(object, "methods")
    if (is.null(entries) || length(entries) == 0L) {
        return(character(0))
    }
    types <- compact(map(entries, .optionsInputType))
    if (length(types) == 0L) {
        return(character(0))
    }
    if (length(reduce(types, intersect)) > 0L) {
        return(character(0))
    }
    detail <- str_flatten(
        map_chr(names(types), .methodsInputTypeLabel, types = types),
        collapse = "; "
    )
    glue(
        "`methods`: these options cannot all be used with one input class ",
        "({detail}). Use `qtlDatasetMethods` / `qtlSumStatsMethods` to say ",
        "which goes where."
    )
}

# One slot's value as a line of text.
#
# A slot holding a LIST -- the per-method entries of a method-selection
# Param, each an Options record or a list of arguments -- is shown by naming
# its entries. format() cannot render those: it falls through to
# as.character() on an S4 object and errors, which is how printing a
# populated FineMappingMethodsParam used to fail.
# @noRd
.paramSlotText <- function(value) {
    # A slot may hold a RECORD rather than a value: an Options bag bound for
    # an engine, or a nested Param. format() prints the S4 class header for
    # those -- "<S4 class 'MethodOptions' with 4 slots>" -- which tells a
    # reader nothing. Name what it is instead.
    if (.isMethodOptions(value)) {
        return(.paramOptionsText(value))
    }
    if (.isMethodParam(value)) {
        return(.paramNestedText(value))
    }
    if (is.list(value) && !is.null(names(value))) {
        return(str_flatten(names(value), ", "))
    }
    if (is.list(value)) {
        return(glue(
            "{length(value)} entr{if (length(value) == 1L) 'y' else 'ies'}"
        ))
    }
    str_flatten(format(value), ", ")
}

# A nested Param by what it CONTAINS, not by its slot names. The
# aggregators hold their per-method entries in one list slot, so reporting
# slot names would say "components" where the reader wants "pca".
# @noRd
.paramNestedText <- function(value) {
    set <- names(value)
    if (length(set) == 0L) {
        return(glue("<{class(value)[[1L]]}: nothing set>"))
    }
    inner <- unlist(map(set, .paramNestedEntryNames, value = value))
    if (length(inner) > 0L) {
        return(str_flatten(unique(inner), ", "))
    }
    str_flatten(set, ", ")
}

# @noRd
.paramNestedEntryNames <- function(nm, value) {
    v <- slot(value, nm)
    if (is.list(v)) names(v) else character(0)
}

# An Options record by the engine it is for, with its settings when it
# carries any.
# @noRd
.paramOptionsText <- function(value) {
    engine <- metadata(value)$engine %||% "options"
    if (length(value) == 0L) {
        return(glue("<{engine}: defaults>"))
    }
    glue("<{engine}: {str_flatten(names(value), ', ')}>")
}

# TRUE for a record one of the *Param() constructors produced. The Param
# counterpart of .isMethodOptions(), and the two are disjoint: a Param is not
# a MethodOptions and never was.
# @noRd
.isMethodParam <- function(x) {
    is(x, "MethodParam")
}

# Demand a record of one specific Param class, naming the constructor that
# builds it. Stricter than the MethodOptions equivalent could be: every
# *Param() has its own class, so handing `ColocPriorParam()` where
# `CredibleSetParam()` was wanted is caught here rather than surfacing as a
# missing field somewhere downstream.
# @noRd
.assertMethodParam <- function(x, ctor, argName, class = NULL) {
    class <- class %||% .paramClassOf(ctor)
    if (is(x, class)) {
        return(invisible(x))
    }
    if (is.list(x) && length(x) == 0L) {
        # An empty list carries no settings, so there is nothing a
        # constructor would have added. Accepting it keeps `list()` working
        # as "no settings", as it does for MethodOptions.
        return(invisible(x))
    }
    got <- if (.isMethodParam(x) || .isMethodOptions(x)) {
        glue("a {class(x)[[1L]]}")
    } else {
        "a bare list"
    }
    abort(glue(
        "`{argName}` must be built with {ctor}(), not {got}. ",
        "A setting pecotmr reads is type-checked when the record is built, ",
        "which a bare list cannot be."
    ))
}

# The engine a Param denotes: "DentistParam" -> "dentist". The inverse of
# .paramClassOf(), for the choices that take either an engine name or that
# engine's record (see .engineOf).
# @noRd
.paramEngineOf <- function(x) {
    base <- str_remove(class(x)[[1L]], "Param$")
    str_c(str_to_lower(str_sub(base, 1L, 1L)), str_sub(base, 2L))
}

# The class a constructor builds: "CredibleSetParam" -> "CredibleSetParam".
# Derived rather than tabulated so the two can never drift apart; every
# *Param() constructor follows the convention without exception.
# @noRd
.paramClassOf <- function(ctor) {
    str_c(str_to_upper(str_sub(ctor, 1L, 1L)), str_sub(ctor, 2L))
}
