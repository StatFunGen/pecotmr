# =============================================================================
# Settings records to and from YAML
# -----------------------------------------------------------------------------
# A config file holds one mapping per record, keyed by the constructor that
# builds it -- `CredibleSetParam:`, `SusieOptions:` -- whose own keys are that
# constructor's arguments. Loading a record is calling its constructor on
# those arguments, so a YAML record gets exactly the checking a record built
# in R gets: a Param's types, an Options record's names against the engine.
#
# Scalars are read under the YAML 1.2 core schema rather than the 1.1 one
# the yaml package implements. Under 1.1 the KEY `n` -- a real susie_rss
# argument -- reads as FALSE and `y` as TRUE, `1e-5` reads as a string
# because a 1.1 float needs a dot, and an integer beyond 32 bits reads as NA.
# =============================================================================

#' @include MethodOptions.R
#' @include MethodParam.R
#' @importFrom methods extends getSlots
#' @importFrom purrr every
#' @importFrom readr read_file write_file
#' @importFrom rlang as_string has_name is_bare_atomic is_scalar_atomic
#' @importFrom rlang caller_env is_string
NULL

# A name is a constructor when it is an exported *Options() / *Param()
# function. The capital excludes qcOptions(), which builds nothing.
# @noRd
.configIsCtor <- function(name) {
    is_string(name) &&
        str_detect(name, "^[A-Z][A-Za-z0-9]*(Options|Param)$") &&
        is_in(name, getNamespaceExports("pecotmr"))
}

# @noRd
.configAssertCtor <- function(name, call = caller_env()) {
    if (.configIsCtor(name)) {
        return(invisible(NULL))
    }
    abort(
        glue("`{name}` is not a pecotmr *Options() or *Param() constructor."),
        call = call
    )
}

# The constructor an argument's default calls, when it calls one:
# `fit = UdFitOptions()` says what a mapping given for `fit` builds.
# @noRd
.configDefaultCtor <- function(default) {
    if (!is.call(default) || !is.symbol(default[[1L]])) {
        return(NULL)
    }
    head <- as_string(default[[1L]])
    if (.configIsCtor(head)) head else NULL
}

# --- reading -----------------------------------------------------------------

# YAML 1.2 core-schema scalars. Only true / false are booleans; anything
# else 1.1 would have taken for one (yes, no, on, off, y, n) stays a string.
# @noRd
.configYamlTrue <- function(x) {
    if (identical(str_to_lower(x), "true")) TRUE else x
}

# @noRd
.configYamlFalse <- function(x) {
    if (identical(str_to_lower(x), "false")) FALSE else x
}

# An exponent float without a dot (`1e-5`) is a number under 1.2, and a hex
# float is how writeConfigToYaml() writes a double no decimal reads back. A
# quoted "1e-5" arrives here too, so such a string cannot be written at all.
# @noRd
.configYamlString <- function(x) {
    float <- "^[-+]?(\\.[0-9]+|[0-9]+(\\.[0-9]*)?)[eE][-+]?[0-9]+$"
    hex <- "^[-+]?0[xX][0-9a-fA-F]+(\\.[0-9a-fA-F]*)?[pP][-+]?[0-9]+$"
    if (str_detect(x, float) || str_detect(x, hex)) as.numeric(x) else x
}

# An integer too wide for R's 32 bits stays a double rather than becoming NA.
# @noRd
.configYamlInteger <- function(x) {
    value <- as.numeric(x)
    if (abs(value) <= .Machine$integer.max) as.integer(value) else value
}

# R's own reader, since the yaml package's turns a subnormal double -- the
# smallPFloor default is one -- into NA. .inf and .nan carry tags of their
# own and never reach here.
# @noRd
.configYamlFloat <- function(x) {
    as.numeric(x)
}

# @noRd
.configReadHandlers <- function() {
    list(
        "bool#yes" = .configYamlTrue,
        "bool#no" = .configYamlFalse,
        str = .configYamlString,
        int = .configYamlInteger,
        "float#fix" = .configYamlFloat,
        "float#exp" = .configYamlFloat
    )
}

# A sequence of scalars is a vector. The yaml package already makes one of a
# sequence whose scalars share a type, but leaves `[1, 2.5]` -- an integer and
# a double -- as a list.
# @noRd
.configSimplify <- function(x) {
    if (!is.list(x)) {
        return(x)
    }
    if (is.null(names(x)) && .configIsScalarRun(x)) {
        return(unlist(x))
    }
    map(x, .configSimplify)
}

# @noRd
.configIsScalarRun <- function(x) {
    if (length(x) == 0L || !every(x, is_scalar_atomic)) {
        return(FALSE)
    }
    types <- unique(map_chr(x, typeof))
    length(types) == 1L || all(is_in(types, c("integer", "double")))
}

# @noRd
.configParse <- function(text) {
    .configSimplify(yaml.load(text, handlers = .configReadHandlers()))
}

# A mapping whose one key is a constructor builds that constructor: the form
# for a field whose default does not say which record it holds.
# @noRd
.configIsTagged <- function(value) {
    is.list(value) &&
        length(value) == 1L &&
        .configIsCtor(names(value)) &&
        (is.list(value[[1L]]) || is.null(value[[1L]]))
}

# An empty sequence takes the type of the slot it fills, so `keepSamples: []`
# is character(0) as written. Without a typed slot it is a list.
# @noRd
.configEmpty <- function(field, slotTypes) {
    if (!has_name(slotTypes, field)) {
        return(list())
    }
    base <- str_remove(slotTypes[[field]], "_OR_NULL$")
    if (!is_in(base, c("character", "numeric", "integer", "logical"))) {
        return(list())
    }
    vector(base, 0L)
}

# One argument's value: a nested record built, an empty sequence typed, and
# any list searched for tagged records -- a method list holds one per method.
# @noRd
.configValue <- function(value, field, defaults, slotTypes) {
    if (.configIsTagged(value)) {
        return(.configBuild(names(value), value[[1L]]))
    }
    if (!is.list(value)) {
        return(value)
    }
    # Before the empty-sequence case: a record with nothing set is written
    # `[]`, and still names the record its default builds.
    nested <- .configDefaultCtor(defaults[[field]])
    if (!is.null(nested) && (length(value) == 0L || is_named(value))) {
        return(.configBuild(nested, value))
    }
    if (length(value) == 0L && is.null(names(value))) {
        return(.configEmpty(field, slotTypes))
    }
    map(value, .configElement)
}

# @noRd
.configElement <- function(value) {
    if (.configIsTagged(value)) {
        return(.configBuild(names(value), value[[1L]]))
    }
    if (is.list(value)) map(value, .configElement) else value
}

# A record from its mapping. Nothing to set -- `SusieOptions:` alone, or an
# empty mapping -- is the constructor's defaults.
# @noRd
.configBuild <- function(name, values) {
    .configAssertCtor(name)
    settings <- values %||% list()
    if (!is.list(settings) || (length(settings) > 0L && !is_named(settings))) {
        abort(glue("`{name}` must be a mapping of settings to values."))
    }
    ctor <- get(name, envir = asNamespace("pecotmr"))
    slotTypes <- if (extends(name, "MethodParam")) getSlots(name) else NULL
    args <- imap(
        settings,
        .configValue,
        defaults = formals(ctor),
        slotTypes = slotTypes
    )
    exec(ctor, !!!args)
}

#' @title Load One Settings Record From A YAML File
#' @description Builds one \code{*Options()} or \code{*Param()} record from a
#'   YAML config file. The file's top-level keys are constructor names; under
#'   each are that constructor's arguments and their values. The record is
#'   built by calling its constructor, so it is checked exactly as one built
#'   in R is.
#' @param filename Path to the YAML file.
#' @param name The constructor whose record to build, such as
#'   \code{"CredibleSetParam"} or \code{"SusieOptions"}. Other records in
#'   the file are ignored.
#' @details A file holds at most one record per constructor:
#'   \preformatted{
#' CredibleSetParam:
#'   coverage: 0.95
#'   perCsColumns: alpha
#' GenotypeFilterParam:
#'   mafCutoff: 0.01
#'   keepSamples: [S1, S2, S3]
#' SusieOptions:
#'   max_iter: 200
#' }
#'   A \code{*Param()} record's keys are its settings; an \code{*Options()}
#'   record's are the engine's own argument names. A sequence gives a vector,
#'   \code{~} gives \code{NULL}, \code{[]} gives an empty vector of the
#'   setting's type, \code{.inf} and \code{.na} give \code{Inf} and
#'   \code{NA}.
#'
#'   Scalars follow the YAML 1.2 core schema: only \code{true} and
#'   \code{false} are logical, so \code{yes}, \code{no}, \code{on} and
#'   \code{off} stay strings, and \code{1e-5} is a number. A string that
#'   looks like an exponent number, even quoted, reads as that number.
#'
#'   A setting that holds another record is a nested mapping. Where the
#'   argument's default already names its constructor, as
#'   \code{GwasFineMappingParam}'s \code{credibleSetParam} does, the mapping
#'   holds that record's settings directly. Elsewhere it is wrapped in a
#'   mapping keyed by the constructor:
#'   \preformatted{
#' MashPriorParam:
#'   engine:
#'     CovUdrOptions:
#'       fit:
#'         maxiter: 50
#' }
#'   Settings that hold R objects rather than values have no YAML form:
#'   \code{BPPARAM}, which belongs to the machine and is set at call time;
#'   matrices such as \code{MashPriorParam}'s \code{priorCovariances};
#'   data frames such as \code{samplePartition}; and engine arguments that
#'   take a fitted model, a function or a per-variant vector.
#' @return The record \code{name} builds.
#' @seealso \code{\link{writeConfigToYaml}}
#' @examples
#' path <- tempfile(fileext = ".yaml")
#' writeLines(c("CredibleSetParam:", "  coverage: 0.9"), path)
#' loadConfigFromYaml(path, "CredibleSetParam")
#' @importFrom checkmate assertFileExists assertString
#' @importFrom yaml yaml.load
#' @export
loadConfigFromYaml <- function(filename, name) {
    assertFileExists(filename, access = "r")
    assertString(name)
    .configAssertCtor(name)
    config <- .configParse(read_file(filename))
    if (!is.list(config) || !is_named(config)) {
        abort(glue(
            "'{filename}' must be a YAML mapping keyed by constructor name."
        ))
    }
    if (!has_name(config, name)) {
        abort(glue(
            "'{filename}' has no `{name}` record. It has: ",
            "{str_flatten(names(config), ', ')}."
        ))
    }
    try_fetch(
        .configBuild(name, config[[name]]),
        error = function(cnd) {
            inner <- .configUnwrap(cnd)
            abort(
                c(
                    glue("Could not build `{name}` from '{filename}'."),
                    .configWhere(inner$path)
                ),
                parent = inner$cause
            )
        }
    )
}

# --- reporting ---------------------------------------------------------------

# The error under the layers purrr wraps around it, one for each map it
# passes through, and the path of elements those maps were on.
# @noRd
.configUnwrap <- function(cnd, path = character()) {
    if (!inherits(cnd, "purrr_error_indexed")) {
        return(list(cause = cnd, path = path))
    }
    .configUnwrap(cnd$parent, c(path, .configPathStep(cnd)))
}

# One map's element: `$name`, or `[[i]]` where the elements are unnamed.
# @noRd
.configPathStep <- function(cnd) {
    if (is_string(cnd$name) && nzchar(cnd$name)) {
        return(str_c("$", cnd$name))
    }
    str_c("[[", cnd$location, "]]")
}

# The setting a failure was reached through, as a message bullet.
# @noRd
.configWhere <- function(path) {
    if (length(path) == 0L) {
        return(NULL)
    }
    setting <- str_remove(str_flatten(path), "^[$]")
    c(i = glue("In setting `{setting}`."))
}

# --- writing -----------------------------------------------------------------

# @noRd
.configIsRecord <- function(x) {
    .isMethodParam(x) || .isMethodOptions(x)
}

# The constructor a record names. A Param's class is its constructor; an
# Options record carries its own in metadata(), NULL for one built some
# other way.
# @noRd
.configRecordCtor <- function(record) {
    if (.isMethodParam(record)) {
        return(class(record)[[1L]])
    }
    metadata(record)$constructor
}

# Only ever the cause of a writer's error, which names the call, so it names
# none itself.
# @noRd
.configCtorName <- function(record) {
    name <- .configRecordCtor(record)
    if (is.null(name)) {
        abort(
            c(
                str_c(
                    "This MethodOptions record does not say which ",
                    "constructor built it."
                ),
                i = "Build it again with its *Options() constructor."
            ),
            call = NULL
        )
    }
    .configAssertCtor(name, call = NULL)
    name
}

# A record's settings as its constructor takes them: an Options record's
# arguments, or a Param's slots, since every Param's constructor takes
# exactly its slots. NULL slots are included: an unset slot is not
# necessarily the default.
# @noRd
.configRecordArgs <- function(record) {
    if (.isMethodOptions(record)) {
        return(as.list(record))
    }
    slots <- slotNames(record)
    set_names(map(slots, slot, object = record), slots)
}

# A record as a mapping. BPPARAM is left out: it names a backend on this
# machine, not a setting of the analysis.
# @noRd
.configDump <- function(record) {
    args <- discard(.configRecordArgs(record), is, "BiocParallelParam")
    defaults <- formals(.configCtorName(record))
    # A named empty list, so a record with nothing set is written `{}`.
    set_names(
        imap(args, .configDumpValue, defaults = defaults),
        names(args) %||% character()
    )
}

# A nested record is written bare where the argument's default names its
# constructor, tagged with the constructor everywhere else.
# @noRd
.configDumpValue <- function(value, field, defaults) {
    if (!.configIsRecord(value)) {
        return(.configDumpElement(value, field))
    }
    name <- .configCtorName(value)
    if (identical(.configDefaultCtor(defaults[[field]]), name)) {
        return(.configDump(value))
    }
    set_names(list(.configDump(value)), name)
}

# @noRd
.configDumpElement <- function(value, field) {
    if (.configIsRecord(value)) {
        return(set_names(list(.configDump(value)), .configCtorName(value)))
    }
    if (is.list(value) && !is.object(value)) {
        return(map(value, .configDumpElement, field = field))
    }
    if (.configIsPlainVector(value)) {
        return(value)
    }
    abort(
        glue(
            "`{field}` holds a {.configKindText(value)}, which has no YAML ",
            "form."
        ),
        call = NULL
    )
}

# What a refused value is: its class, with "named" for a named vector, whose
# class alone would read as writable.
# @noRd
.configKindText <- function(value) {
    kind <- class(value)[[1L]]
    if (is_bare_atomic(value) && !is.null(names(value))) {
        return(str_c("named ", kind))
    }
    kind
}

# @noRd
.configIsPlainVector <- function(value) {
    is.null(value) ||
        (is_bare_atomic(value) && is.null(names(value)) && is.null(dim(value)))
}

# A double written so it reads back the same number: the fewest decimal
# digits that do, else a C99 hex float. R's decimal reader is off by an ulp
# for some 17-digit strings, and a hex float it reads exactly. A decimal
# always carries a point, so a whole number reads back a double.
# @noRd
.configNumberText <- function(x) {
    if (is.nan(x)) {
        return(".nan")
    }
    if (is.na(x)) {
        return(".na.real")
    }
    if (is.infinite(x)) {
        return(if (x > 0) ".inf" else "-.inf")
    }
    decimal <- map_chr(15:17, .configDecimalText, x = x)
    exact <- decimal[as.numeric(decimal) == x]
    if (length(exact) == 0L) {
        return(sprintf("%a", x))
    }
    text <- exact[[1L]]
    if (str_detect(text, "[.]")) {
        return(text)
    }
    if (str_detect(text, "e")) {
        return(str_replace(text, "e", ".0e"))
    }
    str_c(text, ".0")
}

# @noRd
.configDecimalText <- function(digits, x) {
    sprintf(str_c("%.", digits, "g"), x)
}

# @noRd
.configWriteNumber <- function(x) {
    structure(map_chr(x, .configNumberText), class = "verbatim")
}

# true / false rather than the yes / no as.yaml writes by default.
# @noRd
.configWriteLogical <- function(x) {
    text <- if_else(x, "true", "false", missing = ".na")
    structure(text, class = "verbatim")
}

# The text must rebuild the record it was written from. The original's
# BPPARAM is carried onto the rebuilt record first, since it was left out of
# the file on purpose. A failure is only ever the cause of the writer's
# error, which names the call.
# @noRd
.configAssertRoundTrip <- function(record, name, parsed) {
    rebuilt <- .configCarryBackends(
        .configBuild(name, parsed[[name]]),
        record
    )
    if (identical(rebuilt, record)) {
        return(invisible(NULL))
    }
    fields <- .configFieldNames(record)
    kept <- map_lgl(fields, .configSameField, a = rebuilt, b = record)
    if (all(kept)) {
        abort(
            str_c(
                "Only its settings are written, and its metadata() would ",
                "come back as its constructor sets it."
            ),
            call = NULL
        )
    }
    abort(
        glue(
            "These settings would change: ",
            "{str_flatten(fields[!kept], ', ')}."
        ),
        call = NULL
    )
}

# @noRd
.configCarryBackends <- function(rebuilt, record) {
    if (!.isMethodParam(record)) {
        return(rebuilt)
    }
    slots <- slotNames(record)
    backends <- slots[map_lgl(slots, .configHoldsBackend, record = record)]
    reduce(backends, .configCopySlot, from = record, .init = rebuilt)
}

# @noRd
.configHoldsBackend <- function(name, record) {
    is(slot(record, name), "BiocParallelParam")
}

# @noRd
.configCopySlot <- function(to, name, from) {
    slot(to, name) <- slot(from, name)
    to
}

# @noRd
.configFieldNames <- function(record) {
    if (.isMethodParam(record)) slotNames(record) else names(record)
}

# @noRd
.configSameField <- function(field, a, b) {
    identical(a[[field]], b[[field]])
}

# Evaluates `step`, a map over `records`, so that a failure in it is raised
# as the failure of the record it happened on, without purrr's framing.
# @noRd
.configWriteStep <- function(
    step,
    records,
    filename,
    reason = NULL,
    call = caller_env()
) {
    try_fetch(
        step,
        purrr_error_indexed = function(cnd) {
            .configWriteFailed(cnd, records, filename, reason, call)
        }
    )
}

# The outermost layer purrr wrapped is the map over the records, so it says
# which record failed; the layers inside it say which setting.
# @noRd
.configWriteFailed <- function(cnd, records, filename, reason, call) {
    inner <- .configUnwrap(cnd$parent)
    label <- .configRecordLabel(records[[cnd$location]], cnd$location)
    abort(
        c(
            glue("Could not write {label} to '{filename}'."),
            i = reason,
            .configWhere(inner$path)
        ),
        parent = inner$cause,
        call = call
    )
}

# How a failure names a record: by its constructor, or by its position among
# the records written when it names none.
# @noRd
.configRecordLabel <- function(record, i) {
    name <- .configRecordCtor(record)
    if (.configIsCtor(name)) str_c("`", name, "`") else str_c("record ", i)
}

#' @title Write Settings Records To A YAML File
#' @description Writes \code{*Options()} and \code{*Param()} records to a
#'   YAML config file that \code{\link{loadConfigFromYaml}} reads back, one
#'   top-level mapping per record keyed by its constructor.
#' @param obj One record, or a list of records built by different
#'   constructors.
#' @param filename Path of the file to write. An existing file is replaced.
#' @details Every setting a record holds is written, defaults included, so
#'   the file pins the analysis against later changes to a default. Each
#'   record is read back from the text before the file is written, and one
#'   that would not come back identical is refused. The error names the
#'   settings that would change or, when the text does not build the record
#'   at all, gives its constructor's error. A string that reads as a number,
#'   such as \code{"1e5"}, is the usual cause. \code{BPPARAM} is never
#'   written: it names a backend on the machine running the analysis.
#'
#'   Settings with no YAML form -- a matrix, a data frame, a function, a
#'   named vector -- are refused. See \code{\link{loadConfigFromYaml}} for
#'   the format.
#' @return \code{filename}, invisibly.
#' @seealso \code{\link{loadConfigFromYaml}}
#' @examples
#' path <- tempfile(fileext = ".yaml")
#' writeConfigToYaml(
#'     list(CredibleSetParam(coverage = 0.9), SusieOptions(max_iter = 200)),
#'     path
#' )
#' loadConfigFromYaml(path, "SusieOptions")
#' @importFrom checkmate assertString
#' @importFrom yaml as.yaml
#' @export
writeConfigToYaml <- function(obj, filename) {
    assertString(filename)
    records <- if (.configIsRecord(obj)) list(obj) else obj
    if (!is.list(records) || !every(records, .configIsRecord)) {
        abort("`obj` must be a *Options() / *Param() record or a list of them.")
    }
    ctorNames <- .configWriteStep(
        map_chr(records, .configCtorName),
        records,
        filename
    )
    dup <- unique(ctorNames[duplicated(ctorNames)])
    if (length(dup) > 0L) {
        abort(glue(
            "A config file holds one record per constructor; ",
            "given more than one of: {str_flatten(dup, ', ')}."
        ))
    }
    dumped <- .configWriteStep(map(records, .configDump), records, filename)
    text <- as.yaml(
        set_names(dumped, ctorNames),
        handlers = list(
            numeric = .configWriteNumber,
            logical = .configWriteLogical
        )
    )
    parsed <- .configParse(text)
    .configWriteStep(
        walk2(records, ctorNames, .configAssertRoundTrip, parsed = parsed),
        records,
        filename,
        reason = "It would not read back as written."
    )
    write_file(text, filename)
    invisible(filename)
}
