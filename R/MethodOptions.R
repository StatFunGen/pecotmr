# =============================================================================
# MethodOptions S4 class
# -----------------------------------------------------------------------------
# One fitting engine's arguments, as built by a constructor.
#
# Every tunable pecotmr forwards to an external engine is built by a
# constructor: pecotmr's own defaults are the constructor's named formals, and
# anything else the engine accepts goes through `...`. The constructor checks
# the extras against the engine's LIVE formals, so the check can never drift
# from the upstream signature the way a transcribed list of names does.
#
# It IS a SimpleList, which is what makes "constructor required" enforceable.
# A plain list carrying an S3 class loses that class on `[`, `c()` and
# `compact()` -- all of which this package performs on these very records
# (.fmNarrowAfterJoint, .fmQssJointPhase, .ctwasInvoke) -- so the marker would
# evaporate on first use. SimpleList's subsetting is endomorphic and its
# metadata() survives, and it still splices with `!!!`, so no existing
# exec(fn, !!!args) call site changes.
#
# Why constructors exist at all: a bare list cannot be checked.
# `summaryStatsQc`'s imputation bundle advertised an `rcond` in its own
# default, but QC always routes RAISS through the SVD solver, which takes
# `svdTol` -- so `rcond` was read by nobody and silently ignored for as long
# as it existed. A constructor turns that into an error at the call site.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Arguments For One Fitting Engine
#' @description A validated bundle of arguments destined for one external
#'   fitting engine, produced by an engine constructor such as
#'   \code{SusieOptions()} or \code{MrmashOptions()} rather than written as a
#'   bare list. Names are checked against the engine's own formals at
#'   construction, so a misspelled option fails immediately instead of being
#'   silently dropped.
#'
#'   It \strong{is} a \code{\link[S4Vectors]{SimpleList}}: \code{$},
#'   \code{[[}, \code{names()} and \code{length()} behave as on a list, it
#'   splices into a call with \code{!!!}, and subsetting returns a
#'   \code{MethodOptions} rather than degrading to a plain list. The engine name
#'   and the callee it forwards to travel in \code{metadata()}.
#' @section Whether argument names are checked:
#'   A constructor checks the names you give it against the \emph{live}
#'   formals of the function it forwards to, so the check can never drift from
#'   the version of that package you have installed. Printing a
#'   \code{MethodOptions} says which applies to it.
#'
#'   Names cannot be checked in two situations, and in both the arguments are
#'   forwarded unchecked rather than rejected:
#'   \itemize{
#'     \item \strong{The engine takes \code{...}}, so any name is legal as
#'       far as it is concerned. As of writing this applies to
#'       \code{mvsusieR::mvsusie_rss}, \code{mvsusieR::create_mixture_prior},
#'       \code{ncvreg::ncvreg}, \code{ncvreg::cv.ncvreg},
#'       \code{RcppDPR::fit_model}, \code{mashr::cov_ed},
#'       \code{ctwas::finemap_regions} and
#'       \code{ctwas::postprocess_region_merging}. A misspelling here reaches
#'       the engine, which may ignore it silently.
#'     \item \strong{The engine's package is not installed}, so its formals
#'       cannot be read. Building the arguments still works -- it is not
#'       running the engine -- and the engine's own requirement check reports
#'       the missing package when you actually run it.
#'   }
#'
#'   Every other engine enumerates its arguments, so a misspelling is an error
#'   at the call site naming the accepted set. Constructors for pecotmr's own
#'   engines take no \code{...} at all, so R itself rejects an unknown name as
#'   an unused argument.
#' @aliases MethodOptions
#' @importClassesFrom S4Vectors SimpleList
#' @export
setClass("MethodOptions", contains = "SimpleList")

# Resolve "pkg::fn" to the function, or NULL when there is nothing to resolve:
# a package that is not installed, or `callee = NULL` for a bundle pecotmr owns
# outright. A missing Suggests package is not an error here -- building an
# argument list is not running the engine, and the engine's own
# requireNamespace() guard gives a better message where it actually matters.
# @noRd
.engineCallee <- function(callee) {
    if (is.null(callee)) {
        return(NULL)
    }
    if (is.function(callee)) {
        return(callee)
    }
    parts <- str_split_1(callee, fixed("::"))
    if (length(parts) != 2L) {
        abort(glue("engine callee must be 'pkg::fn', got '{callee}'"))
    }
    if (!requireNamespace(parts[[1L]], quietly = TRUE)) {
        return(NULL)
    }
    tryCatch(
        get(parts[[2L]], envir = asNamespace(parts[[1L]])),
        error = function(cnd) NULL
    )
}

# The names an engine will accept, or NULL when that cannot be determined --
# either the package is absent or the engine takes `...`, in which case it
# accepts anything and there is nothing to check against.
#
# Several callees may be given, for a bundle that reaches more than one
# function. What can be checked then depends on HOW the caller forwards it:
#
#   filtered = TRUE  -- the caller intersects the bundle with each callee's
#     explicit formals before calling it (see .ctwasInvoke). A name outside
#     the union of those formals reaches no callee, so it can only be a
#     mistake, and `...` on any one callee changes nothing.
#   filtered = FALSE -- the caller splices the bundle straight in, so a
#     callee's `...` really will accept anything and nothing can be rejected.
# @noRd
.engineAcceptedNames <- function(callee, filtered = FALSE) {
    if (length(callee) > 1L) {
        if (!filtered && any(map_lgl(callee, .engineTakesDots))) {
            return(NULL)
        }
        return(.engineUnionNames(callee))
    }
    fn <- .engineCallee(callee)
    if (is.null(fn)) {
        return(NULL)
    }
    fm <- tryCatch(names(formals(fn)), error = function(cnd) NULL)
    if (is.null(fm) || is_in("...", fm)) {
        return(NULL)
    }
    fm
}

# TRUE when a callee takes `...`, so it accepts any name. An unreadable
# callee is treated as not taking dots; .engineUnionNames drops it anyway.
# @noRd
.engineTakesDots <- function(callee) {
    fn <- .engineCallee(callee)
    if (is.null(fn)) {
        return(FALSE)
    }
    fm <- tryCatch(names(formals(fn)), error = function(cnd) NULL)
    !is.null(fm) && is_in("...", fm)
}

# Union of several callees' explicit formals; NULL when none could be read,
# which keeps a missing package from turning every name into an error.
# @noRd
.engineUnionNames <- function(callees) {
    fns <- map(callees, .engineCallee)
    known <- compact(fns)
    if (length(known) == 0L) {
        return(NULL)
    }
    sort(unique(list_c(map(known, .engineExplicitFormals))))
}

# One callee's formals with `...` removed, NULL when they cannot be read.
# @noRd
.engineExplicitFormals <- function(fn) {
    setdiff(tryCatch(names(formals(fn)), error = function(cnd) NULL), "...")
}

# --- which input classes a set of engine arguments belongs to -------------
#
# A method that runs on both individual-level and summary-statistic data
# reaches a DIFFERENT engine on each path -- often in a different package
# (glmnet vs lassosum, RcppDPR vs SDPR) -- so an Options record is valid for
# one path and meaningless for the other. Recording which lets a pipeline
# reject `SusieRssOptions()` for a QtlDataset by asking the record, instead
# of joining against a registry at the point of use.
#
# Derived from the two dispatch tables, never declared per constructor: a
# hand-written `inputType = "QtlSumStats"` in 24 constructors would be 24
# more things to keep in step with the tables that actually drive dispatch.
#
# NULL for an engine with no input-path dimension at all -- mash prior
# components, p-value tests, ctwas, coloc -- which is most of them.
# @noRd
.engineInputType <- function(callee, label) {
    unique(c(
        .fmCalleeInputType(callee),
        .twasLabelInputType(label)
    ))
}

# Fine-mapping stores the callees themselves, so the lookup is direct. A
# callee may serve several methods (susieR::susie backs susie, susieInf and
# susieAsh), and GwasSumStats is a summary-statistics class that only the
# univariate methods accept -- exactly what `gwasAllowed` already records.
# @noRd
.fmCalleeInputType <- function(callee) {
    if (!is.character(callee) || length(callee) == 0L) {
        return(character(0))
    }
    caps <- .fineMappingMethodCapabilities
    ind <- compact(map(caps, "individualImpl"))
    sum <- compact(map(caps, "sumstatImpl"))
    gwas <- keep(names(sum), function(tk) isTRUE(caps[[tk]]$gwasAllowed))
    c(
        if (any(is_in(callee, unlist(ind)))) "QtlDataset",
        if (any(is_in(callee, unlist(sum)))) "QtlSumStats",
        if (any(is_in(callee, unlist(sum[gwas])))) "GwasSumStats"
    )
}

# TWAS stores pecotmr implementations rather than callees, and each one
# names its own Options constructor through its `methodArgs` default -- so
# the lookup is by constructor name, which is what `label` carries.
# @noRd
.twasLabelInputType <- function(label) {
    if (!is.character(label) || length(label) != 1L) {
        return(character(0))
    }
    caps <- .twasMethodCapabilities
    named <- function(field) {
        unlist(compact(map(
            compact(map(caps, field)),
            .twasImplCtorName
        )))
    }
    c(
        if (is_in(label, named("individualImpl"))) "QtlDataset",
        if (is_in(label, named("sumstatImpl"))) "QtlSumStats"
    )
}

# The input classes a record is valid for, or NULL when it has no path
# dimension. Reads what the constructor derived; never recomputes.
# @noRd
.optionsInputType <- function(x) {
    if (!.isMethodOptions(x)) {
        return(NULL)
    }
    it <- metadata(x)$inputType
    if (length(it) == 0L) NULL else it
}

# Reject extras the engine cannot accept. `accepted = NULL` means the check is
# not possible (absent package, or an engine with `...`), so the extras pass
# through untouched -- deliberately, since such an engine really does take
# arbitrary names.
# @noRd
.engineCheckExtra <- function(extra, accepted, label, calleeName) {
    if (length(extra) == 0L || is.null(accepted)) {
        return(invisible(NULL))
    }
    unknown <- setdiff(names(extra), accepted)
    if (length(unknown) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{label}: unknown argument(s) ",
        "{str_flatten(unknown, ', ')}. ",
        "`{calleeName}` accepts: {str_flatten(accepted, ', ')}."
    ))
}

# The record's entries: pecotmr's defaults, with an explicit value
# replacing one of the same name.
#
# Two names colliding inside `extra` is the caller saying the same thing
# twice -- including a wrapper's own contribution meeting the caller's
# record under .methodOptionsWith() -- and there is no basis for picking a
# winner. A `defaults` entry is only pecotmr's preference, so an explicit
# value replaces it in place, keeping the default's position; an argument
# pecotmr truly owns is refused outright by the constructor rather than
# left to collide here.
# @noRd
.configMerge <- function(defaults, extra, label) {
    dup <- names(extra)[duplicated(names(extra))]
    if (length(dup) > 0L) {
        abort(glue(
            "{label}: argument(s) given twice: {str_flatten(dup, ', ')}"
        ))
    }
    defaults <- discard(defaults, is.null)
    overridden <- intersect(names(defaults), names(extra))
    defaults[overridden] <- extra[overridden]
    c(defaults, extra[!names(extra) %in% overridden])
}

# Build one engine's argument record: pecotmr's defaults, overlaid with the
# caller's extras, validated against the engine where that is possible.
#
# `defaults` is what the constructor's own formals resolved to; `extra` is its
# `...`. A NULL default means "leave it to the engine", so it is dropped rather
# than forwarded as an explicit NULL, which some engines treat differently from
# an absent argument. Only NULL: `compact()` would also drop a zero-length
# value, silently losing `keepSamples = character(0)` and
# `greedy_args = list()` -- both meaningful settings, not absent ones.
# Values are carried as given -- matrices, fits, functions and all -- since
# what an engine accepts is the engine's business; only the argument NAMES
# are checkable.
#
# `callee = NULL` marks a bundle pecotmr owns rather than forwards. Such a
# constructor declares every field as a formal and takes no `...`, so R's own
# argument matching rejects a misspelling before this is ever reached; there
# is no foreign signature to check against.
# @noRd
#' @importFrom rlang is_named
#' @importFrom purrr map2_lgl
#' @importFrom S4Vectors SimpleList metadata metadata<-
.newMethodOptions <- function(
    callee,
    defaults,
    extra,
    label,
    engine = NULL,
    accepted = NULL,
    filtered = FALSE,
    check = TRUE,
    inputType = NULL
) {
    if (length(extra) > 0L && !is_named(extra)) {
        abort(glue("{label}: all arguments passed through `...` must be named"))
    }
    # A callee handed over as a function -- or absent entirely, for a bundle
    # pecotmr owns -- has no spellable name, so the label stands in for it.
    calleeName <- if (is.character(callee)) {
        str_flatten(callee, " / ")
    } else {
        label
    }
    # `accepted` given explicitly is for an engine whose options are not its
    # formals -- udr takes a `control` list, whose valid names come from
    # udr::ud_fit_control_default(). Still a live source, not a transcription.
    # `check = FALSE` is for a bundle that reaches an engine whose accepted
    # names pecotmr has not enumerated. Rejecting a name there would be a
    # guess, and a false rejection is worse than no check.
    accepted <- if (!check) {
        NULL
    } else {
        accepted %||% .engineAcceptedNames(callee, filtered)
    }
    .engineCheckExtra(extra, accepted, label, calleeName)
    res <- new(
        "MethodOptions",
        SimpleList(.configMerge(defaults, extra, label))
    )
    metadata(res) <- list(
        engine = engine %||% label,
        callee = if (is.character(callee)) callee else NA_character_,
        filtered = filtered,
        check = check,
        inputType = inputType %||% .engineInputType(callee, label)
    )
    res
}

# Refuse names the caller cannot usefully set, naming where the setting
# lives instead. Same shape as .ctwasRefusePipelineOwned: a real formal of
# the engine that pecotmr supplies itself, so a value given here would be
# overwritten rather than honoured. `owned` is a named character vector,
# name = the refused argument, value = where the caller sets it instead.
# @noRd
.configRefuseOwned <- function(extra, owned, label) {
    clash <- intersect(names(extra), names(owned))
    if (length(clash) == 0L) {
        return(invisible(NULL))
    }
    where <- str_flatten(
        sprintf("`%s` -> %s", clash, unname(owned[clash])),
        collapse = "; "
    )
    abort(glue(
        "{label}: {str_flatten(clash, \', \')} ",
        "{if (length(clash) == 1L) \'is\' else \'are\'} supplied by ",
        "pecotmr, so a value given here would be overwritten. Use: {where}."
    ))
}

# TRUE for a record a constructor produced. Used by the pipelines to refuse a
# bare list, which is the whole point of requiring constructors.
# @noRd
.isMethodOptions <- function(x) {
    is(x, "MethodOptions")
}

# Demand a constructor-built record, naming the constructor the caller should
# have used so the error is actionable rather than merely correct.
# @noRd
.assertMethodOptions <- function(x, ctor, argName) {
    if (.isMethodOptions(x)) {
        return(invisible(x))
    }
    if (is.list(x) && length(x) == 0L) {
        # An empty list carries no settings, so there is nothing to validate
        # and nothing a constructor would have added. Accepting it keeps
        # `list()` working as "no options".
        return(invisible(x))
    }
    abort(glue(
        "`{argName}` must be built with {ctor}(), not a bare list. ",
        "A bare list cannot be checked against the engine, so a misspelled ",
        "option would be silently ignored."
    ))
}

# A wrapper that contributes an engine option of its own -- bayesCWeights's
# `pi`, dprVbWeights's `n_k` -- cannot simply c() it onto the caller's record:
# c() appends a MethodOptions as one opaque element rather than splicing it, and
# even flattening first would hand the next hop a bare list its
# .assertMethodOptions() then refuses. Rebuild through the constructor instead,
# so what travels on is still a checked record -- and a caller who also set
# that option through `methodArgs` gets .newMethodOptions()'s "given twice"
# error rather than a silent override.
# @noRd
.methodOptionsWith <- function(x, ctor, ...) {
    exec(match.fun(ctor), ..., !!!as.list(x))
}

# --- character-or-constructor selection ------------------------------------
#
# Anywhere pecotmr lets the caller pick an engine, the argument takes either
# the engine's name or that engine's constructor. The constructor already
# carries the identity in metadata(), so `engine = CovUdrOptions(...)` says
# which engine AND how to configure it in one value, and a bare string stays
# the no-options shorthand.

# The engine a selection denotes, whether it arrived as a name or a record.
# @noRd
.engineOf <- function(x) {
    if (.isMethodOptions(x)) {
        return(metadata(x)$engine)
    }
    # A Param carries no metadata() -- its class IS its identity, so the
    # engine comes from the class name. DentistParam() and SlalomParam() are
    # both engine choices for `ldMismatchQcMethod`, so this arm is load
    # bearing, not defensive.
    if (.isMethodParam(x)) {
        return(.paramEngineOf(x))
    }
    if (is.character(x) && length(x) == 1L) {
        return(x)
    }
    NULL
}

# Resolve a character-or-constructor selection into the chosen engine and its
# arguments. `allowed` is the permitted engine vocabulary; a record built by a
# constructor outside it is rejected by name rather than left to fail further
# down.
# @noRd
.resolveEngineChoice <- function(x, allowed, argName) {
    engine <- .engineOf(x)
    if (is.null(engine)) {
        abort(glue(
            "`{argName}` must be one of {str_flatten(allowed, ', ')}, ",
            "or the matching constructor."
        ))
    }
    if (!is_in(engine, allowed)) {
        abort(glue(
            "`{argName}`: unknown engine '{engine}'. ",
            "Known: {str_flatten(allowed, ', ')}."
        ))
    }
    list(
        engine = engine,
        args = if (.isMethodOptions(x) || .isMethodParam(x)) {
            x
        } else {
            .newMethodOptions(
                NULL,
                defaults = list(),
                extra = list(),
                label = argName,
                engine = engine
            )
        }
    )
}

# --- nested per-engine bundles ---------------------------------------------

# Every entry must have been built by a constructor its key accepts.
# Checked before building so the message is the user's, not purrr's "In
# index: 1. With name: foo." framing around it. The acceptable engines come
# from the key's own constructor rather than from the key: several keys may
# share one engine (lasso and enet are both glmnet). `engines` overrides
# that for a key whose paths run more than one -- TWAS's lasso fits with
# glmnet on individual data and with pecotmr's lassosum solver on summary
# statistics, and either constructor is a legitimate way to configure it.
# @noRd
.nestedAssertEngines <- function(elements, ctors, engines, label) {
    declared <- map_chr(elements, .nestedDeclaredEngine)
    allowed <- map(
        names(elements),
        .nestedAllowedEngines,
        ctors = ctors,
        engines = engines
    )
    mismatch <- which(!map2_lgl(declared, allowed, .nestedEngineMatches))
    if (length(mismatch) == 0L) {
        return(invisible(NULL))
    }
    i <- mismatch[[1L]]
    abort(glue(
        "{label}: the `{names(elements)[i]}` entry was built with the ",
        "constructor for '{declared[[i]]}', but `{names(elements)[i]}` ",
        "uses the ",
        "{str_flatten(sprintf(\"'%s'\", allowed[[i]]), ' or ')} ",
        "constructor. ",
        "{if (length(allowed[[i]]) > 1L) 'Pass one of those' else ",
        "'Pass that'}, or a plain list."
    ))
}

# A plain-list entry declares no engine, so nothing to match.
# @noRd
.nestedEngineMatches <- function(declared, allowed) {
    is.na(declared) || is_in(declared, allowed)
}

# The engines a key accepts: the caller's list for that key when given,
# otherwise the single one the key's own constructor produces.
# @noRd
.nestedAllowedEngines <- function(key, ctors, engines) {
    engines[[key]] %||% .nestedExpectedEngine(key, ctors)
}

# The engine a key's own constructor produces, so an entry built by a
# different engine's constructor can be told apart from one built by a
# sibling key that shares the same engine.
# @noRd
.nestedExpectedEngine <- function(key, ctors) {
    rec <- tryCatch(ctors[[key]](), error = function(cnd) NULL)
    if (is.null(rec)) {
        return(key)
    }
    metadata(rec)$engine %||% key
}

# The engine a nested entry was built for, or NA when it is a plain list.
# @noRd
.nestedDeclaredEngine <- function(value) {
    if (.isMethodOptions(value)) metadata(value)$engine else NA_character_
}

# One element of a nested bundle: pass an already-built record through, or
# splice a plain list into the method's constructor.
# @noRd
.nestedElement <- function(value, key, ctors, label) {
    if (.isMethodOptions(value)) {
        return(value)
    }
    if (!is.list(value)) {
        abort(glue(
            "{label}: the `{key}` entry must be a list or a ",
            "constructor result."
        ))
    }
    exec(ctors[[key]], !!!value)
}

# Why a given callee's arguments can or cannot be checked. Three reasons a
# name goes unchecked, and the user is entitled to know which one applies:
# the engine takes `...` so any name is legal, the engine's package is not
# installed so its formals cannot be read, or -- the normal case -- the
# engine enumerates its arguments and the names were checked.
# @noRd
.engineCheckNote <- function(callee, filtered = FALSE) {
    if (length(callee) > 1L) {
        accepted <- .engineAcceptedNames(callee, filtered)
        if (is.null(accepted)) {
            return(glue(
                "argument names NOT checked: one of these functions takes ",
                "`...`, so any name is accepted"
            ))
        }
        return(glue(
            "argument names checked against the {length(accepted)} arguments ",
            "these {length(callee)} functions accept between them"
        ))
    }
    fn <- .engineCallee(callee)
    if (is.null(fn)) {
        return(glue(
            "argument names NOT checked: package ",
            "'{str_split_1(callee, fixed('::'))[[1L]]}' is not installed"
        ))
    }
    fm <- tryCatch(names(formals(fn)), error = function(cnd) NULL)
    if (is.null(fm)) {
        return("argument names NOT checked: the engine's formals are unknown")
    }
    if (is_in("...", fm)) {
        return(glue(
            "argument names NOT checked: {callee}() takes `...`, ",
            "so any name is accepted"
        ))
    }
    glue("argument names checked against {callee}()")
}

#' @rdname MethodOptions-class
#' @param object A \code{MethodOptions} object.
#' @return \code{show} is called for its side effect and returns
#'   \code{invisible(NULL)}.
#' @export
setMethod("show", "MethodOptions", function(object) {
    md <- metadata(object)
    cat(glue("<MethodOptions: {md$engine %||% 'unknown'}>"), "\n", sep = "")
    if (!is.null(md$callee) && !anyNA(md$callee)) {
        cat(
            glue("  forwarded to {str_flatten(md$callee, '(), ')}()"),
            "\n",
            sep = ""
        )
        # Whether names were checkable is a property of the callee, not of
        # this record, so say it here rather than leaving the user to guess
        # why a typo was or was not caught.
        note <- if (isFALSE(md$check)) {
            glue(
                "argument names NOT checked: this method reaches a different ",
                "engine on the summary-statistics path, whose accepted names ",
                "pecotmr has not enumerated"
            )
        } else {
            .engineCheckNote(md$callee, md$filtered %||% FALSE)
        }
        cat("  ", note, "\n", sep = "")
    }
    # Which input classes these arguments belong to, when that is a
    # meaningful question: a method running on both data paths reaches a
    # different engine on each, so a record built for one is not usable on
    # the other.
    it <- .optionsInputType(object)
    if (!is.null(it)) {
        cat(
            "  for input class ",
            str_flatten(sort(it), ", "),
            "\n",
            sep = ""
        )
    }
    if (length(object) == 0L) {
        cat("  (no options set)\n")
        return(invisible(NULL))
    }
    for (nm in names(object)) {
        cat(sprintf(
            "  %-22s %s\n",
            nm,
            str_flatten(format(object[[nm]]), ", ")
        ))
    }
    invisible(NULL)
})
