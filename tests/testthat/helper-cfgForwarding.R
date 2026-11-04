# Shared by the cfg-forwarding tests, which live in
# test_fineMappingPipeline.R and test_mashPipeline.R.
#
# The per-study / per-sumstats recursions hand their cfg record on as
# explicit `name = cfg$name` pairs. R catches the wrong direction for free:
# a name the entry point does not accept is an "unused argument" error, even
# through an S4 generic's `...`. The direction R cannot catch is OMISSION --
# a field added to the record and forgotten at one call site silently takes
# the callee's default, which is how multi-study TWAS once inherited the
# per-study cross-validation default. This pins that direction.

# The field names a cfg builder emits, read off its terminal list().

.cfgFieldNames <- function(fnName) {
    f <- get(fnName, envir = asNamespace("pecotmr"))
    lst <- NULL
    rec <- function(e) {
        if (!is.call(e)) {
            return(invisible(NULL))
        }
        if (identical(e[[1]], as.name("list")) && !is.null(names(e))) {
            nm <- names(e)[nzchar(names(e))]
            if (length(nm) > length(lst)) lst <<- nm
        }
        for (i in seq_along(e)) {
            if (i == 1L || rlang::is_missing(e[[i]])) {
                next
            }
            rec(e[[i]])
        }
    }
    rec(body(f))
    lst
}

# Every `name = ...` the function's body passes to anything.
.argNamesPassed <- function(fnName) {
    f <- get(fnName, envir = asNamespace("pecotmr"))
    out <- character(0)
    rec <- function(e) {
        if (!is.call(e)) {
            return(invisible(NULL))
        }
        nm <- names(e)
        if (!is.null(nm)) {
            out <<- c(out, nm[nzchar(nm)])
        }
        for (i in seq_along(e)) {
            if (i == 1L || rlang::is_missing(e[[i]])) {
                next
            }
            rec(e[[i]])
        }
    }
    rec(body(f))
    unique(out)
}

# The arguments an S4 method takes, minus `data` and `...`. A method whose
# signature differs from its generic's is stored as a wrapper around an inner
# `.local`, whose formals are the ones the caller sees.
.methodFormalNames <- function(generic, cls) {
    f <- methods::getMethod(generic, cls)@.Data
    b <- body(f)
    isLocal <- is.call(b) &&
        identical(b[[1]], as.name("{")) &&
        is.call(b[[2]]) &&
        identical(b[[2]][[2]], as.name(".local"))
    if (isLocal) {
        f <- eval(b[[2]][[3]])
    }
    setdiff(names(formals(f)), c("data", "..."))
}

# Stand in for the joint engine's entry point and record what each joint
# dispatcher hands it: the data form, the dataset, the pipeline config and the
# engine `args`. One record per call in `cap$calls`. It returns NULL (no fits),
# so a run whose joint phase consumed every method then stops with "no joint
# fits produced" -- the tests expect that and read the records afterwards.
.captureJointEngine <- function(cap, env = parent.frame()) {
    cap$calls <- list()
    testthat::local_mocked_bindings(
        .jointSpecsWork = function(
            parsedJointSpec,
            data,
            dataForm,
            pipeline,
            jointMethods,
            contexts,
            traitIds,
            args = list()
        ) {
            cap$calls[[length(cap$calls) + 1L]] <- list(
                dataForm = dataForm,
                data = data,
                config = pecotmr:::.jpConfig(pipeline),
                args = args,
                jointMethods = jointMethods
            )
            pecotmr:::.qtlWorkValue(NULL)
        },
        .package = "pecotmr",
        .env = env
    )
}

# The recorded call for one data form ("individual" / "sumstats").
.capturedCall <- function(cap, dataForm) {
    hits <- Filter(function(r) identical(r$dataForm, dataForm), cap$calls)
    if (length(hits) == 0L) NULL else hits[[1L]]
}
