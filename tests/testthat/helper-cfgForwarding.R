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
