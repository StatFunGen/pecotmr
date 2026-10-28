# The fit-coverage warning an expression raises, or the expression's value
# when it raises none. Every other warning is muffled: a pipeline run on test
# data warns about things these checks do not look at.
.fitCoverageCondition <- function(expr) {
    suppressWarnings(tryCatch(expr, pecotmrFitCoverage = identity))
}

# The message of every fit-coverage warning an expression raises, run to the
# end; other warnings are muffled.
.fitCoverageMessages <- function(expr) {
    msgs <- character(0)
    suppressWarnings(withCallingHandlers(
        expr,
        pecotmrFitCoverage = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    ))
    msgs
}
