# Skips for tests of what needs a newer susieR, mvsusieR or mr.mashr than the
# one installed. Each asks what the package itself asks before it calls the
# engine, so the test runs once the installed engine has the feature.

# Skip unless mvsusieR fits summary statistics across studies
# (.mvsusieFitsStudies()).
.skipIfMvsusieCannotFitStudies <- function() {
    skip_if_not_installed("mvsusieR")
    skip_if_not(.mvsusieFitsStudies(), "mvsusieR predates fits across studies")
}

# Skip unless mr.mashr fits summary statistics across studies
# (.mrmashFitsStudies()).
.skipIfMrmashCannotFitStudies <- function() {
    skip_if_not_installed("mr.mashr")
    skip_if_not(.mrmashFitsStudies(), "mr.mashr predates fits across studies")
}

# Skip unless mr.mashr takes a prior relative to the residual standard
# deviations (.mrmashRelativeS0()).
.skipIfNoMrmashRelativeS0 <- function() {
    skip_if_not_installed("mr.mashr")
    skip_if_not(.mrmashRelativeS0(), "mr.mashr predates relative_S0")
}

# Skip unless susieR has its internal function `fn`, which came with `what`:
# a behaviour of a newer susieR that none of its arguments shows.
.skipIfSusieRPredates <- function(fn, what) {
    skip_if_not(
        exists(fn, envir = asNamespace("susieR"), inherits = FALSE),
        str_c("susieR predates ", what)
    )
}

# Skip unless susieR gives the SuSiEx marginals without enumerating the
# patterns, and enumerates them past 20 conditions (.susiexTakes()).
.skipIfSusiexPredates <- function() {
    skip_if_not(
        .susiexTakes("config_probability"),
        "susieR predates SuSiEx marginals without enumeration"
    )
}
