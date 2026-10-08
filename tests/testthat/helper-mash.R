# Shared fixtures for test_mashPipeline.R.
#
# mashr / flashier fits dominate that file's runtime. These are coverage tests
# for pecotmr's wiring, not numerical checks of mashr, so the fixture is the
# smallest slice that still drives every branch, and the prior is supplied
# outright wherever the prior itself is not what the test is about.

# The multi-context example trimmed to `n` variants per context. Both prior
# estimation and the mash() weight fit scale with the variant count, and the
# full 200 buys nothing a branch test can see.
mashFixture <- function(n = 60L) {
    utils::data(
        "qtlSumStatsMulticontextExample",
        package = "pecotmr",
        envir = environment()
    )
    full <- get("qtlSumStatsMulticontextExample", envir = environment())
    if (n >= min(lengths(full))) {
        return(full)
    }
    endoapply(full, function(g) g[seq_len(n)])
}

# A two-component prior for tests where the prior is not under test. Supplying
# it takes the `priorCovariances` branch and skips the flash/pca/ed chain --
# most of the cost -- while still running the mash() weight fit, so anything
# downstream of the prior is exercised exactly as before.
mashTinyPrior <- function(k = 3L) {
    list(identity = diag(k), shared = matrix(1, k, k))
}

# Shared by test_mashPipeline.R and test_mashWrapper.R.
.mashTestModel <- function(ss) {
    suppressMessages(suppressWarnings(
        mashModelFit(
            list(random = ss),
            alpha = 0,
            priorCovariances = list(
                identity = diag(3),
                effectA = diag(c(1, 0, 0))
            ),
            vhat = diag(3),
            setSeed = 1L
        )
    ))
}
