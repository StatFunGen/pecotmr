# =============================================================================
# Rebuild tests/testthat/test_data/sldscUpstreamZ.rds -- the z-scores the
# upstream-parity tests in test_h2EstimationWrappers.R are compared against.
# -----------------------------------------------------------------------------
# These are stored rather than re-simulated inside the test because the draw
# uses eigen(), and eigen() fixes eigenvectors only up to SIGN. The sign a
# given LAPACK build returns is arbitrary, so regenerating from a seed
# produced different data on linux-64 than on osx-arm64 / linux-aarch64 --
# and the expected numbers came from running upstream (ldsc's
# ldscore.regressions.Hsq, HDL, LDER, gldsc) on one platform's realization.
# That made three of those tests fail on linux-64 only.
#
# Two things must move together, or the tests silently stop meaning anything:
#
#   1. Rerunning this script invalidates the expected values in
#      test_h2EstimationWrappers.R. Re-derive them from upstream on the NEW
#      vectors at the same time.
#   2. These vectors are drawn over ldScoreExample's LD blocks, so rebuilding
#      that object (build_h2_examples.R) requires rebuilding these too.
# =============================================================================

suppressMessages(devtools::load_all(".", quiet = TRUE, export_all = FALSE))
data(ldScoreExample)

# EXACT copy of the in-test generator, so the captured z is the one the
# pinned upstream values were produced from.
simZ <- function(ref, h2, n, seed) {
    M <- length(ref)
    perSnpVar <- rep(h2 / M, M)
    set.seed(seed)
    z <- numeric(M)
    for (b in pecotmr:::getLdMatrixList(ref)) {
        idx <- b$snpIdx
        p <- length(idx)
        sigma <- n * (b$R %*% diag(perSnpVar[idx], p) %*% b$R) + b$R
        e <- eigen((sigma + t(sigma)) / 2, symmetric = TRUE)
        z[idx] <- as.vector(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(p)))
    }
    z
}
# Every (h2, n, seed) the test file asks for. Keep in step with the call
# sites of .sldscUpstreamZ() in test_h2EstimationWrappers.R -- a missing
# entry is a hard error there, not a silent re-simulation.
cases <- c(
    list(
        c(0.4, 1e4, 42), c(0.15, 4e3, 7), c(0.6, 2.5e4, 99),
        c(0.3, 1e4, 5), c(0.6, 5e4, 9), c(0.4, 5e4, 3)
    ),
    # the HDL recovery test sweeps seeds 1:5 at (h2 = 0.4, n = 1e5)
    lapply(1:5, function(s) c(0.4, 1e5, s))
)
out <- list()
for (cs in cases) {
    key <- paste(cs[[1]], cs[[2]], cs[[3]], sep = "_")
    out[[key]] <- simZ(ldScoreExample, cs[[1]], cs[[2]], cs[[3]])
}
saveRDS(out, "tests/testthat/test_data/sldscUpstreamZ.rds", compress = "xz")
cat("keys:", paste(names(out), collapse = " | "), "\n")
cat("length each:", unique(lengths(out)), "\n")
cat("file size:",
    file.size("tests/testthat/test_data/sldscUpstreamZ.rds"), "bytes\n")
