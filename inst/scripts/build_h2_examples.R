# =============================================================================
# Rebuild ldEigenExample / ldScoreExample from a real LD reference.
# -----------------------------------------------------------------------------
# These two objects feed estimateH2(). The versions they replace were 20 SNPs
# in 2 blocks -- enough to demonstrate the accessors, but far too small to run
# a heritability estimator: every method returned a boundary value near zero,
# and two blocks cannot support the delete-one-block jackknife the standard
# errors come from.
#
# They are now built from the xqtl-protocol MWE's chr22 LD reference: real
# genotype-derived LD over the same de-identified `protocol_example` panel the
# package already ships a slice of, using the LD blocks the protocol defines.
# `buildLdEigen()` / `buildLdScore()` do the work, so this script doubles as a
# worked example of them.
#
# Requires the MWE checkout; set MWE_DIR to point at it:
#   MWE_DIR=/path/to/xqtl-protocol Rscript inst/scripts/build_h2_examples.R
# =============================================================================

devtools::load_all(".", quiet = TRUE)

mweDir <- Sys.getenv(
    "MWE_DIR",
    "/Users/danielnachun/Downloads/fungen_xqtl/xqtl-protocol"
)
metaPath <- file.path(
    mweDir,
    "input/ld_reference/protocol_example.ld_meta_file.tsv"
)
if (!file.exists(metaPath)) {
    stop("Set MWE_DIR to an xqtl-protocol checkout; not found: ", metaPath)
}

# How many of chr22's LD blocks to keep. All 20 work, but the eigenvectors of
# a real LD block compress poorly, so the pair costs ~2.4 MB at 20 against
# ~1.2 MB at 14 -- and 14 blocks still leaves a usable jackknife.
N_BLOCKS <- 14L

blocks <- read.delim(metaPath, check.names = FALSE)
names(blocks)[1] <- "chr"
chr22 <- blocks[blocks$chr == "chr22", ]
chr22 <- head(chr22[order(chr22$start), ], N_BLOCKS)
regions <- sprintf("chr22:%d-%d", chr22$start, chr22$end)
message("Loading ", length(regions), " chr22 LD blocks...")

# One LdData per block: buildLd*() takes the list and gives one LD block per
# element. Loading the whole span at once would instead return a single dense
# matrix, i.e. one block, which no jackknife can work with.
ldList <- lapply(regions, function(r) loadLdMatrix(metaPath, region = r))
message(
    "  ",
    sum(vapply(ldList, length, integer(1))),
    " variants across ",
    length(ldList),
    " blocks"
)

ldEigenExample <- buildLdEigen(ldList, genome = "hg38")
ldScoreExample <- buildLdScore(ldList, genome = "hg38")

stopifnot(
    length(ldEigenExample) == length(ldScoreExample),
    length(getEigenList(ldEigenExample)) == N_BLOCKS,
    length(getLdMatrixList(ldScoreExample)) == N_BLOCKS,
    identical(names(ldEigenExample), names(ldScoreExample))
)

# The two routes to an LD score must agree, or one of them has drifted.
stopifnot(all.equal(
    as.vector(getLdScores(ldScoreExample)[, 1]),
    as.vector(computeLdScores(ldEigenExample)[, 1])
))

# -----------------------------------------------------------------------------
# Confirm every estimator recovers a known h2 on this reference before saving.
# z is simulated under the standard model z ~ N(0, N R diag(v) R + R).
# -----------------------------------------------------------------------------
M <- length(ldEigenExample)
N <- 1e5
h2True <- 0.4
perSnpVar <- rep(h2True / M, M)

simZ <- function() {
    z <- numeric(M)
    for (block in getLdMatrixList(ldScoreExample)) {
        idx <- block$snpIdx
        p <- length(idx)
        sigma <- N * (block$R %*% diag(perSnpVar[idx], p) %*% block$R) + block$R
        sigma <- (sigma + t(sigma)) / 2
        e <- eigen(sigma, symmetric = TRUE)
        z[idx] <- as.vector(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(p)))
    }
    z
}

set.seed(1)
recovered <- replicate(5, {
    z <- simZ()
    c(
        sldsc = pecotmr:::sldscUnivariate(z, N, ldScoreExample)$h2,
        gldsc = pecotmr:::gldscUnivariate(z, N, ldScoreExample)$h2,
        lder = pecotmr:::lderUnivariate(z, N, ldEigenExample)$h2
    )
})
message("h2 recovery (true 0.4):")
print(round(rowMeans(recovered), 3))
stopifnot(all(abs(rowMeans(recovered) - h2True) < 0.15))

# HDL is deliberately excluded from that check: it models reference-panel
# noise, and this panel's nRef is 1000, so it correctly shrinks hard. That is
# the estimator working, not failing -- see the heritability vignette.

usethis::use_data(ldEigenExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ldScoreExample, overwrite = TRUE, compress = "xz")
