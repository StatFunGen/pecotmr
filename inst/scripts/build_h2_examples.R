# =============================================================================
# Rebuild ldEigenExample / ldScoreExample: LD references for estimateH2().
# -----------------------------------------------------------------------------
# These two objects feed estimateH2(). An h2 estimator's reference panel has
# to be large: ldsc's unbiased r^2 and HDL's finite-reference term both
# correct for LD estimated from n samples, and at the size of the panel the
# package can ship (the xqtl-protocol MWE's 60 people) those corrections swamp
# the signal. What the examples hold is LD, not samples, so their size does
# not grow with n -- the panel can be large without being shipped.
#
# So the reference is SIMULATED from real LD:
#
#   1. The population LD is the xqtl-protocol MWE's chr22 LD, block by block
#      over the LD blocks the protocol defines (the first N_BLOCKS of them).
#   2. A panel of N_SIM people is drawn from that population: multivariate
#      normal genotypes with that correlation.
#   3. Each block's LD is that panel's sample correlation, wrapped in an
#      LdMatrix that records N_SIM and the variants' real allele frequencies,
#      and buildLdEigen() / buildLdScore() do the rest -- so this script
#      doubles as a worked example of them.
#
# Requires the MWE checkout; set MWE_DIR to point at it:
#   MWE_DIR=/path/to/xqtl-protocol Rscript inst/scripts/build_h2_examples.R
#
# The parity tests in test_h2EstimationWrappers.R do not read these objects
# (they hold their own frozen LD), so rebuilding them moves no pinned value.
# =============================================================================

devtools::load_all(".", quiet = TRUE)
library(purrr)

mweDir <- Sys.getenv(
    "MWE_DIR",
    "/Users/danielnachun/Downloads/fungen_xqtl/xqtl-protocol"
)
refDir <- file.path(mweDir, "input/ld_reference")
metaPath <- file.path(refDir, "protocol_example.ld_meta_file.tsv")
if (!file.exists(metaPath)) {
    stop("Set MWE_DIR to an xqtl-protocol checkout; not found: ", metaPath)
}

# How many of chr22's LD blocks to keep. All 20 work, but the eigenvectors of
# a real LD block compress poorly, so the pair costs ~2.4 MB at 20 against
# ~1.2 MB at 14 -- and 14 blocks still leaves a usable jackknife.
N_BLOCKS <- 14L
# The simulated reference panel's size: large enough that LD-estimation noise
# is small next to the signal, as in a real biobank-scale reference.
N_SIM <- 20000L

blocks <- read.delim(metaPath, check.names = FALSE)
names(blocks)[1] <- "chr"
chr22 <- blocks[blocks$chr == "chr22", ]
chr22 <- head(chr22[order(chr22$start), ], N_BLOCKS)
message("Reading the chr22 population LD over ", nrow(chr22), " blocks...")

# The MWE reference is a stochastic sketch of 60 people; readGenotypes()
# records their real allele frequencies on the panel from its .afreq.
panel <- readGenotypes(
    plink2Prefix = file.path(refDir, "chr22", "protocol_example.LD.chr22")
)
gr <- SummarizedExperiment::rowRanges(panel)

# The panel rows inside one block. Blocks share their boundary coordinate, so
# each is taken half-open, [start, end), to give every variant one block.
blockIdx <- function(start, end) {
    which(GenomicRanges::start(gr) >= start & GenomicRanges::start(gr) < end)
}
idxList <- compact(map2(chr22$start, chr22$end, function(s, e) {
    idx <- blockIdx(s, e)
    if (length(idx) >= 2L) idx
}))

# Draw N_SIM people with correlation R (PSD, possibly rank-deficient: the MWE
# LD comes from 60 people) and return their sample correlation.
simulatedLd <- function(R) {
    e <- eigen((R + t(R)) / 2, symmetric = TRUE)
    keep <- e$values > 1e-10
    loadings <- e$vectors[, keep, drop = FALSE] %*%
        diag(sqrt(e$values[keep]), sum(keep))
    Z <- matrix(rnorm(N_SIM * sum(keep)), N_SIM) %*% t(loadings)
    cor(Z)
}

set.seed(2026)
popLd <- map(idxList, function(idx) unname(computeLd(panel, snpIdx = idx)))
ldList <- map2(idxList, popLd, function(idx, R) {
    LdMatrix(
        simulatedLd(R),
        variants = as(gr[idx], "GRanges"),
        nSamples = N_SIM
    )
})
message(
    "  ",
    sum(map_int(ldList, length)),
    " variants across ",
    length(ldList),
    " blocks, reference n = ",
    N_SIM
)

ldEigenExample <- buildLdEigen(ldList, genome = "hg38")
ldScoreExample <- buildLdScore(ldList, genome = "hg38")

stopifnot(
    length(ldEigenExample) == length(ldScoreExample),
    length(eigenList(ldEigenExample)) == length(ldList),
    length(ldMatrixList(ldScoreExample)) == length(ldList),
    identical(names(ldEigenExample), names(ldScoreExample)),
    ldPanelNSamples(ldScoreExample) == N_SIM
)

# The two routes to an LD score must agree, or one of them has drifted.
stopifnot(all.equal(
    as.vector(mcols(ldScoreExample)$ldScores[, 1]),
    as.vector(computeLdScores(ldEigenExample)[, 1])
))

# -----------------------------------------------------------------------------
# Confirm every estimator recovers a known h2 before saving. The GWAS is drawn
# from the POPULATION LD, z ~ N(0, N R diag(v) R + R), and estimated against
# the simulated reference -- a finite sample of that population, as a real
# reference is.
# -----------------------------------------------------------------------------
M <- length(ldEigenExample)
N <- 1e5
h2True <- 0.4
perSnpVar <- rep(h2True / M, M)
snpIdx <- map(eigenList(ldEigenExample), "snpIdx")

simZ <- function() {
    z <- numeric(M)
    for (b in seq_along(popLd)) {
        R <- popLd[[b]]
        idx <- snpIdx[[b]]
        p <- length(idx)
        sigma <- N * (R %*% diag(perSnpVar[idx], p) %*% R) + R
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
        sldsc = pecotmr:::variantRegressionUnivariate(z, N, ldScoreExample)$h2,
        gldsc = pecotmr:::variantGlsUnivariate(z, N, ldScoreExample)$h2,
        lder = pecotmr:::eigenRegressionUnivariate(z, N, ldEigenExample)$h2,
        hdl = suppressWarnings(
            pecotmr:::eigenLikelihoodUnivariate(z, N, ldEigenExample)$h2
        )
    )
})
message("h2 recovery (true 0.4):")
print(round(rowMeans(recovered), 3))
stopifnot(all(abs(rowMeans(recovered) - h2True) < 0.15))

usethis::use_data(ldEigenExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ldScoreExample, overwrite = TRUE, compress = "xz")
