#!/usr/bin/env Rscript
#
# prepare_paired_finemapping_example.R
#
# Build `qtlFineMappingPairedExample`: a multi-credible-set
# `QtlFineMappingResult` whose fit is over the SAME variants as the shipped LD
# sources `qtlDatasetExample` / `qtlSumStatsExample` (chr22:14.5Mb, 200 SNPs).
#
# The other fine-mapping fixtures (`qtlFineMappingExample`,
# `gwasFineMappingExample`) are fit over a different region (chr22:32.1Mb) than
# any shipped LD source, so cross-object, LD-derived computations such as
# `computeCsCorrelation()` cannot use them. This fixture closes that gap: it is
# self-consistent with `qtlSumStatsExample`, so
# `computeCsCorrelation(getFineMappingResult(qtlFineMappingPairedExample),
# qtlSumStatsExample)` returns a real between-credible-set correlation matrix.
#
# A phenotype with two causal variants in low mutual LD is simulated so the fit
# yields at least two credible sets (a single CS has no between-CS correlation).
#
# Usage:
#   pixi run --environment r45 Rscript \
#     inst/scripts/prepare_paired_finemapping_example.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
library(susieR)

data(qtlDatasetExample)

# Genotypes for the 200 shipped variants (columns match the LD-source panels).
X <- getGenotypes(qtlDatasetExample)
variantIds <- colnames(X)

# Two causal variants far apart and in low mutual LD -> two credible sets.
R <- cor(X)
causalOne <- 20L
lowLd <- which(abs(R[causalOne, ]) < 0.1)
causalTwo <- lowLd[which.min(abs(lowLd - 150L))]

set.seed(42)
Xs <- scale(X)
y <- 3 * Xs[, causalOne] + 3 * Xs[, causalTwo] + rnorm(nrow(X))

fit <- susie(
    X = X,
    y = y,
    L = 10,
    max_iter = 500,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    verbose = FALSE
)
cat(sprintf("SuSiE converged, %d credible sets\n", length(fit$sets$cs)))
if (length(fit$sets$cs) < 2L) {
    stop("expected at least two credible sets for a paired-correlation example")
}

# Trim the fit to the fields the package uses (matches the other fixtures).
trimmed <- list(
    alpha = fit$alpha,
    pip = fit$pip,
    V = fit$V,
    sets = fit$sets
)

# Per-variant top-loci table (one row per fit variant, PIP aligned to the fit).
topLoci <- data.frame(
    variant_id = variantIds,
    pip = fit$pip,
    stringsAsFactors = FALSE
)

entry <- FineMappingEntry(
    variantIds = variantIds,
    susieFit = trimmed,
    topLoci = topLoci
)

qtlFineMappingPairedExample <- QtlFineMappingResult(
    study = "study1",
    context = "context1",
    trait = "gene1",
    method = "susie",
    entry = list(entry),
    ldSketch = NULL
)

save(
    qtlFineMappingPairedExample,
    file = file.path("data", "qtlFineMappingPairedExample.rda"),
    compress = "xz"
)

# Sanity check: the paired correlation resolves against qtlSumStatsExample.
data(qtlSumStatsExample)
cc <- computeCsCorrelation(
    getFineMappingResult(qtlFineMappingPairedExample),
    qtlSumStatsExample
)
cat(sprintf(
    "computeCsCorrelation() -> %d x %d matrix\n",
    nrow(cc),
    ncol(cc)
))
cat("Done.\n")
