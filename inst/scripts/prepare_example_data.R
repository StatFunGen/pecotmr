#!/usr/bin/env Rscript
#
# prepare_example_data.R
#
# Rebuilds the four CLU-locus example objects shipped in data/:
#
#   gwasSumStatsExample    : data.frame of GWAS summary statistics
#   eqtlRegionExample      : list(X, yRes) of individual-level eQTL data
#   qtlFineMappingExample  : QtlFineMappingResult (individual-level SuSiE)
#   gwasFineMappingExample : GwasFineMappingResult (SuSiE-RSS)
#
# The script recovers the original CLU locus data from git history and
# de-identifies it: synthetic sample names, positions shifted onto chr22, and
# no gene or locus identifiers.
#
# Both fine-mapping fixtures are fit over the SAME 2828 variants. That is what
# lets intersectVariants() and the enrichment examples pair them, so the two
# must always be rebuilt together.
#
# The per-variant tables are produced by postprocessFinemappingFits(), the same
# entry point fineMappingPipeline() uses, rather than being assembled by hand.
# Earlier versions of this script hand-built a `cs` column of bare integer
# labels; that predates the cs_95 / cs_70 / cs_50 (+ _purity) schema, and a
# fixture carrying it makes getCs() and getCredibleSetSummary() silently return
# nothing. Going through the real post-processor keeps the fixtures on whatever
# schema the package currently emits.
#
# Prerequisites
#   1. Recover the original files from git history:
#        git show b25ccd7:inst/prototype/CLU_gwas.rds > /tmp/CLU_gwas.rds
#        git show b25ccd7:inst/prototype/pseudo_bulk_CLU.ENSG00000120885.rds \
#          > /tmp/pseudo_bulk_CLU.rds
#      Override the directory with PECOTMR_EXAMPLE_INPUT_DIR if not /tmp.
#   2. susieR must be installed.
#
# Usage
#   pixi run --environment r45 Rscript inst/scripts/prepare_example_data.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
library(susieR)

inputDir <- Sys.getenv("PECOTMR_EXAMPLE_INPUT_DIR", unset = "/tmp")
dataDir <- "data"

# De-identification: every variant is moved onto chr22 by a fixed offset, so
# relative spacing (and therefore LD structure) is preserved exactly.
newChrom <- "chr22"
posOffset <- 5000000L

# The one credible set at this eQTL locus spans 42 variants with
# min.abs.corr = 0.626, so the pipeline default (minAbsCorr = 0.8) would reject
# it and label every variant "susie_0" -- leaving getCs() on the package's
# headline QTL example with nothing to show. Building at 0.5 keeps the set
# labelled while cs_95_purity still records the true 0.626, so a caller who
# asks for getCs(x, minPurity = 0.8) correctly gets nothing back. The GWAS fit
# below needs no such allowance: its credible sets are all above 0.95.
qtlMinAbsCorr <- 0.5

cat("=== Loading recovered data ===\n")
gwasPath <- file.path(inputDir, "CLU_gwas.rds")
eqtlPath <- file.path(inputDir, "pseudo_bulk_CLU.rds")
for (p in c(gwasPath, eqtlPath)) {
    if (!file.exists(p)) {
        stop(
            "missing input '",
            p,
            "'; see the Prerequisites block at the top of this script"
        )
    }
}
gwasOrig <- readRDS(gwasPath)
eqtlOrig <- readRDS(eqtlPath)

cat(sprintf("GWAS: %d variants\n", nrow(gwasOrig)))
cat(sprintf(
    "eQTL: %d samples x %d variants\n",
    nrow(eqtlOrig$X),
    ncol(eqtlOrig$X)
))

# ---------------------------------------------------------------------------
# 1. De-identify variants and samples
# ---------------------------------------------------------------------------
cat("\n=== De-identifying variants and samples ===\n")
newPos <- gwasOrig$POS + posOffset
newVariantIds <- sprintf(
    "%s:%d:%s:%s",
    newChrom,
    newPos,
    gwasOrig$A1,
    gwasOrig$A2
)
newSampleNames <- sprintf("sample_%03d", seq_len(nrow(eqtlOrig$X)))

cat(sprintf(
    "Original position range: %d - %d\n",
    min(gwasOrig$POS),
    max(gwasOrig$POS)
))
cat(sprintf("New position range: %d - %d\n", min(newPos), max(newPos)))

# ---------------------------------------------------------------------------
# 2. gwasSumStatsExample
# ---------------------------------------------------------------------------
cat("\n=== Creating gwasSumStatsExample ===\n")
gwasSumStatsExample <- data.frame(
    variant_id = newVariantIds,
    chrom = newChrom,
    pos = newPos,
    A1 = gwasOrig$A1,
    A2 = gwasOrig$A2,
    beta = gwasOrig$beta,
    se = gwasOrig$standard_error,
    z = gwasOrig$Z,
    stringsAsFactors = FALSE
)
cat(sprintf(
    "  %d variants, %d columns\n",
    nrow(gwasSumStatsExample),
    ncol(gwasSumStatsExample)
))

# ---------------------------------------------------------------------------
# 3. eqtlRegionExample
# ---------------------------------------------------------------------------
cat("\n=== Creating eqtlRegionExample ===\n")
# The recovered fixture names the residualized phenotype `y_res`; accept the
# camelCase spelling too so a re-exported input still works. `[[` rather than
# `$`, which would prefix-match one spelling onto the other.
yRes <- if (!is.null(eqtlOrig[["y_res"]])) {
    eqtlOrig[["y_res"]]
} else {
    eqtlOrig[["yRes"]]
}
if (is.null(yRes)) {
    stop("recovered eQTL input has neither `y_res` nor `yRes`")
}

X <- eqtlOrig$X
colnames(X) <- newVariantIds
rownames(X) <- newSampleNames
names(yRes) <- newSampleNames

eqtlRegionExample <- list(X = X, yRes = yRes)
cat(sprintf(
    "  X: %d x %d, yRes: %d\n",
    nrow(eqtlRegionExample$X),
    ncol(eqtlRegionExample$X),
    length(eqtlRegionExample$yRes)
))

# ---------------------------------------------------------------------------
# 4. qtlFineMappingExample: individual-level SuSiE
# ---------------------------------------------------------------------------
cat("\n=== Running SuSiE for QTL fine-mapping ===\n")
set.seed(42)
qtlFit <- susie(
    X = X,
    y = yRes,
    L = 10,
    max_iter = 500,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    verbose = FALSE
)
if (!isTRUE(qtlFit$converged)) {
    stop("QTL SuSiE did not converge; refusing to ship a non-converged fit")
}
cat(sprintf(
    "  SuSiE converged, %d credible sets (min.abs.corr %s)\n",
    length(qtlFit$sets$cs),
    paste(sprintf("%.3f", qtlFit$sets$purity$min.abs.corr), collapse = ", ")
))

qtlPost <- postprocessFinemappingFits(
    fits = list(susie = qtlFit),
    dataX = X,
    dataY = yRes,
    minAbsCorr = qtlMinAbsCorr
)
qtlRow <- qtlPost$finemappingResults$susie$finemappingEntry

qtlFineMappingExample <- QtlFineMappingResult(
    study = "study_1",
    context = "context_1",
    trait = "gene_1",
    method = "susie",
    entry = list(qtlRow),
    ldSketch = NULL
)

# ---------------------------------------------------------------------------
# 5. gwasFineMappingExample: SuSiE-RSS over the same variants
# ---------------------------------------------------------------------------
cat("\n=== Running SuSiE-RSS for GWAS fine-mapping ===\n")
R <- cor(X)
zGwas <- gwasSumStatsExample$z
names(zGwas) <- newVariantIds

# A realistic GWAS sample size; the exact value is not critical for the example
# but susie_rss needs it to calibrate effect sizes.
gwasN <- 400000L
set.seed(42)
gwasFit <- susie_rss(
    z = zGwas,
    R = R,
    n = gwasN,
    L = 10,
    max_iter = 500,
    estimate_prior_variance = TRUE,
    verbose = FALSE
)
# Unlike the QTL fit, this one is not held to susie_rss's ELBO convergence
# flag. The GWAS z-scores come from a large external study while R is computed
# from these 415 eQTL samples, so the RSS likelihood is misspecified and the
# ELBO does not reach tolerance -- n = 415 fails outright with "estimated prior
# variance is unreasonably large", and omitting n does too. What the fixture is
# for does converge: the same four credible sets, with the same purities and a
# max PIP of 1.0, come back for every n in {5e3, 5e4, 4e5} and every max_iter
# in {500, 2000}. So assert on the credible sets, and report the ELBO flag
# rather than hiding it.
if (length(gwasFit$sets$cs) == 0L) {
    stop("GWAS SuSiE-RSS found no credible set; the fixture would be empty")
}
cat(sprintf(
    "  SuSiE-RSS: %d credible sets (min.abs.corr %s), ELBO converged=%s\n",
    length(gwasFit$sets$cs),
    paste(sprintf("%.3f", gwasFit$sets$purity$min.abs.corr), collapse = ", "),
    isTRUE(gwasFit$converged)
))

# Named `susieRss` so the post-processor dispatches on the RSS method: it takes
# `dataX` as a correlation matrix (csInput = "Xcorr") rather than genotypes.
gwasPost <- postprocessFinemappingFits(
    fits = list(susieRss = gwasFit),
    dataX = R
)
gwasRow <- gwasPost$finemappingResults$susieRss$finemappingEntry

# qtlEnrichmentPipeline() requires the GWAS collection to carry a non-NULL
# ldSketch: that is its "this fit is RSS-derived" marker. It checks the
# handle's presence, and that it agrees with the QTL side when that one is also
# non-NULL (here it is NULL, so the check is skipped) -- the panel itself is
# never read. No genotype file for this CLU region is shipped, and the LD
# actually used above is the in-memory cor(X) from the eQTL genotypes, so the
# bundled toy_canonical reference stands in, as it did in the previously
# shipped fixture. Note it spans chr22:14.5-15.8Mb and so shares no variant
# with this chr22:32.1-33.1Mb fit; it satisfies the marker, nothing more.
ldStem <- file.path("inst", "extdata", "toy_canonical")
gwasLdSketch <- GenotypeHandle(plink1Prefix = ldStem)
# Portable bundled-resource path, so the .rda carries no source-tree path.
gwasLdSketch@path <- paste0(
    "pecotmr://extdata/",
    basename(gwasLdSketch@path)
)

gwasFineMappingExample <- GwasFineMappingResult(
    study = "study_1",
    method = "susie",
    entry = list(gwasRow),
    blockId = "region_1",
    ldSketch = gwasLdSketch
)

# ---------------------------------------------------------------------------
# 6. Save
# ---------------------------------------------------------------------------
cat("\n=== Saving data objects ===\n")
save(
    gwasSumStatsExample,
    file = file.path(dataDir, "gwasSumStatsExample.rda"),
    compress = "xz"
)
save(
    eqtlRegionExample,
    file = file.path(dataDir, "eqtlRegionExample.rda"),
    compress = "xz"
)
save(
    gwasFineMappingExample,
    file = file.path(dataDir, "gwasFineMappingExample.rda"),
    compress = "xz"
)
save(
    qtlFineMappingExample,
    file = file.path(dataDir, "qtlFineMappingExample.rda"),
    compress = "xz"
)

for (nm in c(
    "gwasSumStatsExample",
    "eqtlRegionExample",
    "gwasFineMappingExample",
    "qtlFineMappingExample"
)) {
    f <- file.path(dataDir, paste0(nm, ".rda"))
    cat(sprintf("  %s: %s bytes\n", basename(f), format(file.size(f))))
}

# ---------------------------------------------------------------------------
# 7. Verify what the fixtures are actually for
# ---------------------------------------------------------------------------
cat("\n=== Verifying ===\n")
for (nm in c("qtlFineMappingExample", "gwasFineMappingExample")) {
    obj <- get(nm)
    cs <- getCs(obj)
    summ <- getCredibleSetSummary(obj)
    cat(sprintf(
        "  %-24s %d variants, getCs() %d rows, %d credible sets\n",
        nm,
        sum(lengths(obj)),
        nrow(cs),
        nrow(summ)
    ))
    if (nrow(cs) == 0L) {
        stop(nm, ": getCs() is empty, which is the defect this build fixes")
    }
}

# intersectVariants() returns the two collections restricted to their shared
# variants (a list of x and y), not the variant vector itself.
both <- intersectVariants(qtlFineMappingExample, gwasFineMappingExample)
nShared <- length(getVariantIds(both$x))
cat(sprintf("  intersectVariants() -> %d shared variants\n", nShared))
if (nShared != sum(lengths(qtlFineMappingExample))) {
    stop(
        "the two fixtures must be fit over the same variant set; shared = ",
        nShared
    )
}

# qtlEnrichmentPipeline() is what these two fixtures are documented for, and
# its non-NULL ldSketch requirement bites only at call time: a GWAS fixture
# built without one still constructs, still validates, and still saves. Run it
# here so that gap fails the build instead of the vignette.
invisible(capture.output(
    enrich <- qtlEnrichmentPipeline(
        gwasFineMappingResult = gwasFineMappingExample,
        qtlFineMappingResult = qtlFineMappingExample,
        verbose = FALSE
    )
))
cat(sprintf(
    "  qtlEnrichmentPipeline() -> %d row(s)\n",
    nrow(as.data.frame(enrich))
))

cat("\nDone! All example data objects created.\n")
