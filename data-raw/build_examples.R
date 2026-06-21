# =============================================================================
# Build canonical S4 example objects for vignettes.
#
# Sources `inst/extdata/toy_ref.{bed,bim,fam}` (165 samples x 17,421 chr22
# variants) for a real GenotypeHandle backing, then constructs synthetic
# phenotypes / sumstats so the vignettes can demonstrate every S4 input
# class with a single `data(<name>)` call.
#
# Outputs (under data/):
#   qtl_dataset_example                : QtlDataset (single-context)
#   qtl_sumstats_example               : QtlSumStats (single-context)
#   qtl_sumstats_multicontext_example  : QtlSumStats (3 contexts, mash demo)
#   gwas_sumstats_s4_example           : GwasSumStats
#   multi_study_qtl_dataset_example    : MultiStudyQtlDataset
#
# Re-run:  pixi run Rscript data-raw/build_examples.R
# =============================================================================

devtools::load_all(".")

set.seed(20260620L)

# -----------------------------------------------------------------------------
# 1. Build a GenotypeHandle from the bundled toy PLINK1 reference and pick a
#    small contiguous window for fast vignette builds.
# -----------------------------------------------------------------------------
toyBed <- system.file("extdata", "toy_ref.bed", package = "pecotmr")
toyStem <- sub("\\.bed$", "", toyBed)
toyRef <- GenotypeHandle(plink1Prefix = toyStem)

# Pick the first 200 variants in a contiguous chr22 window for a small
# but realistic LD block.
snpIdx <- seq_len(200L)
windowSnpInfo <- toyRef@snpInfo[snpIdx, ]
region <- sprintf("chr22:%d-%d",
                  min(windowSnpInfo$BP), max(windowSnpInfo$BP))

# Build a small canonical-ID PLINK1 fileset in inst/extdata so the
# bundled S4 objects can name variants in canonical `chr:pos:A2:A1`
# form (what buildTopLoci -> parseVariantId expects). The .bed file
# is byte-packed by position so we copy the original; the .bim gets
# a fresh ID column.
canonicalIds <- formatVariantId(
  chrom = windowSnpInfo$CHR, pos = windowSnpInfo$BP,
  A2    = windowSnpInfo$A2,  A1  = windowSnpInfo$A1)

smallStem <- file.path("inst", "extdata", "toy_canonical")
# Subset the .bed to the first 200 variants via snpStats round-trip.
sm <- snpStats::read.plink(paste0(toyStem, ".bed"),
                            paste0(toyStem, ".bim"),
                            paste0(toyStem, ".fam"),
                            select.snps = snpIdx)
# Apply canonical IDs and write the small canonical-ID PLINK fileset.
colnames(sm$genotypes) <- canonicalIds
snpStats::write.plink(
  file.base    = smallStem,
  snps         = sm$genotypes,
  pedigree     = sm$fam$pedigree,
  id           = sm$fam$member,
  father       = sm$fam$father,
  mother       = sm$fam$mother,
  sex          = sm$fam$sex,
  phenotype    = sm$fam$affected,
  chromosome   = sm$map$chromosome,
  genetic.distance = sm$map$cM,
  position     = sm$map$position,
  allele.1     = sm$map$allele.1,
  allele.2     = sm$map$allele.2,
  snp.major    = TRUE,
  na.code      = 0L)
gh <- GenotypeHandle(plink1Prefix = smallStem)
windowSnpInfo <- gh@snpInfo

# Materialise the dosage matrix once for synthetic-phenotype generation.
block <- extractBlockGenotypes(gh, snpIdx, meanImpute = TRUE)
X <- t(SummarizedExperiment::assay(block, "dosage"))   # samples x variants

# -----------------------------------------------------------------------------
# 2. Synthesise a phenotype with two causal variants.
# -----------------------------------------------------------------------------
nSample <- nrow(X)
sampleIds <- gh@sampleIds
causalIdx <- c(50L, 130L)
beta <- numeric(ncol(X)); beta[causalIdx] <- c(0.7, -0.5)
Y <- as.numeric(X %*% beta + stats::rnorm(nSample, sd = 1))
Yref <- Y  # store for later GWAS sumstats reuse

# Wrap in a SummarizedExperiment indexed by a single trait "ENSG_example".
phenoMat <- matrix(Y, nrow = 1L, dimnames = list("ENSG_example", sampleIds))
phenoRng <- GenomicRanges::GRanges(
  seqnames = "chr22",
  ranges = IRanges::IRanges(start = stats::median(windowSnpInfo$BP),
                            width = 500L))
names(phenoRng) <- "ENSG_example"
phenoSe <- SummarizedExperiment::SummarizedExperiment(
  assays    = list(expression = phenoMat),
  rowRanges = phenoRng,
  colData   = S4Vectors::DataFrame(
    sex = rep(c(0, 1), length.out = nSample),
    age = stats::runif(nSample, 20, 80),
    row.names = sampleIds))

qtl_dataset_example <- QtlDataset(
  study              = "study1",
  genotypes          = gh,
  phenotypes         = list(brain = phenoSe),
  genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))

# -----------------------------------------------------------------------------
# 3. Compute summary statistics from (X, Y) -> QtlSumStats.
# -----------------------------------------------------------------------------
nVar <- ncol(X)
maf <- pmin(colMeans(X) / 2, 1 - colMeans(X) / 2)
ssZ <- numeric(nVar); ssBeta <- numeric(nVar); ssSe <- numeric(nVar)
for (j in seq_len(nVar)) {
  fit <- summary(stats::lm(Y ~ X[, j]))$coefficients
  ssBeta[j] <- fit[2L, 1L]
  ssSe[j]   <- fit[2L, 2L]
  ssZ[j]    <- fit[2L, 3L]
}
qtlGr <- GenomicRanges::GRanges(
  seqnames = "chr22",
  ranges   = IRanges::IRanges(start = windowSnpInfo$BP, width = 1L))
S4Vectors::mcols(qtlGr) <- S4Vectors::DataFrame(
  SNP  = windowSnpInfo$SNP,
  A1   = windowSnpInfo$A1,
  A2   = windowSnpInfo$A2,
  Z    = ssZ,
  N    = rep(nSample, nVar),
  BETA = ssBeta,
  SE   = ssSe,
  MAF  = maf)

qtl_sumstats_example <- QtlSumStats(
  study   = "study1",
  context = "brain",
  trait   = "ENSG_example",
  entry   = list(qtlGr),
  genome  = "hg19",
  ldSketch = gh,
  qcInfo  = list(prebuilt = "synthetic example data; QC bypassed"))

# -----------------------------------------------------------------------------
# 4. Synthesise GWAS sumstats with one shared causal (50) + one GWAS-only
#    causal (75) so coloc demos have a co-localising signal.
# -----------------------------------------------------------------------------
nGwas <- 50000L
betaGwas <- numeric(nVar); betaGwas[c(50L, 75L)] <- c(0.4, 0.3)
ssZg <- (betaGwas + stats::rnorm(nVar, sd = 1 / sqrt(nGwas))) * sqrt(nGwas)
ssBetaG <- ssZg / sqrt(nGwas); ssSeG <- rep(1 / sqrt(nGwas), nVar)

gwasGr <- GenomicRanges::GRanges(
  seqnames = "chr22",
  ranges   = IRanges::IRanges(start = windowSnpInfo$BP, width = 1L))
S4Vectors::mcols(gwasGr) <- S4Vectors::DataFrame(
  SNP  = windowSnpInfo$SNP,
  A1   = windowSnpInfo$A1,
  A2   = windowSnpInfo$A2,
  Z    = ssZg,
  N    = rep(nGwas, nVar),
  BETA = ssBetaG,
  SE   = ssSeG,
  MAF  = maf)

gwas_sumstats_s4_example <- GwasSumStats(
  study    = "trait1",
  entry    = list(gwasGr),
  genome   = "hg19",
  ldSketch = gh,
  qcInfo   = list(prebuilt = "synthetic example data; QC bypassed"))

# -----------------------------------------------------------------------------
# 5. MultiStudyQtlDataset: a second synthetic QtlDataset on the same genotype
#    handle (different causal variants and a noisier signal).
# -----------------------------------------------------------------------------
beta2 <- numeric(nVar); beta2[c(110L, 180L)] <- c(0.5, 0.6)
Y2 <- as.numeric(X %*% beta2 + stats::rnorm(nSample, sd = 1))
phenoMat2 <- matrix(Y2, nrow = 1L, dimnames = list("ENSG_example", sampleIds))
phenoSe2 <- SummarizedExperiment::SummarizedExperiment(
  assays    = list(expression = phenoMat2),
  rowRanges = phenoRng,
  colData   = S4Vectors::DataFrame(row.names = sampleIds))

qd2 <- QtlDataset(
  study              = "study2",
  genotypes          = gh,
  phenotypes         = list(brain = phenoSe2),
  genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))

multi_study_qtl_dataset_example <- MultiStudyQtlDataset(
  qtlDatasets = list(study1 = qtl_dataset_example, study2 = qd2))

# -----------------------------------------------------------------------------
# 6. Multi-context QtlSumStats for mash demos.
#    Same toy genotype panel, same trait `ENSG_example`, three synthetic
#    contexts (brain / blood / muscle) wired with:
#      - a shared causal variant (index 50, present in all 3 contexts)
#      - a brain-specific causal (index 130)
#      - a blood-specific causal (index 75)
#      - muscle has only the shared signal
#    This is the canonical fixture for exercising the mash pipeline's
#    shared / context-unique pattern recovery.
# -----------------------------------------------------------------------------
multiCtxNames <- c("brain", "blood", "muscle")
sharedIdx     <- 50L
brainOnlyIdx  <- 130L
bloodOnlyIdx  <- 75L
multiCtxBetas <- list(
  brain  = { b <- numeric(nVar); b[c(sharedIdx, brainOnlyIdx)] <- c(0.6, -0.4); b },
  blood  = { b <- numeric(nVar); b[c(sharedIdx, bloodOnlyIdx)] <- c(0.6,  0.5); b },
  muscle = { b <- numeric(nVar); b[sharedIdx] <- 0.6; b })

multiCtxEntries <- lapply(multiCtxNames, function(ctx) {
  bj <- multiCtxBetas[[ctx]]
  yCtx <- as.numeric(X %*% bj + stats::rnorm(nSample, sd = 1))
  zc <- numeric(nVar); bc <- numeric(nVar); sc <- numeric(nVar)
  for (j in seq_len(nVar)) {
    fit <- summary(stats::lm(yCtx ~ X[, j]))$coefficients
    bc[j] <- fit[2L, 1L]
    sc[j] <- fit[2L, 2L]
    zc[j] <- fit[2L, 3L]
  }
  gr <- GenomicRanges::GRanges(
    seqnames = "chr22",
    ranges   = IRanges::IRanges(start = windowSnpInfo$BP, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP  = windowSnpInfo$SNP,
    A1   = windowSnpInfo$A1,
    A2   = windowSnpInfo$A2,
    Z    = zc,
    N    = rep(nSample, nVar),
    BETA = bc,
    SE   = sc,
    MAF  = maf)
  gr
})

qtl_sumstats_multicontext_example <- QtlSumStats(
  study    = rep("study1", length(multiCtxNames)),
  context  = multiCtxNames,
  trait    = rep("ENSG_example", length(multiCtxNames)),
  entry    = multiCtxEntries,
  genome   = "hg19",
  ldSketch = gh,
  qcInfo   = list(prebuilt = "synthetic multi-context example; QC bypassed"))

# -----------------------------------------------------------------------------
# 7. Strip the GenotypeHandle paths down to bare basenames so the bundled
#    .rda objects don't carry source-tree-absolute paths. Vignettes /
#    users resolve them at use-time via fixupExampleGenotypePaths().
# -----------------------------------------------------------------------------
qtl_dataset_example@genotypes@path <- basename(
  qtl_dataset_example@genotypes@path)
qtl_sumstats_example@ldSketch@path <- basename(
  qtl_sumstats_example@ldSketch@path)
qtl_sumstats_multicontext_example@ldSketch@path <- basename(
  qtl_sumstats_multicontext_example@ldSketch@path)
gwas_sumstats_s4_example@ldSketch@path <- basename(
  gwas_sumstats_s4_example@ldSketch@path)
for (nm in names(multi_study_qtl_dataset_example@qtlDatasets))
  multi_study_qtl_dataset_example@qtlDatasets[[nm]]@genotypes@path <-
    basename(multi_study_qtl_dataset_example@qtlDatasets[[nm]]@genotypes@path)

# -----------------------------------------------------------------------------
# 8. Save.
# -----------------------------------------------------------------------------
usethis::use_data(qtl_dataset_example, overwrite = TRUE, compress = "xz")
usethis::use_data(qtl_sumstats_example, overwrite = TRUE, compress = "xz")
usethis::use_data(qtl_sumstats_multicontext_example, overwrite = TRUE,
                  compress = "xz")
usethis::use_data(gwas_sumstats_s4_example, overwrite = TRUE, compress = "xz")
usethis::use_data(multi_study_qtl_dataset_example, overwrite = TRUE,
                  compress = "xz")

cat("\nBuilt:\n",
    "  data/qtl_dataset_example.rda\n",
    "  data/qtl_sumstats_example.rda\n",
    "  data/qtl_sumstats_multicontext_example.rda\n",
    "  data/gwas_sumstats_s4_example.rda\n",
    "  data/multi_study_qtl_dataset_example.rda\n", sep = "")
