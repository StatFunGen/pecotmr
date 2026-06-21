context("show methods")

# ===========================================================================
# Shared fixtures
# ===========================================================================

.sh_makeGenotypeHandle <- function(snp_n = 3L) {
  new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("rs", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = 50L,
    sampleIds = paste0("s", seq_len(50)),
    pgenPtr = NULL)
}

.sh_makeFmEntry <- function(n = 3, with_cs = TRUE) {
  tl <- data.frame(
    variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    pip        = seq(0.9, by = -0.1, length.out = n),
    stringsAsFactors = FALSE)
  if (with_cs) tl$cs <- c(1L, 1L, 0L)[seq_len(n)]
  FineMappingEntry(
    variantIds = tl$variant_id,
    trimmedFit = list(),
    topLoci    = tl)
}

.sh_makeTwEntry <- function(p = 4, standardized = FALSE) {
  TwasWeightsEntry(
    variantIds   = paste0("v", seq_len(p)),
    weights      = rep(0.1, p),
    cvPerformance = list(rsq = 0.5),
    standardized = standardized)
}

.sh_makeSe <- function(traits = c("ENSG1", "ENSG2"), n_samples = 6) {
  rng <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(traits)),
    ranges = IRanges::IRanges(
      start = seq(1000L, by = 1000L, length.out = length(traits)),
      width = 500L))
  names(rng) <- traits
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  cd <- S4Vectors::DataFrame(sex = rep(c("M", "F"), length.out = n_samples),
                             row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr),
    rowRanges = rng,
    colData = cd)
}

.sh_makeQtlDataset <- function(study = "study1") {
  QtlDataset(
    study              = study,
    genotypes          = .sh_makeGenotypeHandle(),
    phenotypes         = list(brain = .sh_makeSe()),
    genotypeCovariates = matrix(0, nrow = 50, ncol = 0))
}

.sh_makeQtlSumstatsGr <- function(n = 3) {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = seq(100L, by = 100L, length.out = n), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("rs", seq_len(n)),
    A1 = rep("A", n), A2 = rep("G", n),
    Z = rnorm(n), N = rep(100L, n))
  gr
}

