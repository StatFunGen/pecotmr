context("s4Constructors")

# ===========================================================================
# Shared test helpers
# ===========================================================================

.sc_makeGenotypeHandle <- function(snp_n = 5L) {
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
    nSamples = 100L,
    sampleIds = paste0("s", seq_len(100)),
    pgenPtr = NULL)
}

.sc_makeTopLoci <- function(n = 3) {
  data.frame(
    variant_id     = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    chrom          = rep("1", n),
    pos            = as.integer(100 * seq_len(n)),
    A1             = rep("G", n),
    A2             = rep("A", n),
    N              = rep(1000, n),
    MAF            = rep(0.1, n),
    marginal_beta  = rep(0.1, n),
    marginal_se    = rep(0.05, n),
    marginal_z     = rep(2.0, n),
    marginal_p     = rep(0.05, n),
    pip            = seq(0.9, by = -0.1, length.out = n),
    posterior_mean = rep(0.05, n),
    posterior_sd   = rep(0.02, n),
    cs_95          = paste0("susie_", c(1L, 1L, 0L)[seq_len(n)]),
    stringsAsFactors = FALSE)
}

.sc_makeFineMappingEntry <- function(n = 3) {
  FineMappingEntry(
    variantIds = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    susieFit = list(fake = TRUE),
    topLoci    = .sc_makeTopLoci(n))
}

.sc_makeTwasWeightsEntry <- function(p = 5L, standardized = FALSE,
                                     dataType = "expression") {
  TwasWeightsEntry(
    variantIds   = paste0("v", seq_len(p)),
    weights      = rnorm(p),
    standardized = standardized,
    dataType     = dataType)
}

# ===========================================================================
# FineMappingEntry
# ===========================================================================
.sc_makeSe <- function(traits = c("ENSG1", "ENSG2"), n_samples = 10) {
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
