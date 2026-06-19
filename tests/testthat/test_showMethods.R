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

# ===========================================================================
# show() outputs
# ===========================================================================

test_that("show.QtlFineMappingResult prints entry/study/context/trait/method counts", {
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susieRss"),
    entry   = list(.sh_makeFmEntry(), .sh_makeFmEntry()))
  out <- capture.output(show(res))
  expect_true(any(grepl("QtlFineMappingResult: 2 entries", out)))
  expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
  expect_true(any(grepl("LD sketch: NULL", out)))
})

test_that("show.QtlFineMappingResult reports the ldSketch source when present", {
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(.sh_makeFmEntry()),
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(res))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})

test_that("show.GwasFineMappingResult prints (study, method) summary", {
  res <- GwasFineMappingResult(
    study  = c("g1", "g1"),
    method = c("susie", "susieRss"),
    entry  = list(.sh_makeFmEntry(), .sh_makeFmEntry()))
  out <- capture.output(show(res))
  expect_true(any(grepl("GwasFineMappingResult: 2 entries", out)))
  expect_true(any(grepl("1 studies.*2 methods", out)))
  expect_true(any(grepl("LD sketch: NULL", out)))
})

test_that("show.GwasFineMappingResult reports the ldSketch source when present", {
  res <- GwasFineMappingResult(
    study = "g1", method = "susie",
    entry = list(.sh_makeFmEntry()),
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(res))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})

test_that("show.TwasWeights prints entry/study/context/trait/method counts", {
  e <- .sh_makeTwEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "enet"),
    entry   = list(e, e))
  out <- capture.output(show(tw))
  expect_true(any(grepl("TwasWeights: 2 entries", out)))
  expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
})

test_that("show.TwasWeights reports ldSketch when present", {
  e <- .sh_makeTwEntry()
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "lasso",
    entry = list(e),
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(tw))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})

test_that("show.FineMappingEntry reports variant count and CS count", {
  e_with_cs <- .sh_makeFmEntry(n = 3, with_cs = TRUE)  # 2 distinct cs > 0
  out <- capture.output(show(e_with_cs))
  expect_true(any(grepl("FineMappingEntry: 3 variants.*1 credible sets", out)))

  # No cs column -> 0 credible sets reported.
  tl <- data.frame(variant_id = c("a", "b"), pip = c(0.1, 0.2),
                   stringsAsFactors = FALSE)
  e_no_cs <- FineMappingEntry(variantIds = c("a", "b"),
                              trimmedFit = list(), topLoci = tl)
  out_no <- capture.output(show(e_no_cs))
  expect_true(any(grepl("0 credible sets", out_no)))
})

test_that("show.TwasWeightsEntry reports standardized flag and CV availability", {
  e <- .sh_makeTwEntry(p = 5, standardized = TRUE)
  out <- capture.output(show(e))
  expect_true(any(grepl("TwasWeightsEntry: 5 variants.*standardized=TRUE", out)))
  expect_true(any(grepl("CV performance: TRUE", out)))

  e_no_cv <- TwasWeightsEntry(variantIds = c("v1", "v2"),
                               weights = c(0.1, 0.2))
  out2 <- capture.output(show(e_no_cv))
  expect_true(any(grepl("CV performance: FALSE", out2)))
})

test_that("show.QtlDataset lists context names and trait count", {
  qd <- .sh_makeQtlDataset()
  out <- capture.output(show(qd))
  expect_true(any(grepl("QtlDataset for study 'study1'", out)))
  expect_true(any(grepl("1 context\\(s\\): brain", out)))
  expect_true(any(grepl("2 unique traits", out)))
  expect_true(any(grepl("Genotypes: gds @ /tmp/test.gds", out)))
})

test_that("show.MultiTaskQtlDataset reports per-source study counts", {
  qd1 <- .sh_makeQtlDataset(study = "s1")
  qd2 <- .sh_makeQtlDataset(study = "s2")
  mt <- MultiTaskQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
  out <- capture.output(show(mt))
  expect_true(any(grepl("MultiTaskQtlDataset: 2 individual-level \\+ 0 sumstats", out)))
  expect_true(any(grepl("Individual-level studies: s1, s2", out)))
})

test_that("show.MultiTaskQtlDataset reports sumstats studies when present", {
  qd <- .sh_makeQtlDataset(study = "s1")
  ss <- QtlSumStats(
    study   = "s2",
    context = "c1",
    trait   = "t1",
    entry   = list(.sh_makeQtlSumstatsGr()),
    genome  = "hg19",
    ldSketch = .sh_makeGenotypeHandle())
  mt <- MultiTaskQtlDataset(qtlDatasets = list(s1 = qd), sumStats = ss)
  out <- capture.output(show(mt))
  expect_true(any(grepl("Sumstats studies: s2", out)))
})

test_that("show.GwasSumStats prints nrow and genome build", {
  ss <- GwasSumStats(
    study = c("g1", "g2"),
    entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
    genome = "hg19",
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(ss))
  expect_true(any(grepl("GwasSumStats: 2 studies, genome build hg19", out)))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})
