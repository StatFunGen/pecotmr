# Tests for R/MultiStudyQtlDataset.R

# === Tests migrated from test_s4Constructors.R (MultiStudyQtlDataset) ===

test_that("MultiStudyQtlDataset: combines two QtlDatasets", {
  qd1 <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = .sc_makeSe()))
  qd2 <- QtlDataset(study = "s2", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = .sc_makeSe()))
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
  expect_s4_class(mt, "MultiStudyQtlDataset")
  expect_setequal(getStudy(mt), c("s1", "s2"))
})


test_that("MultiStudyQtlDataset: rejects single dataset with no sumStats", {
  qd <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                   phenotypes = list(brain = .sc_makeSe()))
  expect_error(
    MultiStudyQtlDataset(qtlDatasets = list(s1 = qd)),
    "at least 2 studies"
  )
})


test_that("MultiStudyQtlDataset: rejects unnamed qtlDatasets list", {
  qd <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                   phenotypes = list(brain = .sc_makeSe()))
  expect_error(
    MultiStudyQtlDataset(qtlDatasets = list(qd, qd)),
    "named list"
  )
})


test_that("MultiStudyQtlDataset: rejects trait/position conflicts across studies", {
  se1 <- .sc_makeSe(traits = "ENSG1")
  # Build se2 from scratch (see note in the QtlDataset trait-conflict test).
  rng2 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 9999L, width = 500L))
  names(rng2) <- "ENSG1"
  expr2 <- matrix(rnorm(10), nrow = 1, ncol = 10,
                  dimnames = list("ENSG1", paste0("s", 1:10)))
  cd2 <- S4Vectors::DataFrame(sex = rep(c("M", "F"), 5),
                              row.names = paste0("s", 1:10))
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr2),
    rowRanges = rng2, colData = cd2)
  qd1 <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = se1))
  qd2 <- QtlDataset(study = "s2", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = se2))
  expect_error(
    MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2)),
    "inconsistent rowRanges"
  )
})


test_that("getSumStats(MultiStudyQtlDataset) rejects selection arguments", {
  # Compose one individual-level QtlDataset with a QtlSumStats of
  # summary-statistic-only studies (1 + 1 = 2 studies total).
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = "rs1", A1 = "A", A2 = "G", Z = 1.0, N = 1000L)
  ss <- QtlSumStats(study = "s3", context = "c1", trait = "t1",
                    entry = list(gr), genome = "hg19",
                    ldSketch = .sc_makeGenotypeHandle())
  qd1 <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = .sc_makeSe()))
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1), sumStats = ss)

  # Bare call returns the embedded QtlSumStats collection ...
  expect_s4_class(getSumStats(mt), "QtlSumStats")
  # ... but any selection argument is rejected.
  expect_error(getSumStats(mt, study = "s1"), "does not accept selection")
})

