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

