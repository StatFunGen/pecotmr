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
    variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    pip        = seq(0.9, by = -0.1, length.out = n),
    cs         = c(1L, 1L, 0L)[seq_len(n)],
    stringsAsFactors = FALSE)
}

.sc_makeFineMappingEntry <- function(n = 3) {
  FineMappingEntry(
    variantIds = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    trimmedFit = list(fake = TRUE),
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

test_that("FineMappingEntry: constructor stores slots and accessors return them", {
  tl <- .sc_makeTopLoci(3)
  entry <- FineMappingEntry(
    variantIds = c("a", "b", "c"),
    trimmedFit = list(payload = 1L),
    topLoci    = tl,
    sumstats   = list(z = c(1, 2, 3)))
  expect_s4_class(entry, "FineMappingEntry")
  expect_equal(getVariantIds(entry), c("a", "b", "c"))
  expect_equal(getTrimmedFit(entry), list(payload = 1L))
  expect_equal(getTopLoci(entry), tl)
})

test_that("FineMappingEntry: getPip returns named pip vector keyed by variant_id", {
  entry <- .sc_makeFineMappingEntry(3)
  pip <- getPip(entry)
  expect_equal(length(pip), 3L)
  expect_equal(names(pip),
               paste0("chr1:", 100 * 1:3, ":A:G"))
})

test_that("FineMappingEntry: getPip returns numeric(0) when topLoci is empty", {
  entry <- FineMappingEntry(
    variantIds = character(0),
    trimmedFit = list(),
    topLoci    = data.frame(variant_id = character(0), pip = numeric(0),
                            stringsAsFactors = FALSE))
  expect_equal(getPip(entry), numeric(0))
})

test_that("FineMappingEntry: getCs filters to rows with cs > 0", {
  entry <- .sc_makeFineMappingEntry(3)  # last row has cs = 0
  res <- getCs(entry)
  expect_equal(nrow(res), 2L)
  expect_true(all(res$cs > 0))
})

test_that("FineMappingEntry: validity errors when topLoci is missing required cols", {
  expect_error(
    FineMappingEntry(
      variantIds = "v1",
      trimmedFit = list(),
      topLoci    = data.frame(other = 1, stringsAsFactors = FALSE)),
    "topLoci missing columns"
  )
})

# ===========================================================================
# TwasWeightsEntry
# ===========================================================================

test_that("TwasWeightsEntry: constructor and accessors round-trip", {
  e <- TwasWeightsEntry(
    variantIds    = c("v1", "v2", "v3"),
    weights       = c(0.1, -0.2, 0.05),
    fits          = list(model = "lasso"),
    cvPerformance = list(rsq = 0.4),
    standardized  = TRUE,
    dataType      = "expression")
  expect_s4_class(e, "TwasWeightsEntry")
  expect_equal(getVariantIds(e), c("v1", "v2", "v3"))
  expect_equal(getWeights(e), c(0.1, -0.2, 0.05))
  expect_equal(getFits(e), list(model = "lasso"))
  expect_equal(getCvPerformance(e), list(rsq = 0.4))
  expect_true(getStandardized(e))
  expect_equal(getDataType(e), "expression")
})

test_that("TwasWeightsEntry: standardized is coerced via isTRUE() semantics", {
  # isTRUE() only returns TRUE for a length-1 logical TRUE. Non-TRUE
  # input lands as FALSE (the safe default for the standardized flag).
  e_logical <- TwasWeightsEntry(
    variantIds = "v1", weights = 0.1, standardized = TRUE)
  expect_true(getStandardized(e_logical))

  e_default <- TwasWeightsEntry(
    variantIds = "v1", weights = 0.1, standardized = "yes-please")
  expect_false(getStandardized(e_default))
})

test_that("TwasWeightsEntry: validity rejects matrix weights with wrong nrow", {
  expect_error(
    TwasWeightsEntry(
      variantIds = c("v1", "v2"),
      weights    = matrix(0, nrow = 5, ncol = 1)),
    "nrow\\(weights\\) must equal length\\(variantIds\\)"
  )
})

# ===========================================================================
# QtlFineMappingResult
# ===========================================================================

test_that("QtlFineMappingResult: builds a collection keyed by 4-tuple", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susie"),
    entry   = list(e1, e2))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)
  expect_null(res@ldSketch)
})

test_that("QtlFineMappingResult: stores an LD sketch when supplied", {
  e <- .sc_makeFineMappingEntry(3)
  gh <- .sc_makeGenotypeHandle()
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e), ldSketch = gh)
  expect_identical(getLdSketch(res), gh)
})

test_that("QtlFineMappingResult: errors on length mismatch", {
  e <- .sc_makeFineMappingEntry(3)
  expect_error(
    QtlFineMappingResult(
      study   = c("s1", "s2"),
      context = c("c1"),
      trait   = c("t1"),
      method  = c("susie"),
      entry   = list(e)),
    "same length"
  )
})

test_that("QtlFineMappingResult: validity rejects duplicate 4-tuples", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  expect_error(
    QtlFineMappingResult(
      study   = c("s1", "s1"),
      context = c("c1", "c1"),
      trait   = c("t1", "t1"),
      method  = c("susie", "susie"),
      entry   = list(e1, e2)),
    "uniqueness violated"
  )
})

test_that("QtlFineMappingResult: validity rejects non-FineMappingEntry rows", {
  expect_error(
    QtlFineMappingResult(
      study = "s1", context = "c1", trait = "t1", method = "susie",
      entry = list("not_an_entry")),
    "every element of the `entry` column must be a FineMappingEntry"
  )
})

test_that("QtlFineMappingResult: getFineMappingResult returns selected entry", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susie"),
    entry   = list(e1, e2))
  picked <- getFineMappingResult(res,
                                 study = "s1", context = "c2",
                                 trait = "t1", method = "susie")
  expect_identical(picked, e2)
})

test_that("QtlFineMappingResult: getFineMappingResult errors on missing tuple", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_error(
    getFineMappingResult(res, study = "ghost", context = "c1",
                         trait = "t1", method = "susie"),
    "No entry for"
  )
})

test_that("QtlFineMappingResult: single-row collection allows omitting selectors", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_identical(getFineMappingResult(res), e)
})

test_that("QtlFineMappingResult: show prints summary", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_output(show(res), "QtlFineMappingResult")
})

# ===========================================================================
# GwasFineMappingResult
# ===========================================================================

test_that("GwasFineMappingResult: builds a collection keyed by 2-tuple", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  res <- GwasFineMappingResult(
    study  = c("g1", "g2"),
    method = c("susie", "susie"),
    entry  = list(e1, e2))
  expect_s4_class(res, "GwasFineMappingResult")
  expect_equal(nrow(res), 2L)
})

test_that("GwasFineMappingResult: errors on length mismatch", {
  e <- .sc_makeFineMappingEntry(3)
  expect_error(
    GwasFineMappingResult(
      study  = c("g1", "g2"),
      method = c("susie"),
      entry  = list(e)),
    "same length"
  )
})

test_that("GwasFineMappingResult: rejects duplicate (study, method) tuples", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  expect_error(
    GwasFineMappingResult(
      study  = c("g1", "g1"),
      method = c("susie", "susie"),
      entry  = list(e1, e2)),
    "uniqueness violated"
  )
})

test_that("GwasFineMappingResult: show prints summary", {
  e <- .sc_makeFineMappingEntry(3)
  res <- GwasFineMappingResult(study = "g1", method = "susie",
                               entry = list(e))
  expect_output(show(res), "GwasFineMappingResult")
})

# ===========================================================================
# TwasWeights collection
# ===========================================================================

test_that("TwasWeights: builds a collection keyed by 4-tuple", {
  e1 <- .sc_makeTwasWeightsEntry()
  e2 <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "enet"),
    entry   = list(e1, e2))
  expect_s4_class(tw, "TwasWeights")
  expect_equal(nrow(tw), 2L)
  expect_setequal(getMethodNames(tw), c("lasso", "enet"))
})

test_that("TwasWeights: getStudy / getContexts / getTraits / getMethodNames", {
  e <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s2"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "lasso"),
    entry   = list(e, e))
  expect_setequal(getContexts(tw), c("c1", "c2"))
  expect_equal(getTraits(tw), "t1")
  expect_equal(getMethodNames(tw), "lasso")
})

test_that("TwasWeights: rejects duplicate 4-tuples", {
  e <- .sc_makeTwasWeightsEntry()
  expect_error(
    TwasWeights(
      study   = c("s1", "s1"),
      context = c("c1", "c1"),
      trait   = c("t1", "t1"),
      method  = c("lasso", "lasso"),
      entry   = list(e, e)),
    "uniqueness violated"
  )
})

test_that("TwasWeights: rejects non-TwasWeightsEntry rows", {
  expect_error(
    TwasWeights(
      study = "s1", context = "c1", trait = "t1", method = "lasso",
      entry = list("not_an_entry")),
    "every element of the `entry` column must be a TwasWeightsEntry"
  )
})

test_that("TwasWeights: getTwasWeights extracts the entry for a tuple", {
  e1 <- .sc_makeTwasWeightsEntry()
  e2 <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "enet"),
    entry   = list(e1, e2))
  expect_identical(
    getTwasWeights(tw, study = "s1", context = "c1",
                   trait = "t1", method = "enet"),
    e2)
})

# ===========================================================================
# LdData
# ===========================================================================

test_that("LdData: pre-computed correlation matrix is returned by getCorrelation", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    A1 = rep("A", 3), A2 = rep("G", 3))
  R <- diag(3)
  block_meta <- S4Vectors::DataFrame(region = "chr1:100-300")
  ld <- LdData(correlation = R, genotypeHandle = NULL, snpIdx = NULL,
               variants = gr, blockMetadata = block_meta, nRef = 100L)
  expect_s4_class(ld, "LdData")
  expect_equal(getCorrelation(ld), R)
})

test_that("LdData: validity rejects both correlation AND genotypeHandle being NULL", {
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 100L, width = 1L))
  expect_error(
    LdData(correlation = NULL, genotypeHandle = NULL,
           variants = gr,
           blockMetadata = S4Vectors::DataFrame(x = 1)),
    "At least one of 'correlation' or 'genotypeHandle' must be non-NULL"
  )
})

test_that("LdData: validity rejects empty variants", {
  expect_error(
    LdData(correlation = diag(0), genotypeHandle = NULL,
           variants = GenomicRanges::GRanges(),
           blockMetadata = S4Vectors::DataFrame(x = 1)),
    "'variants' must not be empty"
  )
})

test_that("LdData: mixtureWeights only valid when genotypeHandle is a list", {
  gh <- .sc_makeGenotypeHandle()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 100L, width = 1L))
  expect_error(
    LdData(correlation = NULL, genotypeHandle = gh,
           variants = gr,
           blockMetadata = S4Vectors::DataFrame(x = 1),
           mixtureWeights = c(0.5, 0.5)),
    "'mixtureWeights' may only be set when 'genotypeHandle' is a list"
  )
})

test_that("LdData: mixtureWeights must be non-negative and sum to 1", {
  gh1 <- .sc_makeGenotypeHandle()
  gh2 <- .sc_makeGenotypeHandle()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 100L, width = 1L))
  expect_error(
    LdData(correlation = NULL,
           genotypeHandle = list(gh1, gh2),
           variants = gr,
           blockMetadata = S4Vectors::DataFrame(x = 1),
           mixtureWeights = c(0.4, 0.4)),
    "must be non-negative and sum to 1"
  )
})

# ===========================================================================
# QtlDataset
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

test_that("QtlDataset: builds and validates with a single-context SE", {
  se <- .sc_makeSe()
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = .sc_makeGenotypeHandle(),
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(0, nrow = 10, ncol = 0))
  expect_s4_class(qd, "QtlDataset")
  expect_equal(getStudy(qd), "study1")
  expect_equal(getContexts(qd), "brain")
})

test_that("QtlDataset: rejects empty study name", {
  se <- .sc_makeSe()
  expect_error(
    QtlDataset(study = "", genotypes = .sc_makeGenotypeHandle(),
               phenotypes = list(brain = se)),
    "non-empty character string"
  )
})

test_that("QtlDataset: rejects empty phenotype list", {
  expect_error(
    QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
               phenotypes = list()),
    "must not be empty"
  )
})

test_that("QtlDataset: rejects unnamed phenotype list", {
  se <- .sc_makeSe()
  expect_error(
    QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
               phenotypes = list(se)),
    "named list"
  )
})

test_that("QtlDataset: rejects non-SE elements in phenotype list", {
  expect_error(
    QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
               phenotypes = list(brain = data.frame(x = 1))),
    "must be a SummarizedExperiment"
  )
})

test_that("QtlDataset: rejects negative QC cutoffs", {
  se <- .sc_makeSe()
  expect_error(
    QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
               phenotypes = list(brain = se),
               mafCutoff = -0.1),
    "non-negative numeric"
  )
})

test_that("QtlDataset: rejects shared traits with inconsistent rowRanges", {
  se1 <- .sc_makeSe(traits = c("ENSG1"))
  # Build se2 from scratch with a different start position for ENSG1 so
  # rownames(se) stays in sync with rowRanges (the validity check skips
  # contexts whose rowRanges length mismatches rownames length).
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
  expect_error(
    QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
               phenotypes = list(brain = se1, liver = se2)),
    "inconsistent rowRanges"
  )
})

# ===========================================================================
# MultiTaskQtlDataset
# ===========================================================================

test_that("MultiTaskQtlDataset: combines two QtlDatasets", {
  qd1 <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = .sc_makeSe()))
  qd2 <- QtlDataset(study = "s2", genotypes = .sc_makeGenotypeHandle(),
                    phenotypes = list(brain = .sc_makeSe()))
  mt <- MultiTaskQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
  expect_s4_class(mt, "MultiTaskQtlDataset")
  expect_setequal(getStudy(mt), c("s1", "s2"))
})

test_that("MultiTaskQtlDataset: rejects single dataset with no sumStats", {
  qd <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                   phenotypes = list(brain = .sc_makeSe()))
  expect_error(
    MultiTaskQtlDataset(qtlDatasets = list(s1 = qd)),
    "at least 2 studies"
  )
})

test_that("MultiTaskQtlDataset: rejects unnamed qtlDatasets list", {
  qd <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                   phenotypes = list(brain = .sc_makeSe()))
  expect_error(
    MultiTaskQtlDataset(qtlDatasets = list(qd, qd)),
    "named list"
  )
})

test_that("MultiTaskQtlDataset: rejects non-QtlDataset entries", {
  qd <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                   phenotypes = list(brain = .sc_makeSe()))
  expect_error(
    MultiTaskQtlDataset(qtlDatasets = list(s1 = qd, s2 = "not a dataset")),
    "must be a QtlDataset"
  )
})

test_that("MultiTaskQtlDataset: rejects trait/position conflicts across studies", {
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
    MultiTaskQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2)),
    "inconsistent rowRanges"
  )
})
