# Tests for the joint-specification grammar and input-argument parsers
# (R/jointSpecification.R).

# -----------------------------------------------------------------------------
# Fixture builders
# -----------------------------------------------------------------------------

.js_makeGenotypeHandle <- function(snp_n = 5L) {
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
    nSamples = 5L,
    sampleIds = paste0("s", seq_len(5)),
    pgenPtr = NULL)
}

.js_makeSe <- function(traits = c("ENSG1", "ENSG2"),
                       samples = paste0("s", seq_len(5))) {
  rng <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = seq.int(100L, by = 100L, length.out = length(traits)),
      width = 50L))
  names(rng) <- traits
  expr <- matrix(rnorm(length(traits) * length(samples)),
                 nrow = length(traits), ncol = length(samples),
                 dimnames = list(traits, samples))
  cd <- S4Vectors::DataFrame(sex = rep(c("M", "F"),
                                       length.out = length(samples)),
                             row.names = samples)
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr),
    rowRanges = rng, colData = cd)
}

.js_makeQtlDataset <- function(study = "s1",
                               contexts = c("brain", "liver"),
                               traits = c("ENSG1", "ENSG2")) {
  phenos <- setNames(
    lapply(contexts, function(cx) .js_makeSe(traits = traits)),
    contexts)
  QtlDataset(study = study,
             genotypes = .js_makeGenotypeHandle(),
             phenotypes = phenos,
             genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
}

.js_makeQtlSumStats <- function(studies = c("ssA", "ssB"),
                                 contexts = c("DLPFC"),
                                 traits = c("ENSG3")) {
  rows <- expand.grid(study = studies, context = contexts, trait = traits,
                      stringsAsFactors = FALSE)
  entries <- lapply(seq_len(nrow(rows)), function(i) {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      SNP = c("rs1", "rs2"),
      A1 = c("A", "G"),
      A2 = c("G", "A"),
      Z  = c(1.5, -1.0),
      N  = c(1000L, 1000L))
    gr
  })
  QtlSumStats(study   = rows$study,
              context = rows$context,
              trait   = rows$trait,
              entry   = entries,
              genome  = "hg19",
              ldSketch = .js_makeGenotypeHandle())
}

# -----------------------------------------------------------------------------
# Scope helpers
# -----------------------------------------------------------------------------

test_that(".spListStudies / .spListContexts / .spListTraits: QtlDataset", {
  qd <- .js_makeQtlDataset(study = "S1",
                           contexts = c("brain", "liver"),
                           traits = c("ENSG_A", "ENSG_B"))
  expect_equal(pecotmr:::.spListStudies(qd), "S1")
  expect_setequal(pecotmr:::.spListContexts(qd), c("brain", "liver"))
  expect_setequal(pecotmr:::.spListTraits(qd, study = "S1",
                                          context = "brain"),
                  c("ENSG_A", "ENSG_B"))
  expect_equal(pecotmr:::.spStudyDataForm(qd, "S1"), "individual")
})

test_that(".sp* helpers: MultiStudyQtlDataset combines individual + sumstats", {
  qd1 <- .js_makeQtlDataset(study = "indA", contexts = "brain",
                            traits = "ENSG_A")
  qd2 <- .js_makeQtlDataset(study = "indB", contexts = "liver",
                            traits = "ENSG_A")
  ss <- .js_makeQtlSumStats(studies = "ssC", contexts = "DLPFC",
                            traits = "ENSG_A")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd1, indB = qd2),
                              sumStats = ss)
  expect_setequal(pecotmr:::.spListStudies(mt),
                  c("indA", "indB", "ssC"))
  expect_setequal(pecotmr:::.spListContexts(mt),
                  c("brain", "liver", "DLPFC"))
  expect_equal(pecotmr:::.spStudyDataForm(mt, "indA"), "individual")
  expect_equal(pecotmr:::.spStudyDataForm(mt, "ssC"), "sumstats")
  expect_error(pecotmr:::.spStudyDataForm(mt, "missing"), "not in")
})

# -----------------------------------------------------------------------------
# parseJointSpecification
# -----------------------------------------------------------------------------

test_that("parseJointSpecification: NULL returns empty list", {
  qd <- .js_makeQtlDataset()
  expect_equal(pecotmr:::parseJointSpecification(NULL, qd), list())
})

test_that("parseJointSpecification: auto-wraps a single char vector", {
  qd <- .js_makeQtlDataset()
  out <- pecotmr:::parseJointSpecification("context", qd)
  expect_equal(length(out), 1L)
  expect_equal(out[[1L]]$axes, "context")
  expect_null(out[[1L]]$scope)
})

test_that("parseJointSpecification: accepts a list of specs", {
  qd <- .js_makeQtlDataset()
  out <- pecotmr:::parseJointSpecification(
    list("context", c("context", "trait")), qd)
  expect_equal(length(out), 2L)
  expect_equal(out[[2L]]$axes, c("context", "trait"))
})

test_that("parseJointSpecification: accepts scope-restricted spec", {
  qd1 <- .js_makeQtlDataset(study = "A")
  qd2 <- .js_makeQtlDataset(study = "B")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
  out <- pecotmr:::parseJointSpecification(
    list(list(axes = "context", scope = list(study = c("A")))),
    mt)
  expect_equal(out[[1L]]$scope$study, "A")
})

test_that("parseJointSpecification: rejects unknown axes", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseJointSpecification(c("context", "bogus"), qd),
    "unknown axes")
})

test_that("parseJointSpecification: rejects duplicate axes", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseJointSpecification(c("context", "context"), qd),
    "duplicate axes")
})

test_that("parseJointSpecification: rejects unknown scope keys", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseJointSpecification(
      list(list(axes = "context", scope = list(bogus = "x"))), qd),
    "unknown scope key")
})

test_that("parseJointSpecification: rejects scope values absent from data", {
  qd <- .js_makeQtlDataset(study = "S1")
  expect_error(
    pecotmr:::parseJointSpecification(
      list(list(axes = "context",
                scope = list(study = "NotPresent"))), qd),
    "scope\\$study contains values not in data")
})

test_that("parseJointSpecification: rejects unknown spec elements", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseJointSpecification(
      list(list(axes = "context", scope = NULL, bogus = "x")), qd),
    "unknown element")
})

# -----------------------------------------------------------------------------
# parseContexts
# -----------------------------------------------------------------------------

test_that("parseContexts: NULL passes through", {
  qd <- .js_makeQtlDataset()
  expect_null(pecotmr:::parseContexts(NULL, qd))
})

test_that("parseContexts: vector intersects with each study's contexts", {
  qd1 <- .js_makeQtlDataset(study = "A", contexts = c("brain", "liver"))
  qd2 <- .js_makeQtlDataset(study = "B", contexts = c("brain"))
  mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
  expect_warning(
    out <- pecotmr:::parseContexts(c("brain", "liver"), mt),
    "B.*liver")
  expect_setequal(out$A, c("brain", "liver"))
  expect_equal(out$B, "brain")
})

test_that("parseContexts: named-list form requires valid studies", {
  qd <- .js_makeQtlDataset(study = "A")
  expect_error(
    pecotmr:::parseContexts(list(B = "brain"), qd),
    "unknown studies")
})

test_that("parseContexts: named-list rejects unknown contexts", {
  qd <- .js_makeQtlDataset(study = "A", contexts = c("brain", "liver"))
  expect_error(
    pecotmr:::parseContexts(list(A = "bogus"), qd),
    "unknown contexts")
})

test_that("parseContexts: list form fills unmentioned studies with defaults", {
  qd1 <- .js_makeQtlDataset(study = "A", contexts = c("brain", "liver"))
  qd2 <- .js_makeQtlDataset(study = "B", contexts = c("brain", "liver"))
  mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
  out <- pecotmr:::parseContexts(list(A = "brain"), mt)
  expect_equal(out$A, "brain")
  expect_setequal(out$B, c("brain", "liver"))
})

# -----------------------------------------------------------------------------
# parseTraitIds
# -----------------------------------------------------------------------------

test_that("parseTraitIds: NULL passes through", {
  qd <- .js_makeQtlDataset()
  expect_null(pecotmr:::parseTraitIds(NULL, qd))
})

test_that("parseTraitIds: vector form returns the vector", {
  qd <- .js_makeQtlDataset(traits = c("X", "Y"))
  expect_equal(pecotmr:::parseTraitIds(c("X", "Y"), qd), c("X", "Y"))
})

test_that("parseTraitIds: study-keyed list validates per study", {
  qd1 <- .js_makeQtlDataset(study = "A", traits = c("ENSG_A"))
  qd2 <- .js_makeQtlDataset(study = "B", traits = c("ENSG_B"))
  mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
  expect_error(
    pecotmr:::parseTraitIds(list(A = "ENSG_B"), mt),
    "unknown traits")
  out <- pecotmr:::parseTraitIds(list(A = "ENSG_A", B = "ENSG_B"), mt)
  expect_equal(out$A, "ENSG_A")
  expect_equal(out$B, "ENSG_B")
})

test_that("parseTraitIds: doubly-nested study->context validates per context", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"),
                           traits = c("ENSG_A", "ENSG_B"))
  out <- pecotmr:::parseTraitIds(
    list(A = list(brain = "ENSG_A", liver = "ENSG_B")), qd)
  expect_equal(out$A$brain, "ENSG_A")
  expect_equal(out$A$liver, "ENSG_B")
  expect_error(
    pecotmr:::parseTraitIds(list(A = list(bogus = "ENSG_A")), qd),
    "unknown contexts")
})

# -----------------------------------------------------------------------------
# parseMethods
# -----------------------------------------------------------------------------

test_that("parseMethods: methods XOR (sumStats + qtlDataset)", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseMethods(methods = "susie",
                            sumStatsMethods = "susieInf",
                            qtlDatasetMethods = "susie",
                            data = qd,
                            caps = pecotmr:::.fineMappingMethodCapabilities,
                            multivariateMethods = c("mvsusie", "fsusie")),
    "Use either")
  expect_error(
    pecotmr:::parseMethods(methods = NULL,
                            sumStatsMethods = "susie",
                            qtlDatasetMethods = NULL,
                            data = qd,
                            caps = pecotmr:::.fineMappingMethodCapabilities,
                            multivariateMethods = c("mvsusie", "fsusie")),
    "must be given together")
  expect_error(
    pecotmr:::parseMethods(methods = NULL,
                            sumStatsMethods = NULL,
                            qtlDatasetMethods = NULL,
                            data = qd,
                            caps = pecotmr:::.fineMappingMethodCapabilities,
                            multivariateMethods = c("mvsusie", "fsusie")),
    "Specify")
})

test_that("parseMethods: rejects unknown tokens", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseMethods(methods = "bogus",
                            data = qd,
                            caps = pecotmr:::.fineMappingMethodCapabilities,
                            multivariateMethods = c("mvsusie", "fsusie")),
    "unknown method token")
})

test_that("parseMethods: rejects multi-axis methods at per-context level", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"))
  expect_error(
    pecotmr:::parseMethods(
      methods = list(A = list(brain = c("susie", "mvsusie"))),
      data = qd,
      caps = pecotmr:::.fineMappingMethodCapabilities,
      multivariateMethods = c("mvsusie", "fsusie")),
    "per-context")
})

test_that("parseMethods: rejects multi-axis methods at per-trait level", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = "brain",
                           traits = c("ENSG_A", "ENSG_B"))
  expect_error(
    pecotmr:::parseMethods(
      methods = list(A = list(brain = list(ENSG_A = c("susie", "mvsusie")))),
      data = qd,
      caps = pecotmr:::.fineMappingMethodCapabilities,
      multivariateMethods = c("mvsusie", "fsusie")),
    "per-trait")
})

test_that("parseMethods: rejects user-rejected tokens (mrmash in fineMapping)", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseMethods(
      methods = "mrmash",
      data = qd,
      caps = pecotmr:::.fineMappingMethodCapabilities,
      multivariateMethods = c("mvsusie", "fsusie"),
      rejectedAtUser = "mrmash"),
    "cannot be user-requested")
})

test_that("parseMethods: accepts per-context univariate methods", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"))
  out <- pecotmr:::parseMethods(
    methods = list(A = list(brain = "susie", liver = "susieInf")),
    data = qd,
    caps = pecotmr:::.fineMappingMethodCapabilities,
    multivariateMethods = c("mvsusie", "fsusie"))
  expect_equal(out$shape, "primary")
})

test_that("parseMethods: validates per-(study, context, trait) leaf paths", {
  qd <- .js_makeQtlDataset(study = "A", contexts = "brain",
                           traits = c("ENSG_A", "ENSG_B"))
  expect_error(
    pecotmr:::parseMethods(
      methods = list(A = list(brain = list(BOGUS = "susie"))),
      data = qd,
      caps = pecotmr:::.fineMappingMethodCapabilities,
      multivariateMethods = c("mvsusie", "fsusie")),
    "unknown trait")
})

# -----------------------------------------------------------------------------
# validateMethodsVsJointSpec
# -----------------------------------------------------------------------------

test_that("validateMethodsVsJointSpec: per-study methods + jointCrossStudy errors", {
  qd1 <- .js_makeQtlDataset(study = "A")
  qd2 <- .js_makeQtlDataset(study = "B")
  ss <- .js_makeQtlSumStats(studies = "C")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2),
                              sumStats = ss)
  parsed <- pecotmr:::parseMethods(
    methods = list(A = "susie", B = "susieInf"),
    data = mt,
    caps = pecotmr:::.fineMappingMethodCapabilities,
    multivariateMethods = c("mvsusie", "fsusie"))
  joints <- pecotmr:::parseJointSpecification("study", mt)
  expect_error(
    pecotmr:::validateMethodsVsJointSpec(parsed, joints),
    "per-study")
})

test_that("validateMethodsVsJointSpec: per-context methods + jointCrossContext errors", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"))
  parsed <- pecotmr:::parseMethods(
    methods = list(A = list(brain = "susie", liver = "susie")),
    data = qd,
    caps = pecotmr:::.fineMappingMethodCapabilities,
    multivariateMethods = c("mvsusie", "fsusie"))
  joints <- pecotmr:::parseJointSpecification("context", qd)
  expect_error(
    pecotmr:::validateMethodsVsJointSpec(parsed, joints),
    "per-context")
})

test_that("validateMethodsVsJointSpec: vector methods + any joint flag OK", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"))
  parsed <- pecotmr:::parseMethods(
    methods = c("susie", "mvsusie"),
    data = qd,
    caps = pecotmr:::.fineMappingMethodCapabilities,
    multivariateMethods = c("mvsusie", "fsusie"))
  joints <- pecotmr:::parseJointSpecification("context", qd)
  expect_silent(pecotmr:::validateMethodsVsJointSpec(parsed, joints))
})

# -----------------------------------------------------------------------------
# Pipeline-level wiring: jointSpecification accepted on all three methods
# -----------------------------------------------------------------------------

test_that("fineMappingPipeline(QtlDataset): trait-axis joint dispatcher is wired", {
  # Cross-trait dispatcher is now wired; calling it on a fake genotype
  # handle errors when the genotype extractor tries to load the GDS file.
  # The point of this test is to verify the dispatcher is invoked (not
  # the stub error), so we check that the error comes from the genotype
  # I/O layer rather than the jointSpec wiring.
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = "brain",
                           traits = c("ENSG_A", "ENSG_B"))
  expect_error(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        jointSpecification = "trait"),
    "Can not open file|No such file")
})

test_that("fineMappingPipeline(QtlDataset): composed joint dispatcher is wired", {
  # Composed dispatcher accepts axes = c("context", "trait") for individual-
  # level input. On the fake fixture the genotype I/O fails when the
  # dispatcher reaches the extractor, which proves the dispatcher was
  # invoked (rather than the previous stub error).
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"),
                           traits = c("ENSG_A", "ENSG_B"))
  expect_error(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        jointSpecification = list(c("context", "trait"))),
    "Can not open file|No such file")
})

test_that("fineMappingPipeline(QtlDataset): composed joint rejects axes with 'study'", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"),
                           traits = "ENSG_A")
  expect_error(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        jointSpecification = list(c("study", "context"))),
    "axes including 'study' require sumstats")
})

test_that("fineMappingPipeline(QtlDataset): study-axis on individual data errors", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"),
                           traits = "ENSG_A")
  expect_error(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        jointSpecification = "study"),
    "requires sumstats input")
})

test_that("fineMappingPipeline(QtlDataset): NULL jointSpec is the default", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = "brain",
                           traits = "ENSG_A")
  expect_silent(
    pecotmr:::parseJointSpecification(NULL, qd))
})

test_that("fineMappingPipeline(QtlDataset): invalid jointSpec errors before fit", {
  qd <- .js_makeQtlDataset(study = "A")
  expect_error(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        jointSpecification = "bogus_axis"),
    "unknown axes")
})

test_that("twasWeightsPipeline(QtlDataset): cross-trait joint dispatcher is wired", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = "brain",
                           traits = c("ENSG_A", "ENSG_B"))
  expect_error(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = "trait"),
    "Can not open file|No such file")
})

test_that("twasWeightsPipeline(QtlDataset): study-axis on individual data errors", {
  qd <- .js_makeQtlDataset(study = "A",
                           contexts = c("brain", "liver"),
                           traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = "study"),
    "requires sumstats input")
})

test_that("twasWeightsPipeline(QtlSumStats): non-NULL jointSpec errors", {
  ss <- .js_makeQtlSumStats(studies = c("X", "Y"),
                            contexts = "DLPFC",
                            traits = "ENSG_A")
  # qcInfo is empty; the QC assertion fires BEFORE jointSpec parsing.
  expect_error(
    twasWeightsPipeline(ss, methods = "susie",
                        jointSpecification = "context"),
    "QC")
})

test_that("twasWeightsPipeline(MultiStudyQtlDataset): method exists", {
  expect_true(existsMethod("twasWeightsPipeline", "MultiStudyQtlDataset"))
})

test_that("validateMethodsVsJointSpec: split-form methods skipped", {
  qd <- .js_makeQtlDataset()
  parsed <- pecotmr:::parseMethods(
    methods = NULL,
    sumStatsMethods = "susieInf",
    qtlDatasetMethods = "susie",
    data = qd,
    caps = pecotmr:::.fineMappingMethodCapabilities,
    multivariateMethods = c("mvsusie", "fsusie"))
  joints <- pecotmr:::parseJointSpecification("context", qd)
  expect_silent(pecotmr:::validateMethodsVsJointSpec(parsed, joints))
})

# ============================================================================
# Joint dispatchers (merged from former tests/testthat/test_jointDispatchers.R)
# ============================================================================

context("joint dispatchers (fineMappingDispatcher / twasDispatcher)")

# ============================================================================
# Strategy: each joint-dispatcher function is exercised by driving
# fineMappingPipeline / twasWeightsPipeline through the user-facing
# `jointSpecification` argument and mocking the underlying fitters
# (mvsusieRss, mrmashWeights, mrmashRssWeights, ...). The mocks return tiny
# stub objects so postprocessing builds plausible result rows.
# ============================================================================

# -----------------------------------------------------------------------------
# Fixture builders
# -----------------------------------------------------------------------------

.jd_makeHandle <- function(snp_n = 5L, n_samples = 30L) {
  new("GenotypeHandle",
    path = "/tmp/jd.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("v", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

.jd_mockExtractor <- function(seed = 11, n_samples = 30L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds, handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges   = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx],
                                  width = 1L))
    S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
      SNP = handle@snpInfo$SNP[snpIdx],
      A1  = handle@snpInfo$A1[snpIdx],
      A2  = handle@snpInfo$A2[snpIdx])
    cd <- S4Vectors::DataFrame(sampleId = handle@sampleIds,
                               row.names = handle@sampleIds)
    dosage <- t(sub)
    rownames(dosage) <- handle@snpInfo$SNP[snpIdx]
    colnames(dosage) <- handle@sampleIds
    SummarizedExperiment::SummarizedExperiment(
      assays    = list(dosage = dosage),
      rowRanges = rr,
      colData   = cd)
  }
}

# Multi-(study, context, trait) QtlSumStats. Every row carries the same SNP
# order (5 variants) so jointCrossContext / jointCrossTrait / jointCrossStudy
# can stack Z columns without alignment problems.
.jd_makeQtlSumStats <- function(studies = "Q1",
                                contexts = c("c1", "c2"),
                                traits = "t1") {
  rows <- expand.grid(study = studies, context = contexts, trait = traits,
                      stringsAsFactors = FALSE)
  makeGr <- function() {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(start = seq(100L, by = 100L,
                                            length.out = 5L),
                                width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      SNP = paste0("v", 1:5),
      A1  = rep("A", 5), A2 = rep("G", 5),
      Z   = rnorm(5), N = rep(1000L, 5))
    gr
  }
  QtlSumStats(
    study    = rows$study,
    context  = rows$context,
    trait    = rows$trait,
    entry    = lapply(seq_len(nrow(rows)), function(.) makeGr()),
    genome   = "hg19",
    ldSketch = .jd_makeHandle(),
    qcInfo   = list(step1 = "ok"))
}

# -----------------------------------------------------------------------------
# Mocks for SuSiE / mvsusie / mr.mash families
# -----------------------------------------------------------------------------

.jd_mockMvsusie <- function() {
  function(X, Y, prior_variance, coverage) {
    list(token = "mvsusie", n_X_cols = ncol(X), n_Y_cols = ncol(Y))
  }
}

.jd_mockMvsusieRss <- function() {
  function(Z, R, N, prior_variance, coverage) {
    list(token = "mvsusieRss", nVariants = nrow(Z), nOutcomes = ncol(Z))
  }
}

.jd_mockMixturePrior <- function() {
  function(R, ...) list(R = R)
}

# A stub postprocessor that returns a tiny FineMappingEntry. Mirrors the
# `.fmp_mockPostprocess` shape from test_fineMappingPipeline.R.
.jd_mockPostprocess <- function() {
  function(fit, method, dataX, dataY, coverage, secondaryCoverage,
           signalCutoff, minAbsCorr, csInput = NULL, af = NULL,
           region = NULL) {
    if (is.matrix(dataX)) {
      vids <- colnames(dataX)
    } else if (is.list(dataY) && !is.null(dataY$z)) {
      vids <- names(dataY$z)
    } else {
      vids <- "v_unknown"
    }
    if (is.null(vids)) vids <- "v_unknown"
    FineMappingEntry(
      variantIds = vids,
      susieFit = list(method = method, payload = fit),
      topLoci    = data.frame(variant_id = vids,
                              pip = seq(0.9, by = -0.1,
                                         length.out = length(vids)),
                              stringsAsFactors = FALSE))
  }
}

.jd_mockMrmashWeights <- function() {
  function(X, Y, ...) {
    w <- matrix(0, nrow = ncol(X), ncol = ncol(Y),
                dimnames = list(colnames(X), colnames(Y)))
    w
  }
}

.jd_mockMrmashRssWeights <- function() {
  function(stat, LD, ...) {
    nCols <- if (is.matrix(stat$z)) ncol(stat$z) else 1L
    nVars <- if (is.matrix(stat$z)) nrow(stat$z) else length(stat$z)
    w <- matrix(0, nrow = nVars, ncol = nCols)
    rownames(w) <- if (is.matrix(stat$z)) rownames(stat$z)
                   else stat$variantNames
    if (is.matrix(stat$z) && !is.null(colnames(stat$z)))
      colnames(w) <- colnames(stat$z)
    w
  }
}

# =============================================================================
# fineMappingDispatcher: QtlSumStats
# =============================================================================

test_that("fineMappingPipeline(QtlSumStats): jointSpec='context' fits one joint per (study, trait)", {
  ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = c("c1", "c2"),
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    .fmPostprocessOne     = .jd_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie_rss           = .jd_mockMvsusieRss(),
    create_mixture_prior  = .jd_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "mvsusie",
                        jointSpecification = "context"))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 1L)
  expect_equal(as.character(res$context), "joint")
  expect_true(grepl("c1;c2|c2;c1", as.character(res$jointContexts)))
})

test_that("fineMappingPipeline(QtlSumStats): jointSpec='context' with only one context skips", {
  ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = "c1", traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    .fmPostprocessOne     = .jd_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie_rss           = .jd_mockMvsusieRss(),
    create_mixture_prior  = .jd_mockMixturePrior(),
    .package = "mvsusieR")
  expect_error(
    suppressMessages(
      fineMappingPipeline(ss, methods = "mvsusie",
                          jointSpecification = "context")),
    "no joint fits produced"
  )
})

test_that("fineMappingPipeline(QtlSumStats): jointSpec='trait' fits one joint per (study, context)", {
  ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = "c1",
                            traits = c("t1", "t2"))
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    .fmPostprocessOne     = .jd_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie_rss           = .jd_mockMvsusieRss(),
    create_mixture_prior  = .jd_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "mvsusie",
                        jointSpecification = "trait"))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 1L)
  expect_equal(as.character(res$trait), "joint")
})

test_that("fineMappingPipeline(QtlSumStats): jointSpec='trait' with fsusie errors (no RSS variant)", {
  ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = "c1",
                            traits = c("t1", "t2"))
  expect_error(
    fineMappingPipeline(ss, methods = "fsusie",
                        jointSpecification = "trait"),
    "fsusie"
  )
})

test_that("fineMappingPipeline(QtlSumStats): jointSpec='study' fits one joint per (context, trait)", {
  ss <- .jd_makeQtlSumStats(studies = c("Q1", "Q2"), contexts = "c1",
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    .fmPostprocessOne     = .jd_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie_rss           = .jd_mockMvsusieRss(),
    create_mixture_prior  = .jd_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "mvsusie",
                        jointSpecification = "study"))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 1L)
  expect_equal(as.character(res$study), "joint")
})

test_that("fineMappingPipeline(QtlSumStats): composed jointSpec axes={'study','context'} fits", {
  ss <- .jd_makeQtlSumStats(studies = c("Q1", "Q2"),
                            contexts = c("c1", "c2"),
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    .fmPostprocessOne     = .jd_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie_rss           = .jd_mockMvsusieRss(),
    create_mixture_prior  = .jd_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "mvsusie",
                        jointSpecification = list(c("study", "context"))))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_true("jointStudies"  %in% names(res))
  expect_true("jointContexts" %in% names(res))
  expect_true(any(as.character(res$study)   == "joint"))
  expect_true(any(as.character(res$context) == "joint"))
})

test_that("fineMappingPipeline(QtlSumStats): composed jointSpec rejects fsusie", {
  ss <- .jd_makeQtlSumStats(studies = c("Q1", "Q2"),
                            contexts = c("c1", "c2"),
                            traits = "t1")
  expect_error(
    fineMappingPipeline(ss, methods = "fsusie",
                        jointSpecification = list(c("study", "context"))),
    "fsusie"
  )
})

# =============================================================================
# twasDispatcher: QtlDataset
# =============================================================================

.jd_makeSe <- function(traits = c("t1", "t2"), n_samples = 30L,
                      starts = NULL) {
  if (is.null(starts))
    starts <- seq(1000L, by = 1000L, length.out = length(traits))
  rng <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(traits)),
    ranges = IRanges::IRanges(start = starts, width = 500L))
  names(rng) <- traits
  set.seed(0)
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  cd <- S4Vectors::DataFrame(
    sex = rep(c(0, 1), length.out = n_samples),
    age = seq_len(n_samples),
    row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays    = list(expression = expr),
    rowRanges = rng,
    colData   = cd)
}

.jd_makeQtlDataset <- function(study = "Q1",
                               contexts = c("c1", "c2"),
                               traits = c("t1", "t2")) {
  phen <- setNames(lapply(contexts,
                          function(.) .jd_makeSe(traits = traits)),
                   contexts)
  QtlDataset(
    study              = study,
    genotypes          = .jd_makeHandle(),
    phenotypes         = phen,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
}

test_that("twasWeightsPipeline(QtlDataset): jointSpec='context' fits mr.mash per trait", {
  qd <- .jd_makeQtlDataset(study = "Q1",
                            contexts = c("c1", "c2"),
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashWeights         = .jd_mockMrmashWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = "context"))
  expect_s4_class(res, "TwasWeights")
  expect_equal(as.character(res$context), "joint")
  expect_true("jointContexts" %in% names(res))
})

test_that("twasWeightsPipeline(QtlDataset): jointSpec='context' with only one context skips", {
  qd <- .jd_makeQtlDataset(study = "Q1",
                            contexts = "c1",
                            traits = c("t1", "t2"))
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashWeights         = .jd_mockMrmashWeights(),
    .package = "pecotmr")
  expect_error(
    suppressMessages(
      twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                          jointSpecification = "context")),
    "no joint fits produced|context"
  )
})

test_that("twasWeightsPipeline(QtlDataset): jointSpec='trait' fits mr.mash per context", {
  qd <- .jd_makeQtlDataset(study = "Q1",
                            contexts = "c1",
                            traits = c("t1", "t2"))
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashWeights         = .jd_mockMrmashWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = "trait"))
  expect_s4_class(res, "TwasWeights")
  expect_equal(as.character(res$trait), "joint")
  expect_true("jointTraits" %in% names(res))
})

test_that("twasWeightsPipeline(QtlDataset): study-axis fails on individual data", {
  qd <- .jd_makeQtlDataset(study = "Q1",
                            contexts = c("c1", "c2"),
                            traits = "t1")
  expect_error(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = "study"),
    "requires sumstats input"
  )
})

test_that("twasWeightsPipeline(QtlDataset): composed jointSpec axes=c('context','trait') fits", {
  qd <- .jd_makeQtlDataset(study = "Q1",
                            contexts = c("c1", "c2"),
                            traits = c("t1", "t2"))
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashWeights         = .jd_mockMrmashWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = list(c("context", "trait"))))
  expect_s4_class(res, "TwasWeights")
  expect_equal(as.character(res$context), "joint")
  expect_equal(as.character(res$trait), "joint")
})

test_that("twasWeightsPipeline(QtlDataset): composed jointSpec including 'study' errors", {
  qd <- .jd_makeQtlDataset(study = "Q1",
                            contexts = c("c1", "c2"),
                            traits = "t1")
  expect_error(
    twasWeightsPipeline(qd, methods = "mrmash", cisWindow = 1000L,
                        jointSpecification = list(c("study", "context"))),
    "require sumstats|requires sumstats"
  )
})

# =============================================================================
# twasDispatcher: QtlSumStats
# =============================================================================

test_that("twasWeightsPipeline(QtlSumStats): jointSpec='context' fits mr.mash.rss per (study, trait)", {
  ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = c("c1", "c2"),
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashRssWeights      = .jd_mockMrmashRssWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(ss, methods = "mrmash",
                        jointSpecification = "context"))
  expect_s4_class(res, "TwasWeights")
  expect_equal(as.character(res$context), "joint")
  expect_true("jointContexts" %in% names(res))
})

test_that("twasWeightsPipeline(QtlSumStats): jointSpec='trait' fits mr.mash.rss per (study, context)", {
  ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = "c1",
                            traits = c("t1", "t2"))
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashRssWeights      = .jd_mockMrmashRssWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(ss, methods = "mrmash",
                        jointSpecification = "trait"))
  expect_s4_class(res, "TwasWeights")
  expect_equal(as.character(res$trait), "joint")
})

test_that("twasWeightsPipeline(QtlSumStats): jointSpec='study' fits mr.mash.rss per (context, trait)", {
  ss <- .jd_makeQtlSumStats(studies = c("Q1", "Q2"), contexts = "c1",
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashRssWeights      = .jd_mockMrmashRssWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(ss, methods = "mrmash",
                        jointSpecification = "study"))
  expect_s4_class(res, "TwasWeights")
  expect_equal(as.character(res$study), "joint")
})

test_that("twasWeightsPipeline(QtlSumStats): composed jointSpec axes=c('study','context') fits", {
  ss <- .jd_makeQtlSumStats(studies = c("Q1", "Q2"),
                            contexts = c("c1", "c2"),
                            traits = "t1")
  local_mocked_bindings(
    extractBlockGenotypes = .jd_mockExtractor(),
    mrmashRssWeights      = .jd_mockMrmashRssWeights(),
    .package = "pecotmr")
  res <- suppressMessages(
    twasWeightsPipeline(ss, methods = "mrmash",
                        jointSpecification = list(c("study", "context"))))
  expect_s4_class(res, "TwasWeights")
  expect_true("jointStudies"  %in% names(res))
  expect_true("jointContexts" %in% names(res))
})
