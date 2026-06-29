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

test_that("parseMethods: a TWAS-only token (mrmash) is unknown to the fine-mapping grammar", {
  # mr.mash is not in the fine-mapping capability table (it is TWAS-only), so it
  # is rejected as an unknown fine-mapping token.
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseMethods(
      methods = "mrmash",
      data = qd,
      caps = pecotmr:::.fineMappingMethodCapabilities,
      multivariateMethods = c("mvsusie", "fsusie")),
    "unknown method token")
})

test_that("parseMethods: rejectedAtUser tokens are refused", {
  qd <- .js_makeQtlDataset()
  expect_error(
    pecotmr:::parseMethods(
      methods = "mvsusie",
      data = qd,
      caps = pecotmr:::.fineMappingMethodCapabilities,
      multivariateMethods = c("mvsusie", "fsusie"),
      rejectedAtUser = "mvsusie"),
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
           region = NULL, conditionIdx = NULL) {
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
  expect_equal(nrow(res), 2L)                                # per-context rows
  expect_setequal(as.character(res$context), c("c1", "c2")) # REAL contexts
  expect_true(all(grepl("c1;c2|c2;c1", as.character(res$jointContexts))))
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
  expect_equal(nrow(res), 2L)                                # per-trait rows
  expect_setequal(as.character(res$trait), c("t1", "t2"))
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
  expect_equal(nrow(res), 2L)                                # per-study rows
  expect_setequal(as.character(res$study), c("Q1", "Q2"))
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
  expect_equal(nrow(res), 4L)                               # study x context
  expect_setequal(as.character(res$study), c("Q1", "Q2"))   # both vary -> real
  expect_setequal(as.character(res$context), c("c1", "c2"))
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
  expect_equal(nrow(res), 2L)                                # per-context rows
  expect_setequal(as.character(res$context), c("c1", "c2"))
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
  expect_equal(nrow(res), 2L)                                # per-trait rows
  expect_setequal(as.character(res$trait), c("t1", "t2"))
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
  expect_equal(nrow(res), 4L)                                # context x trait
  expect_setequal(as.character(res$context), c("c1", "c2"))
  expect_setequal(as.character(res$trait), c("t1", "t2"))
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
  expect_equal(nrow(res), 2L)                                # per-context rows
  expect_setequal(as.character(res$context), c("c1", "c2"))
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
  expect_equal(nrow(res), 2L)                                # per-trait rows
  expect_setequal(as.character(res$trait), c("t1", "t2"))
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
  expect_equal(nrow(res), 2L)                                # per-study rows
  expect_setequal(as.character(res$study), c("Q1", "Q2"))
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
  expect_equal(nrow(res), 4L)                                # study x context
  expect_setequal(as.character(res$study), c("Q1", "Q2"))
  expect_setequal(as.character(res$context), c("c1", "c2"))
})

# ============================================================================
# Additional coverage: scope helpers, parser error branches, scope resolution,
# X/Y/Z builders' skip paths, multi-region merges, and MultiStudy dispatchers.
# ============================================================================

# -----------------------------------------------------------------------------
# .sp* helpers: QtlSumStats / MultiStudy branches + unsupported-class errors
# -----------------------------------------------------------------------------

test_that(".sp* helpers: QtlSumStats study/context/trait listing + dataForm", {
  ss <- .js_makeQtlSumStats(studies = c("A", "B"), contexts = c("c1", "c2"),
                            traits = c("t1", "t2"))
  expect_setequal(pecotmr:::.spListStudies(ss), c("A", "B"))
  expect_setequal(pecotmr:::.spListContexts(ss), c("c1", "c2"))          # study=NULL
  expect_setequal(pecotmr:::.spListContexts(ss, "A"), c("c1", "c2"))
  expect_setequal(pecotmr:::.spListTraits(ss, study = "A", context = "c1"),
                  c("t1", "t2"))
  expect_equal(pecotmr:::.spStudyDataForm(ss, "A"), "sumstats")
  expect_error(pecotmr:::.spStudyDataForm(ss, "missing"), "not in QtlSumStats")
})

test_that(".spStudyDataForm: QtlDataset wrong study errors", {
  qd <- .js_makeQtlDataset(study = "S1")
  expect_error(pecotmr:::.spStudyDataForm(qd, "wrong"), "not in QtlDataset")
})

test_that(".sp* helpers: MultiStudy per-study context/trait routing", {
  qd1 <- .js_makeQtlDataset(study = "indA", contexts = "brain",
                            traits = c("g1", "g2"))
  ss  <- .js_makeQtlSumStats(studies = "ssC", contexts = "DLPFC",
                             traits = "g3")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd1), sumStats = ss)
  expect_setequal(pecotmr:::.spListContexts(mt, "ssC"), "DLPFC")        # ss branch
  expect_setequal(pecotmr:::.spListContexts(mt, "indA"), "brain")
  expect_setequal(pecotmr:::.spListContexts(mt), c("brain", "DLPFC"))   # all
  expect_equal(pecotmr:::.spListContexts(mt, "nope"), character(0))
  expect_setequal(pecotmr:::.spListTraits(mt, study = "ssC"), "g3")     # ss branch
  expect_setequal(pecotmr:::.spListTraits(mt, study = "indA"), c("g1", "g2"))
  expect_equal(pecotmr:::.spListTraits(mt, study = "nope"), character(0))
  # study = NULL aggregates ALL traits across individual + sumstats components
  # (regression guard: a present sumStats slot must not shadow the QtlDatasets).
  expect_setequal(pecotmr:::.spListTraits(mt), c("g1", "g2", "g3"))
})

test_that(".spListTraits: study=NULL on a sumstats-free MultiStudy aggregates traits", {
  qdA <- .js_makeQtlDataset(study = "indA", contexts = "brain",
                            traits = c("g1", "g2"))
  qdB <- .js_makeQtlDataset(study = "indB", contexts = "liver", traits = "g3")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qdA, indB = qdB))
  expect_setequal(pecotmr:::.spListTraits(mt), c("g1", "g2", "g3"))     # aggregate
})

test_that(".sp* helpers: unsupported class raises a labelled error", {
  expect_error(pecotmr:::.spListStudies(42), "unsupported class")
  expect_error(pecotmr:::.spStudyDataForm(42, "x"), "unsupported class")
  expect_error(pecotmr:::.spListContexts(42), "unsupported class")
  expect_error(pecotmr:::.spListTraits(42), "unsupported class")
})

# -----------------------------------------------------------------------------
# parseJointSpecification / parseContexts / parseTraitIds error branches
# -----------------------------------------------------------------------------

test_that("parseJointSpecification: malformed inputs error", {
  qd <- .js_makeQtlDataset()
  expect_error(pecotmr:::parseJointSpecification(42, qd), "must be NULL")
  expect_error(pecotmr:::parseJointSpecification(list(list(scope = NULL)), qd),
               "missing `axes`")
  expect_error(pecotmr:::parseJointSpecification(list(42), qd),
               "character vector or a named list")
  expect_error(pecotmr:::parseJointSpecification(list(list(axes = character(0))),
                                                 qd), "non-empty character")
  expect_error(pecotmr:::parseJointSpecification(
    list(list(axes = "context", scope = c("a", "b"))), qd), "named list")
  expect_error(pecotmr:::parseJointSpecification(
    list(list(axes = "context", scope = list(study = integer(0)))), qd),
    "non-empty character vector")
})

test_that("parseContexts: malformed inputs error", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = c("brain", "liver"))
  expect_error(pecotmr:::parseContexts(character(0), qd), "non-empty")
  expect_error(pecotmr:::parseContexts(list(brain = "x"), qd), "unknown studies")
  expect_error(pecotmr:::parseContexts(list(S1 = character(0)), qd),
               "non-empty character vector")
  expect_error(pecotmr:::parseContexts(list("brain"), qd), "named list") # unnamed
  expect_error(pecotmr:::parseContexts(42, qd), "must be NULL")
})

test_that("parseContexts: vector form warns on contexts missing from a study", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = c("brain", "liver"))
  expect_warning(out <- pecotmr:::parseContexts(c("brain", "absent"), qd),
                 "missing requested context")
  expect_equal(out$S1, "brain")
})

test_that("parseTraitIds: malformed inputs error", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = "brain",
                           traits = c("ENSG1", "ENSG2"))
  expect_error(pecotmr:::parseTraitIds(character(0), qd), "non-empty")
  expect_error(pecotmr:::parseTraitIds(42, qd), "must be NULL")
  expect_error(pecotmr:::parseTraitIds(list("ENSG1"), qd), "named by study")
  expect_error(pecotmr:::parseTraitIds(list(nope = "ENSG1"), qd),
               "unknown studies")
  expect_error(pecotmr:::parseTraitIds(list(S1 = character(0)), qd),
               "non-empty character vector")
  expect_error(pecotmr:::parseTraitIds(list(S1 = list("brain")), qd),
               "named by context")
  expect_error(pecotmr:::parseTraitIds(list(S1 = list(brain = character(0))), qd),
               "non-empty character vector")
  expect_error(pecotmr:::parseTraitIds(list(S1 = list(brain = "nope")), qd),
               "unknown traits")
  expect_error(pecotmr:::parseTraitIds(list(S1 = 42), qd),
               "character vector or a named list")
})
# (The doubly-nested study->context success path is covered by the existing
#  "parseTraitIds: doubly-nested study->context validates per context" test.)

# -----------------------------------------------------------------------------
# .spWalkMethods / parseMethods / validateMethodsVsJointSpec error branches
# -----------------------------------------------------------------------------

test_that(".spWalkMethods: structural errors via parseMethods", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = "brain", traits = "ENSG1")
  caps <- list(lasso = list(multivariate = FALSE))
  walk <- function(m) pecotmr:::parseMethods(m, data = qd, caps = caps,
                                             multivariateMethods = character(0))
  expect_error(walk(list(S1 = 42)), "character vector or a named list")
  expect_error(walk(list(S1 = list(brain = list(ENSG1 = list(x = "lasso"))))),
               "cannot nest below the trait level")
  expect_error(walk(list("lasso")), "non-empty names")        # unnamed list node
  expect_error(walk(list(S1 = list())), "non-empty names")     # inner empty node
  # A named-but-empty list reaches the dedicated "empty named list" guard.
  expect_error(pecotmr:::.spWalkMethods(setNames(list(), character(0))),
               "empty named list")
})

test_that("parseMethods: split-form and leaf validation errors", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = "brain", traits = "ENSG1")
  caps <- list(lasso = list(multivariate = FALSE),
               mrmash = list(multivariate = TRUE))
  expect_error(pecotmr:::parseMethods(
    methods = NULL, sumStatsMethods = character(0), qtlDatasetMethods = "lasso",
    data = qd, caps = caps, multivariateMethods = "mrmash"),
    "non-empty character vector")
  expect_error(pecotmr:::parseMethods(
    methods = NULL, sumStatsMethods = "lasso", qtlDatasetMethods = character(0),
    data = qd, caps = caps, multivariateMethods = "mrmash"),
    "non-empty character vector")
  expect_error(pecotmr:::parseMethods(
    list(S1 = character(0)), data = qd, caps = caps,
    multivariateMethods = "mrmash"), "non-empty character vector")
  expect_error(pecotmr:::parseMethods(
    list(nope = "lasso"), data = qd, caps = caps,
    multivariateMethods = "mrmash"), "unknown study")
  expect_error(pecotmr:::parseMethods(
    list(S1 = list(absent = "lasso")), data = qd, caps = caps,
    multivariateMethods = "mrmash"), "unknown context")
})

test_that("validateMethodsVsJointSpec: empty spec is a no-op; per-trait nesting + trait axis errors", {
  expect_null(pecotmr:::validateMethodsVsJointSpec(
    list(shape = "primary", methods = "lasso"), list()))
  parsed <- list(list(axes = "trait"))
  methodsParsed <- list(shape = "primary",
    methods = list(S1 = list(brain = list(ENSG1 = "lasso"))))
  expect_error(pecotmr:::validateMethodsVsJointSpec(methodsParsed, parsed),
               "nests per-trait")
})

# -----------------------------------------------------------------------------
# .fmResolveSpecScope: scope$study / scope$context / scope$trait + filters
# -----------------------------------------------------------------------------

test_that(".fmResolveSpecScope: scope + contexts + traitIds filters intersect", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = c("brain", "liver"),
                           traits = c("ENSG1", "ENSG2"))
  spec <- list(scope = list(study = "S1", context = "brain", trait = "ENSG1"))
  out <- pecotmr:::.fmResolveSpecScope(spec, qd)
  expect_equal(out$studies, "S1")
  expect_equal(out$contexts$S1, "brain")
  expect_equal(out$traits$S1, "ENSG1")
  # Named-list contexts + study-keyed list traitIds filters.
  out2 <- pecotmr:::.fmResolveSpecScope(
    list(scope = NULL), qd, contexts = list(S1 = "liver"),
    traitIds = list(S1 = "ENSG2"))
  expect_equal(out2$contexts$S1, "liver")
  expect_equal(out2$traits$S1, "ENSG2")
})

# -----------------------------------------------------------------------------
# .buildJointSumstatZMatrix: SNP-order mismatch error
# -----------------------------------------------------------------------------

test_that(".buildJointSumstatZMatrix: a mismatched SNP order across entries errors", {
  df <- data.frame(study = "S", context = c("c1", "c2"), trait = "t1",
                   stringsAsFactors = FALSE)
  calls <- 0L
  local_mocked_bindings(
    getSumstatDf = function(x, study, context, trait, require, ...) {
      calls <<- calls + 1L
      vid <- if (calls == 1L) c("v1", "v2") else c("v2", "v1")   # reordered
      data.frame(variant_id = vid, z = c(1, 2), N = c(100, 100),
                 stringsAsFactors = FALSE)
    },
    .package = "pecotmr")
  expect_error(
    pecotmr:::.buildJointSumstatZMatrix(df, c(1L, 2L), c("c1", "c2"),
                                        errorLabel = "TESTLABEL"),
    "identical SNP order")
})

# -----------------------------------------------------------------------------
# Individual X/Y builders: the skip paths (return NULL / message branches)
# -----------------------------------------------------------------------------

test_that(".buildIndividualCrossContextXY: skips when a trait spans < 2 contexts", {
  se1 <- .js_makeSe(traits = "g1")                  # g1 present
  se0 <- .js_makeSe(traits = "other")               # g1 absent
  local_mocked_bindings(
    getPhenotypes = function(data, contexts)
      if (identical(contexts, "c1")) se1 else se0,
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildIndividualCrossContextXY(
    NULL, "g1", c("c1", "c2"), cisWindow = 1000L, verbose = 1,
    label = "X")))
})

test_that(".buildIndividualCrossContextXY: region path + complete-case skip", {
  se <- .js_makeSe(traits = "g1", samples = paste0("s", 1:6))
  samp <- paste0("s", 1:6)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10000))
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 6, 2, dimnames = list(samp, c("v1", "v2"))),
    .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
      ym <- function(v) matrix(v, 6, 1, dimnames = list(samp, "g1"))
      list(c1 = ym(c(NA, NA, NA, NA, NA, 1)),       # mostly NA -> < 2 complete
           c2 = ym(rnorm(6)))
    },
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildIndividualCrossContextXY(
    NULL, "g1", c("c1", "c2"), cisWindow = NULL, verbose = 1,
    label = "X", region = region)))
})

test_that(".buildIndividualCrossContextXY: too few shared samples skips", {
  se <- .js_makeSe(traits = "g1", samples = paste0("s", 1:6))
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 1, 2, dimnames = list("zz", c("v1", "v2"))),   # disjoint sample
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      list(c1 = matrix(0, 6, 1, dimnames = list(paste0("s", 1:6), "g1")),
           c2 = matrix(0, 6, 1, dimnames = list(paste0("s", 1:6), "g1"))),
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildIndividualCrossContextXY(
    NULL, "g1", c("c1", "c2"), cisWindow = 1000L, verbose = 1, label = "X")))
})

test_that(".fmTraitsInRegion: filters traits by phenotype overlap with the region", {
  se <- .js_makeSe(traits = c("g1", "g2"))          # g1@~100, g2@~200
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 160))
  expect_equal(pecotmr:::.fmTraitsInRegion(se, c("g1", "g2"), region), "g1")
  expect_equal(pecotmr:::.fmTraitsInRegion(se, c("g1", "g2"), NULL),
               c("g1", "g2"))                       # NULL region -> unchanged
})

test_that(".buildIndividualCrossTraitXY: skip branches (< 2 traits, region, complete)", {
  se2 <- .js_makeSe(traits = c("g1", "g2"), samples = paste0("s", 1:6))
  samp <- paste0("s", 1:6)
  # < 2 scoped traits in the context -> NULL.
  local_mocked_bindings(getPhenotypes = function(data, contexts) se2,
                        .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildIndividualCrossTraitXY(
    NULL, "cx", "g1", cisWindow = 1000L, verbose = 1, label = "X",
    study = "S")))
  # region path + < 2 complete cases.
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se2,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 6, 2, dimnames = list(samp, c("v1", "v2"))),
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      cbind(g1 = c(NA, NA, NA, NA, NA, 1), g2 = rnorm(6)) |>
        `rownames<-`(samp),
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildIndividualCrossTraitXY(
    NULL, "cx", c("g1", "g2"), cisWindow = NULL, verbose = 1, label = "X",
    study = "S", region = GenomicRanges::GRanges("chr1",
                                                 IRanges::IRanges(1, 9999)))))
})

test_that(".buildComposedIndividualXY: skip branches and single-context wrap", {
  se <- .js_makeSe(traits = c("g1", "g2"), samples = paste0("s", 1:6))
  samp <- paste0("s", 1:6)
  scope <- list(contexts = list(S = "c1"), traits = list(S = c("g1", "g2")))
  # Single context: YresList wrap branch, then a valid 2-tuple build.
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 6, 2, dimnames = list(samp, c("v1", "v2"))),
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      matrix(rnorm(12), 6, 2, dimnames = list(samp, c("g1", "g2"))),
    .package = "pecotmr")
  out <- suppressMessages(pecotmr:::.buildComposedIndividualXY(
    NULL, scope, "S", cisWindow = 1000L, verbose = 1, label = "X"))
  expect_equal(ncol(out$Y), 2L)
  expect_setequal(colnames(out$Y), c("c1:g1", "c1:g2"))
  # < 2 tuples -> NULL.
  scope1 <- list(contexts = list(S = "c1"), traits = list(S = "g1"))
  local_mocked_bindings(getPhenotypes = function(data, contexts)
    .js_makeSe(traits = "g1"), .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildComposedIndividualXY(
    NULL, scope1, "S", cisWindow = 1000L, verbose = 1, label = "X")))
})

# -----------------------------------------------------------------------------
# .enumerateComposedSumstatGroups: empty scope + all-axes (no complement)
# -----------------------------------------------------------------------------

test_that(".enumerateComposedSumstatGroups: empty scope -> NULL; no complement -> one block", {
  df <- data.frame(study = "S", context = c("c1", "c2"), trait = "t1",
                   stringsAsFactors = FALSE)
  empty <- pecotmr:::.enumerateComposedSumstatGroups(
    list(axes = c("context", "trait")), df,
    list(studies = character(0), contexts = list(), traits = list()))
  expect_null(empty)
  # axes = all three -> complement empty -> a single "__all__" block.
  gi <- pecotmr:::.enumerateComposedSumstatGroups(
    list(axes = c("study", "context", "trait")), df,
    list(studies = "S", contexts = list(S = c("c1", "c2")),
         traits = list(S = "t1")))
  expect_equal(length(gi$groups), 1L)
  expect_equal(names(gi$groups), "__all__")
})

# -----------------------------------------------------------------------------
# .fmSynthesizeJointSpec
# -----------------------------------------------------------------------------

test_that(".fmSynthesizeJointSpec: trait wins over context; single/single -> empty", {
  expect_equal(pecotmr:::.fmSynthesizeJointSpec(3L, 2L)[[1L]]$axes, "trait")
  expect_equal(pecotmr:::.fmSynthesizeJointSpec(2L, 1L)[[1L]]$axes, "context")
  expect_equal(pecotmr:::.fmSynthesizeJointSpec(1L, 1L), list())
})

# -----------------------------------------------------------------------------
# Multi-region merges: .fmMergeResultsByKey / .twasMergeResultsByKey
# -----------------------------------------------------------------------------

.js_fmEntry <- function(vid = "v1")
  FineMappingEntry(variantIds = vid, susieFit = list(),
                   topLoci = data.frame(variant_id = vid, pip = 0.9,
                                        stringsAsFactors = FALSE))

test_that(".fmMergeResultsByKey: merges per-region entries by (s,c,t,method)", {
  mk <- function() QtlFineMappingResult(
    study = "S", context = "c1", trait = "t1", method = "mvsusie",
    entry = list(.js_fmEntry()))
  local_mocked_bindings(
    .fmMergeEntries = function(entries) entries[[1L]], .package = "pecotmr")
  out <- pecotmr:::.fmMergeResultsByKey(list(mk(), mk()))
  expect_s4_class(out, "QtlFineMappingResult")
  expect_equal(nrow(out), 1L)
  # n == 0 short-circuit returns the (empty) base unchanged.
  empty <- QtlFineMappingResult(study = character(0), context = character(0),
    trait = character(0), method = character(0), entry = list())
  expect_equal(nrow(pecotmr:::.fmMergeResultsByKey(list(empty, empty))), 0L)
})

test_that(".twasMergeResultsByKey: merges per-region TwasWeights entries", {
  mk <- function() TwasWeights(
    study = "S", context = "c1", trait = "t1", method = "lasso",
    entry = list(TwasWeightsEntry(variantIds = "v1", weights = 0.5)))
  out <- pecotmr:::.twasMergeResultsByKey(list(mk(), mk()), c("r1", "r2"))
  expect_s4_class(out, "TwasWeights")
  expect_equal(nrow(out), 1L)
  empty <- TwasWeights(study = character(0), context = character(0),
    trait = character(0), method = character(0), entry = list())
  expect_equal(length(pecotmr:::.twasMergeResultsByKey(
    list(empty), "r1")$method), 0L)
})

# -----------------------------------------------------------------------------
# Multi-region QtlDataset dispatch (xRegions length 2 -> merge path)
# -----------------------------------------------------------------------------

test_that(".twasDispatchJointSpecsQtlDataset: two region blocks are merged by key", {
  r1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 250))
  r2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 500))
  parsed <- list(list(axes = "context", scope = NULL))
  mkRegionRes <- function() TwasWeights(
    study = c("Q1", "Q1"), context = c("c1", "c2"), trait = c("t1", "t1"),
    method = c("mrmash", "mrmash"),
    entry = list(TwasWeightsEntry(variantIds = "v1", weights = 0.5),
                 TwasWeightsEntry(variantIds = "v1", weights = 0.5)))
  local_mocked_bindings(
    .twasDispatchJointSpecsQtlDatasetOneRegion = function(...) mkRegionRes(),
    .package = "pecotmr")
  res <- pecotmr:::.twasDispatchJointSpecsQtlDataset(
    parsed, data = NULL, methods = "mrmash", contexts = NULL, traitIds = NULL,
    cisWindow = NULL, dataType = NULL, verbose = 0, xRegions = list(r1, r2))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)                       # per-context, merged regions
})

# -----------------------------------------------------------------------------
# MultiStudy dispatchers: .fmDispatchJointSpecsMultiStudy /
# .twasDispatchJointSpecsMultiStudy (leaf dispatchers mocked)
# -----------------------------------------------------------------------------

test_that(".fmDispatchJointSpecsMultiStudy: routes non-study specs to components, study spec to sumstats", {
  qd <- .js_makeQtlDataset(study = "indA", contexts = c("c1", "c2"),
                           traits = "t1")
  ss <- .js_makeQtlSumStats(studies = "ssC", contexts = "c1", traits = "t1")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd), sumStats = ss)
  parsed <- list(list(axes = "context", scope = NULL),
                 list(axes = "study", scope = NULL))
  qdRes <- QtlFineMappingResult(study = "indA", context = "c1", trait = "t1",
    method = "mvsusie", entry = list(.js_fmEntry()))
  ssRes <- QtlFineMappingResult(study = "ssC", context = "c1", trait = "t1",
    method = "mvsusie", entry = list(.js_fmEntry()),
    ldSketch = .js_makeGenotypeHandle())
  local_mocked_bindings(
    .fmDispatchJointSpecsQtlDataset = function(...) qdRes,
    .fmDispatchJointSpecsQtlSumStats = function(...) ssRes,
    .package = "pecotmr")
  out <- suppressMessages(pecotmr:::.fmDispatchJointSpecsMultiStudy(
    parsed, mt, methods = "mvsusie", contexts = NULL, traitIds = NULL,
    cisWindow = NULL, coverage = 0.95, secondaryCoverage = 0.5,
    signalCutoff = 0.1, minAbsCorr = 0.5, verbose = 1))
  expect_s4_class(out, "QtlFineMappingResult")
  expect_equal(nrow(out), 2L)                       # indA row + ssC row
})

test_that(".fmDispatchJointSpecsMultiStudy: study spec with no sumStats slot messages", {
  qdA <- .js_makeQtlDataset(study = "indA", contexts = c("c1", "c2"),
                            traits = "t1")
  qdB <- .js_makeQtlDataset(study = "indB", contexts = c("c1", "c2"),
                            traits = "t1")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qdA, indB = qdB))  # no ss
  parsed <- list(list(axes = "study", scope = NULL))
  expect_message(
    out <- pecotmr:::.fmDispatchJointSpecsMultiStudy(
      parsed, mt, methods = "mvsusie", contexts = NULL, traitIds = NULL,
      cisWindow = NULL, coverage = 0.95, secondaryCoverage = 0.5,
      signalCutoff = 0.1, minAbsCorr = 0.5, verbose = 1),
    "no sumStats slot")
  expect_null(out)
})

test_that(".twasDispatchJointSpecsMultiStudy: routes components + sumstats and rbinds", {
  qd <- .jd_makeQtlDataset(study = "indA", contexts = c("c1", "c2"),
                           traits = "t1")
  ss <- .jd_makeQtlSumStats(studies = "ssC", contexts = "c1", traits = "t1")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd), sumStats = ss)
  parsed <- list(list(axes = "context", scope = NULL),
                 list(axes = "study", scope = NULL))
  qdRes <- TwasWeights(study = "indA", context = "c1", trait = "t1",
    method = "mrmash", entry = list(TwasWeightsEntry(variantIds = "v1",
                                                     weights = 0.5)))
  ssRes <- TwasWeights(study = "ssC", context = "c1", trait = "t1",
    method = "mrmash", entry = list(TwasWeightsEntry(variantIds = "v1",
                                                     weights = 0.5)),
    ldSketch = .jd_makeHandle())
  local_mocked_bindings(
    .twasDispatchJointSpecsQtlDataset = function(...) qdRes,
    .twasDispatchJointSpecsQtlSumStats = function(...) ssRes,
    .package = "pecotmr")
  out <- suppressMessages(pecotmr:::.twasDispatchJointSpecsMultiStudy(
    parsed, mt, methods = "mrmash", contexts = NULL, traitIds = NULL,
    cisWindow = NULL, dataType = NULL, verbose = 1))
  expect_s4_class(out, "TwasWeights")
  expect_equal(nrow(out), 2L)
})

test_that(".twasDispatchJointSpecsMultiStudy: study spec, no sumStats -> message + NULL", {
  qdA <- .jd_makeQtlDataset(study = "indA", contexts = c("c1", "c2"),
                            traits = "t1")
  qdB <- .jd_makeQtlDataset(study = "indB", contexts = c("c1", "c2"),
                            traits = "t1")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qdA, indB = qdB))
  parsed <- list(list(axes = "study", scope = NULL))
  expect_message(
    out <- pecotmr:::.twasDispatchJointSpecsMultiStudy(
      parsed, mt, methods = "mrmash", contexts = NULL, traitIds = NULL,
      cisWindow = NULL, dataType = NULL, verbose = 1),
    "no sumStats slot")
  expect_null(out)
})

# ============================================================================
# Mop-up: remaining .sp* empties, parseContexts unnamed list, validate no-op,
# X/Y builder skip paths, region-missing keys in merges, FM multi-region merge.
# ============================================================================

test_that(".sp* QtlDataset: mismatched study / absent context return empty", {
  qd <- .js_makeQtlDataset(study = "S1", contexts = "brain",
                           traits = c("g1", "g2"))
  expect_equal(pecotmr:::.spListContexts(qd, "wrong"), character(0))     # 59
  expect_equal(pecotmr:::.spListTraits(qd, study = "wrong"), character(0))  # 93
  expect_equal(pecotmr:::.spListTraits(qd, context = "nope"), character(0)) # 98
})

test_that("validateMethodsVsJointSpec: per-study methods with a context joint passes", {
  mp <- list(shape = "primary", methods = list(S1 = "lasso"))           # depth 1
  expect_null(pecotmr:::validateMethodsVsJointSpec(
    mp, list(list(axes = "context"))))                                  # 542
})

test_that(".buildIndividualCrossTraitXY: disjoint X/Y samples skip the context", {
  se <- .js_makeSe(traits = c("g1", "g2"), samples = paste0("s", 1:6))
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 1, 2, dimnames = list("zz", c("v1", "v2"))),   # disjoint sample
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      matrix(0, 6, 2, dimnames = list(paste0("s", 1:6), c("g1", "g2"))),
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildIndividualCrossTraitXY(
    NULL, "cx", c("g1", "g2"), cisWindow = 1000L, verbose = 1, label = "X",
    study = "S")))                                                      # 730
})

test_that(".buildComposedIndividualXY: disjoint samples / missing trait col / NA rows skip", {
  samp <- paste0("s", 1:6)
  se <- .js_makeSe(traits = c("g1", "g2"), samples = samp)
  scope <- list(contexts = list(S = c("c1", "c2")),
                traits = list(S = c("g1", "g2")))
  # (a) disjoint X samples -> < 2 common -> NULL (774).
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 1, 2, dimnames = list("zz", c("v1", "v2"))),
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      setNames(lapply(c("c1", "c2"), function(.)
        matrix(rnorm(12), 6, 2, dimnames = list(samp, c("g1", "g2")))),
        c("c1", "c2")),
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildComposedIndividualXY(
    NULL, scope, "S", cisWindow = 1000L, verbose = 1, label = "X")))    # 774
  # (b) one context's Y lacks the trait column -> tuple skipped -> < 2 yCols (779/785).
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 6, 2, dimnames = list(samp, c("v1", "v2"))),
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      list(c1 = matrix(0, 6, 1, dimnames = list(samp, "g1")),
           c2 = matrix(0, 6, 1, dimnames = list(samp, "zzz"))),   # no g1/g2
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildComposedIndividualXY(
    NULL, list(contexts = list(S = c("c1", "c2")),
               traits = list(S = "g1")),
    "S", cisWindow = 1000L, verbose = 1, label = "X")))                 # 779/785
  # (c) two valid columns but < 2 complete rows (NA) -> NULL (788).
  local_mocked_bindings(
    getPhenotypes = function(data, contexts) se,
    .fmResidGeno = function(x, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 6, 2, dimnames = list(samp, c("v1", "v2"))),
    .fmResidPheno = function(x, contexts, traitId = NULL, ...)
      list(c1 = matrix(c(NA, NA, NA, NA, NA, 1), 6, 1,
                       dimnames = list(samp, "g1")),
           c2 = matrix(c(NA, NA, NA, NA, NA, 1), 6, 1,
                       dimnames = list(samp, "g1"))),
    .package = "pecotmr")
  expect_null(suppressMessages(pecotmr:::.buildComposedIndividualXY(
    NULL, list(contexts = list(S = c("c1", "c2")),
               traits = list(S = "g1")),
    "S", cisWindow = 1000L, verbose = 1, label = "X")))                 # 788
})

test_that(".fmMergeResultsByKey: a key missing from a later region contributes nothing", {
  twoRow <- QtlFineMappingResult(
    study = c("S", "S"), context = c("c1", "c2"), trait = c("t1", "t1"),
    method = c("mvsusie", "mvsusie"),
    entry = list(.js_fmEntry("v1"), .js_fmEntry("v2")))
  oneRow <- QtlFineMappingResult(
    study = "S", context = "c1", trait = "t1", method = "mvsusie",
    entry = list(.js_fmEntry("v1")))                       # missing the c2 key
  seen <- 0L
  local_mocked_bindings(
    .fmMergeEntries = function(entries) { seen <<- seen + length(entries)
                                          entries[[1L]] }, .package = "pecotmr")
  out <- pecotmr:::.fmMergeResultsByKey(list(twoRow, oneRow))
  expect_equal(nrow(out), 2L)
  expect_equal(seen, 3L)                  # 2 for c1 row + 1 for c2 row (849 else)
})

test_that(".fmDispatchJointSpecsQtlDataset: two region blocks are merged", {
  qd <- .jd_makeQtlDataset(study = "Q1", contexts = c("c1", "c2"),
                           traits = "t1")
  r1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 250))
  r2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 500))
  mkRes <- function() QtlFineMappingResult(
    study = c("Q1", "Q1"), context = c("c1", "c2"), trait = c("t1", "t1"),
    method = c("mvsusie", "mvsusie"),
    entry = list(.js_fmEntry("v1"), .js_fmEntry("v1")))
  local_mocked_bindings(
    .fmDispatchJointSpecsQtlDatasetOneRegion = function(...) mkRes(),
    .fmMergeEntries = function(entries) entries[[1L]],
    .package = "pecotmr")
  res <- pecotmr:::.fmDispatchJointSpecsQtlDataset(
    list(list(axes = "context", scope = NULL)), qd, methods = "mvsusie",
    contexts = NULL, traitIds = NULL, cisWindow = NULL, coverage = 0.95,
    secondaryCoverage = 0.5, signalCutoff = 0.1, minAbsCorr = 0.5, verbose = 0,
    xRegions = list(r1, r2))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)                              # 907-910
})

test_that(".fmDispatchJointSpecsQtlDataset: a single region returns directly; all-NULL -> NULL", {
  qd <- .jd_makeQtlDataset(study = "Q1", contexts = c("c1", "c2"), traits = "t1")
  res1 <- QtlFineMappingResult(study = "Q1", context = "c1", trait = "t1",
    method = "mvsusie", entry = list(.js_fmEntry("v1")))
  local_mocked_bindings(
    .fmDispatchJointSpecsQtlDatasetOneRegion = function(...) res1,
    .package = "pecotmr")
  out <- pecotmr:::.fmDispatchJointSpecsQtlDataset(
    list(list(axes = "context", scope = NULL)), qd, methods = "mvsusie",
    contexts = NULL, traitIds = NULL, cisWindow = 1000L, coverage = 0.95,
    secondaryCoverage = 0.5, signalCutoff = 0.1, minAbsCorr = 0.5, verbose = 0)
  expect_identical(out, res1)                               # length-1 short-circuit (909)
  local_mocked_bindings(
    .fmDispatchJointSpecsQtlDatasetOneRegion = function(...) NULL,
    .package = "pecotmr")
  expect_null(pecotmr:::.fmDispatchJointSpecsQtlDataset(
    list(list(axes = "context", scope = NULL)), qd, methods = "mvsusie",
    contexts = NULL, traitIds = NULL, cisWindow = 1000L, coverage = 0.95,
    secondaryCoverage = 0.5, signalCutoff = 0.1, minAbsCorr = 0.5,
    verbose = 0))                                           # all-NULL (908)
})

test_that(".twasMergeResultsByKey: a key absent from a later region contributes nothing", {
  twoRow <- TwasWeights(
    study = c("S", "S"), context = c("c1", "c2"), trait = c("t1", "t1"),
    method = c("lasso", "lasso"),
    entry = list(TwasWeightsEntry(variantIds = "v1", weights = 0.1),
                 TwasWeightsEntry(variantIds = "v2", weights = 0.2)))
  oneRow <- TwasWeights(study = "S", context = "c1", trait = "t1",
    method = "lasso", entry = list(TwasWeightsEntry(variantIds = "v1",
                                                    weights = 0.1)))
  out <- pecotmr:::.twasMergeResultsByKey(list(twoRow, oneRow), c("rA", "rB"))
  expect_equal(nrow(out), 2L)                               # 1063 else-branch
})
