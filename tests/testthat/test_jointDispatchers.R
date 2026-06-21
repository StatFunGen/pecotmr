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
      trimmedFit = list(method = method, payload = fit),
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
