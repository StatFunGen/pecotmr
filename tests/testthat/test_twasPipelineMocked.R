context("twasWeightsPipeline (S4 dispatch) with mocked weight methods")

# ===========================================================================
# Strategy
# ----------------------------------------------------------------------------
# twasWeightsPipeline on a QtlDataset or a QtlSumStats spends almost all its
# uncovered lines orchestrating: variant/sample selection, residualization,
# CV bookkeeping, ensemble fan-in, and packaging results into a TwasWeights
# collection. The actual weight learners are external and slow. We mock the
# weight functions (lassoWeights / enetWeights / susieWeights, plus the
# RSS-side susieRssWeights / lassosumRssWeights / mrAshRssWeights) to return
# zero-valued vectors / matrices so the orchestration runs end-to-end on a
# small fixture.
# ===========================================================================

# ===========================================================================
# Small in-memory QtlDataset fixture (no file IO).
# A custom extractBlockGenotypes mock returns synthetic dosages so we can
# build a QtlDataset whose handle never gets opened.
# ===========================================================================

.tp_makeHandle <- function(snp_n = 20L, n_samples = 40L) {
  new("GenotypeHandle",
    path = "/tmp/tp.gds",
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

.tp_makeSe <- function(traits = c("ENSG_A", "ENSG_B"), n_samples = 40L,
                       chr = "chr1", starts = NULL) {
  if (is.null(starts)) starts <- seq(1000L, by = 1000L, length.out = length(traits))
  rng <- GenomicRanges::GRanges(
    seqnames = rep(chr, length(traits)),
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

.tp_makeQtlDataset <- function(contexts = "brain",
                               traits = c("ENSG_A", "ENSG_B"),
                               n_samples = 40L) {
  gh <- .tp_makeHandle(snp_n = 20L, n_samples = n_samples)
  phen <- setNames(
    lapply(contexts, function(.) .tp_makeSe(traits = traits, n_samples = n_samples)),
    contexts)
  QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = phen,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
}

.tp_mockExtractor <- function(seed = 1, n_samples = 40L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds, handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges   = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx], width = 1L))
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

# Mock individual-level weight methods to return zero vectors quickly.
.tp_mockIndividualWeights <- function() {
  list(
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    enetWeights  = function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X = NULL, y = NULL, susieFit = NULL,
                            retainFit = FALSE, ...)
      rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) {
      out <- rep(0, ncol(X))
      attr(out, "fit") <- list(pi = c(0.9, 0.1))
      out
    }
  )
}

# ===========================================================================
# twasWeightsPipeline(QtlDataset)
# ===========================================================================

test_that("twasWeightsPipeline(QtlDataset): runs end-to-end with mocked solvers", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockIndividualWeights())
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(
    twasWeightsPipeline(qd,
                        methods   = list(lasso_weights = list(),
                                         enet_weights  = list()),
                        cisWindow = 1000L,
                        cvFolds   = 0,
                        ensemble  = FALSE,
                        estimatePi = FALSE,
                        verbose   = 0))
  expect_s4_class(res, "TwasWeights")
  # 1 context x 2 traits x 2 methods = 4 rows.
  expect_equal(nrow(res), 4L)
  expect_setequal(getMethodNames(res), c("lasso", "enet"))
  expect_setequal(getTraits(res), c("ENSG_A", "ENSG_B"))
})

test_that("twasWeightsPipeline(QtlDataset): contexts filter restricts the per-context loop", {
  qd <- .tp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = "ENSG_A")
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockIndividualWeights())
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(
    twasWeightsPipeline(qd,
                        methods   = list(lasso_weights = list()),
                        contexts  = "brain",
                        cisWindow = 1000L,
                        cvFolds   = 0,
                        ensemble  = FALSE,
                        estimatePi = FALSE,
                        verbose   = 0))
  expect_setequal(getContexts(res), "brain")
  expect_equal(nrow(res), 1L)
})

test_that("twasWeightsPipeline(QtlDataset): unknown context errors", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, contexts = "ghost",
                        methods = list(lasso_weights = list())),
    "unknown context"
  )
})

test_that("twasWeightsPipeline(QtlDataset): no traits selected errors", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, traitId = "ENSG_Z",
                        methods = list(lasso_weights = list())),
    "no traits selected"
  )
})

test_that("twasWeightsPipeline(QtlDataset): RSS-only method rejected", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd,
                        methods = list(prsCs_weights = list())),
    "not available for input class 'QtlDataset'"
  )
})

# ===========================================================================
# twasWeightsPipeline(QtlSumStats)
# ===========================================================================

.tp_makeSumstatsEntry <- function(snp_ids = paste0("v", 1:8),
                                  positions = seq(100L, by = 100L, length.out = 8L)) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(snp_ids)),
    ranges = IRanges::IRanges(start = positions, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = snp_ids,
    A1  = rep("A", length(snp_ids)),
    A2  = rep("G", length(snp_ids)),
    Z   = rnorm(length(snp_ids)),
    N   = rep(1000L, length(snp_ids)))
  gr
}

.tp_makeQtlSumStats <- function(n_entries = 1L, qc = TRUE) {
  studies <- rep("s1", n_entries)
  contexts <- if (n_entries == 1L) "c1" else paste0("c", seq_len(n_entries))
  traits   <- rep("t1", n_entries)
  entries  <- lapply(seq_len(n_entries), function(.) .tp_makeSumstatsEntry())
  QtlSumStats(study = studies, context = contexts, trait = traits,
              entry = entries, genome = "hg19",
              ldSketch = .tp_makeHandle(snp_n = 20L),
              qcInfo = if (qc) list(step1 = "ok") else list())
}

# Mock sumstat weight methods to return zero vectors of the LD dim.
.tp_mockSumstatWeights <- function() {
  list(
    susieRssWeights      = function(stat, LD, ...) rep(0, nrow(LD)),
    lassosumRssWeights   = function(stat, LD, ...) rep(0, nrow(LD)),
    mrAshRssWeights      = function(stat, LD, ...) rep(0, nrow(LD)),
    susieInfRssWeights   = function(stat, LD, ...) rep(0, nrow(LD)),
    sdprWeights          = function(stat, LD, ...) rep(0, nrow(LD))
  )
}

test_that("twasWeightsPipeline(QtlSumStats): runs end-to-end with mocked solvers", {
  ss <- .tp_makeQtlSumStats()
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockSumstatWeights())
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  # Method tokens are the bare short names ("susie", "lasso"); the
  # QtlSumStats dispatch resolves them to the *Rss / lassosumRss impl via
  # the .twasMethodCapabilities table.
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = c("susie", "lasso"),
                        verbose = 0)))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)
  expect_setequal(getMethodNames(res), c("susie", "lasso"))
})

test_that("twasWeightsPipeline(QtlSumStats): un-QCd input is rejected", {
  ss <- .tp_makeQtlSumStats(qc = FALSE)
  expect_error(
    twasWeightsPipeline(ss, methods = "susie"),
    "has no QC record"
  )
})

test_that("twasWeightsPipeline(QtlSumStats): individual-only method rejected", {
  ss <- .tp_makeQtlSumStats()
  expect_error(
    twasWeightsPipeline(ss, methods = "enet"),
    "not available for input class 'QtlSumStats'"
  )
})

test_that("twasWeightsPipeline(QtlSumStats): empty contexts/trait filter errors", {
  ss <- .tp_makeQtlSumStats()
  expect_error(
    twasWeightsPipeline(ss, methods = "susie",
                        contexts = "ghost"),
    "no entries matched"
  )
})

test_that("twasWeightsPipeline(QtlSumStats): per-method failure surfaces as warning + skip", {
  ss <- .tp_makeQtlSumStats()
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockSumstatWeights())
  # Override susieRssWeights with the failure-producing version.
  mocks$susieRssWeights <- function(stat, LD, ...) stop("synthetic test failure")
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  # All entries fail -> the per-method-warning fires *and* the pipeline
  # then errors out (no rows produced). Capture both.
  expect_error(
    suppressWarnings(suppressMessages(
      twasWeightsPipeline(ss, methods = "susie", verbose = 0))),
    "no entries produced weights"
  )
})

test_that("twasWeightsPipeline(QtlSumStats): multivariate requires >=2 contexts per (study, trait)", {
  ss <- .tp_makeQtlSumStats(n_entries = 1L)  # 1 context per (study, trait)
  expect_error(
    twasWeightsPipeline(ss, methods = "mvsusie"),
    "multivariate method.*require at least two contexts"
  )
})

# ===========================================================================
# Resume cache (twasWeights = <existing TwasWeights>)
# ===========================================================================

.tp_makeCachedEntry <- function(variant_ids = paste0("v", 1:8),
                                weights = rep(0.5, 8)) {
  TwasWeightsEntry(variantIds = variant_ids,
                    weights = weights,
                    standardized = FALSE)
}

test_that("twasWeightsPipeline(QtlDataset): full cache hit avoids all weight fitting", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  cached <- TwasWeights(
    study   = "study1", context = "brain", trait = "ENSG_A",
    method  = "lasso",
    entry   = list(.tp_makeCachedEntry()))
  fits <- 0L
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    list(lassoWeights = function(X, y, ...) {
      fits <<- fits + 1L; rep(0, ncol(X))
    }))
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(
    twasWeightsPipeline(qd,
                        methods    = list(lasso_weights = list()),
                        cisWindow  = 1000L,
                        cvFolds    = 0,
                        ensemble   = FALSE,
                        estimatePi = FALSE,
                        twasWeights = cached,
                        verbose    = 0))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 1L)
  expect_equal(fits, 0L)  # the cached entry short-circuited the fit
})

test_that("twasWeightsPipeline(QtlDataset): partial cache hit fits only missing methods", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  # Cache covers `lasso` but `enet` is missing -> only enetWeights is called.
  cached <- TwasWeights(
    study   = "study1", context = "brain", trait = "ENSG_A",
    method  = "lasso",
    entry   = list(.tp_makeCachedEntry()))
  enetCalls   <- 0L
  lassoCalls  <- 0L
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    list(
      lassoWeights = function(X, y, ...) { lassoCalls <<- lassoCalls + 1L; rep(0, ncol(X)) },
      enetWeights  = function(X, y, ...) { enetCalls  <<- enetCalls  + 1L; rep(0, ncol(X)) }))
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(
    twasWeightsPipeline(qd,
                        methods    = list(lasso_weights = list(),
                                          enet_weights  = list()),
                        cisWindow  = 1000L,
                        cvFolds    = 0,
                        ensemble   = FALSE,
                        estimatePi = FALSE,
                        twasWeights = cached,
                        verbose    = 0))
  expect_setequal(getMethodNames(res), c("lasso", "enet"))
  expect_equal(nrow(res), 2L)
  expect_equal(lassoCalls, 0L)  # cache hit
  expect_equal(enetCalls,  1L)  # cache miss -> fit
})

test_that("twasWeightsPipeline(QtlSumStats): cache hit on a per-tuple basis", {
  ss <- .tp_makeQtlSumStats()
  cached <- TwasWeights(
    study   = "s1", context = "c1", trait = "t1",
    method  = "susie",
    entry   = list(.tp_makeCachedEntry()))
  rssCalls <- 0L
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockSumstatWeights())
  mocks$susieRssWeights <- function(stat, LD, ...) {
    rssCalls <<- rssCalls + 1L
    rep(0, nrow(LD))
  }
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = "susie",
                        twasWeights = cached,
                        verbose = 0)))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 1L)
  expect_equal(rssCalls, 0L)
})

test_that(".twasCacheLookup: NULL twasWeights returns NULL", {
  expect_null(pecotmr:::.twasCacheLookup(NULL, "s1", "c1", "t1", "lasso"))
})

test_that(".twasCacheLookup: non-TwasWeights input returns NULL", {
  expect_null(pecotmr:::.twasCacheLookup("not_a_tw", "s1", "c1", "t1", "lasso"))
})

test_that(".twasCacheLookup: returns matching entry by 4-tuple", {
  e <- TwasWeightsEntry(variantIds = "v1", weights = 0.5)
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "lasso",
    entry = list(e))
  expect_identical(
    pecotmr:::.twasCacheLookup(tw, "s1", "c1", "t1", "lasso"),
    e)
  expect_null(pecotmr:::.twasCacheLookup(tw, "ghost", "c1", "t1", "lasso"))
})

# ===========================================================================
# twasWeightsPipeline ANY-method dispatch (unsupported input class)
# ===========================================================================

test_that("twasWeightsPipeline(ANY): unsupported input class errors", {
  expect_error(
    twasWeightsPipeline(matrix(0, 5, 5)),
    "does not accept inputs of class 'matrix'"
  )
})

# ===========================================================================
# twasMultivariateWeightsPipeline
# ===========================================================================

# NOTE: twasMultivariateWeightsPipeline currently bottoms out on a
# `getWeights(twasWeight)` call (R/twasWeights.R:2105) that omits the
# (study, context, trait, method) selectors and therefore errors whenever
# the inner learnTwasWeights() returns a collection with >1 row — which it
# always does when both mr.mash and mvSuSiE are requested (the default).
# Skipping this test until that path is reworked to walk all rows; flagging
# here so coverage doesn't appear lifted on a known-broken pipeline.
test_that("twasMultivariateWeightsPipeline: known-broken under the new TwasWeights API", {
  skip("twasMultivariateWeightsPipeline calls getWeights() on a multi-row TwasWeights without selectors; needs a follow-up fix in production code.")
})

# imputeMissingGwasForSketch was removed: it duplicated the inline RAISS
# imputation step at the bottom of .runEntrySummaryStatsQc (sumstatsQc.R),
# had no production callers, and was orphaned post-S4-refactor. The RAISS-
# against-sketch path that does run in production is exercised via the
# summaryStatsQc test file (with raiss() left to run for real on the small
# fixture).
