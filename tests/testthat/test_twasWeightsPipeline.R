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

# Build a minimal FineMappingResult that satisfies .twasCheckFineMappingMethods.
# Each entry stores a stub fit list under the susieFit slot; the per-tuple
# loop pulls these via .twasFineMappingFits() and threads them into the
# corresponding *Weights wrapper as susieFit / susieInfFit / mvsusieFit / etc.
.tp_makeStubFineMappingResult <- function(study   = "study1",
                                          contexts = "brain",
                                          traits   = "ENSG_A",
                                          method   = "susie",
                                          fitPayload = NULL) {
  rows <- expand.grid(study = study, context = contexts, trait = traits,
                      method = method, stringsAsFactors = FALSE)
  entries <- lapply(seq_len(nrow(rows)), function(i) {
    if (is.null(fitPayload))
      fitPayload <- list(method = rows$method[[i]],
                          context = rows$context[[i]],
                          trait = rows$trait[[i]])
    tl <- data.frame(
      variant_id = paste0("v", seq_len(3L)),
      pip        = c(0.9, 0.5, 0.1),
      stringsAsFactors = FALSE)
    FineMappingEntry(
      variantIds = tl$variant_id,
      susieFit   = fitPayload,
      topLoci    = tl)
  })
  if ("brain" %in% rows$context || any(rows$context != "")) {
    QtlFineMappingResult(
      study   = rows$study,
      context = rows$context,
      trait   = rows$trait,
      method  = rows$method,
      entry   = entries)
  } else {
    GwasFineMappingResult(
      study  = rows$study,
      method = rows$method,
      entry  = entries)
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
  # Method tokens are the bare short names; the QtlSumStats dispatch
  # resolves them to the *Rss impl via the .twasMethodCapabilities table.
  # Fine-mapping methods (susie / susieInf / etc.) require a
  # FineMappingResult and are covered by separate tests.
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = c("mrash", "lasso"),
                        verbose = 0)))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)
  expect_setequal(getMethodNames(res), c("mrash", "lasso"))
})

test_that("twasWeightsPipeline(QtlSumStats): un-QCd input is rejected", {
  ss <- .tp_makeQtlSumStats(qc = FALSE)
  expect_error(
    twasWeightsPipeline(ss, methods = "lasso"),
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
    twasWeightsPipeline(ss, methods = "lasso",
                        contexts = "ghost"),
    "no entries matched"
  )
})

test_that("twasWeightsPipeline(QtlSumStats): per-method failure surfaces as warning + skip", {
  ss <- .tp_makeQtlSumStats()
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockSumstatWeights())
  # Override lassosumRssWeights with the failure-producing version.
  mocks$lassosumRssWeights <- function(stat, LD, ...) stop("synthetic test failure")
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  # All entries fail -> the per-method-warning fires *and* the pipeline
  # then errors out (no rows produced). Capture both.
  expect_error(
    suppressWarnings(suppressMessages(
      twasWeightsPipeline(ss, methods = "lasso", verbose = 0))),
    "no entries produced weights"
  )
})

test_that("twasWeightsPipeline(QtlSumStats): multivariate requires >=2 contexts per (study, trait)", {
  ss <- .tp_makeQtlSumStats(n_entries = 1L)  # 1 context per (study, trait)
  # Provide a stub FineMappingResult so the gate passes; the multivariate
  # guard is what we want to exercise here.
  fmr <- .tp_makeStubFineMappingResult(study = "s1", contexts = "c1",
                                        traits = "t1", method = "mvsusie")
  expect_error(
    twasWeightsPipeline(ss, methods = "mvsusie", fineMappingResult = fmr),
    "multivariate method.*require at least two contexts"
  )
})

# ===========================================================================
# Fine-mapping method gate: every fine-mapping method (susie / susieInf /
# susieAsh / mvsusie / fsusie) must be paired with a FineMappingResult.
# twasWeightsPipeline is not allowed to re-fit them from scratch. Input-
# class compatibility is delegated to .fmCheckMethodCapabilities so the
# rule set stays in sync with fineMappingPipeline.
# ===========================================================================

test_that("gate: QtlDataset + susie without fineMappingResult errors", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, methods = "susie"),
    "are fine-mapping methods and may not be re-fit")
})

test_that("gate: QtlSumStats + susieInf without fineMappingResult errors", {
  ss <- .tp_makeQtlSumStats()
  expect_error(
    twasWeightsPipeline(ss, methods = "susieInf"),
    "are fine-mapping methods and may not be re-fit")
})

test_that("gate: QtlDataset + susieAsh without fineMappingResult errors", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, methods = "susieAsh"),
    "are fine-mapping methods and may not be re-fit")
})

test_that("gate: composite (susie + susieInf) without fineMappingResult errors", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  # Mixed: lasso is allowed, but susie/susieInf must still be gated.
  expect_error(
    twasWeightsPipeline(qd, methods = c("susie", "susieInf", "lasso")),
    "are fine-mapping methods and may not be re-fit")
})

test_that("gate: mvsusie without fineMappingResult errors (QtlDataset)", {
  qd <- .tp_makeQtlDataset(contexts = c("brain", "liver"), traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, methods = "mvsusie"),
    "are fine-mapping methods and may not be re-fit")
})

test_that("gate: fsusie has no TWAS-weight extractor (rejected by name)", {
  qd <- .tp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = c("ENSG_A", "ENSG_B"))
  fmr <- .tp_makeStubFineMappingResult(
    study = "study1", contexts = c("brain", "liver"),
    traits = "ENSG_A", method = "fsusie")
  # Even with a fineMappingResult, fsusie can't produce TWAS weights.
  expect_error(
    twasWeightsPipeline(qd, methods = "fsusie", fineMappingResult = fmr),
    "have no TWAS-weight extractor")
})

test_that("gate: fsusie on QtlSumStats delegates to .fmCheckMethodCapabilities", {
  ss <- .tp_makeQtlSumStats()
  fmr <- .tp_makeStubFineMappingResult(study = "s1", contexts = "c1",
                                        traits = "t1", method = "fsusie")
  # fineMappingPipeline rejects fsusie on QtlSumStats ("sumstat-only on this
  # pipeline" / no sumstat impl); twasWeightsPipeline reuses that check.
  expect_error(
    twasWeightsPipeline(ss, methods = "fsusie", fineMappingResult = fmr),
    "individual-only|sumstat-only|not supported")
})

test_that("gate: non-FineMappingResult object passed in errors", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, methods = "susie",
                        fineMappingResult = list()),
    "must be a FineMappingResult")
})

test_that("gate: unknown method tokens still error with full menu", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    twasWeightsPipeline(qd, methods = "totallyMadeUpMethod"),
    "Unknown TWAS method|unknown method")
})

# ---------------------------------------------------------------------------
# Success: when a FineMappingResult is supplied, the fit is threaded into
# the corresponding *Weights wrapper via its *Fit argument. The wrappers
# are mocked to verify the fit arrives and the underlying fitter is NOT
# re-invoked.
# ---------------------------------------------------------------------------

test_that("gate: QtlDataset + susie + fineMappingResult threads the susieFit", {
  qd <- .tp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  fmr <- .tp_makeStubFineMappingResult(
    study = "study1", contexts = "brain", traits = "ENSG_A",
    method = "susie")
  sawFit <- FALSE
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    susieWeights = function(X = NULL, y = NULL, susieFit = NULL, ...) {
      sawFit <<- !is.null(susieFit)
      rep(0, ncol(X))
    })
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(qd,
                        methods           = "susie",
                        fineMappingResult = fmr,
                        cisWindow         = 1000L,
                        cvFolds           = 0,
                        ensemble          = FALSE,
                        estimatePi        = FALSE,
                        verbose           = 0)))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 1L)
  expect_setequal(getMethodNames(res), "susie")
  expect_true(sawFit)
})

test_that("gate: QtlSumStats + susie + fineMappingResult threads the susieRssFit", {
  ss <- .tp_makeQtlSumStats()
  fmr <- .tp_makeStubFineMappingResult(study = "s1", contexts = "c1",
                                        traits = "t1", method = "susie")
  fitsSeen <- 0L
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    susieRssWeights = function(stat, LD, susieRssFit = NULL, ...) {
      if (!is.null(susieRssFit)) fitsSeen <<- fitsSeen + 1L
      rep(0, nrow(LD))
    })
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = "susie",
                        fineMappingResult = fmr, verbose = 0)))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 1L)
  expect_setequal(getMethodNames(res), "susie")
  expect_equal(fitsSeen, 1L)
})

test_that("gate: QtlSumStats + susieAsh + fineMappingResult threads the susieAshRssFit", {
  ss <- .tp_makeQtlSumStats()
  fmr <- .tp_makeStubFineMappingResult(study = "s1", contexts = "c1",
                                        traits = "t1", method = "susieAsh")
  sawFit <- FALSE
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    susieAshRssWeights = function(stat, LD, susieAshRssFit = NULL, ...) {
      sawFit <<- !is.null(susieAshRssFit)
      rep(0, nrow(LD))
    })
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = "susieAsh",
                        fineMappingResult = fmr, verbose = 0)))
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 1L)
  expect_setequal(getMethodNames(res), "susieAsh")
  expect_true(sawFit)
})

test_that("gate: missing matching tuple in fineMappingResult warns and skips", {
  ss <- .tp_makeQtlSumStats()
  # FineMappingResult has a row for a DIFFERENT (study, context, trait):
  # the gate passes (any fine-mapping fit present satisfies the object
  # check), but the per-tuple lookup finds nothing and warns.
  fmr <- .tp_makeStubFineMappingResult(study = "other_study",
                                        contexts = "other_ctx",
                                        traits = "other_trait",
                                        method = "susie")
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockSumstatWeights())
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  expect_error(
    suppressWarnings(suppressMessages(
      twasWeightsPipeline(ss, methods = "susie",
                          fineMappingResult = fmr, verbose = 0))),
    "no entries produced weights")
})

# ===========================================================================
# Multivariate dispatch paths (mvsusie / mr.mash) on individual-level data
# (QtlDataset.runMultivariate) and on summary-statistics input (QtlSumStats
# multivariate-Z dispatch). Both paths build a multi-column Y / Z and invoke
# a single multivariate solver that returns a (variants x conditions) weight
# matrix; we mock that solver to return a zero matrix of the expected shape.
# ===========================================================================

test_that("twasWeightsPipeline(QtlDataset): mvsusie multivariate path returns one row per (context, trait)", {
  qd <- .tp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = "ENSG_A")
  fmr <- .tp_makeStubFineMappingResult(
    study = "study1", contexts = c("brain", "liver"),
    traits = "ENSG_A", method = "mvsusie")
  mvCalls <- 0L
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    mvsusieWeights = function(X, Y, mvsusieFit = NULL, ...) {
      mvCalls <<- mvCalls + 1L
      # Verify the gate threaded the precomputed fit through.
      stopifnot(!is.null(mvsusieFit))
      matrix(0, nrow = ncol(X), ncol = ncol(Y),
             dimnames = list(colnames(X), colnames(Y)))
    })
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(qd,
                        methods           = "mvsusie",
                        fineMappingResult = fmr,
                        cisWindow         = 1000L,
                        cvFolds           = 0,
                        ensemble          = FALSE,
                        estimatePi        = FALSE,
                        verbose           = 0)))
  expect_s4_class(res, "TwasWeights")
  # Joint fit on 2 contexts x 1 trait -> 2 rows back, one per (context, trait).
  expect_equal(nrow(res), 2L)
  expect_setequal(getContexts(res), c("brain", "liver"))
  expect_setequal(getTraits(res), "ENSG_A")
  expect_setequal(getMethodNames(res), "mvsusie")
  expect_gte(mvCalls, 1L)
})

test_that("twasWeightsPipeline(QtlDataset): mr.mash multivariate path with 2 traits x 2 contexts", {
  qd <- .tp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = c("ENSG_A", "ENSG_B"))
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    mrmashWeights = function(X, Y, ...) {
      matrix(0, nrow = ncol(X), ncol = ncol(Y),
             dimnames = list(colnames(X), colnames(Y)))
    })
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(qd,
                        methods    = "mrmash",
                        cisWindow  = 1000L,
                        cvFolds    = 0,
                        ensemble   = FALSE,
                        estimatePi = FALSE,
                        verbose    = 0)))
  expect_s4_class(res, "TwasWeights")
  # 2 contexts x 2 traits = 4 output rows.
  expect_equal(nrow(res), 4L)
  expect_setequal(getContexts(res), c("brain", "liver"))
  expect_setequal(getTraits(res), c("ENSG_A", "ENSG_B"))
  expect_setequal(getMethodNames(res), "mrmash")
})

# ---------------------------------------------------------------------------
# QtlSumStats: multivariate dispatch builds a (variants x contexts) Z matrix
# and invokes the *Rss solver once per (study, trait) group.
# ---------------------------------------------------------------------------

.tp_makeMultiCtxQtlSumStats <- function(contexts = c("c1", "c2"),
                                         snp_ids  = paste0("v", 1:8),
                                         positions = seq(100L, by = 100L, length.out = 8L)) {
  n <- length(contexts)
  entries <- lapply(seq_len(n), function(i) {
    # Same SNP order across contexts -- required by the multivariate path
    # (it errors on any divergence after summaryStatsQc).
    .tp_makeSumstatsEntry(snp_ids = snp_ids, positions = positions)
  })
  QtlSumStats(study   = rep("s1", n),
              context = contexts,
              trait   = rep("t1", n),
              entry   = entries,
              genome  = "hg19",
              ldSketch = .tp_makeHandle(snp_n = 20L),
              qcInfo   = list(step1 = "ok"))
}

test_that("twasWeightsPipeline(QtlSumStats): mvsusie multivariate path returns one row per context", {
  ss <- .tp_makeMultiCtxQtlSumStats(contexts = c("c1", "c2"))
  fmr <- .tp_makeStubFineMappingResult(
    study = "s1", contexts = c("c1", "c2"), traits = "t1",
    method = "mvsusie")
  shapesSeen <- list()
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    mvsusieRssWeights = function(stat, LD, mvsusieRssFit = NULL, ...) {
      shapesSeen[[length(shapesSeen) + 1L]] <<- list(
        Zdim = dim(stat$z), LDdim = dim(LD),
        ctxNames = colnames(stat$z),
        sawFit = !is.null(mvsusieRssFit))
      matrix(0, nrow = nrow(LD), ncol = ncol(stat$z),
             dimnames = list(stat$variantNames, colnames(stat$z)))
    })
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = "mvsusie",
                        fineMappingResult = fmr, verbose = 0)))
  expect_s4_class(res, "TwasWeights")
  # 1 (study, trait) group with 2 contexts -> 2 output rows.
  expect_equal(nrow(res), 2L)
  expect_setequal(getContexts(res), c("c1", "c2"))
  expect_setequal(getTraits(res), "t1")
  expect_setequal(getMethodNames(res), "mvsusie")
  # The solver was called once with a (variants x 2) Z matrix and was
  # handed the precomputed fit from fineMappingResult.
  expect_length(shapesSeen, 1L)
  expect_equal(shapesSeen[[1L]]$Zdim[[2L]], 2L)
  expect_setequal(shapesSeen[[1L]]$ctxNames, c("c1", "c2"))
  expect_true(shapesSeen[[1L]]$sawFit)
})

test_that("twasWeightsPipeline(QtlSumStats): mr.mash multivariate solver failure surfaces as warning + empty result", {
  ss <- .tp_makeMultiCtxQtlSumStats(contexts = c("c1", "c2"))
  mocks <- list(
    extractBlockGenotypes = .tp_mockExtractor(),
    mrmashRssWeights = function(stat, LD, ...) stop("synthetic multivariate failure"))
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  # All multivariate fits fail -> no rows -> the pipeline errors out at the
  # end, surfacing the per-group warning along the way.
  expect_error(
    suppressWarnings(suppressMessages(
      twasWeightsPipeline(ss, methods = "mrmash", verbose = 0))),
    "no entries produced weights"
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
    method  = "lasso",
    entry   = list(.tp_makeCachedEntry()))
  rssCalls <- 0L
  mocks <- c(
    list(extractBlockGenotypes = .tp_mockExtractor()),
    .tp_mockSumstatWeights())
  mocks$lassosumRssWeights <- function(stat, LD, ...) {
    rssCalls <<- rssCalls + 1L
    rep(0, nrow(LD))
  }
  do.call(local_mocked_bindings,
          c(mocks, list(.package = "pecotmr")))
  res <- suppressMessages(suppressWarnings(
    twasWeightsPipeline(ss, methods = "lasso",
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


context("ensembleWeights")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Build a synthetic twasWeightsCv() output with K methods. Each method's
# prediction is a convex combination of the truth + noise, letting us control
# per-method accuracy. Returns a list shaped exactly like twasWeightsCv()'s
# output (with $prediction, $performance, $samplePartition).
make_cv_result <- function(n = 100, K = 4, seed = 1, method_quality = NULL) {
  set.seed(seed)
  y <- rnorm(n)
  sampleNames <- paste0("sample_", seq_len(n))

  if (is.null(method_quality)) {
    # Methods with decreasing quality (noise amounts)
    method_quality <- seq(0.1, 0.9, length.out = K)
  }
  stopifnot(length(method_quality) == K)

  method_names <- paste0("method", seq_len(K))
  pred_names <- paste0(method_names, "_predicted")

  prediction <- setNames(lapply(seq_len(K), function(k) {
    noise_sd <- method_quality[k]
    pred <- y + rnorm(n, sd = noise_sd)
    mat <- matrix(pred, ncol = 1)
    rownames(mat) <- sampleNames
    colnames(mat) <- "outcome_1"
    mat
  }), pred_names)

  # Dummy performance (not used by ensembleWeights)
  performance <- setNames(lapply(seq_len(K), function(k) {
    m <- matrix(NA, nrow = 1, ncol = 6)
    colnames(m) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    m
  }), paste0(method_names, "_performance"))

  list(
    samplePartition = data.frame(Sample = sampleNames,
                                   Fold = rep(1:5, length.out = n),
                                   stringsAsFactors = FALSE),
    prediction = prediction,
    performance = performance,
    timeElapsed = 0,
    .y = y,
    .method_names = method_names
  )
}

# Build synthetic learnTwasWeights() output
make_weight_list <- function(p = 20, method_names, seed = 2) {
  set.seed(seed)
  setNames(lapply(method_names, function(m) {
    w <- matrix(rnorm(p), ncol = 1)
    rownames(w) <- paste0("var_", seq_len(p))
    colnames(w) <- "outcome_1"
    w
  }), paste0(method_names, "_weights"))
}

# ===========================================================================
#  Input validation
# ===========================================================================

test_that("ensembleWeights: NULL cv_results errors", {
  expect_error(ensembleWeights(NULL, Y = rnorm(10)), "cvResults")
})

test_that("ensembleWeights: NULL Y errors", {
  cv <- make_cv_result(n = 20, K = 3)
  expect_error(ensembleWeights(cv, Y = NULL), "'Y' is required")
})

test_that("ensembleWeights: single method errors (need >= 2 for ensemble)", {
  cv <- make_cv_result(n = 20, K = 1)
  expect_error(ensembleWeights(cv, Y = cv$.y),
               "at least 2 methods")
})

test_that("ensembleWeights: invalid context_index errors", {
  cv <- make_cv_result(n = 20, K = 3)
  expect_error(ensembleWeights(cv, Y = cv$.y, contextIndex = 0),
               "contextIndex")
  expect_error(ensembleWeights(cv, Y = cv$.y, contextIndex = "a"),
               "contextIndex")
})

test_that("ensembleWeights: context_index beyond Y columns errors", {
  cv <- make_cv_result(n = 20, K = 3)
  Y_mat <- matrix(cv$.y, ncol = 1)
  expect_error(ensembleWeights(cv, Y = Y_mat, contextIndex = 5),
               "contextIndex")
})

test_that("ensembleWeights: multi-dataset with mismatched lengths errors", {
  cv1 <- make_cv_result(n = 20, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 20, K = 3, seed = 2)
  expect_error(ensembleWeights(list(cv1, cv2), Y = list(cv1$.y)),
               "same length")
})

test_that("ensembleWeights: multi-dataset with different methods errors", {
  cv1 <- make_cv_result(n = 20, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 20, K = 4, seed = 2)
  expect_error(
    ensembleWeights(list(cv1, cv2), Y = list(cv1$.y, cv2$.y)),
    "same method names"
  )
})

# ===========================================================================
#  Core algorithm correctness
# ===========================================================================

test_that("ensembleWeights: coefficients are non-negative and sum to 1", {
  cv <- make_cv_result(n = 100, K = 4, seed = 42)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_true(all(res$methodCoef >= 0))
  expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: best method receives the largest coefficient", {
  # Method 1 is best (lowest noise), method K is worst
  cv <- make_cv_result(n = 200, K = 4, seed = 7,
                        method_quality = c(0.1, 0.5, 0.8, 1.2))
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(names(which.max(res$methodCoef)), "method1")
})

test_that("ensembleWeights: does not return ensemble_performance (in-sample R^2 omitted)", {
  cv <- make_cv_result(n = 300, K = 5, seed = 13)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_null(res$ensemble_performance)
  expect_false("ensemble_performance" %in% names(res))
})

test_that("ensembleWeights: per-method R^2 values are sensible (between 0 and 1)", {
  cv <- make_cv_result(n = 200, K = 4, seed = 21)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_true(all(res$methodPerformance >= 0, na.rm = TRUE))
  expect_true(all(res$methodPerformance <= 1, na.rm = TRUE))
  expect_equal(length(res$methodPerformance), 4)
})

test_that("ensembleWeights: method names are stripped of _predicted suffix", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(names(res$methodCoef),
               c("method1", "method2", "method3"))
  expect_equal(names(res$methodPerformance),
               c("method1", "method2", "method3"))
})

# ===========================================================================
#  Sample name alignment
# ===========================================================================

test_that("ensembleWeights: aligns Y and predictions by sample name", {
  cv <- make_cv_result(n = 50, K = 3, seed = 10)

  # Shuffle Y order relative to predictions
  shuffled_order <- sample(50)
  y_shuffled <- cv$.y[shuffled_order]
  names(y_shuffled) <- paste0("sample_", shuffled_order)

  res_aligned <- ensembleWeights(cv, Y = y_shuffled)
  res_original <- ensembleWeights(cv, Y = cv$.y)

  # Results should be identical regardless of Y order
  expect_equal(res_aligned$methodCoef, res_original$methodCoef, tolerance = 1e-10)
})

test_that("ensembleWeights: aligns Y matrix and predictions by sample name", {
  cv <- make_cv_result(n = 50, K = 3, seed = 10)

  # Create Y as a matrix with shuffled row order
  shuffled_order <- sample(50)
  Y_mat <- matrix(cv$.y[shuffled_order], ncol = 1)
  rownames(Y_mat) <- paste0("sample_", shuffled_order)

  res_aligned <- ensembleWeights(cv, Y = Y_mat)
  res_original <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(res_aligned$methodCoef, res_original$methodCoef, tolerance = 1e-10)
})

test_that("ensembleWeights: errors when no common sample names", {
  cv <- make_cv_result(n = 20, K = 3, seed = 1)
  y_bad <- setNames(rnorm(20), paste0("other_", seq_len(20)))

  expect_error(ensembleWeights(cv, Y = y_bad), "No common sample names")
})

# ===========================================================================
#  Zero-variance / edge cases
# ===========================================================================

test_that("ensembleWeights: zero-variance method gets coefficient 0", {
  cv <- make_cv_result(n = 100, K = 3, seed = 5)
  # Force method 2 to have constant predictions
  cv$prediction$method2_predicted[, 1] <- 0.5
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(res$methodCoef["method2"], c(method2 = 0))
  expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: NA predictions in some samples are dropped", {
  cv <- make_cv_result(n = 100, K = 3, seed = 5)
  cv$prediction$method1_predicted[1:5, 1] <- NA
  expect_message(
    res <- ensembleWeights(cv, Y = cv$.y),
    "Dropping"
  )
  expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: all zero-variance methods errors", {
  cv <- make_cv_result(n = 50, K = 2, seed = 5)
  cv$prediction$method1_predicted[, 1] <- 0
  cv$prediction$method2_predicted[, 1] <- 0
  expect_error(ensembleWeights(cv, Y = cv$.y),
               "zero-variance predictions")
})

# ===========================================================================
#  Weight combination
# ===========================================================================

test_that("ensembleWeights: ensembleTwasWeights is sum of zeta_k * w_k", {
  cv <- make_cv_result(n = 100, K = 3, seed = 42)
  wt <- make_weight_list(p = 10, method_names = cv$.method_names)

  res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = wt)

  expect_false(is.null(res$ensembleTwasWeights))

  # Verify the combination is correct
  expected <- matrix(0, nrow = 10, ncol = 1)
  for (k in seq_along(cv$.method_names)) {
    m <- cv$.method_names[k]
    expected <- expected + res$methodCoef[m] * wt[[paste0(m, "_weights")]]
  }
  expect_equal(as.numeric(res$ensembleTwasWeights),
               as.numeric(expected),
               tolerance = 1e-10)
})

test_that("ensembleWeights: NULL twas_weight_list returns NULL ensembleTwasWeights", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = NULL)
  expect_null(res$ensembleTwasWeights)
})

test_that("ensembleWeights: weights with no matching keys warns and skips", {
  cv <- make_cv_result(n = 50, K = 2, seed = 1)
  wt <- list(unknown_weights = matrix(1, nrow = 10, ncol = 1))

  expect_warning(
    res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = wt),
    "No matching weight keys"
  )
  expect_null(res$ensembleTwasWeights)
})

# ===========================================================================
#  Multi-dataset ensemble
# ===========================================================================

test_that("ensembleWeights: multi-dataset combines predictions correctly", {
  cv1 <- make_cv_result(n = 80, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 80, K = 3, seed = 2)

  res <- ensembleWeights(
    cvResults = list(cv1, cv2),
    Y = list(cv1$.y, cv2$.y)
  )

  expect_true(all(res$methodCoef >= 0))
  expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
  expect_equal(length(res$methodPerformance), 3)
})

test_that("ensembleWeights: Y as matrix with context_index works", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  Y_mat <- matrix(cv$.y, ncol = 1)
  colnames(Y_mat) <- "ctx1"

  res <- ensembleWeights(cv, Y = Y_mat, contextIndex = 1)
  expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
})

# ===========================================================================
#  End-to-end with twasWeightsCv (integration)
# ===========================================================================

test_that("ensembleWeights: end-to-end with twasWeightsCv output", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  cv <- suppressMessages(twasWeightsCv(
    X, y, fold = 3,
    weightMethods = list(
      lassoWeights = list(),
      enetWeights = list()
    )
  ))

  res <- ensembleWeights(cv, Y = y)

  expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
  expect_true(all(res$methodCoef >= 0))
  expect_equal(names(res$methodCoef), c("lasso", "enet"))
  expect_null(res$ensemble_performance)
})

# ===========================================================================
#  twasWeightsPipeline ensemble integration
# ===========================================================================

test_that("pipeline: ensemble=TRUE with only 1 method prints skip message", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list()),
      ensemble = TRUE
    )
  )

  # Should see the skip message
  expect_true(any(grepl("Ensemble model skipped.*only 1 weight method provided", msgs)))

  # No ensemble result should be present
  expect_null(res$ensemble)
  expect_false("ensemble" %in% getMethodNames(res$twasWeights))
})

test_that("pipeline: ensemble=TRUE skips when methods fail R^2 cutoff", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  # Use signal so methods produce non-zero weights, but set threshold very high
  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleR2Threshold = 0.99  # impossibly high threshold
    )
  )

  expect_true(any(grepl("Ensemble TWAS skipped", msgs)))
  expect_null(res$ensemble)
  expect_false("ensemble" %in% getMethodNames(res$twasWeights))
})

test_that("pipeline: ensemble=TRUE succeeds and adds ensembleWeights", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))

  # Ensemble weights added alongside individual methods
  expect_true("ensemble" %in% getMethodNames(res$twasWeights))
  expect_true("lasso" %in% getMethodNames(res$twasWeights))
  expect_true("enet" %in% getMethodNames(res$twasWeights))

  # Ensemble predictions added
  expect_true("ensemble_predicted" %in% names(res$twasPredictions))

  # Ensemble result metadata present
  expect_false(is.null(res$ensemble))
  expect_true(all(res$ensemble$methodCoef >= 0))
  expect_equal(sum(res$ensemble$methodCoef), 1, tolerance = 1e-6)

  # Ensemble weights should have same length as individual weights
  expect_equal(length(getWeights(res$twasWeights,
                                 study = "", context = "", trait = "",
                                 method = "ensemble")),
               length(getWeights(res$twasWeights,
                                 study = "", context = "", trait = "",
                                 method = "lasso")))
})

test_that("pipeline: ensemble=FALSE does not run ensemble", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  res <- suppressMessages(pecotmr:::.twasWeightsPipelineMatrix(
    X, y, cvFolds = 3,
    weightMethods = list(lassoWeights = list(), enetWeights = list()),
    ensemble = FALSE
  ))

  expect_null(res$ensemble)
  expect_false("ensemble" %in% getMethodNames(res$twasWeights))
})

test_that("pipeline: ensemble_r2_threshold filters methods for ensemble", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  # Run with very low threshold - both methods should pass
  msgs_low <- testthat::capture_messages(
    res_low <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleR2Threshold = 0.001
    )
  )
  expect_false(is.null(res_low$ensemble))

  # Run with very high threshold - neither should pass
  msgs_high <- testthat::capture_messages(
    res_high <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleR2Threshold = 0.99
    )
  )
  expect_true(any(grepl("Ensemble TWAS skipped", msgs_high)))
  expect_null(res_high$ensemble)
})

# ===========================================================================
#  Solver alternatives
# ===========================================================================

for (slv in c("quadprog", "nnls", "lbfgsb", "glmnet")) {
  test_that(paste0("ensembleWeights: solver='", slv, "' produces valid coefficients"), {
    if (slv == "quadprog") skip_if_not_installed("quadprog")
    if (slv == "nnls") skip_if_not_installed("nnls")
    if (slv == "glmnet") skip_if_not_installed("glmnet")

    cv <- make_cv_result(n = 100, K = 4, seed = 42)
    res <- ensembleWeights(cv, Y = cv$.y, solver = slv)

    expect_true(all(res$methodCoef >= 0))
    expect_equal(sum(res$methodCoef), 1, tolerance = 1e-6)
    expect_equal(length(res$methodCoef), 4)
  })

  test_that(paste0("ensembleWeights: solver='", slv, "' assigns best method largest coef"), {
    if (slv == "quadprog") skip_if_not_installed("quadprog")
    if (slv == "nnls") skip_if_not_installed("nnls")
    if (slv == "glmnet") skip_if_not_installed("glmnet")

    cv <- make_cv_result(n = 200, K = 4, seed = 7,
                          method_quality = c(0.1, 0.5, 0.8, 1.2))
    res <- ensembleWeights(cv, Y = cv$.y, solver = slv)

    expect_equal(names(which.max(res$methodCoef)), "method1")
  })

  test_that(paste0("ensembleWeights: solver='", slv, "' combines weights correctly"), {
    if (slv == "quadprog") skip_if_not_installed("quadprog")
    if (slv == "nnls") skip_if_not_installed("nnls")
    if (slv == "glmnet") skip_if_not_installed("glmnet")

    cv <- make_cv_result(n = 100, K = 3, seed = 42)
    wt <- make_weight_list(p = 10, method_names = cv$.method_names)
    res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = wt, solver = slv)

    expect_false(is.null(res$ensembleTwasWeights))

    expected <- matrix(0, nrow = 10, ncol = 1)
    for (k in seq_along(cv$.method_names)) {
      m <- cv$.method_names[k]
      expected <- expected + res$methodCoef[m] * wt[[paste0(m, "_weights")]]
    }
    expect_equal(as.numeric(res$ensembleTwasWeights),
                 as.numeric(expected),
                 tolerance = 1e-10)
  })
}

test_that("ensembleWeights: invalid solver errors", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  expect_error(ensembleWeights(cv, Y = cv$.y, solver = "bogus"),
               "arg")
})

test_that("pipeline: ensemble_solver='nnls' works end-to-end", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("nnls")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleSolver = "nnls"
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))
  expect_true("ensemble" %in% getMethodNames(res$twasWeights))
  expect_true(all(res$ensemble$methodCoef >= 0))
  expect_equal(sum(res$ensemble$methodCoef), 1, tolerance = 1e-6)
})

test_that("pipeline: ensemble_solver='lbfgsb' works end-to-end", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleSolver = "lbfgsb"
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))
  expect_true("ensemble" %in% getMethodNames(res$twasWeights))
  expect_true(all(res$ensemble$methodCoef >= 0))
  expect_equal(sum(res$ensemble$methodCoef), 1, tolerance = 1e-6)
})

test_that("pipeline: ensemble_solver='glmnet' works end-to-end", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- pecotmr:::.twasWeightsPipelineMatrix(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleSolver = "glmnet"
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))
  expect_true("ensemble" %in% getMethodNames(res$twasWeights))
  expect_true(all(res$ensemble$methodCoef >= 0))
  expect_equal(sum(res$ensemble$methodCoef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: solver='glmnet' respects alpha parameter", {
  skip_if_not_installed("glmnet")

  cv <- make_cv_result(n = 200, K = 4, seed = 42)

  res_lasso <- ensembleWeights(cv, Y = cv$.y, solver = "glmnet", alpha = 1)
  res_ridge <- ensembleWeights(cv, Y = cv$.y, solver = "glmnet", alpha = 0)

  # Both should be valid
  expect_true(all(res_lasso$methodCoef >= 0))
  expect_equal(sum(res_lasso$methodCoef), 1, tolerance = 1e-6)
  expect_true(all(res_ridge$methodCoef >= 0))
  expect_equal(sum(res_ridge$methodCoef), 1, tolerance = 1e-6)

  # Lasso should be at least as sparse as ridge (fewer or equal non-zero coefs)
  n_nonzero_lasso <- sum(res_lasso$methodCoef > 1e-8)
  n_nonzero_ridge <- sum(res_ridge$methodCoef > 1e-8)
  expect_true(n_nonzero_lasso <= n_nonzero_ridge)
})


context("multivariate_pipeline")

# multivariateAnalysisPipeline and the legacy regionDataToMvsusieRssInput /
# regionDataToSusieRssInput entry points have been removed in the S4
# refactor. Joint multi-trait/multi-context analyses now run through
# fineMappingPipeline(qtlDataset, methods = "mvsusie") (individual-level)
# or fineMappingPipeline(..., methods = "mvsusieRSS") (RSS-based).
# Coverage for those pipelines lives alongside their implementations.

test_that("multivariateAnalysisPipeline is a deprecated no-op", {
  expect_warning(
    res <- multivariateAnalysisPipeline(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})


context("twasWeights internal helpers")

# ===========================================================================
# .twasNormalizeMethods
# ===========================================================================

test_that(".twasNormalizeMethods: NULL falls through to the 'default' preset", {
  res <- pecotmr:::.twasNormalizeMethods(NULL)
  default_names <- names(pecotmr:::.twasMethodLookup("default"))
  expect_equal(sort(names(res$methodList)), sort(default_names))
})

test_that(".twasNormalizeMethods: character preset string forwards to .twasMethodLookup", {
  res <- pecotmr:::.twasNormalizeMethods("fast_default")
  fast_names <- names(pecotmr:::.twasMethodLookup("fast_default"))
  expect_equal(sort(names(res$methodList)), sort(fast_names))
})

test_that(".twasNormalizeMethods: character vector of short names forwards to .twasMethodLookup", {
  res <- pecotmr:::.twasNormalizeMethods(c("lasso", "enet"))
  expect_equal(sort(names(res$methodList)),
               sort(c("lasso_weights", "enet_weights")))
})

test_that(".twasNormalizeMethods: named list passes through unchanged", {
  ml <- list(lassoWeights = list(), enetWeights = list(alpha = 0.5))
  res <- pecotmr:::.twasNormalizeMethods(ml)
  expect_identical(res$methodList, ml)
})

test_that(".twasNormalizeMethods: tokens strip both _weights and Weights suffixes", {
  res_snake <- pecotmr:::.twasNormalizeMethods(list(lasso_weights = list(),
                                                    enet_weights  = list()))
  res_camel <- pecotmr:::.twasNormalizeMethods(list(lassoWeights = list(),
                                                    enetWeights  = list()))
  expect_equal(res_snake$tokens, c("lasso", "enet"))
  expect_equal(res_camel$tokens, c("lasso", "enet"))
})

test_that(".twasNormalizeMethods: unrecognised input type errors", {
  expect_error(
    pecotmr:::.twasNormalizeMethods(42L),
    "must be a character vector, preset string, or named list"
  )
})

# ===========================================================================
# .twasCheckMethodCapabilities
# ===========================================================================

test_that(".twasCheckMethodCapabilities: empty token list is a no-op", {
  expect_silent(pecotmr:::.twasCheckMethodCapabilities(character(0), "QtlDataset"))
})

test_that(".twasCheckMethodCapabilities: individual-only token + QtlDataset is fine", {
  expect_silent(pecotmr:::.twasCheckMethodCapabilities(c("lasso", "enet"),
                                                       "QtlDataset"))
  expect_silent(pecotmr:::.twasCheckMethodCapabilities(c("lasso", "enet"),
                                                       "MultiStudyQtlDataset"))
})

test_that(".twasCheckMethodCapabilities: sumstat-only token + QtlDataset errors", {
  # prsCs has sumstatImpl only -- not legal for QtlDataset.
  expect_error(
    pecotmr:::.twasCheckMethodCapabilities("prsCs", "QtlDataset"),
    "is sumstat-only"
  )
})

test_that(".twasCheckMethodCapabilities: individual-only token + QtlSumStats errors", {
  # enet has individualImpl only -- not legal for QtlSumStats.
  expect_error(
    pecotmr:::.twasCheckMethodCapabilities("enet", "QtlSumStats"),
    "is individual-only"
  )
})

test_that(".twasCheckMethodCapabilities: unknown token errors with the full menu", {
  expect_error(
    pecotmr:::.twasCheckMethodCapabilities("bogus", "QtlDataset"),
    "unknown method token"
  )
})

# ===========================================================================
# .twasCheckMultivariateY
# ===========================================================================

test_that(".twasCheckMultivariateY: no multivariate tokens is a no-op", {
  expect_silent(pecotmr:::.twasCheckMultivariateY(c("lasso", "enet"),
                                                   nTraits = 1L, nContexts = 1L))
})

test_that(".twasCheckMultivariateY: mvsusie passes when contexts >= 2", {
  expect_silent(pecotmr:::.twasCheckMultivariateY("mvsusie",
                                                   nTraits = 1L, nContexts = 2L))
})

test_that(".twasCheckMultivariateY: mvsusie passes when traits >= 2", {
  expect_silent(pecotmr:::.twasCheckMultivariateY("mvsusie",
                                                   nTraits = 2L, nContexts = 1L))
})

test_that(".twasCheckMultivariateY: mvsusie/mrmash error with single trait, single context", {
  expect_error(
    pecotmr:::.twasCheckMultivariateY(c("mvsusie", "mrmash"),
                                       nTraits = 1L, nContexts = 1L),
    "require multi-trait or multi-context input"
  )
})

# ===========================================================================
# .twasAssertQcd
# ===========================================================================

.tw_makeSumStatsBare <- function() {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 200L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = c("rs1", "rs2"), A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1, 2), N = c(100L, 100L))
  gh <- new("GenotypeHandle",
            path = "/tmp/x.gds", format = "gds",
            snpInfo = data.frame(), nSamples = 0L,
            sampleIds = character(), pgenPtr = NULL)
  GwasSumStats(study = "g1", entry = list(gr),
                genome = "hg19", ldSketch = gh)
}

test_that(".twasAssertQcd: errors when qcInfo is empty", {
  ss <- .tw_makeSumStatsBare()
  expect_error(
    pecotmr:::.twasAssertQcd(ss),
    "has no QC record"
  )
})

test_that(".twasAssertQcd: passes when qcInfo is populated", {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 200L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = c("rs1", "rs2"), A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1, 2), N = c(100L, 100L))
  gh <- new("GenotypeHandle",
            path = "/tmp/x.gds", format = "gds",
            snpInfo = data.frame(), nSamples = 0L,
            sampleIds = character(), pgenPtr = NULL)
  ss <- GwasSumStats(study = "g1", entry = list(gr),
                      genome = "hg19", ldSketch = gh,
                      qcInfo = list(step1 = "ok"))
  expect_silent(pecotmr:::.twasAssertQcd(ss))
})

# ===========================================================================
# getSumstatDf (public method on GwasSumStats / QtlSumStats; replaces
# the now-deleted `.twasSumstatsEntryToDf` shim)
# ===========================================================================

.gsd_makeHandle <- function(snp_n = 2L) {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(
      SNP = paste0("rs", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n), A2 = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = 10L, sampleIds = paste0("s", seq_len(10L)),
    pgenPtr = NULL)
}

.gsd_makeGwasSumStats <- function(mcolsDf, snp_n = nrow(mcolsDf)) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", snp_n),
    ranges = IRanges::IRanges(
      start = seq(100L, by = 100L, length.out = snp_n), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcolsDf)
  GwasSumStats(study = "G1", entry = list(gr), genome = "hg19",
               ldSketch = .gsd_makeHandle(snp_n),
               qcInfo = list(prebuilt = "synthetic"))
}

test_that("getSumstatDf: returns the canonical column layout", {
  ss <- .gsd_makeGwasSumStats(data.frame(
    SNP = c("rs1", "rs2"), A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1.0, -2.5), N = c(1000L, 1500L), MAF = c(0.1, 0.3)))
  df <- getSumstatDf(ss)
  expect_s3_class(df, "data.frame")
  expect_equal(df$variant_id, c("rs1", "rs2"))
  expect_equal(df$chrom, c("chr1", "chr1"))
  expect_equal(df$pos, c(100L, 200L))
  expect_equal(df$z, c(1.0, -2.5))
  expect_equal(df$N, c(1000, 1500))
  expect_equal(df$maf, c(0.1, 0.3))
})

test_that("getSumstatDf: derives z from beta/se when derive='zFromBetaSe'", {
  ss <- .gsd_makeGwasSumStats(data.frame(
    SNP = "rs1", A1 = "A", A2 = "G",
    BETA = 0.5, SE = 0.1, N = 1000L), snp_n = 1L)
  df <- getSumstatDf(ss, derive = "zFromBetaSe")
  expect_equal(df$beta, 0.5)
  expect_equal(df$se, 0.1)
  expect_equal(df$z, 0.5 / 0.1)
})

test_that("getSumstatDf: omits optional columns when absent", {
  ss <- .gsd_makeGwasSumStats(data.frame(
    SNP = "rs1", A1 = "A", A2 = "G"), snp_n = 1L)
  df <- getSumstatDf(ss)
  expect_false("z"    %in% names(df))
  expect_false("beta" %in% names(df))
  expect_false("se"   %in% names(df))
  expect_false("N"    %in% names(df))
  expect_false("maf"  %in% names(df))
})

test_that("getSumstatDf: require=c('Z','N') errors when columns missing", {
  ss <- .gsd_makeGwasSumStats(data.frame(
    SNP = "rs1", A1 = "A", A2 = "G", Z = 1.5), snp_n = 1L)
  expect_error(getSumstatDf(ss, require = c("Z", "N")),
                "no N mcol")
})

# ===========================================================================
# estimateSparsity
# ===========================================================================

test_that("estimateSparsity: legacy list input reads attr(.,'fit')$pi", {
  # Build a weight result that mimics learnTwasWeights(retainFits = TRUE):
  # element name carries the `_weights` suffix; attr 'fit' carries an mr.ash
  # object whose pi[1] is the spike weight.
  fake_w <- structure(c(0.1, 0, 0.3),
                      fit = list(pi = c(0.6, 0.2, 0.2)))
  weightResults <- list(mrash_weights = fake_w)
  expect_equal(estimateSparsity(weightResults), 1 - 0.6,
               tolerance = 1e-12)
})

test_that("estimateSparsity: TwasWeights collection input reads from the mrash entry", {
  entry <- TwasWeightsEntry(
    variantIds    = c("v1", "v2", "v3"),
    weights       = c(0.1, 0, 0.3),
    fits          = list(pi = c(0.7, 0.2, 0.1)),
    standardized  = FALSE)
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "mrash",
    entry = list(entry))
  expect_equal(estimateSparsity(tw), 1 - 0.7, tolerance = 1e-12)
})

test_that("estimateSparsity: TwasWeights without mrash entry errors", {
  entry <- TwasWeightsEntry(variantIds = c("v1"), weights = c(0.1))
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "lasso",
    entry = list(entry))
  expect_error(
    estimateSparsity(tw),
    "mr.ash entry not found in TwasWeights"
  )
})

test_that("estimateSparsity: TwasWeights with mrash entry but no fit$pi errors", {
  entry <- TwasWeightsEntry(variantIds = c("v1"), weights = c(0.1),
                            fits = list(other = 1))
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "mrash",
    entry = list(entry))
  expect_error(
    estimateSparsity(tw),
    "mr.ash fit object not found"
  )
})

test_that("estimateSparsity: legacy list input without mrash_weights errors", {
  expect_error(
    estimateSparsity(list(lasso_weights = c(0.1, 0.2))),
    "'mrash_weights'.*not found"
  )
})

test_that("estimateSparsity: legacy list input without fit attr errors", {
  expect_error(
    estimateSparsity(list(mrash_weights = c(0.1, 0.2))),
    "mr.ash fit object not found"
  )
})


context("twasWeights internal helpers (extra)")

# ===========================================================================
# .twasFineMappingFits
# ===========================================================================

.tw_makeFmEntry <- function(method_tag = "susie", n = 3) {
  FineMappingEntry(
    variantIds = paste0("v", seq_len(n)),
    susieFit = list(payload = method_tag),
    topLoci    = data.frame(variant_id = paste0("v", seq_len(n)),
                            pip = seq(0.9, by = -0.1, length.out = n),
                            stringsAsFactors = FALSE))
}

test_that(".twasFineMappingFits: NULL input returns an empty list", {
  out <- pecotmr:::.twasFineMappingFits(NULL,
                                        study = "s1", context = "c1",
                                        trait = "t1")
  expect_equal(out, list())
})

test_that(".twasFineMappingFits: non-FineMappingResult input errors", {
  expect_error(
    pecotmr:::.twasFineMappingFits("not_a_result",
                                    study = "s1", context = "c1", trait = "t1"),
    "must be a FineMappingResult or NULL"
  )
})

test_that(".twasFineMappingFits: pulls matching (study,context,trait) susie fit", {
  e_susie    <- .tw_makeFmEntry("susie")
  e_susieInf <- .tw_makeFmEntry("susieInf")
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susieInf"),
    entry   = list(e_susie, e_susieInf))
  out <- pecotmr:::.twasFineMappingFits(res,
                                         study = "s1", context = "c1",
                                         trait = "t1")
  expect_true("susie" %in% names(out))
  expect_true("susieInf" %in% names(out))
  expect_equal(out$susie$payload, "susie")
  expect_equal(out$susieInf$payload, "susieInf")
})

test_that(".twasFineMappingFits: returns empty list when no tuple matches", {
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(.tw_makeFmEntry("susie")))
  out <- pecotmr:::.twasFineMappingFits(res,
                                         study = "ghost", context = "c1",
                                         trait = "t1")
  expect_equal(out, list())
})

test_that(".twasFineMappingFits: ignores non-susie methods (e.g. lasso)", {
  e <- .tw_makeFmEntry("susie")
  e_lasso <- .tw_makeFmEntry("lasso")
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("susie", "lasso"),
    entry   = list(e, e_lasso))
  out <- pecotmr:::.twasFineMappingFits(res,
                                         study = "s1", context = "c1",
                                         trait = "t1")
  # Only susie should be picked up; lasso is not in the susie/susieInf/susieAsh
  # canonical menu.
  expect_equal(names(out), "susie")
})

# ===========================================================================
# .twasLdFromSketch (mocked extractBlockGenotypes)
# ===========================================================================

.tw_makeSketchHandle <- function(snp_n = 6L, n_samples = 30L) {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds",
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

.tw_mockExtractor <- function(seed = 7, n_samples = 30L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds,
                                    handle@snpInfo$SNP))
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

test_that(".twasLdFromSketch: rejects non-GenotypeHandle ldSketch", {
  expect_error(
    pecotmr:::.twasLdFromSketch("not_a_handle", c("v1", "v2")),
    "ldSketch must be a GenotypeHandle"
  )
})

test_that(".twasLdFromSketch: variants not in the panel error", {
  h <- .tw_makeSketchHandle()
  expect_error(
    pecotmr:::.twasLdFromSketch(h, c("v1", "ghost")),
    "variant id.*not present in the LD sketch"
  )
})

test_that(".twasLdFromSketch: returns a square LD matrix named by variantIds", {
  h <- .tw_makeSketchHandle()
  local_mocked_bindings(extractBlockGenotypes = .tw_mockExtractor(),
                        .package = "pecotmr")
  ids <- c("v2", "v4", "v5")
  R <- pecotmr:::.twasLdFromSketch(h, ids)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(rownames(R), ids)
  expect_equal(colnames(R), ids)
  # Symmetric, diagonal == 1 (sample correlation).
  expect_equal(unname(diag(R)), c(1, 1, 1), tolerance = 1e-12)
  expect_equal(R, t(R), tolerance = 1e-12)
})

# ===========================================================================
# .twasWeightsPipelineMatrix: susieFit pre-fit pass-through
# ===========================================================================

test_that(".twasWeightsPipelineMatrix: susieFit pre-fit is recorded in res", {
  set.seed(0)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:p)))
  y <- as.numeric(X %*% c(1.0, -0.5, 0, 0, 0) + rnorm(n, sd = 0.2))

  # Build a stub susie fit shape; the pipeline records its intermediate.
  fake_susie <- list(
    alpha = matrix(1/p, nrow = 2, ncol = p),
    mu    = matrix(0,    nrow = 2, ncol = p),
    X_column_scale_factors = rep(1, p),
    pip   = rep(0.1, p))

  # The intermediate-recording branch keys on snake_case `susie_weights`.
  res <- suppressMessages(
    pecotmr:::.twasWeightsPipelineMatrix(
      X = X, y = y,
      susieFit = fake_susie,
      cvFolds = 0,
      weightMethods = list(susie_weights = list()),
      estimatePi = FALSE,
      verbose = 0))
  expect_true("susieWeightsIntermediate" %in% names(res))
  expect_true("twasWeights" %in% names(res))
})

# ===========================================================================
# .twasWeightsPipelineMatrix: empirical pi path via mr.ash (mocked)
# ===========================================================================

test_that(".twasWeightsPipelineMatrix: empirical pi from mr.ash gets propagated", {
  set.seed(1)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:p)))
  y <- as.numeric(X %*% c(0.5, 0, 0, 0, 0) + rnorm(n, sd = 0.2))

  # Mock mrashWeights to return a fake matrix carrying a fit$pi attribute.
  local_mocked_bindings(
    mrashWeights = function(X, y, ...) {
      out <- matrix(rep(0.05, ncol(X)), ncol = 1)
      attr(out, "fit") <- list(pi = c(0.8, 0.1, 0.1))
      rownames(out) <- colnames(X)
      out
    },
    bayesCWeights = function(X, y, pi, ...) {
      # Capture the pi the pipeline injected.
      out <- matrix(pi, nrow = ncol(X), ncol = 1)
      rownames(out) <- colnames(X)
      out
    },
    .package = "pecotmr"
  )

  res <- suppressMessages(
    pecotmr:::.twasWeightsPipelineMatrix(
      X = X, y = y,
      cvFolds = 0,
      weightMethods = list(mrash_weights = list(),
                           bayes_c_weights = list()),
      estimatePi = TRUE,
      verbose = 0))
  expect_true("empiricalPi" %in% names(res))
  expect_equal(as.numeric(res$empiricalPi), 1 - 0.8, tolerance = 1e-12)
})
