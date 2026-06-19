context("fineMappingPipeline")

# ===========================================================================
# Strategy
# ----------------------------------------------------------------------------
# Mock the SuSiE fitters (.fmFitSusieIndiv / .fmFitSusieRss) and the post-
# processor (.fmPostprocessOne) so the pipeline orchestration runs end-to-
# end without firing real susieR / susie_rss / postprocess_finemapping_fits
# calls. The fixture uses a small in-memory QtlDataset (mocked
# extractBlockGenotypes) and small QtlSumStats / GwasSumStats collections
# (with mocked extractBlockGenotypes for the LD sketch).
# ===========================================================================

.fmp_makeHandle <- function(snp_n = 6L, n_samples = 40L) {
  new("GenotypeHandle",
    path = "/tmp/fmsketch.gds",
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

.fmp_mockExtractor <- function(seed = 3, n_samples = 40L) {
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

.fmp_makeSe <- function(traits = c("ENSG_A", "ENSG_B"), n_samples = 40L,
                        starts = NULL) {
  if (is.null(starts)) starts <- seq(1000L, by = 1000L, length.out = length(traits))
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

.fmp_makeQtlDataset <- function(contexts = "brain",
                                traits = c("ENSG_A", "ENSG_B")) {
  gh <- .fmp_makeHandle()
  phen <- setNames(lapply(contexts, function(.) .fmp_makeSe(traits = traits)),
                   contexts)
  QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = phen,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
}

.fmp_makeSumstatsGr <- function(snp_ids = paste0("v", 1:5)) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L,
                                          length.out = length(snp_ids)),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = snp_ids, A1 = rep("A", length(snp_ids)),
    A2 = rep("G", length(snp_ids)),
    Z = rnorm(length(snp_ids)), N = rep(1000L, length(snp_ids)))
  gr
}

.fmp_makeQtlSumStats <- function(qc = TRUE) {
  QtlSumStats(
    study    = "Q1", context = "c1", trait = "t1",
    entry    = list(.fmp_makeSumstatsGr()),
    genome   = "hg19",
    ldSketch = .fmp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.fmp_makeGwasSumStats <- function(qc = TRUE, study = "G1") {
  GwasSumStats(
    study    = study,
    entry    = list(.fmp_makeSumstatsGr()),
    genome   = "hg19",
    ldSketch = .fmp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

# Mocks for the SuSiE fitters + post-processor. Return tiny payloads keyed
# only by the token so post-process knows what to wrap.
.fmp_mockFitIndiv <- function() {
  function(X, y, token, chainFromInf = NULL, coverage = 0.95) {
    list(token = token, X_cols = ncol(X))
  }
}

.fmp_mockFitRss <- function() {
  function(z, R, n, token, chainFromInf = NULL, coverage = 0.95) {
    list(token = token, n_variants = length(z))
  }
}

.fmp_mockPostprocess <- function() {
  function(fit, method, dataX, dataY, coverage, secondaryCoverage,
           signalCutoff, minAbsCorr, csInput = NULL, af = NULL,
           region = NULL) {
    # Capture the requesting method on the FineMappingEntry so the test can
    # verify the right dispatch happened.
    if (is.matrix(dataX)) {
      vids <- colnames(dataX)
    } else {
      vids <- names(dataY)
      if (is.null(vids) && is.list(dataY) && !is.null(dataY$z))
        vids <- names(dataY$z)
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

# ===========================================================================
# .fmNormalizeMethods
# ===========================================================================

test_that(".fmNormalizeMethods: rejects NULL / empty / non-character", {
  expect_error(pecotmr:::.fmNormalizeMethods(NULL),
               "non-empty character vector")
  expect_error(pecotmr:::.fmNormalizeMethods(character(0)),
               "non-empty character vector")
  expect_error(pecotmr:::.fmNormalizeMethods(42L),
               "must be a character vector")
})

test_that(".fmNormalizeMethods: deduplicates", {
  expect_equal(pecotmr:::.fmNormalizeMethods(c("susie", "susie", "susieInf")),
               c("susie", "susieInf"))
})

# ===========================================================================
# .fmCheckMethodCapabilities
# ===========================================================================

test_that(".fmCheckMethodCapabilities: unknown token errors with full menu", {
  expect_error(
    pecotmr:::.fmCheckMethodCapabilities("bogus", "QtlDataset"),
    "unknown method token"
  )
})

test_that(".fmCheckMethodCapabilities: mrmash always rejected", {
  expect_error(
    pecotmr:::.fmCheckMethodCapabilities("mrmash", "QtlDataset"),
    "TWAS-weight-oriented"
  )
})

test_that(".fmCheckMethodCapabilities: fsusie on QtlSumStats rejected (no sumstatImpl)", {
  expect_error(
    pecotmr:::.fmCheckMethodCapabilities("fsusie", "QtlSumStats"),
    "individual-only"
  )
})

test_that(".fmCheckMethodCapabilities: mvsusie on GwasSumStats rejected", {
  expect_error(
    pecotmr:::.fmCheckMethodCapabilities("mvsusie", "GwasSumStats"),
    "not supported on GwasSumStats"
  )
})

# ===========================================================================
# .fmResolveSusieChain
# ===========================================================================

test_that(".fmResolveSusieChain: chains susie from susieInf when both are requested", {
  res <- pecotmr:::.fmResolveSusieChain(c("susieInf", "susie"),
                                         addSusieInf = TRUE)
  expect_true(res$chainSusie)
  expect_true(res$runInf)
  expect_true(res$keepInf)
})

test_that(".fmResolveSusieChain: keeps susieInf when explicitly requested", {
  res <- pecotmr:::.fmResolveSusieChain(c("susieInf", "susie"), addSusieInf = FALSE)
  expect_true(res$runInf)
  expect_true(res$keepInf)
})

test_that(".fmResolveSusieChain: no chain when addSusieInf=FALSE", {
  res <- pecotmr:::.fmResolveSusieChain(c("susie"), addSusieInf = FALSE)
  expect_false(res$chainSusie)
  expect_false(res$runInf)
})

# ===========================================================================
# .fmCacheLookup / .fmCacheLookupGwas
# ===========================================================================

test_that(".fmCacheLookup: NULL fineMappingResult returns NULL", {
  expect_null(pecotmr:::.fmCacheLookup(NULL, "s1", "c1", "t1", "susie"))
})

test_that(".fmCacheLookup: returns matching entry by 4-tuple", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  fmr <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  hit <- pecotmr:::.fmCacheLookup(fmr, "s1", "c1", "t1", "susie")
  expect_identical(hit, e)
  expect_null(pecotmr:::.fmCacheLookup(fmr, "ghost", "c1", "t1", "susie"))
})

test_that(".fmCacheLookupGwas: returns matching entry by (study, method)", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  fmr <- GwasFineMappingResult(study = "g1", method = "susie",
                                entry = list(e))
  expect_identical(
    pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie"),
    e)
  expect_null(pecotmr:::.fmCacheLookupGwas(fmr, "ghost", "susie"))
})

test_that(".fmCacheLookup: non-QtlFineMappingResult input returns NULL", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  gwasFmr <- GwasFineMappingResult(study = "g1", method = "susie",
                                    entry = list(e))
  expect_null(pecotmr:::.fmCacheLookup(gwasFmr, "g1", "c1", "t1", "susie"))
})

test_that(".fmCacheLookupGwas: non-GwasFineMappingResult input returns NULL", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  qtlFmr <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_null(pecotmr:::.fmCacheLookupGwas(qtlFmr, "s1", "susie"))
})

# ===========================================================================
# .fmBuildQtlResult / .fmBuildGwasResult — empty-entries errors
# ===========================================================================

test_that(".fmBuildQtlResult: empty entries errors", {
  expect_error(
    pecotmr:::.fmBuildQtlResult(character(0), character(0), character(0),
                                 character(0), list()),
    "no \\(study, context, trait, method\\) tuples"
  )
})

test_that(".fmBuildGwasResult: empty entries errors", {
  expect_error(
    pecotmr:::.fmBuildGwasResult(character(0), character(0), list()),
    "no \\(study, method\\) tuples"
  )
})

# ===========================================================================
# .rbindFineMappingResult — class-check branches
# ===========================================================================

test_that(".rbindFineMappingResult: rejects non-FineMappingResultBase input", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  fmr <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_error(
    pecotmr:::.rbindFineMappingResult(fmr, "not_an_fmr"),
    "expects two FineMappingResultBase inputs"
  )
  expect_error(
    pecotmr:::.rbindFineMappingResult("not_an_fmr", fmr),
    "expects two FineMappingResultBase inputs"
  )
})

test_that(".rbindFineMappingResult: rejects mixed Qtl/Gwas inputs", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  qtlFmr <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  gwasFmr <- GwasFineMappingResult(
    study = "g1", method = "susie", entry = list(e))
  expect_error(
    pecotmr:::.rbindFineMappingResult(qtlFmr, gwasFmr),
    "inputs must be the same concrete class"
  )
})

test_that(".rbindFineMappingResult: concatenates two GwasFineMappingResult collections", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  a <- GwasFineMappingResult(study = "g1", method = "susie", entry = list(e))
  b <- GwasFineMappingResult(study = "g2", method = "susie", entry = list(e))
  out <- pecotmr:::.rbindFineMappingResult(a, b)
  expect_s4_class(out, "GwasFineMappingResult")
  expect_equal(nrow(out), 2L)
})

# ===========================================================================
# .fmExtractZN
# ===========================================================================

test_that(".fmExtractZN: errors on missing SNP / Z / N columns", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100))
  expect_error(pecotmr:::.fmExtractZN(gr, "x"), "no SNP mcol")
  S4Vectors::mcols(gr)$SNP <- "v1"
  expect_error(pecotmr:::.fmExtractZN(gr, "x"), "no Z mcol")
  S4Vectors::mcols(gr)$Z <- 1.0
  expect_error(pecotmr:::.fmExtractZN(gr, "x"), "no N mcol")
})

# ===========================================================================
# .fmLdFromSketch
# ===========================================================================

test_that(".fmLdFromSketch: returns named LD matrix; missing variants error", {
  h <- .fmp_makeHandle()
  local_mocked_bindings(extractBlockGenotypes = .fmp_mockExtractor(),
                        .package = "pecotmr")
  R <- pecotmr:::.fmLdFromSketch(h, c("v1", "v3"))
  expect_equal(dim(R), c(2L, 2L))
  expect_equal(rownames(R), c("v1", "v3"))
  expect_error(pecotmr:::.fmLdFromSketch(h, c("v1", "ghost")),
               "not present in the LD sketch")
})

# ===========================================================================
# fineMappingPipeline(QtlDataset)
# ===========================================================================

test_that("fineMappingPipeline(QtlDataset): runs univariate dispatch with mocked fitters", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie",
                        cisWindow = 1000L,
                        addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  # 1 context x 2 traits x 1 method = 2 rows.
  expect_equal(nrow(res), 2L)
  expect_setequal(getMethodNames(res), "susie")
})

test_that("fineMappingPipeline(QtlDataset): RSS-only method rejected by capability check", {
  qd <- .fmp_makeQtlDataset()
  expect_error(
    fineMappingPipeline(qd, methods = "mrmash"),
    "TWAS-weight-oriented"
  )
})

test_that("fineMappingPipeline(QtlDataset): unknown context errors", {
  qd <- .fmp_makeQtlDataset()
  expect_error(
    fineMappingPipeline(qd, methods = "susie", contexts = "ghost"),
    "unknown context"
  )
})

test_that("fineMappingPipeline(QtlDataset): empty traitId filter errors", {
  qd <- .fmp_makeQtlDataset(traits = "ENSG_A")
  expect_error(
    fineMappingPipeline(qd, methods = "susie", traitId = "ENSG_Z"),
    "no traits selected"
  )
})

test_that("fineMappingPipeline(QtlDataset): mvsusie with single trait/context rejected", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    fineMappingPipeline(qd, methods = "mvsusie"),
    "mvsusie requires multi-trait or multi-context"
  )
})

test_that("fineMappingPipeline(QtlDataset): fsusie with single trait rejected", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  expect_error(
    fineMappingPipeline(qd, methods = "fsusie"),
    "fsusie requires multi-trait"
  )
})

# ===========================================================================
# fineMappingPipeline(QtlDataset): mvsusie + fsusie dispatch (mocked)
# ===========================================================================

# Mock mvsusieR::mvsusie / create_mixture_prior so the joint-fit branches run
# without actually fitting. Returns a stub fit object tagged with the input
# shape so the test can assert it was constructed as expected.
.fmp_mockMvsusie <- function() {
  function(X, Y, prior_variance, coverage) {
    list(token = "mvsusie",
         n_X_cols = ncol(X),
         n_Y_cols = ncol(Y))
  }
}
.fmp_mockMixturePrior <- function() {
  function(R, ...) list(R = R)
}
.fmp_mockSusiF <- function() {
  function(X, Y, pos) {
    list(token = "fsusie",
         n_X_cols = ncol(X),
         n_Y_cols = ncol(Y),
         pos = pos)
  }
}

test_that("fineMappingPipeline(QtlDataset): mvsusie multi-trait single-context dispatch", {
  qd <- .fmp_makeQtlDataset(contexts = "brain",
                            traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie               = .fmp_mockMvsusie(),
    create_mixture_prior  = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L))
  expect_s4_class(res, "QtlFineMappingResult")
  # mvsusie multi-trait fans the joint fit out across both traits.
  expect_equal(nrow(res), 2L)
  expect_setequal(getTraits(res), c("ENSG_A", "ENSG_B"))
  expect_setequal(getMethodNames(res), "mvsusie")
})

test_that("fineMappingPipeline(QtlDataset): mvsusie multi-context single-trait dispatch", {
  qd <- .fmp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie               = .fmp_mockMvsusie(),
    create_mixture_prior  = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L))
  # Multi-context fan-out: one row per context for the shared trait.
  expect_equal(nrow(res), 2L)
  expect_setequal(getContexts(res), c("brain", "liver"))
  expect_setequal(getTraits(res), "ENSG_A")
})

test_that("fineMappingPipeline(QtlDataset): mvsusie both multi falls back to per-context multi-trait", {
  qd <- .fmp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie               = .fmp_mockMvsusie(),
    create_mixture_prior  = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L))
  # 2 contexts * 2 traits = 4 rows (joint fit reused per context).
  expect_equal(nrow(res), 4L)
})

test_that("fineMappingPipeline(QtlDataset): fsusie multi-trait per context dispatch", {
  qd <- .fmp_makeQtlDataset(contexts = "brain",
                            traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    susiF                 = .fmp_mockSusiF(),
    .package = "fsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "fsusie", cisWindow = 1000L))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)
  expect_setequal(getTraits(res), c("ENSG_A", "ENSG_B"))
  expect_setequal(getMethodNames(res), "fsusie")
})

# ===========================================================================
# fineMappingPipeline(MultiTaskQtlDataset)
# ===========================================================================

test_that("fineMappingPipeline(MultiTaskQtlDataset): aggregates results across constituent QtlDatasets", {
  qd1 <- QtlDataset(
    study              = "s1",
    genotypes          = .fmp_makeHandle(),
    phenotypes         = list(brain = .fmp_makeSe(traits = "ENSG_A")),
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
  qd2 <- QtlDataset(
    study              = "s2",
    genotypes          = .fmp_makeHandle(),
    phenotypes         = list(brain = .fmp_makeSe(traits = "ENSG_A")),
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
  mt <- MultiTaskQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(mt, methods = "susie",
                        cisWindow = 1000L, addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  # One row per (study, context, trait, method) tuple.
  expect_equal(nrow(res), 2L)
  expect_setequal(getStudy(res), c("s1", "s2"))
  # Pure individual-level -> ldSketch should be NULL.
  expect_null(getLdSketch(res))
})

test_that("fineMappingPipeline(MultiTaskQtlDataset): with embedded QtlSumStats stamps the ldSketch", {
  qd <- QtlDataset(
    study              = "s1",
    genotypes          = .fmp_makeHandle(),
    phenotypes         = list(brain = .fmp_makeSe(traits = "ENSG_A")),
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
  ss <- .fmp_makeQtlSumStats()
  mt <- MultiTaskQtlDataset(qtlDatasets = list(s1 = qd), sumStats = ss)
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmFitSusieRss        = .fmp_mockFitRss(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(mt, methods = "susie",
                        cisWindow = 1000L, addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_true(nrow(res) >= 2L)
  # Embedded sumstats has an LD sketch -> the merged result carries it.
  expect_s4_class(getLdSketch(res), "GenotypeHandle")
})

# ===========================================================================
# fineMappingPipeline(QtlSumStats)
# ===========================================================================

test_that("fineMappingPipeline(QtlSumStats): runs end-to-end with mocked RSS fitters", {
  ss <- .fmp_makeQtlSumStats()
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss        = .fmp_mockFitRss(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "susie", addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 1L)
})

test_that("fineMappingPipeline(QtlSumStats): un-QCd input rejected", {
  ss <- .fmp_makeQtlSumStats(qc = FALSE)
  expect_error(
    fineMappingPipeline(ss, methods = "susie"),
    "has no QC record"
  )
})

test_that("fineMappingPipeline(QtlSumStats): empty selection rejected", {
  ss <- .fmp_makeQtlSumStats()
  expect_error(
    fineMappingPipeline(ss, methods = "susie", contexts = "ghost"),
    "no entries matched"
  )
})

# ===========================================================================
# fineMappingPipeline(GwasSumStats)
# ===========================================================================

test_that("fineMappingPipeline(GwasSumStats): runs end-to-end with mocked RSS fitters", {
  gss <- .fmp_makeGwasSumStats()
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss        = .fmp_mockFitRss(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(gss, methods = "susie", addSusieInf = FALSE))
  expect_s4_class(res, "GwasFineMappingResult")
  expect_equal(nrow(res), 1L)
  expect_setequal(getMethodNames(res), "susie")
})

test_that("fineMappingPipeline(GwasSumStats): un-QCd input rejected", {
  gss <- .fmp_makeGwasSumStats(qc = FALSE)
  expect_error(
    fineMappingPipeline(gss, methods = "susie"),
    "has no QC record"
  )
})

test_that("fineMappingPipeline(GwasSumStats): non-RSS family rejected by capability check", {
  gss <- .fmp_makeGwasSumStats()
  expect_error(
    fineMappingPipeline(gss, methods = "fsusie"),
    "not supported on GwasSumStats"
  )
})

# ===========================================================================
# fineMappingPipeline(ANY)
# ===========================================================================

test_that("fineMappingPipeline(ANY): unsupported input class errors", {
  expect_error(
    fineMappingPipeline(matrix(0, 5, 5), methods = "susie"),
    "does not accept inputs of class 'matrix'"
  )
})

# ===========================================================================
# Cache hit short-circuits the fit
# ===========================================================================

# ===========================================================================
# .fmFitSusieIndiv / .fmFitSusieRss — chained-init and branch coverage
# ===========================================================================

# Capture the args passed to susieR::susie / susie_rss by mocking each to
# stash its first invocation's args into a global. The captured args let
# us assert which code path was taken.
.fmp_capturingSusie <- function(captured) {
  function(X, y, ...) {
    captured$lastArgs <<- list(X = X, y = y, ...)
    # Return a minimal "fit" shape downstream cares about; .setFinemappingFitClass
    # only attaches an S3 class, so any list works.
    list(token = "test", V = 0.1)
  }
}

.fmp_capturingSusieRss <- function(captured) {
  function(z, R, n, ...) {
    captured$lastArgs <<- list(z = z, R = R, n = n, ...)
    list(token = "test_rss", V = 0.1)
  }
}

test_that(".fmFitSusieIndiv: susieInf branch passes convergence_method='pip', refine=FALSE, model_init=NULL", {
  captured <- new.env(parent = emptyenv())
  X <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  local_mocked_bindings(susie = .fmp_capturingSusie(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieIndiv(X, y, "susieInf")
  expect_equal(captured$lastArgs$convergence_method, "pip")
  expect_false(captured$lastArgs$refine)
  expect_null(captured$lastArgs$model_init)
  expect_equal(captured$lastArgs$unmappable_effects, "inf")
})

test_that(".fmFitSusieIndiv: chained branch (chainFromInf) propagates susieInf fit as model_init", {
  captured <- new.env(parent = emptyenv())
  X <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  # Build a stub susieInf fit with a V slot so prepareSusieFromInfArgs can read L.
  infFit <- list(V = c(0.1, 0.2))
  local_mocked_bindings(susie = .fmp_capturingSusie(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieIndiv(X, y, "susie", chainFromInf = infFit)
  # prepareSusieFromInfArgs writes the susieInf fit into model_init and
  # sets unmappable_effects to "none" for the `susie` token.
  expect_identical(captured$lastArgs$model_init, infFit)
  expect_equal(captured$lastArgs$unmappable_effects, "none")
})

test_that(".fmFitSusieIndiv: chained susieAsh branch sets unmappable_effects='ash'", {
  captured <- new.env(parent = emptyenv())
  X <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  infFit <- list(V = c(0.1, 0.2))
  local_mocked_bindings(susie = .fmp_capturingSusie(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieIndiv(X, y, "susieAsh", chainFromInf = infFit)
  expect_equal(captured$lastArgs$unmappable_effects, "ash")
  expect_identical(captured$lastArgs$model_init, infFit)
})

test_that(".fmFitSusieIndiv: unchained susieAsh branch sets convergence_method='pip'", {
  captured <- new.env(parent = emptyenv())
  X <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  local_mocked_bindings(susie = .fmp_capturingSusie(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieIndiv(X, y, "susieAsh")
  expect_equal(captured$lastArgs$convergence_method, "pip")
  expect_equal(captured$lastArgs$unmappable_effects, "ash")
})

test_that(".fmFitSusieIndiv: rejects non-SuSiE-family token", {
  expect_error(
    pecotmr:::.fmFitSusieIndiv(matrix(0, 2, 2), c(0, 0), "mvsusie"),
    "not a SuSiE-family method"
  )
  expect_error(
    pecotmr:::.fmFitSusieIndiv(matrix(0, 2, 2), c(0, 0), "ghost"),
    "not a SuSiE-family method"
  )
})

test_that(".fmFitSusieRss: susieInf branch passes convergence_method='pip', refine=FALSE, model_init=NULL", {
  captured <- new.env(parent = emptyenv())
  z <- rnorm(3); R <- diag(3); n <- 1000
  local_mocked_bindings(susie_rss = .fmp_capturingSusieRss(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieRss(z, R, n, "susieInf")
  expect_equal(captured$lastArgs$convergence_method, "pip")
  expect_false(captured$lastArgs$refine)
  expect_null(captured$lastArgs$model_init)
  expect_equal(captured$lastArgs$unmappable_effects, "inf")
})

test_that(".fmFitSusieRss: chained branch (chainFromInf) propagates susieInf fit as model_init", {
  captured <- new.env(parent = emptyenv())
  z <- rnorm(3); R <- diag(3); n <- 1000
  infFit <- list(V = c(0.1, 0.2))
  local_mocked_bindings(susie_rss = .fmp_capturingSusieRss(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieRss(z, R, n, "susie", chainFromInf = infFit)
  expect_identical(captured$lastArgs$model_init, infFit)
  expect_equal(captured$lastArgs$unmappable_effects, "none")
})

test_that(".fmFitSusieRss: chained susieAsh branch sets unmappable_effects='ash'", {
  captured <- new.env(parent = emptyenv())
  z <- rnorm(3); R <- diag(3); n <- 1000
  infFit <- list(V = c(0.1, 0.2))
  local_mocked_bindings(susie_rss = .fmp_capturingSusieRss(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieRss(z, R, n, "susieAsh", chainFromInf = infFit)
  expect_equal(captured$lastArgs$unmappable_effects, "ash")
  expect_identical(captured$lastArgs$model_init, infFit)
})

test_that(".fmFitSusieRss: unchained susieAsh branch sets convergence_method='pip'", {
  captured <- new.env(parent = emptyenv())
  z <- rnorm(3); R <- diag(3); n <- 1000
  local_mocked_bindings(susie_rss = .fmp_capturingSusieRss(captured),
                        .package = "susieR")
  pecotmr:::.fmFitSusieRss(z, R, n, "susieAsh")
  expect_equal(captured$lastArgs$convergence_method, "pip")
  expect_equal(captured$lastArgs$unmappable_effects, "ash")
})

test_that(".fmFitSusieRss: rejects non-SuSiE-family token", {
  expect_error(
    pecotmr:::.fmFitSusieRss(c(0, 0), diag(2), 1000, "mvsusie"),
    "not a SuSiE-family method"
  )
})

# ===========================================================================
# QtlSumStats: empty selection, mvsusie single-context rejection, cache hits
# ===========================================================================

.fmp_makeMultiCtxQtlSumStats <- function() {
  # Multi-row QtlSumStats: 2 contexts x 1 trait so the test can filter to
  # a single context and exercise both selRows filters.
  e1 <- .fmp_makeSumstatsGr()
  e2 <- .fmp_makeSumstatsGr()
  QtlSumStats(
    study    = c("Q1", "Q1"),
    context  = c("c1", "c2"),
    trait    = c("t1", "t1"),
    entry    = list(e1, e2),
    genome   = "hg19",
    ldSketch = .fmp_makeHandle(),
    qcInfo   = list(step1 = "ok"))
}

test_that("fineMappingPipeline(QtlSumStats): traitId filter that selects no rows errors", {
  ss <- .fmp_makeMultiCtxQtlSumStats()
  expect_error(
    fineMappingPipeline(ss, methods = "susie", traitId = "ghost"),
    "no entries matched"
  )
})

test_that("fineMappingPipeline(QtlSumStats): mvsusie rejected when every (study, trait) has only one context", {
  # Build a single-context-per-(study, trait) collection. Mvsusie requires
  # at least two contexts per (study, trait) group.
  ss <- QtlSumStats(
    study    = c("Q1", "Q1"),
    context  = c("c1", "c1"),
    trait    = c("t1", "t2"),
    entry    = list(.fmp_makeSumstatsGr(), .fmp_makeSumstatsGr()),
    genome   = "hg19",
    ldSketch = .fmp_makeHandle(),
    qcInfo   = list(step1 = "ok"))
  expect_error(
    fineMappingPipeline(ss, methods = "mvsusie"),
    "mvsusie requires at least two"
  )
})

test_that("fineMappingPipeline(QtlSumStats): cache hit short-circuits the RSS fitter", {
  ss <- .fmp_makeQtlSumStats()
  cachedEntry <- FineMappingEntry(
    variantIds = paste0("v", 1:5),
    trimmedFit = list(token = "susie_cached"),
    topLoci    = data.frame(variant_id = paste0("v", 1:5),
                             pip = seq(0.9, 0.1, length.out = 5),
                             stringsAsFactors = FALSE))
  cache <- QtlFineMappingResult(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(cachedEntry),
    ldSketch = .fmp_makeHandle())
  rss_calls <- 0
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss = function(...) {
      rss_calls <<- rss_calls + 1L
      .fmp_mockFitRss()(...)
    },
    .fmPostprocessOne = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "susie",
                        addSusieInf = FALSE,
                        fineMappingResult = cache))
  expect_equal(rss_calls, 0L)
  expect_equal(nrow(res), 1L)
})

# ===========================================================================
# GwasSumStats: cache hit by (study, method)
# ===========================================================================

test_that("fineMappingPipeline(GwasSumStats): cache hit short-circuits the RSS fitter", {
  gss <- .fmp_makeGwasSumStats()
  cachedEntry <- FineMappingEntry(
    variantIds = paste0("v", 1:5),
    trimmedFit = list(token = "susie_cached"),
    topLoci    = data.frame(variant_id = paste0("v", 1:5),
                             pip = seq(0.9, 0.1, length.out = 5),
                             stringsAsFactors = FALSE))
  cache <- GwasFineMappingResult(
    study = "G1", method = "susie",
    entry = list(cachedEntry),
    ldSketch = .fmp_makeHandle())
  rss_calls <- 0
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss = function(...) {
      rss_calls <<- rss_calls + 1L
      .fmp_mockFitRss()(...)
    },
    .fmPostprocessOne = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(gss, methods = "susie",
                        addSusieInf = FALSE,
                        fineMappingResult = cache))
  expect_equal(rss_calls, 0L)
  expect_equal(nrow(res), 1L)
})

test_that("fineMappingPipeline(GwasSumStats): wrong-shape cache (QtlFineMappingResult) is ignored", {
  # When `fineMappingResult` is a QtlFineMappingResult the GwasSumStats
  # method's cache-lookup branch should treat it as a cache miss and
  # still invoke the RSS fitter.
  gss <- .fmp_makeGwasSumStats()
  cachedEntry <- FineMappingEntry(
    variantIds = paste0("v", 1:5),
    trimmedFit = list(token = "susie_cached"),
    topLoci    = data.frame(variant_id = paste0("v", 1:5),
                             pip = rep(0.5, 5),
                             stringsAsFactors = FALSE))
  wrongCache <- QtlFineMappingResult(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(cachedEntry))
  rss_calls <- 0
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss = function(...) {
      rss_calls <<- rss_calls + 1L
      .fmp_mockFitRss()(...)
    },
    .fmPostprocessOne = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(gss, methods = "susie",
                        addSusieInf = FALSE,
                        fineMappingResult = wrongCache))
  expect_equal(rss_calls, 1L)
  expect_equal(nrow(res), 1L)
})

test_that("fineMappingPipeline(QtlDataset): cache hit avoids the fitter", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  # Build a cache that already has the (study1, brain, ENSG_A, susie) row.
  cachedEntry <- FineMappingEntry(
    variantIds = paste0("v", 1:3),
    trimmedFit = list(token = "susie_cached"),
    topLoci    = data.frame(variant_id = paste0("v", 1:3),
                             pip = c(0.9, 0.5, 0.1),
                             stringsAsFactors = FALSE))
  cache <- QtlFineMappingResult(
    study = "study1", context = "brain", trait = "ENSG_A", method = "susie",
    entry = list(cachedEntry))
  fitter_calls <- 0
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = function(...) {
      fitter_calls <<- fitter_calls + 1L
      .fmp_mockFitIndiv()(...)
    },
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        addSusieInf = FALSE,
                        fineMappingResult = cache))
  expect_equal(fitter_calls, 0L)
  expect_equal(nrow(res), 1L)
})
