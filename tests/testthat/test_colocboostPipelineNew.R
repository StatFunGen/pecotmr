context("colocboostPipeline (S4 dispatch)")

# ===========================================================================
# Strategy
# ----------------------------------------------------------------------------
# colocboost::colocboost is the heavy compute; we mock it to return a stub
# result so the pipeline orchestration runs end-to-end. The helpers
# (.cbBuildLdArgs, .cbMergeSumstatBundles, .cbApplyEventFilters,
# .cbRequireMatchingLdSketches, .cbEmptyResult) are exercised directly.
# ===========================================================================

.cbp_makeHandle <- function(snp_n = 6L, n_samples = 30L,
                            sample_prefix = "s") {
  new("GenotypeHandle",
    path = "/tmp/cb.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("v", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0(sample_prefix, seq_len(n_samples)),
    pgenPtr = NULL)
}

.cbp_mockExtractor <- function(seed = 11, n_samples = 30L) {
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

.cbp_makeSe <- function(traits = c("ENSG_A", "ENSG_B"), n_samples = 30L) {
  rng <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(traits)),
    ranges = IRanges::IRanges(start = seq(1000L, by = 1000L,
                                          length.out = length(traits)),
                              width = 500L))
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

.cbp_makeQtlDataset <- function(contexts = "brain",
                                traits = c("ENSG_A", "ENSG_B")) {
  gh <- .cbp_makeHandle()
  phen <- setNames(lapply(contexts, function(.) .cbp_makeSe(traits = traits)),
                   contexts)
  QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = phen,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
}

.cbp_makeQtlSumStats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = rnorm(5), N = rep(1000L, 5))
  QtlSumStats(
    study    = "Q1", context = "c1", trait = "t1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .cbp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.cbp_makeGwasSumStats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = rnorm(5), N = rep(1000L, 5))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .cbp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

# ===========================================================================
# Internal helpers (run-anywhere)
# ===========================================================================

test_that(".cbBuildLdArgs: square matrices route to LD", {
  R1 <- diag(4); R2 <- diag(4)
  res <- pecotmr:::.cbBuildLdArgs(list(R1, R2))
  expect_true("LD" %in% names(res))
  expect_false("X_ref" %in% names(res))
})

test_that(".cbBuildLdArgs: non-square matrices route to X_ref", {
  X1 <- matrix(0, 10, 4)
  res <- pecotmr:::.cbBuildLdArgs(list(X1))
  expect_true("X_ref" %in% names(res))
})

test_that(".cbBuildLdArgs: empty list returns empty list", {
  expect_equal(pecotmr:::.cbBuildLdArgs(list()), list())
})

test_that(".cbRequireSumStatsQc: un-QCd input errors", {
  ss <- .cbp_makeQtlSumStats(qc = FALSE)
  expect_error(
    pecotmr:::.cbRequireSumStatsQc(ss, "qtlData"),
    "summaryStatsQc"
  )
})

test_that(".cbRequireSumStatsQc: NULL input is a no-op", {
  expect_silent(pecotmr:::.cbRequireSumStatsQc(NULL, "x"))
})

test_that(".cbRequireSumStatsQc: QCd input passes", {
  ss <- .cbp_makeQtlSumStats(qc = TRUE)
  expect_silent(pecotmr:::.cbRequireSumStatsQc(ss, "x"))
})

test_that(".cbApplyEventFilters: NULL/empty filters return Y unchanged", {
  Y <- matrix(1, 3, 2, dimnames = list(NULL, c("a", "b")))
  expect_identical(pecotmr:::.cbApplyEventFilters(Y, NULL, "c1"), Y)
  expect_identical(pecotmr:::.cbApplyEventFilters(Y, list(), "c1"), Y)
})

test_that(".cbApplyEventFilters: missing type_pattern errors", {
  Y <- matrix(1, 3, 2, dimnames = list(NULL, c("a", "b")))
  expect_error(
    pecotmr:::.cbApplyEventFilters(
      Y, list(list(valid_pattern = "x")), "c1"),
    "Each event filter must specify type_pattern"
  )
})

test_that(".cbApplyEventFilters: exclude_pattern drops matching events", {
  Y <- matrix(1, 3, 4, dimnames = list(NULL,
              c("eQTL_geneA", "eQTL_geneB", "sQTL_intronA", "sQTL_intronB")))
  filters <- list(list(type_pattern = "^eQTL_",
                       exclude_pattern = "geneB$"))
  out <- suppressMessages(
    pecotmr:::.cbApplyEventFilters(Y, filters, "c1"))
  expect_setequal(colnames(out),
                  c("eQTL_geneA", "sQTL_intronA", "sQTL_intronB"))
})

test_that(".cbApplyEventFilters: filter that drops everything returns NULL", {
  Y <- matrix(1, 3, 2, dimnames = list(NULL, c("eQTL_a", "eQTL_b")))
  filters <- list(list(type_pattern = "^eQTL_",
                       exclude_pattern = ".*"))
  expect_null(suppressMessages(
    pecotmr:::.cbApplyEventFilters(Y, filters, "c1")))
})

test_that(".cbRequireMatchingLdSketches: NULL sides are allowed", {
  expect_silent(pecotmr:::.cbRequireMatchingLdSketches(
    NULL, .cbp_makeHandle()))
  expect_silent(pecotmr:::.cbRequireMatchingLdSketches(
    .cbp_makeHandle(), NULL))
})

test_that(".cbRequireMatchingLdSketches: variant-count mismatch errors", {
  expect_error(
    pecotmr:::.cbRequireMatchingLdSketches(
      .cbp_makeHandle(snp_n = 4L),
      .cbp_makeHandle(snp_n = 5L)),
    "differ in size"
  )
})

test_that(".cbMergeSumstatBundles: empty input gives empty dict", {
  res <- pecotmr:::.cbMergeSumstatBundles(list())
  expect_equal(length(res$sumstat), 0L)
  expect_equal(length(res$LD), 0L)
  expect_equal(nrow(res$dict_sumstatLD), 0L)
})

test_that(".cbMergeSumstatBundles: identical LD matrices are deduplicated", {
  R <- diag(3)
  bundles <- list(
    a = list(sumstat = data.frame(z = 1:3), LD = R),
    b = list(sumstat = data.frame(z = 4:6), LD = R))
  res <- pecotmr:::.cbMergeSumstatBundles(bundles)
  expect_equal(length(res$LD), 1L)
  expect_equal(unique(res$dict_sumstatLD[, 2L]), 1L)
})

test_that(".cbEmptyResult: matches the documented schema", {
  res <- pecotmr:::.cbEmptyResult()
  expect_true(all(c("xqtl_coloc", "joint_gwas", "separate_gwas",
                    "computing_time") %in% names(res)))
})

# ===========================================================================
# colocboostPipeline(QtlDataset)
# ===========================================================================

test_that("colocboostPipeline(QtlDataset): runs xqtl-only ColocBoost with mocked engine", {
  qd <- .cbp_makeQtlDataset()
  capturedArgs <- NULL
  local_mocked_bindings(
    extractBlockGenotypes = .cbp_mockExtractor(),
    .package = "pecotmr")
  local_mocked_bindings(
    colocboost = function(...) {
      capturedArgs <<- list(...)
      list(stub = TRUE)
    },
    .package = "colocboost")
  out <- suppressMessages(
    colocboostPipeline(qd, xqtlColoc = TRUE, jointGwas = FALSE,
                       separateGwas = FALSE))
  expect_true(!is.null(out$xqtl_coloc))
  expect_equal(out$xqtl_coloc$stub, TRUE)
  expect_true("X" %in% names(capturedArgs))
})

test_that("colocboostPipeline(QtlDataset): no QTL bundle and no GWAS returns empty result", {
  qd <- .cbp_makeQtlDataset()
  local_mocked_bindings(extractBlockGenotypes = .cbp_mockExtractor(),
                        .package = "pecotmr")
  # An impossible traitId narrows the phenotype matrix to zero columns;
  # the bundle ends up NULL and the driver short-circuits.
  out <- suppressMessages(suppressWarnings(
    colocboostPipeline(qd, contexts = "brain",
                       traitId = "ENSG_DOES_NOT_EXIST",
                       xqtlColoc = TRUE)))
  expect_null(out$xqtl_coloc)
})

# ===========================================================================
# colocboostPipeline(QtlSumStats)
# ===========================================================================

test_that("colocboostPipeline(QtlSumStats): un-QCd input rejected", {
  ss <- .cbp_makeQtlSumStats(qc = FALSE)
  expect_error(
    colocboostPipeline(ss, xqtlColoc = TRUE),
    "summaryStatsQc"
  )
})

test_that("colocboostPipeline(QtlSumStats): eventFilters trigger a warning (no per-trait Y)", {
  ss <- .cbp_makeQtlSumStats()
  # The warning fires up-front, before .cbQtlSumStatsBundle's extractor
  # call. Mock the extractor anyway so the rest of the pipeline can complete.
  local_mocked_bindings(extractBlockGenotypes = .cbp_mockExtractor(),
                        .package = "pecotmr")
  expect_warning(
    suppressMessages(
      colocboostPipeline(ss, eventFilters = list(list(type_pattern = "x")),
                         xqtlColoc = FALSE)),
    "eventFilters is ignored for QtlSumStats"
  )
})

# ===========================================================================
# colocboostPipeline(ANY)
# ===========================================================================

test_that("colocboostPipeline(ANY): unsupported input class errors", {
  expect_error(
    colocboostPipeline(matrix(0, 3, 3)),
    "does not accept inputs of class"
  )
})

# ===========================================================================
# Driver: jointGwas / separateGwas paths via QtlSumStats + GwasSumStats
# ===========================================================================

test_that("colocboostPipeline: jointGwas merges qtl + gwas sumstats and runs once", {
  ss <- .cbp_makeQtlSumStats()
  gs <- .cbp_makeGwasSumStats()
  capturedArgs <- NULL
  local_mocked_bindings(extractBlockGenotypes = .cbp_mockExtractor(),
                        .package = "pecotmr")
  local_mocked_bindings(
    colocboost = function(...) {
      capturedArgs <<- list(...)
      list(jointly_run = TRUE)
    },
    .package = "colocboost")
  out <- suppressMessages(
    colocboostPipeline(ss, gwasSumStats = gs,
                       xqtlColoc = FALSE, jointGwas = TRUE,
                       separateGwas = FALSE))
  expect_true(!is.null(out$joint_gwas))
})

test_that("colocboostPipeline: separateGwas runs once per merged sumstat study", {
  ss <- .cbp_makeQtlSumStats()
  gs <- .cbp_makeGwasSumStats()
  callCount <- 0
  local_mocked_bindings(extractBlockGenotypes = .cbp_mockExtractor(),
                        .package = "pecotmr")
  local_mocked_bindings(
    colocboost = function(...) {
      callCount <<- callCount + 1L
      list(round = callCount)
    },
    .package = "colocboost")
  out <- suppressMessages(
    colocboostPipeline(ss, gwasSumStats = gs,
                       xqtlColoc = FALSE, jointGwas = FALSE,
                       separateGwas = TRUE))
  # Driver merges QTL + GWAS sumstats into a single bundle and the
  # separate-loop iterates over every merged study label (Q1:c1:t1 + G1).
  expect_equal(callCount, 2L)
  expect_true(!is.null(out$separate_gwas))
})

test_that("colocboostPipeline: no analysis flag set emits a message and returns empty", {
  ss <- .cbp_makeQtlSumStats()
  local_mocked_bindings(extractBlockGenotypes = .cbp_mockExtractor(),
                        .package = "pecotmr")
  out <- suppressMessages(
    colocboostPipeline(ss, xqtlColoc = FALSE, jointGwas = FALSE,
                       separateGwas = FALSE))
  expect_null(out$xqtl_coloc)
  expect_null(out$joint_gwas)
  expect_null(out$separate_gwas)
})

test_that("colocboostPipeline: GWAS ldSketch mismatch errors during the driver", {
  ss <- .cbp_makeQtlSumStats()
  # Build a GwasSumStats whose ldSketch has a different sample set.
  gh_diff <- .cbp_makeHandle(sample_prefix = "z")
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = rnorm(5), N = rep(1000L, 5))
  gs <- GwasSumStats(
    study = "G1", entry = list(gr), genome = "hg19",
    ldSketch = gh_diff, qcInfo = list(step1 = "ok"))
  # .cbQtlSumStatsBundle reads the qtl sketch via extractBlockGenotypes
  # before the LD-sketch mismatch check fires further down the driver, so
  # mock the extractor here too.
  local_mocked_bindings(extractBlockGenotypes = .cbp_mockExtractor(),
                        .package = "pecotmr")
  expect_error(
    suppressMessages(
      colocboostPipeline(ss, gwasSumStats = gs,
                         xqtlColoc = FALSE, jointGwas = TRUE)),
    "different sample sets"
  )
})
