context("QtlDataset residualization methods")

# ===========================================================================
# Fixture: mirrors test_qtlDatasetHelpers.R so we can mock extractBlockGenotypes
# the same way.
# ===========================================================================

.qr_makeHandle <- function(snp_n = 6L, n_samples = 12L) {
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
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

.qr_makeSe <- function(traits = c("ENSG1", "ENSG2"), n_samples = 12,
                       starts = NULL, chr = "chr1") {
  if (is.null(starts)) starts <- seq(1000L, by = 1000L, length.out = length(traits))
  rng <- GenomicRanges::GRanges(
    seqnames = rep(chr, length(traits)),
    ranges = IRanges::IRanges(start = starts, width = 500L))
  names(rng) <- traits
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  # Use numeric covariates only — .qtlBuildResidualizationDesign coerces the
  # full colData via as.matrix(as.data.frame(...)), so character columns
  # would coerce to NA and break lm.fit downstream.
  cd <- S4Vectors::DataFrame(sex = rep(c(0, 1), length.out = n_samples),
                             age = seq_len(n_samples),
                             row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr),
    rowRanges = rng,
    colData = cd)
}

.qr_makeDataset <- function(contexts = c("brain", "liver"),
                            n_samples = 12L, geno_cov = NULL,
                            scaleResiduals = TRUE) {
  gh <- .qr_makeHandle(n_samples = n_samples)
  pheno <- setNames(lapply(contexts, function(.) .qr_makeSe(n_samples = n_samples)),
                    contexts)
  if (is.null(geno_cov)) {
    geno_cov <- matrix(numeric(0), nrow = 0, ncol = 0)
  }
  QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = pheno,
    genotypeCovariates = geno_cov,
    scaleResiduals     = scaleResiduals)
}

.qr_mockExtractor <- function(seed = 42, n_samples = 12L, n_snp = 6L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * n_snp, 2, 0.3),
                    nrow = n_samples, ncol = n_snp,
                    dimnames = list(handle@sampleIds,
                                    handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx], width = 1L))
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

# ===========================================================================
# .qtlResolveResidualizationFlag (pure helper used by both methods)
# ===========================================================================

test_that(".qtlResolveResidualizationFlag: both missing returns TRUE", {
  res <- pecotmr:::.qtlResolveResidualizationFlag(
    conveniencePassed = NA, convenienceMissing = TRUE,
    precisePassed     = NA, preciseMissing     = TRUE,
    convenienceName = "conv", preciseName = "prec")
  expect_true(res)
})

test_that(".qtlResolveResidualizationFlag: only convenience set returns that value", {
  expect_true(pecotmr:::.qtlResolveResidualizationFlag(
    TRUE, FALSE, NA, TRUE, "conv", "prec"))
  expect_false(pecotmr:::.qtlResolveResidualizationFlag(
    FALSE, FALSE, NA, TRUE, "conv", "prec"))
})

test_that(".qtlResolveResidualizationFlag: only precise set returns that value", {
  expect_true(pecotmr:::.qtlResolveResidualizationFlag(
    NA, TRUE, TRUE, FALSE, "conv", "prec"))
  expect_false(pecotmr:::.qtlResolveResidualizationFlag(
    NA, TRUE, FALSE, FALSE, "conv", "prec"))
})

test_that(".qtlResolveResidualizationFlag: both set + agreeing returns the shared value", {
  expect_true(pecotmr:::.qtlResolveResidualizationFlag(
    TRUE, FALSE, TRUE, FALSE, "conv", "prec"))
})

test_that(".qtlResolveResidualizationFlag: both set + conflicting errors", {
  expect_error(
    pecotmr:::.qtlResolveResidualizationFlag(
      TRUE, FALSE, FALSE, FALSE, "conv", "prec"),
    "Conflicting values: `conv`"
  )
})

# ===========================================================================
# getResidualizedGenotypes (QtlDataset)
# ===========================================================================

test_that("getResidualizedGenotypes: requires contexts", {
  qd <- .qr_makeDataset()
  expect_error(getResidualizedGenotypes(qd),
               "`contexts` is required")
  expect_error(getResidualizedGenotypes(qd, contexts = NULL),
               "`contexts` is required")
  expect_error(getResidualizedGenotypes(qd, contexts = character(0)),
               "`contexts` is required")
})

test_that("getResidualizedGenotypes: unknown context errors", {
  qd <- .qr_makeDataset()
  expect_error(getResidualizedGenotypes(qd, contexts = "ghost"),
               "Unknown context")
})

test_that("getResidualizedGenotypes: empty genotype block short-circuits to G", {
  qd <- .qr_makeDataset()
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  # region with no SNPs in the panel.
  region <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1, 1000))
  G <- getResidualizedGenotypes(qd, contexts = "brain", region = region)
  expect_equal(ncol(G), 0L)
})

test_that("getResidualizedGenotypes: produces residualized matrix shape", {
  qd <- .qr_makeDataset(contexts = "brain")
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  G <- getResidualizedGenotypes(qd, contexts = "brain")
  expect_equal(nrow(G), 12L)
  expect_equal(ncol(G), 6L)
  # When scaleResiduals = TRUE (the default), kept columns should have unit sd
  # (constant columns are clamped to zero in .qtlResidualizeQR).
  sds <- apply(G, 2L, sd)
  nonZero <- sds > 1e-6
  expect_true(all(abs(sds[nonZero] - 1) < 1e-6))
})

test_that("getResidualizedGenotypes: residualizes only against selected pheno covariate", {
  qd <- .qr_makeDataset(contexts = "brain")
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  G <- getResidualizedGenotypes(
    qd, contexts = "brain",
    phenotypeCovariatesToResidualize = "age")
  expect_equal(nrow(G), 12L)
  # Resulting columns should be uncorrelated with 'age'.
  age <- seq_len(12)
  for (j in seq_len(ncol(G))) {
    expect_lt(abs(cor(G[, j], age)), 1e-6)
  }
})

test_that("getResidualizedGenotypes: respects residualizePhenotypeCovariates = FALSE", {
  qd <- .qr_makeDataset(contexts = "brain")
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  # When pheno is disabled, the design becomes intercept-only, so the result
  # is just the centered (and scaled) raw block.
  G1 <- getResidualizedGenotypes(qd, contexts = "brain",
                                  residualizePhenotypeCovariates = FALSE)
  expect_equal(nrow(G1), 12L)
  expect_equal(ncol(G1), 6L)
})

test_that("getResidualizedGenotypes: precise-name kwarg routes correctly", {
  qd <- .qr_makeDataset(contexts = "brain")
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  G_precise <- getResidualizedGenotypes(
    qd, contexts = "brain",
    residualizePhenotypeCovariatesFromGenotypes = FALSE)
  G_conv <- getResidualizedGenotypes(
    qd, contexts = "brain",
    residualizePhenotypeCovariates = FALSE)
  expect_equal(G_precise, G_conv)
})

test_that("getResidualizedGenotypes: conflict between convenience and precise errors", {
  qd <- .qr_makeDataset(contexts = "brain")
  expect_error(
    getResidualizedGenotypes(
      qd, contexts = "brain",
      residualizePhenotypeCovariates = TRUE,
      residualizePhenotypeCovariatesFromGenotypes = FALSE),
    "Conflicting values"
  )
})

test_that("getResidualizedGenotypes: joint-context mode intersects samples", {
  qd <- .qr_makeDataset(contexts = c("brain", "liver"))
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  G <- getResidualizedGenotypes(qd, contexts = c("brain", "liver"))
  expect_equal(nrow(G), 12L)
  expect_setequal(rownames(G), paste0("s", 1:12))
})

test_that("getResidualizedGenotypes: includes genotype covariates when supplied", {
  gc <- matrix(rnorm(12 * 2), nrow = 12, ncol = 2,
               dimnames = list(paste0("s", 1:12), c("pc1", "pc2")))
  qd <- .qr_makeDataset(contexts = "brain", geno_cov = gc)
  local_mocked_bindings(extractBlockGenotypes = .qr_mockExtractor(),
                        .package = "pecotmr")
  G <- getResidualizedGenotypes(qd, contexts = "brain",
                                 genotypeCovariatesToResidualize = c("pc1", "pc2"))
  # Columns should be uncorrelated with the included PCs.
  for (j in seq_len(ncol(G))) {
    expect_lt(abs(cor(G[, j], gc[rownames(G), 1])), 1e-6)
    expect_lt(abs(cor(G[, j], gc[rownames(G), 2])), 1e-6)
  }
})

# ===========================================================================
# getResidualizedPhenotypes (QtlDataset)
# ===========================================================================

test_that("getResidualizedPhenotypes: requires contexts", {
  qd <- .qr_makeDataset()
  expect_error(getResidualizedPhenotypes(qd),
               "`contexts` is required")
})

test_that("getResidualizedPhenotypes: unknown context errors", {
  qd <- .qr_makeDataset()
  expect_error(getResidualizedPhenotypes(qd, contexts = "ghost"),
               "Unknown context")
})

test_that("getResidualizedPhenotypes: returns one matrix per context", {
  qd <- .qr_makeDataset(contexts = c("brain", "liver"))
  res <- getResidualizedPhenotypes(qd, contexts = c("brain", "liver"))
  expect_equal(names(res), c("brain", "liver"))
  expect_equal(nrow(res$brain), 12L)
  expect_equal(ncol(res$brain), 2L)
  expect_equal(nrow(res$liver), 12L)
})

test_that("getResidualizedPhenotypes: residualizes against age covariate", {
  qd <- .qr_makeDataset(contexts = "brain")
  res <- getResidualizedPhenotypes(qd, contexts = "brain",
                                    phenotypeCovariatesToResidualize = "age")
  Y <- res$brain
  age <- seq_len(12)
  for (j in seq_len(ncol(Y))) {
    expect_lt(abs(cor(Y[, j], age)), 1e-6)
  }
})

test_that("getResidualizedPhenotypes: respects residualizePhenotypeCovariates = FALSE", {
  qd <- .qr_makeDataset(contexts = "brain")
  res <- getResidualizedPhenotypes(qd, contexts = "brain",
                                    residualizePhenotypeCovariates = FALSE)
  Y <- res$brain
  expect_equal(nrow(Y), 12L)
  expect_equal(ncol(Y), 2L)
})

test_that("getResidualizedPhenotypes: precise-name kwarg routes correctly", {
  qd <- .qr_makeDataset(contexts = "brain")
  Y_precise <- getResidualizedPhenotypes(
    qd, contexts = "brain",
    residualizePhenotypeCovariatesFromPhenotypes = FALSE)
  Y_conv <- getResidualizedPhenotypes(
    qd, contexts = "brain",
    residualizePhenotypeCovariates = FALSE)
  expect_equal(Y_precise$brain, Y_conv$brain)
})

test_that("getResidualizedPhenotypes: traitId subsets to requested traits", {
  qd <- .qr_makeDataset(contexts = "brain")
  res <- getResidualizedPhenotypes(qd, contexts = "brain", traitId = "ENSG1")
  expect_equal(ncol(res$brain), 1L)
  expect_equal(colnames(res$brain), "ENSG1")
})

test_that("getResidualizedPhenotypes: scaleResiduals = FALSE skips the rescale step", {
  qd <- .qr_makeDataset(contexts = "brain", scaleResiduals = FALSE)
  res <- getResidualizedPhenotypes(qd, contexts = "brain")
  # Without scaling the residual columns generally won't have sd = 1.
  Y <- res$brain
  sds <- apply(Y, 2L, sd)
  expect_false(any(abs(sds - 1) < 1e-6))
})
