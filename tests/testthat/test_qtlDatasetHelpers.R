context("QtlDataset internal helpers")

# ===========================================================================
# Fixture builder: a QtlDataset whose GenotypeHandle is a stub. extraction
# functions are stubbed in individual tests via local_mocked_bindings.
# ===========================================================================

.qh_makeHandle <- function(snp_n = 6L, n_samples = 12L) {
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

.qh_makeSe <- function(traits = c("ENSG1", "ENSG2"), n_samples = 12,
                       starts = NULL,
                       chr = "chr1",
                       extra_cov = NULL) {
  if (is.null(starts)) starts <- seq(1000L, by = 1000L, length.out = length(traits))
  rng <- GenomicRanges::GRanges(
    seqnames = rep(chr, length(traits)),
    ranges = IRanges::IRanges(start = starts, width = 500L))
  names(rng) <- traits
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  cd_list <- list(sex = rep(c("M", "F"), length.out = n_samples),
                  age = seq_len(n_samples))
  if (!is.null(extra_cov)) cd_list <- c(cd_list, extra_cov)
  cd <- S4Vectors::DataFrame(cd_list,
                             row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr),
    rowRanges = rng,
    colData = cd)
}

.qh_makeDataset <- function(contexts = c("brain", "liver"),
                            n_samples = 12L,
                            geno_cov = NULL) {
  gh <- .qh_makeHandle(n_samples = n_samples)
  pheno <- setNames(lapply(contexts, function(.) .qh_makeSe(n_samples = n_samples)),
                    contexts)
  if (is.null(geno_cov)) {
    geno_cov <- matrix(numeric(0), nrow = 0, ncol = 0)
  }
  QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = pheno,
    genotypeCovariates = geno_cov)
}

# ===========================================================================
# .qtlResidualizeQR — pure linear algebra
# ===========================================================================

test_that(".qtlResidualizeQR: intercept-only residualization centers Y", {
  set.seed(0)
  Y <- matrix(rnorm(20) + 5, nrow = 10, ncol = 2)
  res <- pecotmr:::.qtlResidualizeQR(Y, C = NULL, scaleResiduals = FALSE)
  # After removing the intercept, columns should have zero mean.
  expect_equal(unname(colMeans(res)), c(0, 0), tolerance = 1e-10)
})

test_that(".qtlResidualizeQR: covariate residualization removes the covariate signal", {
  set.seed(1)
  n <- 50
  C <- matrix(rnorm(n * 2), nrow = n, ncol = 2,
              dimnames = list(NULL, c("c1", "c2")))
  # Y = 0.5 * c1 - 0.3 * c2 + noise
  Y <- matrix(0.5 * C[, 1] - 0.3 * C[, 2] + rnorm(n, sd = 0.1),
              nrow = n, ncol = 1)
  res <- pecotmr:::.qtlResidualizeQR(Y, C = C, scaleResiduals = FALSE)
  # Residuals should be near-zero (only contain the noise).
  expect_lt(max(abs(res)), 0.5)
  # And uncorrelated with the covariates.
  expect_lt(abs(cor(res[, 1], C[, 1])), 1e-8)
  expect_lt(abs(cor(res[, 1], C[, 2])), 1e-8)
})

test_that(".qtlResidualizeQR: scaleResiduals = TRUE gives unit variance per column", {
  set.seed(2)
  Y <- matrix(rnorm(30), nrow = 10, ncol = 3)
  res <- pecotmr:::.qtlResidualizeQR(Y, C = NULL, scaleResiduals = TRUE)
  sds <- apply(res, 2, sd)
  expect_equal(sds, c(1, 1, 1), tolerance = 1e-10)
})

test_that(".qtlResidualizeQR: constant residual columns survive the rescale step", {
  # Y is exactly its own mean -> residuals are 0, sd is 0 (and clamped to 1).
  Y <- matrix(5, nrow = 5, ncol = 1)
  res <- pecotmr:::.qtlResidualizeQR(Y, C = NULL, scaleResiduals = TRUE)
  expect_true(all(abs(res) < 1e-10))
})

test_that(".qtlResidualizeQR: rank-deficient covariates are dropped by pivoted QR", {
  set.seed(3)
  n <- 30
  c1 <- rnorm(n)
  C <- matrix(cbind(c1, 2 * c1, rnorm(n)), nrow = n,
              dimnames = list(NULL, c("a", "b", "c")))
  Y <- matrix(rnorm(n), nrow = n, ncol = 1)
  # Even though `a` and `b` are collinear, the QR should not error.
  expect_no_error(pecotmr:::.qtlResidualizeQR(Y, C = C, scaleResiduals = FALSE))
})

# ===========================================================================
# .qtlResolveVariantRegion
# ===========================================================================

test_that(".qtlResolveVariantRegion: both traitId and region errors", {
  qd <- .qh_makeDataset()
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, traitId = "ENSG1", region = region),
    "Specify either `traitId` or `region`, not both"
  )
})

test_that(".qtlResolveVariantRegion: neither argument returns NULL", {
  qd <- .qh_makeDataset()
  expect_null(pecotmr:::.qtlResolveVariantRegion(qd))
})

test_that(".qtlResolveVariantRegion: traitId requires cisWindow", {
  qd <- .qh_makeDataset()
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, traitId = "ENSG1"),
    "`cisWindow` is required"
  )
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, traitId = "ENSG1", cisWindow = -1),
    "non-negative"
  )
})

test_that(".qtlResolveVariantRegion: traitId expands by cisWindow on each side", {
  qd <- .qh_makeDataset()
  gr <- pecotmr:::.qtlResolveVariantRegion(qd, traitId = "ENSG1",
                                            cisWindow = 200L)
  # ENSG1 spans 1000-1499. With cisWindow=200, span is 800-1699.
  expect_equal(as.character(GenomicRanges::seqnames(gr)), "chr1")
  expect_equal(GenomicRanges::start(gr), 800L)
  expect_equal(GenomicRanges::end(gr), 1699L)
})

test_that(".qtlResolveVariantRegion: traitId span is clipped at 1", {
  qd <- .qh_makeDataset()
  # ENSG1 starts at 1000; a 5000-bp window would push us below 1.
  gr <- pecotmr:::.qtlResolveVariantRegion(qd, traitId = "ENSG1",
                                            cisWindow = 5000L)
  expect_equal(GenomicRanges::start(gr), 1L)
})

test_that(".qtlResolveVariantRegion: unknown trait errors", {
  qd <- .qh_makeDataset()
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, traitId = "GHOST", cisWindow = 0L),
    "None of the requested traitId values were found"
  )
})

test_that(".qtlResolveVariantRegion: traits across chromosomes error", {
  # Build a dataset where two contexts hold traits on different chromosomes
  # but share a name (validity allows this — names differ).
  gh <- .qh_makeHandle()
  se1 <- .qh_makeSe(traits = "ENSG_A", chr = "chr1")
  se2 <- .qh_makeSe(traits = "ENSG_B", chr = "chr2")
  qd <- QtlDataset(study = "s1", genotypes = gh,
                   phenotypes = list(brain = se1, liver = se2),
                   genotypeCovariates = matrix(0, nrow = 12, ncol = 0))
  # Combining single-row GRanges from chr1 and chr2 emits a Bioconductor
  # warning about disjoint seqlevels — that's exactly the cross-chromosome
  # case we are exercising, so suppress it.
  expect_error(
    suppressWarnings(pecotmr:::.qtlResolveVariantRegion(
      qd, traitId = c("ENSG_A", "ENSG_B"), cisWindow = 0L)),
    "share a chromosome"
  )
})

test_that(".qtlResolveVariantRegion: region path requires single-range GRanges", {
  qd <- .qh_makeDataset()
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, region = "chr1:100-200"),
    "must be a GRanges object"
  )
  multi <- GenomicRanges::GRanges(c("chr1", "chr1"),
                                  IRanges::IRanges(c(1, 100), c(50, 200)))
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, region = multi),
    "single range"
  )
})

test_that(".qtlResolveVariantRegion: region path expands by cisWindow", {
  qd <- .qh_makeDataset()
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(500, 1000))
  gr <- pecotmr:::.qtlResolveVariantRegion(qd, region = region,
                                            cisWindow = 250L)
  expect_equal(GenomicRanges::start(gr), 250L)
  expect_equal(GenomicRanges::end(gr), 1250L)
})

# ===========================================================================
# .qtlVariantIndices
# ===========================================================================

test_that(".qtlVariantIndices: NULL region returns all SNP indices", {
  qd <- .qh_makeDataset()
  idx <- pecotmr:::.qtlVariantIndices(qd)
  expect_equal(idx, seq_len(nrow(qd@genotypes@snpInfo)))
})

test_that(".qtlVariantIndices: filters by chromosome and BP range", {
  qd <- .qh_makeDataset()
  # The handle has SNPs at chr1:100, 200, ..., 600.
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 350))
  idx <- pecotmr:::.qtlVariantIndices(qd, region = region)
  expect_equal(idx, c(2L, 3L))
})

test_that(".qtlVariantIndices: accepts chr-prefixed and bare chromosome names", {
  qd <- .qh_makeDataset()
  r1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 250))
  r2 <- GenomicRanges::GRanges("1",    IRanges::IRanges(50, 250))
  expect_equal(pecotmr:::.qtlVariantIndices(qd, r1),
               pecotmr:::.qtlVariantIndices(qd, r2))
})

test_that(".qtlVariantIndices: returns integer(0) when no overlap", {
  qd <- .qh_makeDataset()
  region <- GenomicRanges::GRanges("chr2", IRanges::IRanges(100, 200))
  expect_equal(pecotmr:::.qtlVariantIndices(qd, region), integer(0))
})

# ===========================================================================
# .qtlResolvePhenoSelection
# ===========================================================================

test_that(".qtlResolvePhenoSelection: NULL returns all colData columns per context", {
  qd <- .qh_makeDataset(contexts = c("brain", "liver"))
  out <- pecotmr:::.qtlResolvePhenoSelection(qd,
                                              contexts = c("brain", "liver"),
                                              toResidualize = NULL)
  expect_equal(names(out), c("brain", "liver"))
  expect_setequal(out$brain, c("sex", "age"))
  expect_setequal(out$liver, c("sex", "age"))
})

test_that(".qtlResolvePhenoSelection: character vector applies to all contexts", {
  qd <- .qh_makeDataset(contexts = c("brain", "liver"))
  out <- pecotmr:::.qtlResolvePhenoSelection(qd,
                                              contexts = c("brain", "liver"),
                                              toResidualize = "age")
  expect_equal(out$brain, "age")
  expect_equal(out$liver, "age")
})

test_that(".qtlResolvePhenoSelection: character vector with unknown name errors", {
  qd <- .qh_makeDataset(contexts = "brain")
  expect_error(
    pecotmr:::.qtlResolvePhenoSelection(qd, contexts = "brain",
                                         toResidualize = "ghost"),
    "no covariate.*ghost"
  )
})

test_that(".qtlResolvePhenoSelection: named list dispatches per context", {
  qd <- .qh_makeDataset(contexts = c("brain", "liver"))
  out <- pecotmr:::.qtlResolvePhenoSelection(qd,
                                              contexts = c("brain", "liver"),
                                              toResidualize = list(brain = "age",
                                                                   liver = "sex"))
  expect_equal(out$brain, "age")
  expect_equal(out$liver, "sex")
})

test_that(".qtlResolvePhenoSelection: list missing keys errors", {
  qd <- .qh_makeDataset(contexts = c("brain", "liver"))
  expect_error(
    pecotmr:::.qtlResolvePhenoSelection(qd,
                                         contexts = c("brain", "liver"),
                                         toResidualize = list(brain = "age")),
    "list does not cover all"
  )
})

test_that(".qtlResolvePhenoSelection: list with extra keys errors", {
  qd <- .qh_makeDataset(contexts = "brain")
  expect_error(
    pecotmr:::.qtlResolvePhenoSelection(qd,
                                         contexts = "brain",
                                         toResidualize = list(brain = "age",
                                                              ghost = "sex")),
    "list key.*not in `contexts`"
  )
})

test_that(".qtlResolvePhenoSelection: unnamed list errors", {
  qd <- .qh_makeDataset(contexts = "brain")
  expect_error(
    pecotmr:::.qtlResolvePhenoSelection(qd,
                                         contexts = "brain",
                                         toResidualize = list("age")),
    "must be named"
  )
})

test_that(".qtlResolvePhenoSelection: unsupported type errors", {
  qd <- .qh_makeDataset(contexts = "brain")
  expect_error(
    pecotmr:::.qtlResolvePhenoSelection(qd,
                                         contexts = "brain",
                                         toResidualize = 42L),
    "must be NULL, a character vector, or a named list"
  )
})

# ===========================================================================
# .qtlResolveGenoSelection
# ===========================================================================

test_that(".qtlResolveGenoSelection: NULL returns all genotypeCovariates columns", {
  gc <- matrix(rnorm(12 * 3), nrow = 12, ncol = 3,
               dimnames = list(paste0("s", 1:12), c("pc1", "pc2", "pc3")))
  qd <- .qh_makeDataset(contexts = "brain", geno_cov = gc)
  expect_setequal(pecotmr:::.qtlResolveGenoSelection(qd, toResidualize = NULL),
                  c("pc1", "pc2", "pc3"))
})

test_that(".qtlResolveGenoSelection: empty genotypeCovariates returns character(0)", {
  qd <- .qh_makeDataset(contexts = "brain")
  expect_equal(pecotmr:::.qtlResolveGenoSelection(qd, toResidualize = NULL),
               character(0))
})

test_that(".qtlResolveGenoSelection: subset selection works", {
  gc <- matrix(0, nrow = 12, ncol = 3,
               dimnames = list(paste0("s", 1:12), c("pc1", "pc2", "pc3")))
  qd <- .qh_makeDataset(contexts = "brain", geno_cov = gc)
  expect_equal(pecotmr:::.qtlResolveGenoSelection(qd,
                                                   toResidualize = c("pc1", "pc3")),
               c("pc1", "pc3"))
})

test_that(".qtlResolveGenoSelection: unknown name errors", {
  gc <- matrix(0, nrow = 12, ncol = 2,
               dimnames = list(paste0("s", 1:12), c("pc1", "pc2")))
  qd <- .qh_makeDataset(contexts = "brain", geno_cov = gc)
  expect_error(
    pecotmr:::.qtlResolveGenoSelection(qd, toResidualize = "pc99"),
    "no covariate.*pc99"
  )
})

# ===========================================================================
# .qtlBuildResidualizationDesign
# ===========================================================================

test_that(".qtlBuildResidualizationDesign: pheno-only single-context builds the colData matrix", {
  qd <- .qh_makeDataset(contexts = "brain")
  phenoSel <- list(brain = c("age", "sex"))
  D <- pecotmr:::.qtlBuildResidualizationDesign(
    qd, contexts = "brain",
    phenoSelection = phenoSel,
    genoSelection  = character(0),
    includePheno = TRUE, includeGeno = FALSE)
  expect_true(is.matrix(D))
  expect_equal(nrow(D), 12L)
  expect_setequal(colnames(D), c("brain.age", "brain.sex"))
})

test_that(".qtlBuildResidualizationDesign: pheno-only multi-context concatenates per context", {
  qd <- .qh_makeDataset(contexts = c("brain", "liver"))
  phenoSel <- list(brain = "age", liver = "sex")
  D <- pecotmr:::.qtlBuildResidualizationDesign(
    qd, contexts = c("brain", "liver"),
    phenoSelection = phenoSel,
    genoSelection  = character(0),
    includePheno = TRUE, includeGeno = FALSE)
  expect_equal(ncol(D), 2L)
  expect_setequal(colnames(D), c("brain.age", "liver.sex"))
})

test_that(".qtlBuildResidualizationDesign: includeGeno-only returns genotype covariates", {
  gc <- matrix(rnorm(12 * 2), nrow = 12, ncol = 2,
               dimnames = list(paste0("s", 1:12), c("pc1", "pc2")))
  qd <- .qh_makeDataset(contexts = "brain", geno_cov = gc)
  D <- pecotmr:::.qtlBuildResidualizationDesign(
    qd, contexts = "brain",
    phenoSelection = list(brain = character(0)),
    genoSelection  = c("pc1", "pc2"),
    includePheno = FALSE, includeGeno = TRUE)
  expect_equal(ncol(D), 2L)
  expect_setequal(colnames(D), c("pc1", "pc2"))
})

test_that(".qtlBuildResidualizationDesign: pheno + geno concatenates both blocks", {
  gc <- matrix(rnorm(12 * 2), nrow = 12, ncol = 2,
               dimnames = list(paste0("s", 1:12), c("pc1", "pc2")))
  qd <- .qh_makeDataset(contexts = "brain", geno_cov = gc)
  D <- pecotmr:::.qtlBuildResidualizationDesign(
    qd, contexts = "brain",
    phenoSelection = list(brain = "age"),
    genoSelection  = "pc1",
    includePheno = TRUE, includeGeno = TRUE)
  expect_equal(ncol(D), 2L)
  expect_setequal(colnames(D), c("brain.age", "pc1"))
})

test_that(".qtlBuildResidualizationDesign: returns NULL when nothing to include", {
  qd <- .qh_makeDataset(contexts = "brain")
  D <- pecotmr:::.qtlBuildResidualizationDesign(
    qd, contexts = "brain",
    phenoSelection = list(brain = character(0)),
    genoSelection  = character(0),
    includePheno = FALSE, includeGeno = FALSE)
  expect_null(D)
})

test_that(".qtlBuildResidualizationDesign: intersects sample sets across blocks", {
  # Build a dataset where the genotype covariates only cover samples s1..s6.
  gc <- matrix(rnorm(6 * 1), nrow = 6, ncol = 1,
               dimnames = list(paste0("s", 1:6), "pc1"))
  qd <- .qh_makeDataset(contexts = "brain", geno_cov = gc)
  D <- pecotmr:::.qtlBuildResidualizationDesign(
    qd, contexts = "brain",
    phenoSelection = list(brain = "age"),
    genoSelection  = "pc1",
    includePheno = TRUE, includeGeno = TRUE)
  expect_equal(nrow(D), 6L)
  expect_setequal(rownames(D), paste0("s", 1:6))
})

# ===========================================================================
# .qtlExtractBlock (uses mocked extractBlockGenotypes)
# ===========================================================================

# Build a stub extractBlockGenotypes that returns a synthetic SE for the
# requested snpIdx, drawing dosages from a per-handle-deterministic seed so
# the same indices always give the same numbers.
.qh_mockExtractor <- function(seed = 42, n_samples = 12L, n_snp = 6L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    # Build a (n_samples x n_snp) dosage matrix for the whole panel; the
    # caller's snpIdx subsets columns.
    panel <- matrix(rbinom(n_samples * n_snp, 2, 0.3),
                    nrow = n_samples, ncol = n_snp,
                    dimnames = list(handle@sampleIds,
                                    handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    # Build the SE in variants x samples orientation (matches the real impl).
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

test_that(".qtlExtractBlock: returns dosage matrix with kept variants and samples", {
  qd <- .qh_makeDataset()
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd)
  expect_equal(nrow(blk$geno), 12L)
  expect_equal(ncol(blk$geno), 6L)
  expect_equal(blk$variantIds, paste0("rs", 1:6))
  expect_equal(blk$sampleIds, paste0("s", 1:12))
  expect_equal(length(blk$maf), 6L)
})

test_that(".qtlExtractBlock: empty snpIdx returns a zero-column block", {
  qd <- .qh_makeDataset()
  region <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1, 1000))
  blk <- pecotmr:::.qtlExtractBlock(qd, region = region)
  expect_equal(ncol(blk$geno), 0L)
  expect_equal(blk$variantIds, character(0))
})

test_that(".qtlExtractBlock: keepVariants restriction narrows the returned set", {
  qd <- .qh_makeDataset()
  qd@keepVariants <- c("rs2", "rs4")
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd)
  expect_setequal(blk$variantIds, c("rs2", "rs4"))
})

test_that(".qtlExtractBlock: keepSamples restriction narrows the returned set", {
  qd <- .qh_makeDataset()
  qd@keepSamples <- paste0("s", 1:6)
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd)
  expect_setequal(blk$sampleIds, paste0("s", 1:6))
})

test_that(".qtlExtractBlock: per-call samples arg further narrows the sample set", {
  qd <- .qh_makeDataset()
  qd@keepSamples <- paste0("s", 1:6)
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd, samples = c("s1", "s3", "s5"))
  expect_setequal(blk$sampleIds, c("s1", "s3", "s5"))
})

test_that(".qtlExtractBlock: keepVariants with empty intersection returns empty block", {
  qd <- .qh_makeDataset()
  qd@keepVariants <- c("rsGHOST")
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd)
  expect_equal(ncol(blk$geno), 0L)
  expect_equal(nrow(blk$geno), 0L)
})

test_that(".qtlExtractBlock: mafCutoff drops low-MAF variants", {
  qd <- .qh_makeDataset()
  # The mocked extractor returns binomial(0.3) dosages: realized MAFs hover
  # around 0.3-0.5 (small sample noise). Cutoff above the realized maximum
  # drops everything.
  qd@mafCutoff <- 0.51
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd)
  expect_equal(ncol(blk$geno), 0L)
})

test_that(".qtlExtractBlock: mafCutoff retains variants above the threshold", {
  qd <- .qh_makeDataset()
  qd@mafCutoff <- 0.4  # realized MAFs include 0.458 and 0.5
  local_mocked_bindings(
    extractBlockGenotypes = .qh_mockExtractor(),
    .package = "pecotmr")
  blk <- pecotmr:::.qtlExtractBlock(qd)
  expect_true(ncol(blk$geno) >= 1L)
  expect_true(all(blk$maf >= 0.4))
})
