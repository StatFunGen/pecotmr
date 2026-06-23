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

test_that(".qtlResolveVariantRegion: region must be a GRanges; multi-range is allowed", {
  qd <- .qh_makeDataset()
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, region = "chr1:100-200"),
    "must be a GRanges object"
  )
  expect_error(
    pecotmr:::.qtlResolveVariantRegion(qd, region = GenomicRanges::GRanges()),
    "at least one range"
  )
  # A multi-range region is now taken literally (joint multi-region extraction).
  multi <- GenomicRanges::GRanges(c("chr1", "chr1"),
                                  IRanges::IRanges(c(1, 100), c(50, 200)))
  gr <- pecotmr:::.qtlResolveVariantRegion(qd, region = multi)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2L)
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
  Y <- getResidualizedPhenotypes(qd, contexts = "brain",
                                  phenotypeCovariatesToResidualize = "age")
  age <- seq_len(12)
  for (j in seq_len(ncol(Y))) {
    expect_lt(abs(cor(Y[, j], age)), 1e-6)
  }
})

test_that("getResidualizedPhenotypes: respects residualizePhenotypeCovariates = FALSE", {
  qd <- .qr_makeDataset(contexts = "brain")
  Y <- getResidualizedPhenotypes(qd, contexts = "brain",
                                  residualizePhenotypeCovariates = FALSE)
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
  expect_equal(Y_precise, Y_conv)
})

test_that("getResidualizedPhenotypes: traitId subsets to requested traits", {
  qd <- .qr_makeDataset(contexts = "brain")
  Y <- getResidualizedPhenotypes(qd, contexts = "brain", traitId = "ENSG1")
  expect_equal(ncol(Y), 1L)
  expect_equal(colnames(Y), "ENSG1")
})

test_that("getResidualizedPhenotypes: scaleResiduals = FALSE skips the rescale step", {
  qd <- .qr_makeDataset(contexts = "brain", scaleResiduals = FALSE)
  Y <- getResidualizedPhenotypes(qd, contexts = "brain")
  # Without scaling the residual columns generally won't have sd = 1.
  sds <- apply(Y, 2L, sd)
  expect_false(any(abs(sds - 1) < 1e-6))
})

# ===========================================================================
# getPhenotypes/getResidualizedPhenotypes naAction
# ===========================================================================

# Helper to build a single-context SE with controlled NA placement.
.qr_makeSeWithNa <- function(n_samples = 8L, na_idx = c(2L, 5L)) {
  traits <- c("ENSG1", "ENSG2")
  rng <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1000L, 2000L), width = 500L))
  names(rng) <- traits
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  # Sprinkle NAs in the first trait at na_idx samples
  expr[1L, na_idx] <- NA_real_
  cd <- S4Vectors::DataFrame(sex = rep(c(0, 1), length.out = n_samples),
                             age = seq_len(n_samples),
                             row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays    = list(expression = expr),
    rowRanges = rng,
    colData   = cd)
}

test_that("getPhenotypes naAction='drop' drops samples with any NA in selected traits", {
  gh <- .qr_makeHandle(n_samples = 8L)
  se <- .qr_makeSeWithNa(n_samples = 8L, na_idx = c(2L, 5L))
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(numeric(0), 0L, 0L))
  out <- getPhenotypes(qd, contexts = "brain", naAction = "drop")
  expect_s4_class(out, "SummarizedExperiment")
  expect_equal(ncol(out), 6L)
  expect_false(any(is.na(SummarizedExperiment::assay(out))))
})

test_that("getPhenotypes naAction='impute' mean-imputes NAs per trait", {
  gh <- .qr_makeHandle(n_samples = 8L)
  se <- .qr_makeSeWithNa(n_samples = 8L, na_idx = c(2L, 5L))
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(numeric(0), 0L, 0L))
  out <- getPhenotypes(qd, contexts = "brain", naAction = "impute")
  Y <- SummarizedExperiment::assay(out)
  expect_equal(ncol(Y), 8L)
  expect_false(any(is.na(Y)))
  # Imputed values equal the mean of the non-missing entries of the same trait.
  orig <- SummarizedExperiment::assay(se)
  obsMean <- mean(orig[1L, !is.na(orig[1L, ])])
  expect_equal(Y[1L, 2L], obsMean)
  expect_equal(Y[1L, 5L], obsMean)
})

test_that("getResidualizedPhenotypes naAction='impute' yields NA-free residuals", {
  gh <- .qr_makeHandle(n_samples = 8L)
  se <- .qr_makeSeWithNa(n_samples = 8L, na_idx = c(2L, 5L))
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(numeric(0), 0L, 0L))
  Y <- getResidualizedPhenotypes(qd, contexts = "brain", naAction = "impute")
  expect_false(any(is.na(Y)))
  expect_equal(nrow(Y), 8L)
})

# ===========================================================================
# getPhenotypes/getResidualizedPhenotypes outlierAction
# ===========================================================================

# Build a single-context SE with a clear multivariate outlier at sample s10.
.qr_makeSeWithOutlier <- function(n_samples = 30L, outlier_idx = 10L,
                                  outlier_z = 20) {
  traits <- c("ENSG1", "ENSG2", "ENSG3")
  rng <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1000L, 2000L, 3000L), width = 500L))
  names(rng) <- traits
  set.seed(7L)
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  # Slam a large value across all traits for the chosen sample.
  expr[, outlier_idx] <- outlier_z
  cd <- S4Vectors::DataFrame(sex = rep(c(0, 1), length.out = n_samples),
                             age = seq_len(n_samples),
                             row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays    = list(expression = expr),
    rowRanges = rng,
    colData   = cd)
}

test_that("getPhenotypes outlierAction='drop' drops a clear multivariate outlier", {
  skip_if_not_installed("robustbase")
  gh <- .qr_makeHandle(n_samples = 30L)
  se <- .qr_makeSeWithOutlier(n_samples = 30L, outlier_idx = 10L)
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(numeric(0), 0L, 0L))
  out <- getPhenotypes(qd, contexts = "brain", outlierAction = "drop")
  expect_s4_class(out, "SummarizedExperiment")
  expect_lt(ncol(out), 30L)
  expect_false("s10" %in% colnames(out))
})

test_that("getPhenotypes outlierAction='keep' is the default and a no-op", {
  gh <- .qr_makeHandle(n_samples = 30L)
  se <- .qr_makeSeWithOutlier(n_samples = 30L, outlier_idx = 10L)
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(numeric(0), 0L, 0L))
  out <- getPhenotypes(qd, contexts = "brain")
  expect_equal(ncol(out), 30L)
})

test_that("getResidualizedPhenotypes outlierAction='drop' drops residualized outliers", {
  skip_if_not_installed("robustbase")
  gh <- .qr_makeHandle(n_samples = 30L)
  se <- .qr_makeSeWithOutlier(n_samples = 30L, outlier_idx = 10L)
  qd <- QtlDataset(
    study              = "study1",
    genotypes          = gh,
    phenotypes         = list(brain = se),
    genotypeCovariates = matrix(numeric(0), 0L, 0L))
  Y <- getResidualizedPhenotypes(qd, contexts = "brain",
                                  outlierAction = "drop")
  expect_lt(nrow(Y), 30L)
  expect_false("s10" %in% rownames(Y))
})

test_that(".qtlOutlierKeepMask: single-trait reduces to z-test (drops huge value)", {
  set.seed(11L)
  Y <- matrix(c(rnorm(29), 50), ncol = 1L,
              dimnames = list(paste0("s", 1:30), "ENSG1"))
  keep <- pecotmr:::.qtlOutlierKeepMask(Y, pvalThreshold = 1e-3)
  expect_equal(length(keep), 30L)
  expect_false(keep[[30L]])
})

test_that(".qtlOutlierKeepMask: returns all-TRUE when n < p + 2", {
  Y <- matrix(rnorm(6), nrow = 2L, ncol = 3L)
  expect_warning(keep <- pecotmr:::.qtlOutlierKeepMask(Y, 1e-3),
                 "covariance estimate")
  expect_true(all(keep))
})


# Tests for functions that consume genotype matrices or compute LD:
#   computeLd, checkLd, ldPruneByCorrelation, ldClumpByScore,
#   enforceDesignFullRank, filterVariantsByLdReference,
#   resolveLdInput, dentistSingleWindow, dentist

# Fixtures: 100 samples x 100 biallelic polymorphic SNPs on chr21
test_data_dir <- test_path("test_data")
plink_prefix <- file.path(test_data_dir, "test_variants")

# Load genotype matrix once for reuse across tests
load_test_genotype <- function() {
  loadGenotypeRegion(plink_prefix, returnVariantInfo = TRUE)
}

# --- computeLd --------------------------------------------------------------

test_that("computeLd produces valid sample correlation matrix", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  R <- computeLd(geno$X, method = "sample")
  expect_true(is.matrix(R))
  expect_equal(nrow(R), ncol(geno$X))
  expect_equal(ncol(R), ncol(geno$X))
  expect_true(isSymmetric(R))
  expect_true(all(abs(diag(R) - 1) < 1e-10))
  expect_false(any(is.nan(R)))
  expect_true(all(R >= -1 - 1e-10 & R <= 1 + 1e-10))
})

test_that("computeLd population method produces valid matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "population")
  expect_true(isSymmetric(R))
  expect_true(all(abs(diag(R) - 1) < 1e-10))
  expect_false(any(is.nan(R)))
})

test_that("computeLd sample and population methods are similar", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R_s <- computeLd(X, method = "sample")
  R_p <- computeLd(X, method = "population")
  # Should be close but not identical (N-1 vs N denominator)
  expect_true(max(abs(R_s - R_p)) < 0.05)
})

test_that("computeLd errors on NULL input", {
  expect_error(computeLd(NULL), "X must be provided")
})

# --- checkLd ----------------------------------------------------------------

test_that("checkLd diagnoses real LD matrix correctly", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  result <- checkLd(R)
  expect_true(is.list(result))
  expect_true(result$isPsd)
  expect_equal(result$methodApplied, "none")
  # min eigenvalue may be near-zero (numerically PSD, not strictly PD)
  expect_true(result$minEigenvalue > -1e-7)
  expect_true(result$nNegative == 0)
  expect_true(is.finite(result$conditionNumber))
})

test_that("checkLd eigenfix improves non-PSD matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  # Make a non-PSD matrix by negating a small block of off-diagonal entries
  R_bad <- R
  R_bad[1:3, 4:6] <- -abs(R_bad[1:3, 4:6]) - 0.5
  R_bad[4:6, 1:3] <- t(R_bad[1:3, 4:6])
  diag(R_bad) <- 1

  result_check <- checkLd(R_bad, method = "check")
  expect_false(result_check$isPsd)
  expect_true(result_check$nNegative > 0)

  result_fix <- checkLd(R_bad, method = "eigenfix")
  expect_equal(result_fix$methodApplied, "eigenfix")
  # Eigenfix should improve (raise) minimum eigenvalue
  fixed_check <- checkLd(result_fix$R)
  expect_true(fixed_check$minEigenvalue > result_check$minEigenvalue)
})

test_that("checkLd shrink repairs perturbed LD matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  R_bad <- R
  R_bad[1, 2] <- R_bad[2, 1] <- 1.5
  diag(R_bad) <- 1
  result <- checkLd(R_bad, method = "shrink", shrinkage = 0.1)
  expect_equal(result$methodApplied, "shrink")
})

# --- ldPruneByCorrelation -------------------------------------------------

test_that("ldPruneByCorrelation prunes correlated variants", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- ldPruneByCorrelation(X, corThres = 0.8)
  expect_true(is.list(result))
  expect_true(is.matrix(result$X.new))
  expect_true(ncol(result$X.new) <= ncol(X))
  expect_true(ncol(result$X.new) > 0)
  expect_equal(length(result$filter.id), ncol(result$X.new))
  # Retained columns are a subset of original
  expect_true(all(result$filter.id %in% seq_len(ncol(X))))
})

test_that("ldPruneByCorrelation with strict threshold prunes more", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  loose <- ldPruneByCorrelation(X, corThres = 0.95)
  strict <- ldPruneByCorrelation(X, corThres = 0.5)
  expect_true(ncol(strict$X.new) <= ncol(loose$X.new))
})

test_that("ldPruneByCorrelation with high threshold keeps most columns", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- ldPruneByCorrelation(X, corThres = 0.999)
  # At threshold near 1, only near-duplicates are pruned; real data may have many
  expect_true(ncol(result$X.new) >= ncol(X) * 0.4)
})

# --- ldClumpByScore -------------------------------------------------------

test_that("ldClumpByScore returns valid indices", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  geno <- load_test_genotype()
  set.seed(42)
  score <- runif(ncol(geno$X))
  chr <- as.integer(geno$variant_info$chrom)
  pos <- geno$variant_info$pos
  keep <- ldClumpByScore(geno$X, score = score, chr = chr, pos = pos, r2 = 0.2)
  expect_true(is.integer(keep))
  expect_true(length(keep) > 0)
  expect_true(length(keep) <= ncol(geno$X))
  expect_true(all(keep %in% seq_len(ncol(geno$X))))
})

# --- enforceDesignFullRank ------------------------------------------------

test_that("enforceDesignFullRank handles genotype matrix with covariates", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  # Create a simple covariate matrix (e.g., first 2 PCs of X)
  pca <- prcomp(X, rank. = 2)
  C <- pca$x
  result <- enforceDesignFullRank(X, C, strategy = "correlation")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(X))
  # Should produce full-rank design
  full_design <- cbind(1, result, C)
  expect_equal(qr(full_design)$rank, ncol(full_design))
})

# --- filterVariantsByLdReference -----------------------------------------

test_that("filterVariantsByLdReference filters against PLINK reference via metadata", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  vi <- geno$variant_info
  variant_ids <- paste0(vi$chrom, ":", vi$pos, ":", vi$A2, ":", vi$A1)
  fake_ids <- c("21:999999:A:G", "21:888888:C:T")
  all_ids <- c(variant_ids, fake_ids)

  # Create a metadata TSV in the same directory as the PLINK files
  # so the relative path resolves correctly
  meta_file <- file.path(test_data_dir, "ld_meta_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  result <- suppressMessages(
    filterVariantsByLdReference(all_ids, meta_file, keepIndel = TRUE)
  )
  expect_true(is.list(result))
  expect_true(length(result$data) <= length(all_ids))
  expect_true(length(result$idx) == length(result$data))
  # Fake variants should be filtered out
  expect_true(length(result$data) <= length(variant_ids))
})

# --- resolveLdInput (internal) ---------------------------------------------

test_that("resolveLdInput computes LD from genotype matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- pecotmr:::resolveLdInput(X = X, needNSample = TRUE)
  expect_true(is.list(result))
  expect_true(is.matrix(result$R))
  expect_equal(nrow(result$R), ncol(X))
  expect_equal(result$nSample, nrow(X))
  expect_true(isSymmetric(result$R))
})

test_that("resolveLdInput passes through pre-computed R", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  result <- pecotmr:::resolveLdInput(R = R, nSample = 100L, needNSample = TRUE)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100L)
})

test_that("resolveLdInput errors when neither R nor X provided", {
  expect_error(pecotmr:::resolveLdInput(), "Either R .* or X .* must be provided")
})

test_that("resolveLdInput errors when both R and X provided", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  expect_error(pecotmr:::resolveLdInput(R = R, X = X), "Provide either R or X, not both")
})

test_that("resolveLdInput errors when R given without nSample and needed", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  expect_error(pecotmr:::resolveLdInput(R = R, needNSample = TRUE),
               "nSample is required")
})

# --- dentistSingleWindow ---------------------------------------------------

test_that("dentistSingleWindow works with genotype matrix X", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  set.seed(42)
  z <- rnorm(ncol(X))
  result <- suppressWarnings(dentistSingleWindow(z, X = X))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(z))
  expect_true("original_z" %in% names(result))
  expect_true("imputed_z" %in% names(result))
  expect_true("outlier" %in% names(result))
  expect_true(is.logical(result$outlier))
})

test_that("dentistSingleWindow works with pre-computed R", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  set.seed(42)
  z <- rnorm(ncol(X))
  result <- suppressWarnings(dentistSingleWindow(z, R = R, nSample = nrow(X)))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(z))
})

test_that("dentistSingleWindow detects injected outliers", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  set.seed(42)
  z <- rnorm(ncol(X))
  # Inject extreme outliers
  z[1] <- 50
  z[2] <- -50
  result <- suppressWarnings(dentistSingleWindow(z, R = R, nSample = nrow(X)))
  # At least one of the injected values should be flagged
  expect_true(any(result$outlier))
})

# --- dentist (multi-window) -------------------------------------------------

test_that("dentist works with genotype matrix and sum_stat data frame", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  set.seed(42)
  sum_stat <- data.frame(
    pos = geno$variant_info$pos,
    z = rnorm(ncol(geno$X))
  )
  # Use count mode with small window since we only have 100 variants
  result <- suppressWarnings(
    dentist(sum_stat, X = geno$X, windowMode = "count", minDim = 50)
  )
  expect_true(is.data.frame(result))
  # Window merging may add overlap rows; result should be >= input size
  expect_true(nrow(result) >= nrow(sum_stat))
  expect_true(all(c("original_z", "imputed_z", "outlier") %in% names(result)))
})

test_that("dentist accepts zscore column name variant", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  R <- computeLd(geno$X, method = "sample")
  set.seed(42)
  sum_stat <- data.frame(
    position = geno$variant_info$pos,
    zscore = rnorm(ncol(geno$X))
  )
  result <- suppressWarnings(
    dentist(sum_stat, R = R, nSample = nrow(geno$X), windowMode = "count", minDim = 50)
  )
  expect_true(nrow(result) >= nrow(sum_stat))
})

test_that("dentist errors when sum_stat missing required columns", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  bad_stat <- data.frame(x = 1:ncol(X), y = rnorm(ncol(X)))
  expect_error(dentist(bad_stat, X = X), "missing either")
})

# === Tests migrated from test_s4Constructors.R (QtlDataset) ===

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
# MultiStudyQtlDataset
# ===========================================================================


test_that("MultiStudyQtlDataset: rejects non-QtlDataset entries", {
  qd <- QtlDataset(study = "s1", genotypes = .sc_makeGenotypeHandle(),
                   phenotypes = list(brain = .sc_makeSe()))
  expect_error(
    MultiStudyQtlDataset(qtlDatasets = list(s1 = qd, s2 = "not a dataset")),
    "must be a QtlDataset"
  )
})



# === Tests migrated from test_showMethods.R (QtlDataset) ===

test_that("show.QtlDataset lists context names and trait count", {
  qd <- .sh_makeQtlDataset()
  out <- capture.output(show(qd))
  expect_true(any(grepl("QtlDataset for study 'study1'", out)))
  expect_true(any(grepl("1 context\\(s\\): brain", out)))
  expect_true(any(grepl("2 unique traits", out)))
  expect_true(any(grepl("Genotypes: gds @ /tmp/test.gds", out)))
})


test_that("show.MultiStudyQtlDataset reports per-source study counts", {
  qd1 <- .sh_makeQtlDataset(study = "s1")
  qd2 <- .sh_makeQtlDataset(study = "s2")
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
  out <- capture.output(show(mt))
  expect_true(any(grepl("MultiStudyQtlDataset: 2 individual-level \\+ 0 sumstats", out)))
  expect_true(any(grepl("Individual-level studies: s1, s2", out)))
})


test_that("show.MultiStudyQtlDataset reports sumstats studies when present", {
  qd <- .sh_makeQtlDataset(study = "s1")
  ss <- QtlSumStats(
    study   = "s2",
    context = "c1",
    trait   = "t1",
    entry   = list(.sh_makeQtlSumstatsGr()),
    genome  = "hg19",
    ldSketch = .sh_makeGenotypeHandle())
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd), sumStats = ss)
  out <- capture.output(show(mt))
  expect_true(any(grepl("Sumstats studies: s2", out)))
})

# ===========================================================================
# Multi-region variant extraction: getGenotypes() with a multi-range region
# (the mechanism behind jointRegions). Single-file (chr21) and sharded
# (chr21+chr22) handles wrapped in a QtlDataset with a minimal phenotype SE.
# ===========================================================================
.mr_makeSE <- function(samples, chrom = "chr21", traits = c("g1", "g2")) {
  rng <- GenomicRanges::GRanges(
    chrom, IRanges::IRanges(
      start = seq(1e6L, by = 1e5L, length.out = length(traits)), width = 1000L))
  names(rng) <- traits
  expr <- matrix(0, nrow = length(traits), ncol = length(samples),
                 dimnames = list(traits, samples))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr), rowRanges = rng,
    colData = S4Vectors::DataFrame(row.names = samples))
}
.mr_ncol <- function(qd, region) ncol(getGenotypes(qd, region = region))
.mr_vids <- function(qd, region) colnames(getGenotypes(qd, region = region))

test_that("multi-range region unions disjoint sub-ranges on one chromosome", {
  skip_if_not_installed("snpStats")
  h  <- GenotypeHandle(plink1Prefix = file.path(test_data_dir, "test_variants"))
  qd <- QtlDataset(study = "S", genotypes = h,
                   phenotypes = list(ctx = .mr_makeSE(h@sampleIds)))
  bp <- h@snpInfo$BP
  lo <- min(bp); hi <- max(bp); mid <- lo + (hi - lo) %/% 2L
  rA <- GenomicRanges::GRanges("chr21", IRanges::IRanges(lo, mid))
  rB <- GenomicRanges::GRanges("chr21", IRanges::IRanges(mid + 1L, hi))
  rAB <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(lo, mid + 1L), c(mid, hi)))

  nA <- .mr_ncol(qd, rA); nB <- .mr_ncol(qd, rB)
  expect_gt(nA, 0L); expect_gt(nB, 0L)
  expect_equal(.mr_ncol(qd, rAB), nA + nB)
  expect_equal(.mr_vids(qd, rAB), c(.mr_vids(qd, rA), .mr_vids(qd, rB)))
})

test_that("multi-range region spans chromosomes on a sharded handle", {
  skip_if_not_installed("snpStats")
  hs <- GenotypeHandle(genoMeta = c(
    "21" = file.path(test_data_dir, "test_variants"),
    "22" = file.path(test_data_dir, "test_variants_chr22")))
  qd <- QtlDataset(study = "S", genotypes = hs,
                   phenotypes = list(ctx = .mr_makeSE(hs@sampleIds)))
  bp <- GenotypeHandle(plink1Prefix = file.path(test_data_dir, "test_variants"))@snpInfo$BP
  lo <- min(bp); hi <- max(bp)
  r21 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(lo, hi))
  r22 <- GenomicRanges::GRanges("chr22", IRanges::IRanges(lo, hi))
  rBoth <- GenomicRanges::GRanges(c("chr21", "chr22"),
                                  IRanges::IRanges(c(lo, lo), c(hi, hi)))

  n21 <- .mr_ncol(qd, r21); n22 <- .mr_ncol(qd, r22)
  expect_gt(n21, 0L); expect_gt(n22, 0L)
  expect_equal(.mr_ncol(qd, rBoth), n21 + n22)

  gBoth <- getGenotypes(qd, region = rBoth)
  g21   <- getGenotypes(qd, region = r21)
  expect_equal(unname(gBoth[, seq_len(n21)]), unname(g21))
  expect_true(all(grepl("_c22$", colnames(gBoth)[(n21 + 1):(n21 + n22)])))
})

test_that("multi-region: single-range extraction is unchanged (regression)", {
  skip_if_not_installed("snpStats")
  h  <- GenotypeHandle(plink1Prefix = file.path(test_data_dir, "test_variants"))
  qd <- QtlDataset(study = "S", genotypes = h,
                   phenotypes = list(ctx = .mr_makeSE(h@sampleIds)))
  bp <- h@snpInfo$BP
  r <- GenomicRanges::GRanges("chr21", IRanges::IRanges(min(bp), max(bp)))
  expect_equal(.mr_ncol(qd, r), nrow(h@snpInfo))
})

test_that(".qtlResolveVariantRegion rejects a non-GRanges / empty region", {
  skip_if_not_installed("snpStats")
  h  <- GenotypeHandle(plink1Prefix = file.path(test_data_dir, "test_variants"))
  qd <- QtlDataset(study = "S", genotypes = h,
                   phenotypes = list(ctx = .mr_makeSE(h@sampleIds)))
  expect_error(getGenotypes(qd, region = "chr21:1-2"), "must be a GRanges")
  expect_error(getGenotypes(qd, region = GenomicRanges::GRanges()),
               "at least one range")
})


