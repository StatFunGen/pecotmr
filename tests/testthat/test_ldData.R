context("LdData accessors")

# ===========================================================================
# Fixture helpers
# ===========================================================================

.ld_makeHandle <- function(snp_n = 4L, n_samples = 30L, path = "/tmp/h.gds") {
  new("GenotypeHandle",
    path = path,
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

.ld_makeVariants <- function(snp_n = 4L) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", snp_n),
    ranges = IRanges::IRanges(start = seq(100L, by = 100L,
                                           length.out = snp_n),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    A1 = rep("A", snp_n),
    A2 = rep("G", snp_n),
    variant_id = paste0("v", seq_len(snp_n)))
  gr
}

.ld_mockExtractor <- function(seed, n_samples = 30L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds,
                                    handle@snpInfo$SNP))
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

# ===========================================================================
# getCorrelation
# ===========================================================================

test_that("getCorrelation: returns the stored correlation matrix when set", {
  R <- diag(4)
  ld <- LdData(correlation = R, variants = .ld_makeVariants(),
               blockMetadata = S4Vectors::DataFrame(x = 1))
  expect_identical(getCorrelation(ld), R)
})

test_that("getCorrelation: computes from a single GenotypeHandle via extractBlockGenotypes", {
  ld <- LdData(correlation = NULL,
               genotypeHandle = .ld_makeHandle(),
               snpIdx         = 1:4,
               variants       = .ld_makeVariants(),
               blockMetadata  = S4Vectors::DataFrame(x = 1))
  local_mocked_bindings(extractBlockGenotypes = .ld_mockExtractor(seed = 7),
                        .package = "pecotmr")
  R <- getCorrelation(ld)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(4L, 4L))
  expect_equal(unname(diag(R)), c(1, 1, 1, 1), tolerance = 1e-12)
})

# NB: The "neither correlation nor genotypeHandle" branch in getCorrelation
# is defensive — LdData validity rejects that state at construction time, so
# the only way to hit the runtime stop() is to mutate the slot post-hoc.
# Skipping that test in favor of paths the validity actually permits.

test_that("getCorrelation: mixture handles produce a weighted-average R", {
  gh1 <- .ld_makeHandle(path = "/tmp/h1.gds")
  gh2 <- .ld_makeHandle(path = "/tmp/h2.gds")
  ld <- LdData(correlation = NULL,
               genotypeHandle = list(gh1, gh2),
               snpIdx         = 1:4,
               variants       = .ld_makeVariants(),
               blockMetadata  = S4Vectors::DataFrame(x = 1),
               mixtureWeights = c(0.3, 0.7))
  # Different seeds per panel so the underlying LD matrices differ.
  panelSeed <- list(11, 23)
  call_n <- 0
  local_mocked_bindings(
    extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
      call_n <<- call_n + 1L
      .ld_mockExtractor(seed = panelSeed[[call_n]])(handle, snpIdx,
                                                    meanImpute = meanImpute)
    },
    .package = "pecotmr")
  R_mix <- getCorrelation(ld)
  expect_equal(dim(R_mix), c(4L, 4L))
  # Recompute per-panel R independently and check weighted-average property.
  call_n <- 0
  local_mocked_bindings(
    extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
      call_n <<- call_n + 1L
      .ld_mockExtractor(seed = panelSeed[[call_n]])(handle, snpIdx,
                                                    meanImpute = meanImpute)
    },
    .package = "pecotmr")
  R_each <- lapply(list(gh1, gh2), function(h) {
    geno <- extractBlockGenotypes(h, 1:4)
    Xt <- t(SummarizedExperiment::assay(geno, "dosage"))
    computeLd(Xt, method = "sample")
  })
  expected <- 0.3 * R_each[[1L]] + 0.7 * R_each[[2L]]
  expect_equal(R_mix, expected, tolerance = 1e-12)
})

test_that("getCorrelation: mixture handles without mixtureWeights errors", {
  gh <- .ld_makeHandle()
  # Build via new() to skip the constructor's mixtureWeights validity check.
  ld <- new("LdData",
            correlation    = NULL,
            genotypeHandle = list(gh, gh),
            snpIdx         = 1:4,
            variants       = .ld_makeVariants(),
            blockMetadata  = S4Vectors::DataFrame(x = 1),
            nRef           = 0L,
            mixtureWeights = NULL)
  expect_error(getCorrelation(ld),
               "Cannot compute mixture LD: `mixtureWeights` is NULL")
})

test_that("getCorrelation: mixture panels of differing dim error", {
  gh_small <- .ld_makeHandle(snp_n = 3L)
  gh <- .ld_makeHandle(snp_n = 4L)
  ld <- new("LdData",
            correlation    = NULL,
            genotypeHandle = list(gh_small, gh),
            snpIdx         = 1:3,
            variants       = .ld_makeVariants(),
            blockMetadata  = S4Vectors::DataFrame(x = 1),
            nRef           = 0L,
            mixtureWeights = c(0.5, 0.5))
  # The first call sees snpIdx 1:3 against gh_small (3 variants); the second
  # against gh (4 variants in its panel) — but we still pass snpIdx 1:3, so
  # both return 3x3. To force a dim mismatch we mock the second call to
  # return a 4x4 panel.
  call_n <- 0
  local_mocked_bindings(
    extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
      call_n <<- call_n + 1L
      idx <- if (call_n == 1L) snpIdx else seq_len(4L)
      .ld_mockExtractor(seed = call_n)(handle, idx, meanImpute = meanImpute)
    },
    .package = "pecotmr")
  expect_error(getCorrelation(ld),
               "panels yielded LD matrices of differing dimensions")
})

# ===========================================================================
# getGenotypes
# ===========================================================================

test_that("getGenotypes: NULL handle returns NULL", {
  ld <- LdData(correlation = diag(4),
               variants = .ld_makeVariants(),
               blockMetadata = S4Vectors::DataFrame(x = 1))
  expect_null(getGenotypes(ld))
})

test_that("getGenotypes: matrix handle is returned unchanged", {
  X <- matrix(0, nrow = 10, ncol = 4,
              dimnames = list(paste0("s", 1:10), paste0("v", 1:4)))
  ld <- new("LdData",
            correlation    = NULL,
            genotypeHandle = X,
            snpIdx         = NULL,
            variants       = .ld_makeVariants(),
            blockMetadata  = S4Vectors::DataFrame(x = 1),
            nRef           = 0L,
            mixtureWeights = NULL)
  expect_identical(getGenotypes(ld), X)
})

test_that("getGenotypes: single handle returns samples x variants dosage", {
  ld <- LdData(correlation = NULL,
               genotypeHandle = .ld_makeHandle(),
               snpIdx         = 1:4,
               variants       = .ld_makeVariants(),
               blockMetadata  = S4Vectors::DataFrame(x = 1))
  local_mocked_bindings(extractBlockGenotypes = .ld_mockExtractor(seed = 7),
                        .package = "pecotmr")
  G <- getGenotypes(ld)
  expect_equal(dim(G), c(30L, 4L))
  expect_equal(colnames(G), paste0("v", 1:4))
})

test_that("getGenotypes: list of handles returns a list of dosage matrices", {
  gh1 <- .ld_makeHandle(path = "/tmp/h1.gds")
  gh2 <- .ld_makeHandle(path = "/tmp/h2.gds")
  ld <- LdData(correlation = NULL,
               genotypeHandle = list(gh1, gh2),
               snpIdx         = 1:4,
               variants       = .ld_makeVariants(),
               blockMetadata  = S4Vectors::DataFrame(x = 1),
               mixtureWeights = c(0.5, 0.5))
  call_n <- 0
  local_mocked_bindings(
    extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
      call_n <<- call_n + 1L
      .ld_mockExtractor(seed = 10 + call_n)(handle, snpIdx, meanImpute = meanImpute)
    },
    .package = "pecotmr")
  G <- getGenotypes(ld)
  expect_true(is.list(G))
  expect_equal(length(G), 2L)
  expect_equal(dim(G[[1L]]), c(30L, 4L))
})

# ===========================================================================
# Other accessors
# ===========================================================================

test_that("hasGenotypes: TRUE when handle present, FALSE otherwise", {
  ld_R <- LdData(correlation = diag(4), variants = .ld_makeVariants(),
                 blockMetadata = S4Vectors::DataFrame(x = 1))
  expect_false(hasGenotypes(ld_R))

  ld_gh <- LdData(correlation = NULL,
                  genotypeHandle = .ld_makeHandle(),
                  snpIdx = 1:4,
                  variants = .ld_makeVariants(),
                  blockMetadata = S4Vectors::DataFrame(x = 1))
  expect_true(hasGenotypes(ld_gh))
})

test_that("getVariantIds returns the variant_id mcol", {
  ld <- LdData(correlation = diag(4), variants = .ld_makeVariants(),
               blockMetadata = S4Vectors::DataFrame(x = 1))
  expect_equal(getVariantIds(ld), paste0("v", 1:4))
})

test_that("getVariantInfo / getBlockMetadata return slots verbatim", {
  vars <- .ld_makeVariants()
  bm <- S4Vectors::DataFrame(region = "chr1:100-400")
  ld <- LdData(correlation = diag(4), variants = vars, blockMetadata = bm)
  expect_identical(getVariantInfo(ld), vars)
  expect_identical(getBlockMetadata(ld), bm)
})

test_that("getRefPanel: assembles the chrom/pos/A1/A2/variant_id data.frame", {
  ld <- LdData(correlation = diag(4), variants = .ld_makeVariants(),
               blockMetadata = S4Vectors::DataFrame(x = 1))
  rp <- getRefPanel(ld)
  expect_s3_class(rp, "data.frame")
  expect_setequal(colnames(rp), c("A1", "A2", "variant_id", "chrom", "pos"))
  expect_equal(rp$variant_id, paste0("v", 1:4))
  expect_equal(rp$pos, c(100L, 200L, 300L, 400L))
})
