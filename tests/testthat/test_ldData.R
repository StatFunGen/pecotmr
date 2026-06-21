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

# ===========================================================================
# Tests migrated from test_dataStructures.R (LdData class)
# ===========================================================================

test_that("LdData constructor works with correlation matrix", {
  R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(R) <- colnames(R) <- c("chr1:100:A:G", "chr1:200:C:T")
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
    A1 = c("G", "T"), A2 = c("A", "C")
  )
  bm <- data.frame(blockId = 1L, chrom = "1", blockStart = 100L,
                    blockEnd = 200L, size = 2L, startIdx = 1L, endIdx = 2L)

  ld <- LdData(correlation = R, variants = gr, blockMetadata = bm)
  expect_s4_class(ld, "LdData")
  expect_false(hasGenotypes(ld))
  expect_true(is.matrix(getCorrelation(ld)))
  expect_equal(getVariantIds(ld), c("chr1:100:A:G", "chr1:200:C:T"))
  expect_equal(nrow(getBlockMetadata(ld)), 1L)
  expect_null(getGenotypes(ld))
})


test_that("LdData validation rejects empty variants", {
  R <- diag(2)
  gr <- GenomicRanges::GRanges()
  expect_error(
    LdData(correlation = R, variants = gr,
           blockMetadata = data.frame()),
    "must not be empty"
  )
})


test_that("LdData validation rejects NULL correlation and handle", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100L, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = "chr1:100:A:G", A1 = "G", A2 = "A"
  )
  expect_error(
    LdData(correlation = NULL, genotypeHandle = NULL,
           variants = gr, blockMetadata = data.frame()),
    "At least one"
  )
})


test_that("LdData show method works", {
  R <- diag(3)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("v", 1:3), A1 = rep("A", 3), A2 = rep("G", 3)
  )
  ld <- LdData(correlation = R, variants = gr,
               blockMetadata = data.frame())
  expect_output(show(ld), "LdData: 3 variants")
})


test_that("LdData supports block-diagonal correlation", {
  R1 <- diag(2)
  R2 <- diag(3)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(start = seq(100L, 500L, 100L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("v", 1:5), A1 = rep("A", 5), A2 = rep("G", 5)
  )
  ld <- LdData(correlation = list(R1, R2), variants = gr,
               blockMetadata = data.frame())
  corr <- getCorrelation(ld)
  expect_true(is.list(corr))
  expect_equal(length(corr), 2)
})


test_that("LdData S4 accessors return correct data", {
  R <- diag(2)
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = c("v1", "v2"), A1 = c("A", "C"), A2 = c("G", "T")
  )
  ld <- LdData(correlation = R, variants = gr,
               blockMetadata = data.frame(blockId = 1L))
  expect_equal(getCorrelation(ld), R)
  expect_equal(getVariantIds(ld), c("v1", "v2"))
  expect_false(hasGenotypes(ld))
  rp <- getRefPanel(ld)
  expect_true(is.data.frame(rp))
  expect_true("variant_id" %in% names(rp))
  expect_equal(rp$variant_id, c("v1", "v2"))
})


test_that(".refPanelToGranges builds GRanges from data.frame", {
  rp <- data.frame(
    chrom = c("1", "1"), pos = c(100L, 200L),
    A1 = c("G", "T"), A2 = c("A", "C"),
    variant_id = c("v1", "v2"),
    allele_freq = c(0.3, 0.7),
    stringsAsFactors = FALSE
  )
  gr <- pecotmr:::.refPanelToGranges(rp)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(S4Vectors::mcols(gr)$variant_id, c("v1", "v2"))
  expect_equal(S4Vectors::mcols(gr)$allele_freq, c(0.3, 0.7))
})
# =============================================================================
# getTopLoci(type = "GRanges")
# =============================================================================



# === Tests migrated from test_s4Constructors.R (LdData) ===

test_that("LdData: pre-computed correlation matrix is returned by getCorrelation", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    A1 = rep("A", 3), A2 = rep("G", 3))
  R <- diag(3)
  block_meta <- S4Vectors::DataFrame(region = "chr1:100-300")
  ld <- LdData(correlation = R, genotypeHandle = NULL, snpIdx = NULL,
               variants = gr, blockMetadata = block_meta, nRef = 100L)
  expect_s4_class(ld, "LdData")
  expect_equal(getCorrelation(ld), R)
})


test_that("LdData: validity rejects both correlation AND genotypeHandle being NULL", {
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 100L, width = 1L))
  expect_error(
    LdData(correlation = NULL, genotypeHandle = NULL,
           variants = gr,
           blockMetadata = S4Vectors::DataFrame(x = 1)),
    "At least one of 'correlation' or 'genotypeHandle' must be non-NULL"
  )
})


test_that("LdData: validity rejects empty variants", {
  expect_error(
    LdData(correlation = diag(0), genotypeHandle = NULL,
           variants = GenomicRanges::GRanges(),
           blockMetadata = S4Vectors::DataFrame(x = 1)),
    "'variants' must not be empty"
  )
})


test_that("LdData: mixtureWeights only valid when genotypeHandle is a list", {
  gh <- .sc_makeGenotypeHandle()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 100L, width = 1L))
  expect_error(
    LdData(correlation = NULL, genotypeHandle = gh,
           variants = gr,
           blockMetadata = S4Vectors::DataFrame(x = 1),
           mixtureWeights = c(0.5, 0.5)),
    "'mixtureWeights' may only be set when 'genotypeHandle' is a list"
  )
})


test_that("LdData: mixtureWeights must be non-negative and sum to 1", {
  gh1 <- .sc_makeGenotypeHandle()
  gh2 <- .sc_makeGenotypeHandle()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 100L, width = 1L))
  expect_error(
    LdData(correlation = NULL,
           genotypeHandle = list(gh1, gh2),
           variants = gr,
           blockMetadata = S4Vectors::DataFrame(x = 1),
           mixtureWeights = c(0.4, 0.4)),
    "must be non-negative and sum to 1"
  )
})

# ===========================================================================
# QtlDataset
# ===========================================================================

.sc_makeSe <- function(traits = c("ENSG1", "ENSG2"), n_samples = 10) {
  rng <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(traits)),
    ranges = IRanges::IRanges(
      start = seq(1000L, by = 1000L, length.out = length(traits)),
      width = 500L))
  names(rng) <- traits
  expr <- matrix(rnorm(length(traits) * n_samples),
                 nrow = length(traits), ncol = n_samples,
                 dimnames = list(traits, paste0("s", seq_len(n_samples))))
  cd <- S4Vectors::DataFrame(sex = rep(c("M", "F"), length.out = n_samples),
                             row.names = paste0("s", seq_len(n_samples)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expr),
    rowRanges = rng,
    colData = cd)
}


