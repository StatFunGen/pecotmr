context("ctwasPipeline")

# ===========================================================================
# Strategy: ctwas::ctwas_sumstats does the heavy work. We mock it to a
# function that just returns its inputs back, so we can verify how the
# pipeline assembles z_snp / weights / region_info / LD loader inputs.
# ===========================================================================

.ctp_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
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

.ctp_mockExtractor <- function(seed = 5, n_samples = 30L) {
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

.ctp_makeGwasSumstats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 6L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:6),
    A1  = rep("A", 6), A2  = rep("G", 6),
    Z   = rnorm(6), N = rep(1000L, 6))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .ctp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.ctp_makeTwasWeights <- function() {
  e <- TwasWeightsEntry(
    variantIds = paste0("v", 1:5),
    weights    = c(0.1, 0.05, -0.2, 0.3, 0.0))
  TwasWeights(
    study    = "Q1", context = "c1", trait = "t1", method = "susie",
    entry    = list(e),
    ldSketch = .ctp_makeHandle())
}

# ===========================================================================
# Input validation
# ----------------------------------------------------------------------------
# The top of ctwasPipeline requires the (non-CRAN) `ctwas` package; without
# it the function errors out before any input-validation branch fires. Skip
# entry-point tests when ctwas isn't installed, but exercise the input-
# building helpers directly (they don't gate on ctwas).
# ===========================================================================

test_that("ctwasPipeline: rejects non-GwasSumStats gwasSumStats", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = "no",
                  twasWeights  = .ctp_makeTwasWeights()),
    "must be a GwasSumStats"
  )
})

test_that("ctwasPipeline: rejects un-QCd GwasSumStats", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(qc = FALSE),
                  twasWeights  = .ctp_makeTwasWeights()),
    "has no QC record"
  )
})

test_that("ctwasPipeline: rejects missing twasWeights", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats()),
    "must be a TwasWeights"
  )
})

test_that("ctwasPipeline: rejects non-GRanges twasZ", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights(),
                  twasZ        = "not a GRanges"),
    "must be a GRanges"
  )
})

test_that("ctwasPipeline: rejects bad regionId", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights(),
                  regionId     = ""),
    "non-empty character"
  )
})

test_that("ctwasPipeline: rejects unknown groupPriorVarStructure value", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights(),
                  groupPriorVarStructure = "bogus"),
    "'arg'"
  )
})

# ===========================================================================
# .ctwasRequireMatchingLdSketches
# ===========================================================================

test_that(".ctwasRequireMatchingLdSketches: NULL twas-side handle is allowed", {
  twNoLd <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(variantIds = paste0("v", 1:5),
                                   weights = rep(0.1, 5))),
    ldSketch = NULL)
  expect_silent(pecotmr:::.ctwasRequireMatchingLdSketches(
    twLd = NULL, gwasLd = .ctp_makeHandle()))
})

test_that(".ctwasRequireMatchingLdSketches: panel-size mismatch errors", {
  twLd  <- .ctp_makeHandle(snp_n = 5L)
  gwasLd <- .ctp_makeHandle(snp_n = 6L)
  expect_error(
    pecotmr:::.ctwasRequireMatchingLdSketches(twLd, gwasLd),
    "ldSketch panels differ in size"
  )
})

# ===========================================================================
# Input-building helpers
# ===========================================================================

test_that(".ctwasBuildZSnp: produces a flat data.frame keyed by SNP/study", {
  ss <- .ctp_makeGwasSumstats()
  df <- pecotmr:::.ctwasBuildZSnp(ss)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 6L)
  expect_setequal(colnames(df),
                  c("id", "chrom", "pos", "A1", "A2", "z", "study"))
  expect_setequal(df$id, paste0("v", 1:6))
  expect_setequal(unique(df$study), "G1")
})

test_that(".ctwasBuildSingleRegionInfo: pulls chrom + bp span from the ldSketch", {
  ri <- pecotmr:::.ctwasBuildSingleRegionInfo("block1", .ctp_makeHandle())
  expect_equal(ri$region_id, "block1")
  expect_equal(ri$chrom, 1L)
  expect_equal(ri$start, 100L)
  expect_equal(ri$stop, 600L)
})

test_that(".ctwasBuildSingleRegionInfo: multi-chromosome sketch errors", {
  h <- .ctp_makeHandle()
  h@snpInfo$CHR[1:3] <- "2"
  expect_error(
    pecotmr:::.ctwasBuildSingleRegionInfo("block1", h),
    "spans multiple chromosomes"
  )
})

test_that(".ctwasSnpInfoForBlock: returns ctwas-required columns chrom/id/pos/alt/ref", {
  df <- pecotmr:::.ctwasSnpInfoForBlock(.ctp_makeHandle())
  # ctwas's read_snp_info_files asserts these exact column names
  expect_setequal(colnames(df),
                  c("chrom", "id", "pos", "alt", "ref"))
})

# Build an ldPanel fixture (matches .ctwasComputeFullPanelLd's return
# shape) for the 6-SNP toy panel from .ctp_makeHandle().
.ctp_makeLdPanel <- function(snp_n = 6L) {
  h <- .ctp_makeHandle(snp_n = snp_n)
  snpInfo <- pecotmr:::.ctwasSnpInfoForBlock(h)
  R <- diag(1, snp_n)
  dimnames(R) <- list(snpInfo$id, snpInfo$id)
  list(R = R, snpInfo = snpInfo)
}

test_that(".ctwasBuildWeights: keys per-tuple weights and stamps gene metadata", {
  tw <- .ctp_makeTwasWeights()
  panel <- .ctp_makeLdPanel()
  wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
  expect_equal(length(wl), 1L)
  expect_equal(names(wl), "Q1|c1|t1|susie")
  expect_equal(wl[[1L]]$study, "Q1")
  expect_equal(wl[[1L]]$context, "c1")
  expect_equal(wl[[1L]]$gene_name, "t1")
  # wgt is a variants x 1 matrix with rownames = SNP IDs
  expect_true(is.matrix(wl[[1L]]$wgt))
  expect_equal(dim(wl[[1L]]$wgt), c(5L, 1L))
  expect_equal(rownames(wl[[1L]]$wgt), paste0("v", 1:5))
  # R_wgt is a 5x5 slice of the cached panel R
  expect_true(is.matrix(wl[[1L]]$R_wgt))
  expect_equal(dim(wl[[1L]]$R_wgt), c(5L, 5L))
  expect_equal(rownames(wl[[1L]]$R_wgt), paste0("v", 1:5))
  expect_equal(wl[[1L]]$n_wgt, 5L)
  # And it is literally a slice of the panel R (no recompute path).
  expect_equal(wl[[1L]]$R_wgt, panel$R[paste0("v", 1:5), paste0("v", 1:5)])
})

test_that(".ctwasBuildWeights: drops variants not present in the LD panel", {
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(
      variantIds = c(paste0("v", 1:3), "missing1", "missing2"),
      weights    = c(0.1, 0.2, 0.3, 0.4, 0.5))),
    ldSketch = .ctp_makeHandle())
  panel <- .ctp_makeLdPanel()
  wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
  expect_equal(nrow(wl[[1L]]$wgt), 3L)
  expect_equal(rownames(wl[[1L]]$wgt), paste0("v", 1:3))
  expect_equal(wl[[1L]]$n_wgt, 3L)
})

test_that(".ctwasComputeFullPanelLd: extracts once + returns cached R + snpInfo", {
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  out <- pecotmr:::.ctwasComputeFullPanelLd(.ctp_makeHandle())
  expect_named(out, c("R", "snpInfo"))
  expect_true(is.matrix(out$R))
  expect_equal(dim(out$R), c(6L, 6L))
  expect_equal(rownames(out$R), paste0("v", 1:6))
  expect_setequal(colnames(out$snpInfo),
                   c("chrom", "id", "pos", "alt", "ref"))
})

test_that(".ctwasBuildZGene: builds z_gene from a TWAS-Z GRanges", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    qtlStudy = "Q1", context = "c1", trait = "t1", method = "susie",
    twasZ = 1.5)
  df <- pecotmr:::.ctwasBuildZGene(gr)
  expect_equal(nrow(df), 1L)
  expect_setequal(colnames(df),
                  c("id", "z", "type", "context", "gene_name",
                    "study", "method"))
  expect_equal(df$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# LD loader / SNP-info loader closures
# ===========================================================================

test_that(".ctwasSingleBlockLdLoader: returns the cached R unchanged on every call", {
  R0 <- matrix(runif(36), 6, 6, dimnames = list(paste0("v", 1:6),
                                                  paste0("v", 1:6)))
  loader <- pecotmr:::.ctwasSingleBlockLdLoader(R0)
  expect_identical(loader("any_token"), R0)
  expect_identical(loader("different_token"), R0)
})

test_that(".ctwasSingleBlockSnpInfoLoader: returns the cached snpInfo unchanged on every call", {
  info0 <- data.frame(chrom = 1L, id = paste0("v", 1:3),
                       pos = c(100L, 200L, 300L),
                       alt = "A", ref = "G", stringsAsFactors = FALSE)
  loader <- pecotmr:::.ctwasSingleBlockSnpInfoLoader(info0)
  expect_identical(loader("any_token"), info0)
  expect_identical(loader("different_token"), info0)
})

# ===========================================================================
# End-to-end with mocked ctwas::ctwas_sumstats
# ===========================================================================

test_that("ctwasPipeline: assembles the documented input shape for ctwas_sumstats", {
  skip_if_not_installed("ctwas")
  ss <- .ctp_makeGwasSumstats()
  tw <- .ctp_makeTwasWeights()
  capturedArgs <- NULL
  local_mocked_bindings(
    ctwas_sumstats = function(...) {
      capturedArgs <<- list(...)
      list(susie_alpha_res = "mocked")
    },
    .package = "ctwas")
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  out <- ctwasPipeline(gwasSumStats = ss, twasWeights = tw,
                       regionId = "myBlock")
  expect_equal(out$susie_alpha_res, "mocked")
  expect_equal(capturedArgs$region_info$region_id, "myBlock")
  expect_equal(capturedArgs$z_snp$id, paste0("v", 1:6))
  expect_equal(length(capturedArgs$weights), 1L)
  expect_equal(capturedArgs$L, 5L)
})

test_that("ctwasPipeline: forwards a twasZ argument as z_gene", {
  skip_if_not_installed("ctwas")
  ss <- .ctp_makeGwasSumstats()
  tw <- .ctp_makeTwasWeights()
  twasZ <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  S4Vectors::mcols(twasZ) <- S4Vectors::DataFrame(
    qtlStudy = "Q1", context = "c1", trait = "t1", method = "susie",
    twasZ = 1.5)
  capturedArgs <- NULL
  local_mocked_bindings(
    ctwas_sumstats = function(...) {
      capturedArgs <<- list(...)
      list(ok = TRUE)
    },
    .package = "ctwas")
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  ctwasPipeline(gwasSumStats = ss, twasWeights = tw, twasZ = twasZ,
                regionId = "block1")
  expect_equal(capturedArgs$z_gene$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# Real-engine end-to-end: drives ctwas::ctwas_sumstats with the bundled
# example PLINK panel + synthetic TwasWeights. Exercises the LD-loader
# and snp-info-loader closure bodies as actually invoked by ctwas
# (mocked tests only construct the closures, never invoke them).
# ===========================================================================

test_that("ctwasPipeline: real-engine end-to-end on the bundled example panel", {
  skip_if_not_installed("ctwas")
  data(gwas_sumstats_s4_example)
  data(qtl_dataset_example)
  gss <- fixupExampleGenotypePaths(gwas_sumstats_s4_example)
  qd  <- fixupExampleGenotypePaths(qtl_dataset_example)
  gh  <- qd@genotypes

  # 5-variant synthetic gene from the bundled panel.
  vids <- gh@snpInfo$SNP[1:5]
  ent  <- TwasWeightsEntry(
    variantIds = vids,
    weights    = c(0.1, -0.2, 0.05, 0.0, 0.3))
  tw <- TwasWeights(
    study   = "study1", context = "brain",
    trait   = "ENSG_example", method = "susie",
    entry   = list(ent),
    ldSketch = gh)

  res <- suppressMessages(suppressWarnings(
    ctwasPipeline(gwasSumStats = gss, twasWeights = tw,
                  regionId    = "block1",
                  niter       = 5L,
                  niterPrefit = 2L,
                  # Single-gene single-block toy: relax the production
                  # filters that gate out tiny inputs.
                  min_group_size       = 1L,
                  min_p_single_effect  = 0,
                  filter_L             = FALSE)))

  # ctwas_sumstats returns these 7 elements on success.
  expect_named(res, c("z_gene", "param", "finemap_res", "susie_alpha_res",
                       "region_data", "boundary_genes", "screen_res"),
                ignore.order = TRUE)
  # The gene we passed in came through.
  expect_equal(nrow(res$z_gene), 1L)
  expect_equal(res$z_gene$id, "study1|brain|ENSG_example|susie")
  # Per-SNP susie_alpha output covers our 5 weight SNPs.
  expect_equal(nrow(res$susie_alpha_res), 5L)
  expect_true(all(c("susie_pip", "susie_alpha", "region_id")
                   %in% colnames(res$susie_alpha_res)))
})

# ===========================================================================
# .ctwasFilterVariants — ported from R/ctwasWrapper.R::trimCtwasVariants
# ===========================================================================
# The filter has four knobs:
#   1. twasWeightCutoff — drop |w| < cutoff
#   2. csMinCor         — high-purity CS rescue (must-keep)
#   3. minPipCutoff     — high-PIP rescue (must-keep)
#   4. maxNumVariants   — per-gene cap, prioritized by PIP then |w|

test_that(".ctwasFilterVariants: twasWeightCutoff drops low-magnitude variants", {
  vids <- paste0("v", 1:6)
  w    <- c(0.5, 0.001, 0.3, 0.0005, -0.4, 0)
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = NULL,
    twasWeightCutoff = 0.01, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = Inf)
  # Survivors: v1 (0.5), v3 (0.3), v5 (-0.4) — three with |w| >= 0.01
  expect_setequal(out$vids, c("v1", "v3", "v5"))
})

test_that(".ctwasFilterVariants: maxNumVariants caps by |w| when no PIP", {
  vids <- paste0("v", 1:5)
  w    <- c(0.1, 0.5, 0.2, 0.4, 0.05)
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = NULL,
    twasWeightCutoff = 0, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = 3)
  # Top 3 by |w|: v2 (0.5), v4 (0.4), v3 (0.2)
  expect_setequal(out$vids, c("v2", "v4", "v3"))
})

test_that(".ctwasFilterVariants: minPipCutoff rescues high-PIP variants from cap", {
  vids <- paste0("v", 1:5)
  w    <- c(0.5, 0.4, 0.3, 0.2, 0.1)
  finemapAux <- list(
    pip = setNames(c(0.01, 0.02, 0.8, 0.01, 0.95), vids),
    csMembers = list(),
    csPurity  = numeric(0))
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = finemapAux,
    twasWeightCutoff = 0, csMinCor = 0.8,
    minPipCutoff = 0.5, maxNumVariants = 2)
  # Must-keep (PIP > 0.5): v3, v5. Cap is 2 → both kept.
  expect_setequal(out$vids, c("v3", "v5"))
})

test_that(".ctwasFilterVariants: csMinCor rescues high-purity CS members from cap", {
  vids <- paste0("v", 1:6)
  w    <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
  finemapAux <- list(
    pip = setNames(rep(0, length(vids)), vids),
    csMembers = list(c("v3", "v6"), c("v2", "v4")),
    csPurity  = c(0.9, 0.5))  # CS 1 (v3, v6) is high-purity
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = finemapAux,
    twasWeightCutoff = 0, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = 3)
  # Must-keep from high-purity CS: v3, v6. Remaining slot filled by
  # next-highest |w| that isn't must-keep: v1 (0.5).
  expect_setequal(out$vids, c("v3", "v6", "v1"))
})

test_that(".ctwasFilterVariants: returns NULL when no variants survive", {
  vids <- paste0("v", 1:3)
  w    <- c(0.001, 0.0005, 0.002)
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = NULL,
    twasWeightCutoff = 0.5, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = Inf)
  expect_null(out)
})

test_that(".ctwasBuildWeights: maxNumVariants caps the per-gene weight matrix", {
  data(qtl_dataset_example)
  qd <- fixupExampleGenotypePaths(qtl_dataset_example)
  gh <- qd@genotypes
  vids <- gh@snpInfo$SNP[1:5]
  ent <- TwasWeightsEntry(
    variantIds = vids,
    weights    = c(0.1, -0.2, 0.05, 0.3, 0.15))
  tw <- TwasWeights(
    study = "study1", context = "brain",
    trait = "ENSG_example", method = "susie",
    entry = list(ent), ldSketch = gh)
  ldPanel <- pecotmr:::.ctwasComputeFullPanelLd(gh)
  wl <- pecotmr:::.ctwasBuildWeights(tw, ldPanel, maxNumVariants = 3L)
  expect_equal(wl[[1L]]$n_wgt, 3L)
  expect_equal(nrow(wl[[1L]]$wgt), 3L)
  # Top 3 by |w| from c(0.1, -0.2, 0.05, 0.3, 0.15): 0.3, -0.2, 0.15
  expect_setequal(rownames(wl[[1L]]$wgt), vids[c(4L, 2L, 5L)])
})

test_that(".ctwasBuildWeights: twasWeightCutoff drops low-magnitude variants", {
  data(qtl_dataset_example)
  qd <- fixupExampleGenotypePaths(qtl_dataset_example)
  gh <- qd@genotypes
  vids <- gh@snpInfo$SNP[1:5]
  ent <- TwasWeightsEntry(
    variantIds = vids,
    # v1 (0.005) and v3 (0.001) will be dropped at cutoff 0.01
    weights    = c(0.005, 0.2, 0.001, 0.3, 0.1))
  tw <- TwasWeights(
    study = "study1", context = "brain",
    trait = "ENSG_example", method = "susie",
    entry = list(ent), ldSketch = gh)
  ldPanel <- pecotmr:::.ctwasComputeFullPanelLd(gh)
  wl <- pecotmr:::.ctwasBuildWeights(tw, ldPanel, twasWeightCutoff = 0.01)
  expect_equal(wl[[1L]]$n_wgt, 3L)
  expect_setequal(rownames(wl[[1L]]$wgt), vids[c(2L, 4L, 5L)])
})
