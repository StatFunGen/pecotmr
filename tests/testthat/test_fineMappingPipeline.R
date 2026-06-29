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
# only by the token so post-process knows what to wrap. The `userArgs`
# parameter (per-method kwargs merged in by .fmMergeUserArgs) is accepted
# but ignored — the mocks don't simulate downstream susie behaviour.
.fmp_mockFitIndiv <- function() {
  function(X, y, token, chainFromInf = NULL, coverage = 0.95,
           userArgs = NULL) {
    list(token = token, X_cols = ncol(X))
  }
}

.fmp_mockFitRss <- function() {
  function(z, R, n, token, chainFromInf = NULL, coverage = 0.95,
           userArgs = NULL) {
    list(token = token, n_variants = length(z))
  }
}

.fmp_mockPostprocess <- function() {
  function(fit, method, dataX, dataY, coverage, secondaryCoverage,
           signalCutoff, minAbsCorr, csInput = NULL, af = NULL,
           region = NULL, conditionIdx = NULL) {
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
      susieFit = list(method = method, payload = fit),
      topLoci    = data.frame(variant_id = vids,
                              pip = seq(0.9, by = -0.1,
                                        length.out = length(vids)),
                              stringsAsFactors = FALSE))
  }
}

# ===========================================================================
# .fmNormalizeMethods
# ===========================================================================

test_that(".fmNormalizeMethods: rejects NULL / empty / non-character/list", {
  expect_error(pecotmr:::.fmNormalizeMethods(NULL),
               "non-empty character")
  expect_error(pecotmr:::.fmNormalizeMethods(character(0)),
               "non-empty character")
  expect_error(pecotmr:::.fmNormalizeMethods(42L),
               "character vector or")
})

test_that(".fmNormalizeMethods: char-vector form deduplicates + seeds susie L defaults", {
  res <- pecotmr:::.fmNormalizeMethods(c("susie", "susie", "susieInf"))
  expect_equal(res$tokens, c("susie", "susieInf"))
  expect_equal(names(res$methodArgs), c("susie", "susieInf"))
  # SuSiE-family tokens get the pipeline L / L_greedy defaults (pecotmr owns
  # these, not the CLI wrappers).
  expect_equal(res$methodArgs$susie$L, 20L)
  expect_equal(res$methodArgs$susie$L_greedy, 5L)
  expect_equal(res$methodArgs$susieInf$L, 20L)
  # Non-susie-family tokens are left untouched.
  expect_length(pecotmr:::.fmNormalizeMethods("mvsusie")$methodArgs$mvsusie, 0L)
})

test_that(".fmNormalizeMethods: named-list keeps kwargs + fills missing susie L", {
  res <- pecotmr:::.fmNormalizeMethods(
    list(susie    = list(L = 1, refine = FALSE),
         susieInf = list()))
  expect_equal(res$tokens, c("susie", "susieInf"))
  expect_equal(res$methodArgs$susie$L, 1)           # explicit kwarg wins
  expect_false(res$methodArgs$susie$refine)
  expect_equal(res$methodArgs$susie$L_greedy, 5L)   # filled-in default
  expect_equal(res$methodArgs$susieInf$L, 20L)      # both filled
})

test_that(".fmNormalizeMethods: L / Lgreedy args override the susie defaults", {
  res <- pecotmr:::.fmNormalizeMethods(c("susie"), L = 30L, Lgreedy = 7L)
  expect_equal(res$methodArgs$susie$L, 30L)
  expect_equal(res$methodArgs$susie$L_greedy, 7L)
})

test_that(".fmNormalizeMethods: list without names errors", {
  expect_error(pecotmr:::.fmNormalizeMethods(list(list(L = 1), list())),
               "must be named")
})

test_that(".fmNormalizeMethods: list with non-list child errors", {
  expect_error(
    pecotmr:::.fmNormalizeMethods(list(susie = 42, susieInf = list())),
    "list of named kwargs")
})

test_that(".fmMergeUserArgs: user overrides win over capability defaults + base", {
  # base sets convergence_method = "pip"; user overrides to "objective"
  out <- pecotmr:::.fmMergeUserArgs(
    list(z = 1:3, convergence_method = "pip"),
    token    = "susie",
    userArgs = list(convergence_method = "objective", L = 1))
  expect_equal(out$convergence_method, "objective")
  expect_equal(out$L, 1)
  expect_identical(out$z, 1:3)
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

test_that(".fmTraitsInRegion: keeps genes overlapping the region; NULL keeps all", {
  # g1: 1000-1500, g2: 2000-2500, g3: 3000-3500 (width 500).
  se <- .fmp_makeSe(traits = c("g1", "g2", "g3"),
                    starts = c(1000L, 2000L, 3000L))
  r1  <- GenomicRanges::GRanges("chr1", IRanges::IRanges(900L, 1600L))  # g1 only
  r12 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(900L, 2600L))  # g1 + g2
  expect_equal(pecotmr:::.fmTraitsInRegion(se, c("g1", "g2", "g3"), r1), "g1")
  expect_setequal(pecotmr:::.fmTraitsInRegion(se, c("g1", "g2", "g3"), r12),
                  c("g1", "g2"))
  # NULL region (gene/cisWindow mode) leaves the set unchanged.
  expect_setequal(pecotmr:::.fmTraitsInRegion(se, c("g1", "g2", "g3"), NULL),
                  c("g1", "g2", "g3"))
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
    susieFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  fmr <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  hit <- pecotmr:::.fmCacheLookup(fmr, "s1", "c1", "t1", "susie")
  expect_identical(hit, e)
  expect_null(pecotmr:::.fmCacheLookup(fmr, "ghost", "c1", "t1", "susie"))
})

test_that(".fmCacheLookupGwas: returns matching entry by (study, method, region_id)", {
  e <- FineMappingEntry(
    variantIds = "v1",
    susieFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  # GwasFineMappingResult assigns the synthetic region_id "region_1"
  # when none is supplied; the lookup must include it in the 3-tuple
  # key (multi-block FMRs disambiguate per-block fits by region_id).
  fmr <- GwasFineMappingResult(study = "g1", method = "susie",
                                entry = list(e))
  expect_identical(
    pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", "region_1"),
    e)
  expect_null(pecotmr:::.fmCacheLookupGwas(fmr, "ghost", "susie", "region_1"))
  # Wrong region_id is also a miss.
  expect_null(pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", "other_region"))
})

test_that(".fmCacheLookup: non-QtlFineMappingResult input returns NULL", {
  e <- FineMappingEntry(
    variantIds = "v1",
    susieFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  gwasFmr <- GwasFineMappingResult(study = "g1", method = "susie",
                                    entry = list(e))
  expect_null(pecotmr:::.fmCacheLookup(gwasFmr, "g1", "c1", "t1", "susie"))
})

test_that(".fmCacheLookupGwas: non-GwasFineMappingResult input returns NULL", {
  e <- FineMappingEntry(
    variantIds = "v1",
    susieFit = list(token = "susie"),
    topLoci    = data.frame(variant_id = "v1", pip = 0.5,
                            stringsAsFactors = FALSE))
  qtlFmr <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_null(pecotmr:::.fmCacheLookupGwas(qtlFmr, "s1", "susie", "region_1"))
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
    "no \\(study, method, region_id\\) tuples"
  )
})

# ===========================================================================
# .rbindFineMappingResult — class-check branches
# ===========================================================================

test_that(".rbindFineMappingResult: rejects non-FineMappingResultBase input", {
  e <- FineMappingEntry(
    variantIds = "v1",
    susieFit = list(token = "susie"),
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
    susieFit = list(token = "susie"),
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
    susieFit = list(token = "susie"),
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

test_that(".fmAfForX: returns directional effect-allele af aligned to colnames(X)", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .package = "pecotmr")
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 100000L))
  # The helper aligns af to dimnames; values come from getAf over the same
  # selection. Columns deliberately reordered to test name-based alignment.
  X <- matrix(0, nrow = 5L, ncol = 3L,
              dimnames = list(paste0("s", 1:5), c("v3", "v1", "v2")))
  af <- pecotmr:::.fmAfForX(qd, X, region = region)
  expect_length(af, 3L)
  expect_false(anyNA(af))  # region matched -> every fitted variant has an af
  expect_equal(
    af,
    unname(getAf(qd, region = region, samples = rownames(X))[colnames(X)]))
})

test_that(".fmAfForX: returns NULL for an empty block or a non-QtlDataset source", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  emptyX <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  expect_null(pecotmr:::.fmAfForX(qd, emptyX))
  X <- matrix(0, nrow = 2L, ncol = 1L,
              dimnames = list(c("s1", "s2"), "v1"))
  expect_null(pecotmr:::.fmAfForX(list(not = "a dataset"), X))
})

test_that("fineMappingPipeline(QtlDataset): threads directional af into postprocess", {
  # Regression for af = NA in getCs: the individual-level univariate path must
  # forward a non-NULL, directional effect-allele frequency to the
  # post-processor (which writes it into the topLoci `af` column).
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  captured <- new.env(parent = emptyenv())
  captured$af   <- "UNSET"
  captured$cols <- NULL
  recordingPostprocess <- function(fit, method, dataX, dataY, coverage,
                                   secondaryCoverage, signalCutoff, minAbsCorr,
                                   csInput = NULL, af = NULL, region = NULL,
                                   conditionIdx = NULL) {
    captured$af   <- af
    captured$cols <- colnames(dataX)
    vids <- colnames(dataX)
    FineMappingEntry(
      variantIds = vids,
      susieFit   = list(method = method),
      topLoci    = data.frame(variant_id = vids,
                              pip = seq(0.9, by = -0.1,
                                        length.out = length(vids)),
                              af  = af,
                              stringsAsFactors = FALSE))
  }
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = recordingPostprocess,
    .package = "pecotmr")
  suppressMessages(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        addSusieInf = FALSE))
  # af forwarded (not left at the NULL default), one value per fitted variant.
  expect_false(identical(captured$af, "UNSET"))
  expect_false(is.null(captured$af))
  expect_length(captured$af, length(captured$cols))
  expect_true(any(!is.na(captured$af)))
  expect_true(all(captured$af >= 0 & captured$af <= 1, na.rm = TRUE))
})

test_that("fineMappingPipeline(QtlDataset): seed argument is accepted and runs", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        addSusieInf = FALSE, seed = 42L))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 1L)
})

test_that(".fmSerScreen: disables on 0, skips no-signal, keeps signal + adaptive", {
  skip_if_not_installed("susieR")
  set.seed(1)
  n <- 150L; p <- 25L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("v", seq_len(p))
  yNull <- rnorm(n)                       # no association
  ySig  <- X[, 1] * 2 + rnorm(n, sd = 0.3)  # strong single effect at v1
  fn <- function(...) suppressMessages(pecotmr:::.fmSerScreen(...))
  expect_true(fn(X, yNull, 0))            # cutoff 0 disables -> always keep
  expect_false(fn(X, yNull, 0.5))         # no PIP that high -> skip
  expect_true(fn(X, ySig, 0.5))           # strong signal clears 0.5 -> keep
  expect_true(fn(X, ySig, -1))            # adaptive 3/p: signal keeps
  expect_false(fn(X, yNull, -1))          # adaptive 3/p: null skips
  expect_true(fn(X, yNull, NA))           # malformed cutoff -> advisory keep
})

test_that(".fmTopPcScores: clean matrix -> samples x min(nPCs, traits) topPC scores", {
  set.seed(7)
  n <- 30L
  Y <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L,
              dimnames = list(paste0("s", seq_len(n)), c("ta", "tb", "tc")))
  fn <- function(...) pecotmr:::.fmTopPcScores(...)
  # (a) 3 traits, nPCs >= traits -> 3 columns named topPC1..topPC3, rows = samples.
  sc <- fn(Y, 10L)
  expect_true(is.matrix(sc))
  expect_equal(ncol(sc), 3L)
  expect_equal(colnames(sc), c("topPC1", "topPC2", "topPC3"))
  expect_equal(nrow(sc), n)
  expect_equal(rownames(sc), rownames(Y))
  # (b) nPCs caps the number of returned columns.
  sc2 <- fn(Y, 2L)
  expect_equal(ncol(sc2), 2L)
  expect_equal(colnames(sc2), c("topPC1", "topPC2"))
  # (c) single-column Y -> NULL (PCA undefined for a single trait).
  expect_null(fn(Y[, 1L, drop = FALSE], 10L))
  # (d) a zero-variance trait is dropped; k reflects only the usable traits.
  Yzv <- cbind(Y, td = rep(1, n))
  scz <- fn(Yzv, 10L)
  expect_equal(ncol(scz), 3L)            # td dropped -> still 3 usable traits
  expect_equal(colnames(scz), c("topPC1", "topPC2", "topPC3"))
  # (e) rows with any NA are dropped before PCA.
  Yna <- Y
  Yna[c(1L, 2L), 1L] <- NA
  scn <- fn(Yna, 10L)
  expect_equal(nrow(scn), n - 2L)
  expect_false(any(c("s1", "s2") %in% rownames(scn)))
})

test_that("fineMappingPipeline(QtlDataset, usePCA): top-PC susie rows keyed topPC{i}", {
  qd <- .fmp_makeQtlDataset(contexts = "brain",
                            traits = c("ENSG_A", "ENSG_B", "ENSG_C"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie",
                        cisWindow = 1000L, addSusieInf = FALSE,
                        usePCA = TRUE, nPCs = 2L))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_setequal(getMethodNames(res), "susie")
  pcRows <- as.character(res$trait) %in% c("topPC1", "topPC2")
  # 3 per-trait univariate susie rows + 2 top-PC rows = 5.
  expect_equal(sum(pcRows), 2L)
  expect_setequal(as.character(res$trait)[pcRows], c("topPC1", "topPC2"))
  expect_setequal(as.character(res$method)[pcRows], "susie")
})

test_that(".buildMvsusieReweightedPrior: canonical fallback when no usable fit", {
  bp <- function(...) pecotmr:::.buildMvsusieReweightedPrior(...)
  # No fit at all -> canonical prior, residualVariance NULL.
  p1 <- bp(NULL, c("c1", "c2"))
  expect_false(is.null(p1$priorVariance))
  expect_null(p1$residualVariance)
  # Fit with no data-driven matrices -> canonical prior, but V carried through.
  p2 <- bp(list(dataDrivenPriorMatrices = NULL, V = diag(2)), c("c1", "c2"))
  expect_equal(p2$residualVariance, diag(2))
})

test_that(".buildMvsusieReweightedPrior: reweights matrices by rescaleCovW0(w0)", {
  ddpm <- list(U = list(compA = diag(2), compB = diag(2) * 2),
               w = c(compA = 0.5, compB = 0.5))
  fit  <- list(dataDrivenPriorMatrices = ddpm,
               w0 = c(compA_grid1 = 0.3, compB_grid1 = 0.7),
               V  = diag(2) * 3)
  captured <- NULL
  # rescaleCovW0 collapses expanded w0 onto the original matrix names; mock it
  # so the test asserts the wiring, not rescaleCovW0's internals.
  local_mocked_bindings(
    rescaleCovW0 = function(w0) c(compA = 0.4, compB = 0.6),
    .package = "pecotmr")
  local_mocked_bindings(
    create_mixture_prior = function(...) { captured <<- list(...); "PRIOR" },
    .package = "mvsusieR")
  res <- pecotmr:::.buildMvsusieReweightedPrior(fit, c("c1", "c2"),
                                                weightsTol = 1e-8)
  expect_identical(res$priorVariance, "PRIOR")
  expect_equal(res$residualVariance, diag(2) * 3)
  expect_equal(captured$mixture_prior$weights, c(compA = 0.4, compB = 0.6))
  expect_equal(names(captured$mixture_prior$matrices), c("compA", "compB"))
  expect_equal(captured$include_indices, c("c1", "c2"))
  expect_equal(captured$weights_tol, 1e-8)
})

test_that(".fmLookupMrmashFit: finds the mr.mash fit by (study, trait)", {
  mkEntry <- function(fits) TwasWeightsEntry(
    variantIds = c("v1", "v2"), weights = c(0.1, 0.2), fits = fits)
  payload <- list(dataDrivenPriorMatrices = list(U = list(a = diag(2))),
                  w0 = c(a = 1), V = diag(2))
  # The joint fit lives on the first mrmash row of the (study, trait) group;
  # the other context row carries fits = NULL. A non-mrmash row is ignored.
  tw <- TwasWeights(
    study   = c("S", "S", "S"),
    context = c("c1", "c2", "c1"),
    trait   = c("G", "G", "G"),
    method  = c("mrmash", "mrmash", "enet"),
    entry   = list(mkEntry(payload), mkEntry(NULL), mkEntry(payload)))
  lk <- function(...) pecotmr:::.fmLookupMrmashFit(...)
  expect_identical(lk(tw, "S", "G"), payload)   # first non-NULL mrmash row
  expect_null(lk(tw, "S", "OTHER"))             # no such trait
  expect_null(lk(tw, "OTHER", "G"))             # no such study
  expect_null(lk(NULL, "S", "G"))               # no TwasWeights supplied
})

test_that(".fmLookupMrmashCv: finds the per-fold CV payload by (study, trait)", {
  mkEntry <- function(cv) TwasWeightsEntry(
    variantIds = c("v1", "v2"), weights = c(0.1, 0.2), cvResult = cv)
  cv <- list(samplePartition = data.frame(Sample = "s1", Fold = 1L),
             foldFits = list(fold_1 = list(w0 = 1)))
  tw <- TwasWeights(
    study = c("S", "S"), context = c("c1", "c2"), trait = c("G", "G"),
    method = c("mrmash", "mrmash"), entry = list(mkEntry(cv), mkEntry(NULL)))
  lk <- function(...) pecotmr:::.fmLookupMrmashCv(...)
  expect_identical(lk(tw, "S", "G"), cv)
  expect_null(lk(tw, "S", "OTHER"))
  expect_null(lk(NULL, "S", "G"))
  # A cvResult without foldFits is not a per-fold prior payload -> NULL.
  tw2 <- TwasWeights(study = "S", context = "c1", trait = "G", method = "mrmash",
                     entry = list(mkEntry(list(predictions = 1))))
  expect_null(lk(tw2, "S", "G"))
})

test_that(".buildMvsusieReweightedPrior: overrideU swaps matrices, keeps fit w0/V", {
  fit <- list(dataDrivenPriorMatrices = list(U = list(K = diag(2)), w = c(K = 1)),
              w0 = c(K_grid1 = 1), V = diag(2) * 7)
  override <- list(U = list(K = diag(2) * 5))
  captured <- NULL
  local_mocked_bindings(rescaleCovW0 = function(w0) c(K = 1), .package = "pecotmr")
  local_mocked_bindings(
    create_mixture_prior = function(...) { captured <<- list(...); "PRIOR" },
    .package = "mvsusieR")
  res <- pecotmr:::.buildMvsusieReweightedPrior(fit, c("c1", "c2"),
                                                overrideU = override)
  expect_equal(captured$mixture_prior$matrices$K, diag(2) * 5)  # the override U
  expect_equal(res$residualVariance, diag(2) * 7)               # the fit's own V
})

test_that(".fmBuildMvsusiePriorCv: mode B reweights each fold's own fit", {
  sp <- data.frame(Sample = paste0("s", 1:6), Fold = rep(1:3, each = 2),
                   stringsAsFactors = FALSE)
  mkFold <- function(uname, v) list(
    dataDrivenPriorMatrices = list(U = setNames(list(diag(2)), uname),
                                   w = setNames(0.5, uname)),
    w0 = setNames(0.5, paste0(uname, "_grid1")),
    V  = diag(2) * v)
  mvCv <- list(samplePartition = sp, foldFits = list(
    fold_1 = mkFold("A", 10), fold_2 = mkFold("B", 20), fold_3 = mkFold("C", 30)))
  local_mocked_bindings(
    rescaleCovW0 = function(w0) setNames(1, sub("_grid1$", "", names(w0))),
    .package = "pecotmr")
  local_mocked_bindings(
    create_mixture_prior = function(...) "PRIOR", .package = "mvsusieR")
  out <- pecotmr:::.fmBuildMvsusiePriorCv(mvCv, fullFitParts = NULL,
                                          conditionNames = c("c1", "c2"))
  expect_equal(names(out), c("1", "2", "3"))
  expect_equal(out[["1"]]$residualVariance, diag(2) * 10)  # fold 1's own V
  expect_equal(out[["3"]]$residualVariance, diag(2) * 30)  # fold 3's own V
})

test_that(".fmBuildMvsusiePriorCv: mode C reuses full-fit w0/V with per-fold U", {
  sp <- data.frame(Sample = paste0("s", 1:4), Fold = rep(1:2, each = 2),
                   stringsAsFactors = FALSE)
  full <- list(dataDrivenPriorMatrices = list(U = list(Z = diag(2)), w = c(Z = 1)),
               w0 = c(Z_grid1 = 1), V = diag(2) * 99)
  # Fold stubs carry only U (no w0) -> mode C: override the full fit's U.
  mvCv <- list(samplePartition = sp, foldFits = list(
    fold_1 = list(dataDrivenPriorMatrices = list(U = list(Z = diag(2) * 2))),
    fold_2 = list(dataDrivenPriorMatrices = list(U = list(Z = diag(2) * 3)))))
  captured <- list()
  local_mocked_bindings(
    rescaleCovW0 = function(w0) setNames(1, sub("_grid1$", "", names(w0))),
    .package = "pecotmr")
  local_mocked_bindings(
    create_mixture_prior = function(...) {
      captured[[length(captured) + 1L]] <<- list(...); "PRIOR" },
    .package = "mvsusieR")
  out <- pecotmr:::.fmBuildMvsusiePriorCv(mvCv, fullFitParts = full,
                                          conditionNames = c("c1", "c2"))
  expect_equal(names(out), c("1", "2"))
  expect_equal(out[["1"]]$residualVariance, diag(2) * 99)        # full fit's V
  expect_equal(captured[[1]]$mixture_prior$matrices$Z, diag(2) * 2)  # fold 1's U
  expect_equal(captured[[2]]$mixture_prior$matrices$Z, diag(2) * 3)  # fold 2's U
})

test_that("fineMappingPipeline(QtlDataset): pipCutoffToSkip skips no-signal univariate traits", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  # Stateful screen: reject the first block (ENSG_A), keep the rest (ENSG_B).
  seen <- 0L
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .fmSerScreen          = function(X, y, cutoff) { seen <<- seen + 1L; seen > 1L },
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        addSusieInf = FALSE, pipCutoffToSkip = -1))
  # ENSG_A screened out, ENSG_B kept -> a single row.
  expect_equal(nrow(res), 1L)
  expect_setequal(getTraits(res), "ENSG_B")
})

test_that("fineMappingPipeline(QtlDataset, cvFolds>1): attaches cvResult end to end", {
  skip_if_not_installed("susieR")
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  # Real susie fits drive CV; the genotype extraction and the full-data
  # post-processor are mocked (the mock SNP ids are not chr:pos:a1:a2, which the
  # real buildTopLoci requires). CV runs its own real .fmFitSusieIndiv, so the
  # cvResult is genuine.
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        addSusieInf = FALSE, cvFolds = 3, verbose = 0))
  cv <- getCvResult(res, study = "study1", context = "brain",
                    trait = "ENSG_A", method = "susie")
  expect_false(is.null(cv))
  expect_setequal(colnames(cv$samplePartition), c("Sample", "Fold"))
  expect_setequal(sort(unique(cv$samplePartition$Fold)), 1:3)
  expect_true("susie_predicted" %in% names(cv$prediction))
  expect_true("susie_performance" %in% names(cv$performance))
  # One out-of-fold prediction per sample (no fold leaves a sample unscored).
  expect_false(anyNA(cv$prediction[["susie_predicted"]]))
})

test_that("fineMappingPipeline(QtlDataset, cvFolds=0): leaves cvResult NULL", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", cisWindow = 1000L,
                        addSusieInf = FALSE))
  expect_null(getCvResult(res, study = "study1", context = "brain",
                          trait = "ENSG_A", method = "susie"))
})

# ===========================================================================
# Cross-validation internals: .fmMakeSamplePartition / .fmCrossValidate /
# .fmSliceCv / .fmAttachCv (unit-level counterparts to the cvFolds end-to-end
# tests above).
# ===========================================================================

test_that(".fmMakeSamplePartition partitions every sample into the requested folds", {
  part <- pecotmr:::.fmMakeSamplePartition(paste0("s", 1:20), fold = 4L)
  expect_setequal(part$Sample, paste0("s", 1:20))
  expect_setequal(sort(unique(part$Fold)), 1:4)
  expect_equal(nrow(part), 20L)
})

test_that(".fmCrossValidate returns twasWeightsCv-shaped output keyed by snake method", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 60L; p <- 12L
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("v", seq_len(p))))
  y <- X[, 2] * 1.5 + rnorm(n, sd = 0.5)
  names(y) <- rownames(X)
  cv <- pecotmr:::.fmCrossValidate(
    X, y, tokens = "susie",
    methodArgs = list(susie = list()), fold = 3L,
    coverage = 0.95, verbose = 0)
  expect_named(cv, c("samplePartition", "prediction", "performance"))
  expect_setequal(colnames(cv$samplePartition), c("Sample", "Fold"))
  # Keyed by the TWAS snake method name (adapter methodKey base).
  expect_true("susie_predicted" %in% names(cv$prediction))
  expect_true("susie_performance" %in% names(cv$performance))
  pred <- cv$prediction[["susie_predicted"]]
  expect_equal(dim(pred), c(n, 1L))
  # Every sample is held out exactly once => no missing out-of-fold predictions.
  expect_false(anyNA(pred))
  perf <- cv$performance[["susie_performance"]]
  expect_equal(colnames(perf), c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
  # A real causal signal should yield positive out-of-fold correlation.
  expect_gt(perf[1, "corr"], 0)
})

test_that(".fmCrossValidate reuses a supplied samplePartition verbatim", {
  skip_if_not_installed("susieR")
  set.seed(7)
  n <- 40L; p <- 8L
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("v", seq_len(p))))
  y <- X[, 1] + rnorm(n, sd = 0.5); names(y) <- rownames(X)
  part <- pecotmr:::.fmMakeSamplePartition(rownames(X), fold = 4L)
  cv <- pecotmr:::.fmCrossValidate(X, y, tokens = "susie",
                                   methodArgs = list(susie = list()),
                                   fold = 4L, samplePartition = part,
                                   coverage = 0.95, verbose = 0)
  expect_identical(cv$samplePartition, part)
})

test_that(".fmSliceCv / .fmAttachCv slice one method and round-trip onto an entry", {
  full <- list(
    samplePartition = data.frame(Sample = c("s1", "s2"), Fold = c(1L, 2L)),
    prediction = list(susie_predicted = matrix(1, 2, 1),
                      susie_inf_predicted = matrix(2, 2, 1)),
    performance = list(susie_performance = matrix(0, 1, 6),
                       susie_inf_performance = matrix(0, 1, 6)))
  sl <- pecotmr:::.fmSliceCv(full, "susieInf")
  expect_identical(names(sl$prediction), "susie_inf_predicted")
  expect_identical(names(sl$performance), "susie_inf_performance")
  expect_identical(sl$samplePartition, full$samplePartition)

  tl <- data.frame(variant_id = "v1", pip = 0.5, stringsAsFactors = FALSE)
  e <- FineMappingEntry("v1", list(), tl)
  e2 <- pecotmr:::.fmAttachCv(e, sl)
  expect_identical(getCvResult(e2), sl)
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

test_that("fineMappingPipeline(QtlDataset): pipCutoffToSkip drops null contexts before joint mvsusie", {
  qd <- .fmp_makeQtlDataset(contexts = c("brain", "liver", "heart"),
                            traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    # Drop the middle context (liver); keep brain + heart.
    .fmSerScreenColumns   = function(X, Y, cutoff) c(TRUE, FALSE, TRUE),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie               = .fmp_mockMvsusie(),
    create_mixture_prior  = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        pipCutoffToSkip = -1))
  # liver screened out -> the joint fit runs on brain + heart only.
  expect_equal(nrow(res), 2L)
  expect_setequal(getContexts(res), c("brain", "heart"))
})

test_that("fineMappingPipeline(QtlDataset): pipCutoffToSkip skips mvsusie when < 2 contexts survive", {
  # susie runs alongside so the result still has rows after mvsusie is skipped.
  qd <- .fmp_makeQtlDataset(contexts = c("brain", "liver"), traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .fmSerScreen          = function(X, y, cutoff) TRUE,   # keep univariate susie
    .fmSerScreenColumns   = function(X, Y, cutoff) c(TRUE, FALSE),  # only brain
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie               = .fmp_mockMvsusie(),
    create_mixture_prior  = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = c("susie", "mvsusie"),
                        cisWindow = 1000L, addSusieInf = FALSE,
                        pipCutoffToSkip = -1))
  # mvsusie skipped (only 1 context survives); susie still produced per-context.
  expect_setequal(getMethodNames(res), "susie")
  expect_false("mvsusie" %in% getMethodNames(res))
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
  # Auto-detection now routes through the joint engine: each per-context group
  # is cross-trait, so every row tags its co-fit trait membership.
  expect_true("jointTraits" %in% names(res))
  expect_true(all(grepl("ENSG_A;ENSG_B|ENSG_B;ENSG_A",
                        as.character(res$jointTraits))))
})

test_that("fineMappingPipeline(QtlDataset): multi-trait auto-detection USES the data-driven mr.mash prior", {
  # Regression for the original bug: the old multi-trait path hardcoded the
  # canonical create_mixture_prior(R=ncol). Routed through the engine, it now
  # looks up a prior cross-trait mr.mash fit (study, context fixed; trait
  # match-any) and reweights from it.
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  fitParts <- list(dataDrivenPriorMatrices = list(U = list(K = diag(2)),
                                                  w = c(K = 1)),
                   w0 = c(K_grid1 = 1), V = diag(2))
  mkE <- function() TwasWeightsEntry(variantIds = c("v1", "v2"),
                                     weights = c(0.1, 0.2), fits = fitParts)
  # What twasWeightsPipeline(jointSpecification='trait') emits: per-trait rows
  # each carrying the SHARED joint fit for (study1, brain).
  tw <- TwasWeights(study = c("study1", "study1"), context = c("brain", "brain"),
                    trait = c("ENSG_A", "ENSG_B"), method = c("mrmash", "mrmash"),
                    entry = list(mkE(), mkE()))
  sawMixturePrior <- FALSE
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    rescaleCovW0          = function(w0) c(K = 1),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie              = function(X, Y, prior_variance, coverage, ...)
      list(token = "mvsusie"),
    create_mixture_prior = function(...) {
      if (!is.null(list(...)$mixture_prior)) sawMixturePrior <<- TRUE
      "PRIOR"
    },
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        twasWeights = tw))
  expect_equal(nrow(res), 2L)
  expect_true(sawMixturePrior)   # data-driven prior built, not canonical
})

test_that("fineMappingPipeline(QtlDataset): multi-trait without twasWeights keeps the canonical prior", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  sawMixturePrior <- FALSE; sawCanonical <- FALSE
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie              = .fmp_mockMvsusie(),
    create_mixture_prior = function(...) {
      a <- list(...)
      if (!is.null(a$mixture_prior)) sawMixturePrior <<- TRUE
      if (!is.null(a$R))             sawCanonical    <<- TRUE
      "PRIOR"
    },
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L))
  expect_false(sawMixturePrior)
  expect_true(sawCanonical)
})

test_that("fineMappingPipeline(QtlDataset): mvsusie resume cache short-circuits the joint fitter", {
  # All conditions of the cross-trait group are already in the prior partial
  # result -> the engine reuses the cached entries and never calls mvsusie.
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  cachedEntry <- function() FineMappingEntry(
    variantIds = paste0("v", 1:3),
    susieFit = list(token = "mvsusie_cached"),
    topLoci  = data.frame(variant_id = paste0("v", 1:3),
                          pip = c(0.9, 0.5, 0.1), stringsAsFactors = FALSE))
  cache <- QtlFineMappingResult(
    study = c("study1", "study1"), context = c("brain", "brain"),
    trait = c("ENSG_A", "ENSG_B"), method = c("mvsusie", "mvsusie"),
    entry = list(cachedEntry(), cachedEntry()))
  mv_calls <- 0
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie              = function(...) { mv_calls <<- mv_calls + 1L; list() },
    create_mixture_prior = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        fineMappingResult = cache))
  expect_equal(mv_calls, 0L)        # cache hit -> fitter never called
  expect_equal(nrow(res), 2L)
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='context' produces one joint row per trait", {
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
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        jointSpecification = "context"))
  expect_s4_class(res, "QtlFineMappingResult")
  # Per-context rows: each trait's 2-context joint emits 2 rows (4 total),
  # sharing the joint fit; jointContexts tags each with the co-fit membership.
  expect_equal(nrow(res), 4L)
  expect_true("jointContexts" %in% names(res))
  expect_setequal(as.character(res$context), c("brain", "liver"))
  expect_setequal(getTraits(res), c("ENSG_A", "ENSG_B"))
  expect_true(all(grepl("brain;liver|liver;brain",
                        as.character(res$jointContexts))))
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='context' + univariate compose without double-fitting mvsusie", {
  qd <- .fmp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie               = .fmp_mockMvsusie(),
    create_mixture_prior  = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = c("susie", "mvsusie"), cisWindow = 1000L,
                        jointSpecification = "context",
                        addSusieInf = FALSE))
  # susie -> 2 univariate rows (one per context); mvsusie joint over 2 contexts
  # -> 2 per-context rows sharing the joint fit. 4 rows total.
  expect_equal(nrow(res), 4L)
  expect_equal(sum(as.character(res$method) == "mvsusie"), 2L)
  expect_equal(sum(as.character(res$method) == "susie"), 2L)
  # Univariate rows have NA in jointContexts; both mvsusie rows carry membership.
  jc <- as.character(res$jointContexts)
  expect_equal(sum(is.na(jc)), 2L)
  expect_equal(sum(!is.na(jc)), 2L)
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='context' with only one context skips with message", {
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
  expect_error(
    suppressMessages(
      fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                          jointSpecification = "context")),
    "no joint fits produced")
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='trait' produces one joint row per context", {
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
    fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L,
                        jointSpecification = "trait"))
  expect_s4_class(res, "QtlFineMappingResult")
  # Per-trait rows: each context's 2-trait joint emits 2 rows (4 total).
  expect_equal(nrow(res), 4L)
  expect_true("jointTraits" %in% names(res))
  expect_setequal(as.character(res$trait), c("ENSG_A", "ENSG_B"))
  expect_setequal(as.character(res$context), c("brain", "liver"))
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='trait' with fsusie wires fsusieR::susiF", {
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
    fineMappingPipeline(qd, methods = "fsusie", cisWindow = 1000L,
                        jointSpecification = "trait"))
  expect_s4_class(res, "QtlFineMappingResult")
  # fsusie over 2 traits emits one per-trait row each (shared functional fit).
  expect_equal(nrow(res), 2L)
  expect_true(all(as.character(res$method) == "fsusie"))
  expect_setequal(as.character(res$trait), c("ENSG_A", "ENSG_B"))
  expect_true("jointTraits" %in% names(res))
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
# fineMappingPipeline(MultiStudyQtlDataset)
# ===========================================================================

test_that("fineMappingPipeline(MultiStudyQtlDataset): aggregates results across constituent QtlDatasets", {
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
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
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

test_that("fineMappingPipeline(MultiStudyQtlDataset): jointRegions=FALSE merges per region in each study", {
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
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  regions <- GenomicRanges::GRanges("chr1",
               IRanges::IRanges(c(50L, 350L), c(350L, 650L)))
  res <- suppressMessages(fineMappingPipeline(
    mt, methods = "susie", traitId = "ENSG_A",
    region = regions, jointRegions = FALSE, addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  # one merged row per study (region collapsed into the entry).
  expect_equal(nrow(res), 2L)
  expect_setequal(getStudy(res), c("s1", "s2"))
  fit <- getSusieFit(res, study = "s1", context = "brain",
                     trait = "ENSG_A", method = "susie")
  expect_equal(names(fit), c("region1", "region2"))
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): with embedded QtlSumStats stamps the ldSketch", {
  qd <- QtlDataset(
    study              = "s1",
    genotypes          = .fmp_makeHandle(),
    phenotypes         = list(brain = .fmp_makeSe(traits = "ENSG_A")),
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0))
  ss <- .fmp_makeQtlSumStats()
  mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd), sumStats = ss)
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

test_that("fineMappingPipeline(QtlSumStats): mvsusie auto-detection fits cross-context per (study, trait)", {
  # Multi-context (c1, c2) single-trait collection: mvsusie without an explicit
  # jointSpecification routes through the joint engine -> per-context rows that
  # share the cross-context RSS fit, each tagged with the co-fit membership.
  ss <- .fmp_makeMultiCtxQtlSumStats()
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie_rss          = function(Z, R, N, prior_variance, coverage, ...)
      list(token = "mvsusie_rss"),
    create_mixture_prior = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  res <- suppressMessages(fineMappingPipeline(ss, methods = "mvsusie"))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)
  expect_setequal(as.character(res$context), c("c1", "c2"))
  expect_true("jointContexts" %in% names(res))
  expect_true(all(grepl("c1;c2|c2;c1", as.character(res$jointContexts))))
})

test_that("fineMappingPipeline(QtlSumStats): cache hit short-circuits the RSS fitter", {
  ss <- .fmp_makeQtlSumStats()
  cachedEntry <- FineMappingEntry(
    variantIds = paste0("v", 1:5),
    susieFit = list(token = "susie_cached"),
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
    susieFit = list(token = "susie_cached"),
    topLoci    = data.frame(variant_id = paste0("v", 1:5),
                             pip = seq(0.9, 0.1, length.out = 5),
                             stringsAsFactors = FALSE))
  # The GwasSumStats branch keys its cache by (study, method, region_id)
  # where region_id is "{seqname}_{minPos}_{maxPos}" derived from the
  # entry's GRanges. .fmp_makeSumstatsGr() yields chr1 positions 100..500,
  # so the cache row must use region_id = "chr1_100_500" to hit.
  cache <- GwasFineMappingResult(
    study = "G1", method = "susie",
    region_id = "chr1_100_500",
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
    susieFit = list(token = "susie_cached"),
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
    susieFit = list(token = "susie_cached"),
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


context("univariate_pipeline")

# ===========================================================================
# Post-S4-refactor cleanup
# ---------------------------------------------------------------------------
# The S4 refactor removed `univariateAnalysisPipeline()`,
# `rssAnalysisPipeline()`, `loadStudyLd()`, `loadRssData()`, and
# `rssBasicQc()` as functional pipelines. They are now `.Deprecated()`
# no-ops that return `NULL` invisibly. The `QcResult` /
# `FineMappingResult` classes and `regionDataToSusieRssInput()` helper
# were also removed. Their replacements live behind a different S4 API
# (`fineMappingPipeline()` dispatched on `GwasSumStats` / `QtlSumStats` /
# `QtlDataset`, `summaryStatsQc()`, `FineMappingEntry()`) with different
# signatures and contracts, so the legacy mocked pipeline tests cannot
# be ported in a meaningful way and have been removed.
#
# What survives in this file:
#   * Deprecation no-op checks for the removed public wrappers.
#   * Tests for `resolveLdInput()`, which is still a live internal helper
#     in `dentistQc.R` with an unchanged signature.
# ===========================================================================

# ===========================================================================
# resolveLdInput (still-live internal helper in dentistQc.R)
# ===========================================================================

test_that("resolveLdInput errors when both R and X are NULL", {
  expect_error(
    pecotmr:::resolveLdInput(R = NULL, X = NULL),
    "Either R .* or X .* must be provided"
  )
})

test_that("resolveLdInput errors when both R and X are provided", {
  R <- diag(5)
  X <- matrix(rnorm(50), 10, 5)
  expect_error(
    pecotmr:::resolveLdInput(R = R, X = X),
    "Provide either R or X, not both"
  )
})

test_that("resolveLdInput with R returns R unchanged", {
  R <- diag(5)
  result <- pecotmr:::resolveLdInput(R = R, nSample = 100)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100)
})

test_that("resolveLdInput with X computes LD and infers nSample", {
  set.seed(42)
  X <- matrix(rbinom(200, 2, 0.3), 20, 10)
  result <- pecotmr:::resolveLdInput(X = X)
  expect_equal(nrow(result$R), 10)
  expect_equal(ncol(result$R), 10)
  expect_equal(result$nSample, 20)
})

test_that("resolveLdInput errors when nSample required but missing", {
  R <- diag(5)
  expect_error(
    pecotmr:::resolveLdInput(R = R, needNSample = TRUE),
    "nSample is required"
  )
})

test_that("resolveLdInput does not error when nSample not needed", {
  R <- diag(5)
  result <- pecotmr:::resolveLdInput(R = R, needNSample = FALSE)
  expect_equal(result$R, R)
  expect_null(result$nSample)
})

# ===========================================================================
# Deprecation no-op checks for removed public wrappers
# ===========================================================================

test_that("univariateAnalysisPipeline is a deprecated no-op", {
  expect_warning(res <- univariateAnalysisPipeline(), "deprecated|removed")
  expect_null(res)
})

test_that("rssAnalysisPipeline is a deprecated no-op", {
  expect_warning(res <- rssAnalysisPipeline(), "deprecated|removed")
  expect_null(res)
})

test_that("loadStudyLd is a deprecated no-op", {
  expect_warning(res <- loadStudyLd(), "deprecated|removed")
  expect_null(res)
})

test_that("loadRssData is a deprecated no-op", {
  expect_warning(res <- loadRssData(), "deprecated|removed")
  expect_null(res)
})

# ===========================================================================
# Residualization flag propagation
# ===========================================================================
# .resPickFlags() walks up the call stack and harvests the four
# residualization flags from whichever frame defines them. The
# fineMappingPipeline / twasWeightsPipeline setMethod signatures
# define these flags so they reach `getResidualized{Phenotypes,
# Genotypes}` via the .fmResid* wrappers without per-call-site
# threading.

test_that(".resPickFlags picks up flags from the enclosing frame", {
  outerFn <- function() {
    # Mirror the QtlDataset setMethod's residualization signature.
    residualizePhenotypeCovariates <- FALSE
    residualizeGenotypeCovariates  <- TRUE
    phenotypeCovariatesToResidualize <- c("age", "sex")
    genotypeCovariatesToResidualize  <- NULL
    innerFn <- function() {
      pecotmr:::.resPickFlags()
    }
    innerFn()
  }
  flags <- outerFn()
  expect_false(flags$residualizePhenotypeCovariates)
  expect_true(flags$residualizeGenotypeCovariates)
  expect_equal(flags$phenotypeCovariatesToResidualize, c("age", "sex"))
  expect_null(flags$genotypeCovariatesToResidualize)
})

test_that(".resPickFlags returns an empty list when nothing is in scope", {
  flags <- pecotmr:::.resPickFlags()
  # Top-level call should not pick up any of the flags (the names are
  # not defined here).
  expect_false(any(c("residualizePhenotypeCovariates",
                     "residualizeGenotypeCovariates") %in% names(flags)))
})

test_that(".fmResidGeno / .fmResidPheno forward picked-up flags to the real accessors", {
  capturedGeno <- NULL
  capturedPheno <- NULL
  fakeGeno <- function(x, ...) {
    capturedGeno <<- list(...); matrix(0, 0, 0)
  }
  fakePheno <- function(x, ...) {
    capturedPheno <<- list(...); matrix(0, 0, 0)
  }
  local_mocked_bindings(
    getResidualizedGenotypes  = fakeGeno,
    getResidualizedPhenotypes = fakePheno,
    .package = "pecotmr")

  # Emulate the setMethod frame: define the four flags then call the
  # wrappers.
  outerFn <- function() {
    residualizePhenotypeCovariates <- FALSE
    residualizeGenotypeCovariates  <- TRUE
    phenotypeCovariatesToResidualize <- "age"
    genotypeCovariatesToResidualize  <- NULL
    pecotmr:::.fmResidGeno(NULL, contexts = "c1")
    pecotmr:::.fmResidPheno(NULL, contexts = "c1")
  }
  outerFn()
  expect_false(capturedGeno$residualizePhenotypeCovariates)
  expect_true(capturedGeno$residualizeGenotypeCovariates)
  expect_equal(capturedGeno$phenotypeCovariatesToResidualize, "age")
  expect_false(capturedPheno$residualizePhenotypeCovariates)
  expect_true(capturedPheno$residualizeGenotypeCovariates)
})

# ===========================================================================
# Removed during the post-S4-refactor cleanup (for traceability)
# ---------------------------------------------------------------------------
# univariateAnalysisPipeline input-validation tests (12):
#   * X non-matrix / non-numeric / Y non-numeric / multi-column Y /
#     single-column Y / X-Y row mismatch / maf length mismatch /
#     maf out of bounds / invalid xScalar / wrong-length xScalar /
#     non-numeric yScalar / non-positive L / non-positive lGreedy.
# univariateAnalysisPipeline functional tests (4):
#   * minimal valid input; twasWeights = TRUE; respects xScalar;
#     cvFolds = 0.
# univariateAnalysisPipeline mocked-pipeline tests (16):
#   * pip_cutoff_to_skip branch (3); LD reference filtering (2);
#     filter_X branch (2); main susie + post-processing (4);
#     TWAS weights (4); coverage forwarding (2); combined LD +
#     filter_X (1); null imiss/maf cutoffs (1).
# univariateAnalysisPipeline af/maf single-source tests (6):
#   * af supplied (no maf); maf supplied (no af); both disagreeing;
#     neither with mafCutoff errors; neither without mafCutoff runs;
#     af-derived MAF parity with explicit maf.
# rssAnalysisPipeline tests (~21):
#   * input-file validation; empty sumstats branches; pip_cutoff_to_skip
#     branches; full QC + imputation + fine-mapping; method-name
#     conventions (NO_QC / SLALOM / DENTIST / RAISS); outlierNumber
#     plumbing; finemappingMethod = NULL skip; zMismatchQc = NULL;
#     impute = FALSE; diagnostics branches (BCR + SER re-analysis);
#     finemappingOpts passthrough; mafCutoff (af single source);
#     finemapping-defaults passthrough; X-as-LD genotype path;
#     mixture-LD passthrough.
# regionDataToSusieRssInput tests (3):
#   * correlation-backed input; genotype-backed input; rejects
#     default rownames as variant IDs.
# loadStudyLd tests (2):
#   * single path returns loadLdMatrix output unchanged; comma-separated
#     paths build a mixture LdData with a list of genotype handles.
#
# Replacement APIs (do NOT speculatively port — that is out of scope for
# this cleanup): `fineMappingPipeline()` dispatched on `GwasSumStats` /
# `QtlSumStats` / `QtlDataset`; `summaryStatsQc()` (returns a SumStats
# with `getQcInfo()` populated); `FineMappingEntry(variantIds, ...)`;
# `result$finemappingEntry` (was `result$finemappingResult`).
# ===========================================================================

# ===========================================================================
# Multi-region / jointRegions (P3)
# ===========================================================================

test_that(".fmCsIdx / .fmRelabelCs parse and renumber <method>_<idx> labels", {
  expect_equal(pecotmr:::.fmCsIdx(c("susie_0", "susie_1", "susie_2")),
               c(0L, 1L, 2L))
  expect_equal(pecotmr:::.fmRelabelCs(c("susie_0", "susie_1", "susie_2"), 3L),
               c("susie_0", "susie_4", "susie_5"))
  # offset 0 is a no-op; the "_0" not-in-CS sentinel is always preserved.
  expect_equal(pecotmr:::.fmRelabelCs(c("susie_0", "susie_1"), 0L),
               c("susie_0", "susie_1"))
})

test_that(".fmMergeEntries concatenates variants, renumbers CS, lists susieFit", {
  mk <- function(vids, cs95, fit) FineMappingEntry(
    variantIds = vids, susieFit = fit,
    topLoci = data.frame(variant_id = vids, pip = rep(0.5, length(vids)),
                         cs_95 = cs95, stringsAsFactors = FALSE))
  e1 <- mk(c("a", "b"), c("susie_1", "susie_0"), list(tag = "f1"))
  e2 <- mk(c("c", "d"), c("susie_1", "susie_1"), list(tag = "f2"))
  m <- pecotmr:::.fmMergeEntries(list(e1, e2))
  expect_s4_class(m, "FineMappingEntry")
  expect_equal(m@variantIds, c("a", "b", "c", "d"))
  expect_equal(as.character(m@topLoci$variant_id), c("a", "b", "c", "d"))
  # e1's max CS index is 1, so e2's "susie_1" is renumbered to "susie_2".
  expect_equal(m@topLoci$cs_95, c("susie_1", "susie_0", "susie_2", "susie_2"))
  expect_true(is.list(m@susieFit))
  expect_equal(names(m@susieFit), c("region1", "region2"))
})

test_that(".fmMergeEntries returns a single entry unchanged", {
  e <- FineMappingEntry(variantIds = "a", susieFit = list(),
                        topLoci = data.frame(variant_id = "a", pip = 0.5))
  expect_identical(pecotmr:::.fmMergeEntries(list(e)), e)
})

test_that("fineMappingPipeline(QtlDataset): region + cisWindow is rejected", {
  qd <- .fmp_makeQtlDataset()
  expect_error(
    fineMappingPipeline(
      qd, methods = "susie",
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 600)),
      cisWindow = 1000L),
    "either `region` or `cisWindow`")
})

test_that("fineMappingPipeline(QtlDataset): jointRegions=FALSE merges per-region fits", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  regions <- GenomicRanges::GRanges("chr1",
               IRanges::IRanges(c(50L, 350L), c(350L, 650L)))
  res <- suppressMessages(fineMappingPipeline(
    qd, methods = "susie", traitId = "ENSG_A",
    region = regions, jointRegions = FALSE, addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  # 1 ctx x 1 trait x 1 method -> a single merged row, not one per region.
  expect_equal(nrow(res), 1L)
  fit <- getSusieFit(res, study = "study1", context = "brain",
                     trait = "ENSG_A", method = "susie")
  expect_equal(names(fit), c("region1", "region2"))  # per-region fit list
})

test_that("fineMappingPipeline(QtlDataset): jointRegions=TRUE fits one concatenated block", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv      = .fmp_mockFitIndiv(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  regions <- GenomicRanges::GRanges("chr1",
               IRanges::IRanges(c(50L, 350L), c(350L, 650L)))
  res <- suppressMessages(fineMappingPipeline(
    qd, methods = "susie", traitId = "ENSG_A",
    region = regions, jointRegions = TRUE, addSusieInf = FALSE))
  expect_equal(nrow(res), 1L)
  fit <- getSusieFit(res, study = "study1", context = "brain",
                     trait = "ENSG_A", method = "susie")
  # one concatenated fit -> the mock fit object, not a per-region list.
  expect_false(identical(names(fit), c("region1", "region2")))
})

test_that("fineMappingPipeline(QtlDataset): mvsusie jointRegions=FALSE merges per-region fits", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie              = .fmp_mockMvsusie(),
    create_mixture_prior = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  regions <- GenomicRanges::GRanges("chr1",
               IRanges::IRanges(c(50L, 350L), c(350L, 650L)))
  res <- suppressMessages(fineMappingPipeline(
    qd, methods = "mvsusie", traitId = c("ENSG_A", "ENSG_B"),
    region = regions, jointRegions = FALSE))
  expect_equal(nrow(res), 2L)  # joint fit fanned out to both traits
  fit <- getSusieFit(res, study = "study1", context = "brain",
                     trait = "ENSG_A", method = "mvsusie")
  expect_equal(names(fit), c("region1", "region2"))  # per-region merged fit list
})

test_that("fineMappingPipeline(QtlDataset): fsusie jointRegions=FALSE merges per-region fits", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(susiF = .fmp_mockSusiF(), .package = "fsusieR")
  regions <- GenomicRanges::GRanges("chr1",
               IRanges::IRanges(c(50L, 350L), c(350L, 650L)))
  res <- suppressMessages(fineMappingPipeline(
    qd, methods = "fsusie", traitId = c("ENSG_A", "ENSG_B"),
    region = regions, jointRegions = FALSE))
  expect_equal(nrow(res), 2L)
  fit <- getSusieFit(res, study = "study1", context = "brain",
                     trait = "ENSG_A", method = "fsusie")
  expect_equal(names(fit), c("region1", "region2"))
})

test_that("fineMappingPipeline(QtlDataset): jointSpec + jointRegions=FALSE merges per region", {
  qd <- .fmp_makeQtlDataset(contexts = c("brain", "liver"),
                            traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusie              = .fmp_mockMvsusie(),
    create_mixture_prior = .fmp_mockMixturePrior(),
    .package = "mvsusieR")
  regions <- GenomicRanges::GRanges("chr1",
               IRanges::IRanges(c(50L, 350L), c(350L, 650L)))
  res <- suppressMessages(fineMappingPipeline(
    qd, methods = "mvsusie", traitId = c("ENSG_A", "ENSG_B"),
    region = regions, jointRegions = FALSE, jointSpecification = "context"))
  expect_s4_class(res, "QtlFineMappingResult")
  # cross-context joint -> per-context rows (2 contexts x 2 traits = 4), each
  # merged across the 2 regions.
  expect_equal(nrow(res), 4L)
  expect_setequal(as.character(res$context), c("brain", "liver"))
  fit <- getSusieFit(res, study = "study1", context = "brain",
                     trait = "ENSG_A", method = "mvsusie")
  expect_equal(names(fit), c("region1", "region2"))  # merged across regions
})

# ===========================================================================
# Additional coverage: tractable pure helpers. (The real-fit internals --
# .fmPostprocessOne, mvsusie/fsusie fold weights, susieInf chaining, the
# top-PC cross-trait path -- are driven by the mocked-fit pipeline tests and
# the actual fitting is left to the external solvers.)
# ===========================================================================

test_that(".fmCheckMethodCapabilities: empty token list is a no-op", {
  expect_null(pecotmr:::.fmCheckMethodCapabilities(character(0), "QtlDataset"))
})

test_that(".fmCacheLookupGwas: NULL / non-GwasFineMappingResult -> NULL", {
  expect_null(pecotmr:::.fmCacheLookupGwas(NULL, "G1", "susie", "chr1:1-100"))
  fmr <- QtlFineMappingResult(study = "S", context = "c1", trait = "t1",
    method = "susie", entry = list(FineMappingEntry("v1", list(),
      data.frame(variant_id = "v1", pip = 0.9))))
  expect_null(pecotmr:::.fmCacheLookupGwas(fmr, "S", "susie", "chr1:1-100"))
})

test_that(".buildMvsusieReweightedPrior: empty reweighted w0 -> canonical(V)", {
  local_mocked_bindings(rescaleCovW0 = function(w0) c(zzz = 1), .package = "pecotmr")
  local_mocked_bindings(create_mixture_prior = function(...) "PV",
                        .package = "mvsusieR")
  fp <- list(dataDrivenPriorMatrices = list(U = list(compA = diag(2))),
             w0 = c(compA_grid1 = 1), V = diag(2))
  res <- pecotmr:::.buildMvsusieReweightedPrior(fp, c("c1", "c2"))
  expect_equal(res$residualVariance, diag(2))          # w0 names disjoint from U (749)
})

test_that(".fmBuildMvsusiePriorCv: NULL CV -> NULL; NULL fold fits are skipped", {
  expect_null(pecotmr:::.fmBuildMvsusiePriorCv(NULL, NULL, c("c1", "c2")))
  local_mocked_bindings(.buildMvsusieReweightedPrior = function(...) "PRIOR",
                        .package = "pecotmr")
  mvCv <- list(samplePartition = data.frame(Sample = paste0("s", 1:4),
                                            Fold = c(1, 1, 2, 2)),
               foldFits = list(fold_1 = list(w0 = 1), fold_2 = NULL))
  out <- pecotmr:::.fmBuildMvsusiePriorCv(mvCv, list(w0 = 1, V = diag(2)),
                                          c("c1", "c2"))
  expect_equal(out[["1"]], "PRIOR")
  expect_null(out[["2"]])                              # NULL fold -> next (831)
})

test_that(".fmTopPcScores: PCA scores for multi-trait Y; degenerate inputs -> NULL", {
  expect_null(pecotmr:::.fmTopPcScores(matrix(1, 5, 1), 2L))         # < 2 traits
  expect_null(pecotmr:::.fmTopPcScores(matrix(c(1, NA), 2, 2), 2L))  # < 2 complete
  expect_null(pecotmr:::.fmTopPcScores(
    cbind(a = rep(1, 4), b = rnorm(4)), 2L))                         # < 2 nonzero-var
  set.seed(1)
  Y <- matrix(rnorm(20), 10, 2, dimnames = list(paste0("s", 1:10), c("t1", "t2")))
  sc <- pecotmr:::.fmTopPcScores(Y, 2L)
  expect_equal(dim(sc), c(10L, 2L))
  expect_equal(colnames(sc), c("topPC1", "topPC2"))
})

test_that(".fmSerScreen / .fmScreenActive / .fmSerScreenColumns", {
  set.seed(2)
  X <- matrix(rnorm(40), 20, 2, dimnames = list(paste0("s", 1:20), c("v1", "v2")))
  y <- rnorm(20)
  expect_true(pecotmr:::.fmSerScreen(X, y, cutoff = 0))              # disabled
  expect_true(pecotmr:::.fmSerScreen(X, c(1, rep(NA, 19)), 0.5))     # < 2 obs (880)
  expect_type(pecotmr:::.fmSerScreen(X, y, 0.5), "logical")          # real susie fit
  local_mocked_bindings(susie = function(...) stop("boom"), .package = "susieR")
  expect_true(pecotmr:::.fmSerScreen(X, y, 0.5))                     # fit fails -> keep (886)
  expect_false(pecotmr:::.fmScreenActive(0))
  expect_true(pecotmr:::.fmScreenActive(0.5))
  expect_length(pecotmr:::.fmSerScreenColumns(X, matrix(rnorm(40), 20, 2), 0), 2L)
})

test_that(".fmMergeEntries: empty -> NULL; merges per-region entries + relabels CS", {
  expect_null(pecotmr:::.fmMergeEntries(list(NULL, NULL)))
  e1 <- FineMappingEntry("v1", list(a = 1),
    data.frame(variant_id = "v1", pip = 0.9, cs_95 = "susie_1"))
  e2 <- FineMappingEntry("v2", list(b = 2),
    data.frame(variant_id = "v2", pip = 0.8, cs_95 = "susie_1"))
  m <- pecotmr:::.fmMergeEntries(list(e1, e2))
  expect_s4_class(m, "FineMappingEntry")
  expect_equal(m@variantIds, c("v1", "v2"))
  expect_equal(m@topLoci$cs_95, c("susie_1", "susie_2"))   # region2 CS relabelled
  expect_equal(names(m@susieFit), c("region1", "region2"))
})

test_that(".fmJointBlocks: all-NULL -> NULL; single -> unchanged; many -> merged", {
  mkE <- function(v) FineMappingEntry(v, list(),
    data.frame(variant_id = v, pip = 0.9))
  expect_null(pecotmr:::.fmJointBlocks(list(1, 2), function(rg) NULL))
  expect_equal(pecotmr:::.fmJointBlocks(list(1),
    function(rg) mkE("v1"))@variantIds, "v1")
  expect_equal(pecotmr:::.fmJointBlocks(list(1, 2),
    function(rg) mkE(paste0("v", rg)))@variantIds, c("v1", "v2"))
})

test_that(".fmTwasMethodKey: bare token without adapter returned unchanged", {
  expect_equal(pecotmr:::.fmTwasMethodKey("lasso"), "lasso")        # no adapter (1170)
  expect_equal(pecotmr:::.fmTwasMethodKey("susie"), "susie")        # adapter -> stripped
})

test_that(".fmCvMetricRow: < 3 usable predictions -> all-NA row", {
  expect_true(all(is.na(pecotmr:::.fmCvMetricRow(c(1, 2), c(1, 2)))))  # < 3 (1182)
  ok <- pecotmr:::.fmCvMetricRow(c(1, 2, 3, 4, 5), c(1.1, 2, 2.9, 4, 5))
  expect_false(is.na(ok[["rsq"]]))
})

test_that(".fmSliceCv: NULL cv or missing predicted key -> NULL", {
  expect_null(pecotmr:::.fmSliceCv(NULL, "susie"))                  # 1323
  cv <- list(samplePartition = NULL,
             prediction = list(enet_predicted = matrix(0, 1, 1)),
             performance = list(enet_performance = matrix(0, 1, 6)))
  expect_null(pecotmr:::.fmSliceCv(cv, "susie"))                    # pk absent (1327)
  expect_true("enet_predicted" %in% names(pecotmr:::.fmSliceCv(cv, "enet")$prediction))
})

test_that(".fmAttachCv: NULL entry or NULL cvResult returns the entry unchanged", {
  e <- FineMappingEntry("v1", list(), data.frame(variant_id = "v1", pip = 0.9))
  expect_identical(pecotmr:::.fmAttachCv(e, NULL), e)               # 1336
  expect_null(pecotmr:::.fmAttachCv(NULL, list(x = 1)))
  expect_equal(getCvResult(pecotmr:::.fmAttachCv(e, list(samplePartition = 1))),
               list(samplePartition = 1))
})

# ---- jointSpec dispatch branches in the QtlSumStats / MultiStudy methods -----
# (dispatchers mocked; the real joint fitting is covered in test_jointSpecification.R)

test_that(".fmTopPcScores: nPCs = 0 -> k < 1 -> NULL", {
  set.seed(3)
  Y <- matrix(rnorm(20), 10, 2, dimnames = list(paste0("s", 1:10), c("t1", "t2")))
  expect_null(pecotmr:::.fmTopPcScores(Y, 0L))                      # k < 1 (858)
})

test_that("fineMappingPipeline(QtlSumStats): mvsusie-only jointSpec returns the joint result", {
  ss <- .fmp_makeQtlSumStats()
  jr <- QtlFineMappingResult(study = "Q1", context = "c1", trait = "t1",
    method = "mvsusie", entry = list(FineMappingEntry("v1", list(),
      data.frame(variant_id = "v1", pip = 0.9))))
  local_mocked_bindings(.fmDispatchJointSpecsQtlSumStats = function(...) jr,
                        .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "mvsusie", jointSpecification = "context"))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_setequal(as.character(res$method), "mvsusie")             # 1843-1857
})

test_that("fineMappingPipeline(QtlSumStats): mvsusie-only jointSpec with no fits -> error", {
  ss <- .fmp_makeQtlSumStats()
  local_mocked_bindings(.fmDispatchJointSpecsQtlSumStats = function(...) NULL,
                        .package = "pecotmr")
  expect_error(
    suppressMessages(fineMappingPipeline(ss, methods = "mvsusie",
                                         jointSpecification = "context")),
    "no joint fits produced")
})

.fmp_makeMultiStudy <- function() MultiStudyQtlDataset(
  qtlDatasets = list(study1 = .fmp_makeQtlDataset(contexts = "brain",
                                                  traits = "ENSG_A")),
  sumStats = .fmp_makeQtlSumStats())

test_that("fineMappingPipeline(MultiStudyQtlDataset): region + cisWindow is rejected", {
  expect_error(
    fineMappingPipeline(.fmp_makeMultiStudy(), methods = "mvsusie",
                        region = GenomicRanges::GRanges("chr1",
                                   IRanges::IRanges(1, 100)),
                        cisWindow = 1000L),
    "specify either")                                              # 1694
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): mvsusie-only jointSpec returns the joint result", {
  mt <- .fmp_makeMultiStudy()
  jr <- QtlFineMappingResult(study = "study1", context = "brain", trait = "ENSG_A",
    method = "mvsusie", entry = list(FineMappingEntry("v1", list(),
      data.frame(variant_id = "v1", pip = 0.9))))
  local_mocked_bindings(.fmDispatchJointSpecsMultiStudy = function(...) jr,
                        .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(mt, methods = "mvsusie", jointSpecification = "context"))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_setequal(as.character(res$method), "mvsusie")             # 1709-1726
})

# =============================================================================
# Coverage top-ups (mock the fitters, exercise the orchestration)
# =============================================================================

test_that(".fmRelabelCs returns non-matching labels unchanged", {
  fn <- pecotmr:::.fmRelabelCs
  # "nomatch" doesn't match ^(.*)_([0-9]+)$ (length(parts) != 3) -> returned
  # as-is; "susie_0" keeps the not-in-CS sentinel; "susie_1" shifts by offset.
  out <- fn(c("susie_1", "nomatch", "susie_0"), offset = 2L)
  expect_equal(out, c("susie_3", "nomatch", "susie_0"))
})

test_that(".fmCrossValidate + .fmFoldWeights cover the mvSuSiE CV path (mocked fitter)", {
  # Mock one level below the orchestration: fitMvsusie/mvsusieWeights return
  # canned outputs (sized to the per-fold training columns), so the real
  # .fmFoldWeights mvsusie branch + .fmCrossValidate fold loop run at ~no cost.
  local_mocked_bindings(
    fitMvsusie = function(X, Y, ...) list(vn = colnames(X), R = ncol(as.matrix(Y))),
    mvsusieWeights = function(mvsusieFit = NULL, ...)
      matrix(0.01, length(mvsusieFit$vn), mvsusieFit$R,
             dimnames = list(mvsusieFit$vn, NULL)),
    .package = "pecotmr")
  set.seed(1)
  n <- 30L; p <- 5L; R <- 2L
  X <- matrix(rbinom(n * p, 2, 0.4), n, p,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:p)))
  Y <- matrix(rnorm(n * R), n, R,
              dimnames = list(rownames(X), c("c1", "c2")))
  cv <- pecotmr:::.fmCrossValidate(
    X, Y, tokens = "mvsusie", methodArgs = list(mvsusie = list()),
    fold = 3L, coverage = 0.95, verbose = 0)
  expect_named(cv, c("samplePartition", "prediction", "performance"))
  expect_true("mvsusie_performance" %in% names(cv$performance))
  expect_equal(dim(cv$prediction[["mvsusie_predicted"]]), c(n, R))
})

test_that(".fmFoldWeights covers the fSuSiE branch (mocked fitter)", {
  local_mocked_bindings(
    fitFsusie = function(...) list(),
    fsusieWeights = function(fsusieFit = NULL, variantIds = NULL, ...)
      matrix(0.02, length(variantIds), 1L, dimnames = list(variantIds, NULL)),
    .package = "pecotmr")
  set.seed(2)
  n <- 30L; p <- 5L
  X <- matrix(rbinom(n * p, 2, 0.4), n, p,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:p)))
  Y <- matrix(rnorm(n * 4L), n, 4L, dimnames = list(rownames(X), NULL))
  W <- pecotmr:::.fmFoldWeights("fsusie", X, Y, coverage = 0.95,
                                userArgs = list(), pos = seq_len(p))
  expect_true(is.matrix(W))
  expect_equal(rownames(W), colnames(X))
})

test_that(".fmFitXBlock fits the susieInf indiv chain + cross-validates (mocked)", {
  local_mocked_bindings(
    .fmFitSusieIndiv = function(...) list(),
    .fmPostprocessOne = function(fit, method, dataX, dataY, ...)
      FineMappingEntry(colnames(dataX), list(),
                       data.frame(variant_id = colnames(dataX), pip = 0.5)),
    .fmFoldWeights = function(token, Xtr, Ytr, ...)
      matrix(0.01, ncol(Xtr), 1L, dimnames = list(colnames(Xtr), NULL)),
    .package = "pecotmr")
  set.seed(1)
  X <- matrix(rbinom(60, 2, 0.4), 20, 3,
              dimnames = list(paste0("s", 1:20), c("v1", "v2", "v3")))
  y <- rnorm(20)
  out <- pecotmr:::.fmFitXBlock(
    X, y, toRun = "susieInf", addSusieInf = FALSE, coverage = 0.95,
    secondaryCoverage = 0.7, signalCutoff = 0.1, minAbsCorr = 0.5,
    methodArgs = list(susieInf = list()), verbose = 1,
    ctx = "brain", tid = "ENSG_A", cvFolds = 3L)
  expect_named(out, "susieInf")
  expect_s4_class(out$susieInf, "FineMappingEntry")
})

test_that(".fmPostprocessOne wraps a fit into a FineMappingEntry", {
  local_mocked_bindings(
    postprocessFinemappingFits = function(...) list(x = 1),
    formatFinemappingOutput = function(post, primaryMethod, ...)
      list(finemappingEntry = FineMappingEntry("v1", list(),
             data.frame(variant_id = "v1", pip = 0.5))),
    .package = "pecotmr")
  ent <- pecotmr:::.fmPostprocessOne(
    fit = list(), method = "susie",
    dataX = matrix(0, 2, 1, dimnames = list(NULL, "v1")), dataY = c(1, 2),
    coverage = 0.95, secondaryCoverage = 0.7, signalCutoff = 0.1,
    minAbsCorr = 0.5)
  expect_s4_class(ent, "FineMappingEntry")
})

test_that(".fmPostprocessOne errors when output carries no FineMappingEntry", {
  local_mocked_bindings(
    postprocessFinemappingFits = function(...) list(),
    formatFinemappingOutput = function(...) list(finemappingEntry = "nope"),
    .package = "pecotmr")
  expect_error(
    pecotmr:::.fmPostprocessOne(list(), "susie",
      matrix(0, 2, 1), c(1, 2), 0.95, 0.7, 0.1, 0.5),
    "FineMappingEntry payload")
})

test_that(".fmCrossValidate covers per-fold prior, NULL-weights, and no-overlap branches", {
  # mvPriorCv supplies a prior for fold "1" only: fold 1 takes the per-fold
  # prior (else-branch) and returns weights whose rownames don't overlap the
  # test columns (no-common `next`); fold 2 has no prior, so .fmFoldWeights
  # returns NULL (NULL-weights `next`).
  local_mocked_bindings(
    .fmFoldWeights = function(token, Xtr, Ytr, coverage, userArgs, pos, mvPrior) {
      if (is.null(mvPrior)) return(NULL)
      matrix(0.5, 1L, 1L, dimnames = list("not_a_variant", NULL))
    },
    .package = "pecotmr")
  set.seed(3)
  X <- matrix(rbinom(40, 2, 0.4), 20, 2,
              dimnames = list(paste0("s", 1:20), c("v1", "v2")))
  Y <- matrix(rnorm(40), 20, 2, dimnames = list(rownames(X), c("c1", "c2")))
  cv <- pecotmr:::.fmCrossValidate(
    X, Y, tokens = "mvsusie", methodArgs = list(mvsusie = list()),
    fold = 2L, coverage = 0.95, verbose = 0,
    mvPriorCv = list("1" = list(priorVariance = diag(2))))
  expect_named(cv, c("samplePartition", "prediction", "performance"))
})

test_that("fineMappingPipeline(QtlSumStats): susieInf RSS chain (mocked)", {
  ss <- .fmp_makeQtlSumStats()
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss        = .fmp_mockFitRss(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(ss, methods = "susieInf", addSusieInf = FALSE,
                        verbose = 1))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_setequal(as.character(res$method), "susieInf")
})

test_that("fineMappingPipeline(GwasSumStats): susieInf RSS chain (mocked)", {
  gss <- .fmp_makeGwasSumStats()
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieRss        = .fmp_mockFitRss(),
    .fmPostprocessOne     = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(gss, methods = "susieInf", addSusieInf = FALSE,
                        verbose = 1))
  expect_s4_class(res, "GwasFineMappingResult")
  expect_setequal(getMethodNames(res), "susieInf")
})

test_that(".fmFoldWeights covers mvPrior residual var, missing rownames, unknown token", {
  local_mocked_bindings(
    fitMvsusie = function(X, Y, ...) list(vn = colnames(X), R = ncol(as.matrix(Y))),
    # weights WITHOUT rownames -> .fmFoldWeights sets them from Xtr (1231)
    mvsusieWeights = function(mvsusieFit = NULL, ...)
      matrix(0.01, length(mvsusieFit$vn), mvsusieFit$R),
    .package = "pecotmr")
  X <- matrix(rbinom(40, 2, 0.4), 20, 2,
              dimnames = list(paste0("s", 1:20), c("v1", "v2")))
  Y <- matrix(rnorm(40), 20, 2)
  # mvPrior carrying a residualVariance exercises line 1227
  W <- pecotmr:::.fmFoldWeights(
    "mvsusie", X, Y, coverage = 0.95, userArgs = list(), pos = NULL,
    mvPrior = list(priorVariance = diag(2), residualVariance = diag(2)))
  expect_equal(rownames(W), colnames(X))
  # an unknown token falls through to NULL (1241)
  expect_null(pecotmr:::.fmFoldWeights("bogus", X, Y, 0.95, list(), NULL))
})

test_that(".fmCrossValidate returns NULL for empty tokens", {
  X <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("v1", "v2")))
  expect_null(pecotmr:::.fmCrossValidate(
    X, matrix(0, 10, 1), tokens = character(0), methodArgs = list(), fold = 2L))
})

test_that(".fmCrossValidate fills Y rownames and reports per-fold fit failures", {
  local_mocked_bindings(.fmFoldWeights = function(...) stop("boom"),
                        .package = "pecotmr")
  X <- matrix(rbinom(40, 2, 0.4), 20, 2,
              dimnames = list(paste0("s", 1:20), c("v1", "v2")))
  Y <- matrix(rnorm(20), 20, 1)   # no rownames -> filled from X (1258)
  expect_message(
    cv <- suppressWarnings(pecotmr:::.fmCrossValidate(
      X, Y, tokens = "susie", methodArgs = list(susie = list()),
      fold = 2L, verbose = 1)),
    "CV fold .* failed")            # 1291-1294
})

test_that(".fmCrossValidate skips a fold that holds out every sample", {
  X <- matrix(rbinom(40, 2, 0.4), 20, 2,
              dimnames = list(paste0("s", 1:20), c("v1", "v2")))
  Y <- matrix(rnorm(20), 20, 1, dimnames = list(rownames(X), NULL))
  sp <- data.frame(Sample = rownames(X), Fold = 1L)   # single fold = all test -> 1273
  cv <- pecotmr:::.fmCrossValidate(
    X, Y, tokens = "susie", methodArgs = list(susie = list()), fold = 1L,
    samplePartition = sp, verbose = 0)
  expect_true(all(is.na(cv$prediction[["susie_predicted"]])))
})

# --- method-level branches (drive the pipeline methods with mocked fitters) ---

test_that(".fmAfForX returns NULL when getAf yields nothing", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(getAf = function(...) NULL, .package = "pecotmr")
  X <- matrix(0, 3, 2, dimnames = list(paste0("s", 1:3), c("v1", "v2")))
  expect_null(pecotmr:::.fmAfForX(qd, X, traitId = "ENSG_A"))
})

test_that("fineMappingPipeline(QtlDataset): explicit valid contexts arg is honored", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", contexts = "brain",
                        cisWindow = 1000L, addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
})

test_that("fineMappingPipeline(QtlDataset): region selects traits by rowRanges overlap", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .package = "pecotmr")
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3000L))
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", region = region,
                        addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
})

test_that("fineMappingPipeline(QtlDataset): too few shared samples errors", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    # residualized X carries a sample absent from Y -> < 2 shared samples
    .fmResidGeno = function(x, ...)
      matrix(0, 1L, 2L, dimnames = list("ghost_sample", c("v1", "v2"))),
    .package = "pecotmr")
  expect_error(
    suppressMessages(fineMappingPipeline(qd, methods = "susie",
                     cisWindow = 1000L, addSusieInf = FALSE)),
    "too few shared samples")
})

test_that("fineMappingPipeline(QtlDataset): errors when no tuple produces a result", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmSerScreen = function(...) FALSE,   # screen out every block
    .package = "pecotmr")
  expect_error(
    suppressMessages(fineMappingPipeline(qd, methods = "susie",
                     cisWindow = 1000L, addSusieInf = FALSE)),
    "no .*tuples")
})

test_that("fineMappingPipeline(QtlDataset): usePCA fine-maps top PCs of a multi-trait context", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .fmSerScreen = function(...) TRUE,
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", usePCA = TRUE, nPCs = 1L,
                        cisWindow = 1000L, addSusieInf = FALSE))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_true(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(QtlDataset): usePCA skips single-trait contexts", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .fmSerScreen = function(...) TRUE,
    .package = "pecotmr")
  res <- suppressMessages(
    fineMappingPipeline(qd, methods = "susie", usePCA = TRUE, nPCs = 1L,
                        cisWindow = 1000L, addSusieInf = FALSE))
  # single-trait context -> PC loop hits `length(traits) < 2L` next; only the
  # univariate susie row survives, no topPC pseudo-trait.
  expect_false(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): jointSpec with no intersecting scope errors", {
  mt <- .fmp_makeMultiStudy()
  # The joint engine yields nothing for this scope -> the no-joint-fits stop.
  local_mocked_bindings(.fmDispatchJointSpecsMultiStudy = function(...) NULL,
                        .package = "pecotmr")
  expect_error(
    suppressMessages(fineMappingPipeline(mt, methods = "mvsusie",
                     jointSpecification = "context")),
    "no joint fits produced")
})

# --- usePCA sub-branches. Use methods="mvsusie" so the univariate dispatch is
# skipped (no 1531 stop / SER entanglement); the PCA path still runs susie, and
# the mvsusie joint dispatch is mocked to keep a non-empty result (no 1647).
.fmp_jr <- function()
  QtlFineMappingResult(study = "study1", context = "brain", trait = "ENSG_A",
    method = "mvsusie",
    entry = list(FineMappingEntry("v1", list(),
                                  data.frame(variant_id = "v1", pip = 0.9))))

test_that("fineMappingPipeline(QtlDataset): usePCA skips a context whose PCA yields no scores", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .fmTopPcScores = function(...) NULL,                       # 1577 next
    .package = "pecotmr")
  res <- suppressMessages(fineMappingPipeline(qd, methods = "susie", usePCA = TRUE,
    nPCs = 1L, cisWindow = 1000L, addSusieInf = FALSE))
  expect_false(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(QtlDataset): usePCA reuses a cached PC entry", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  cachedFMR <- QtlFineMappingResult(study = "study1", context = "brain",
    trait = "topPC1", method = "susie",
    entry = list(FineMappingEntry("v1", list(),
                                  data.frame(variant_id = "v1", pip = 0.9))))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .fmTopPcScores = function(Y, nPCs) matrix(rnorm(nrow(Y)), nrow(Y), 1L,
      dimnames = list(rownames(Y), "topPC1")),
    .package = "pecotmr")
  res <- suppressMessages(fineMappingPipeline(qd, methods = "susie", usePCA = TRUE,
    nPCs = 1L, cisWindow = 1000L, addSusieInf = FALSE, fineMappingResult = cachedFMR))
  expect_true(any(as.character(res$trait) == "topPC1"))       # 1584 cache hit
})

test_that("fineMappingPipeline(QtlDataset): usePCA + region uses the region genotype block", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmFitSusieIndiv = .fmp_mockFitIndiv(), .fmPostprocessOne = .fmp_mockPostprocess(),
    .fmSerScreen = function(...) TRUE,
    .fmDispatchJointSpecsQtlDataset = function(...) NULL,
    .package = "pecotmr")
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3000L))
  res <- suppressMessages(fineMappingPipeline(qd, methods = "mvsusie", usePCA = TRUE,
    nPCs = 1L, region = region))
  expect_true(any(grepl("PC", as.character(res$trait))))      # 1591-1592 region block
})

test_that("fineMappingPipeline(QtlDataset): usePCA skips a PC block with too few shared samples", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmResidGeno = function(x, ...) matrix(0, 1L, 2L,
      dimnames = list("ghost_sample", c("v1", "v2"))),        # PC common < 2 -> 1595
    .fmDispatchJointSpecsQtlDataset = function(...) .fmp_jr(),
    .package = "pecotmr")
  res <- suppressMessages(fineMappingPipeline(qd, methods = "mvsusie", usePCA = TRUE,
    nPCs = 1L, cisWindow = 1000L))
  expect_false(any(grepl("PC", as.character(res$trait))))     # 1595 + 1607
})

test_that("fineMappingPipeline(QtlDataset): usePCA skips a PC block screened out by SER", {
  qd <- .fmp_makeQtlDataset(contexts = "brain", traits = c("ENSG_A", "ENSG_B"))
  local_mocked_bindings(
    extractBlockGenotypes = .fmp_mockExtractor(),
    .fmSerScreen = function(...) FALSE,                        # PC screened -> 1597
    .fmDispatchJointSpecsQtlDataset = function(...) .fmp_jr(),
    .package = "pecotmr")
  res <- suppressMessages(fineMappingPipeline(qd, methods = "mvsusie", usePCA = TRUE,
    nPCs = 1L, cisWindow = 1000L))
  expect_false(any(grepl("PC", as.character(res$trait))))     # 1597 + 1607
})
