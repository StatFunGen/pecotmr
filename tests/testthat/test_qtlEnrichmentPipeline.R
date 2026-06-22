context("qtlEnrichmentPipeline")

# ===========================================================================
# Strategy: mock qtlEnrichment so the pipeline runs end-to-end on a
# small fixture, but the heavy mixture-of-enrichment estimator never fires.
# ===========================================================================

.qep_makeHandle <- function(snp_n = 6L, n_samples = 30L,
                            path = "/tmp/sketch.gds") {
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

.qep_makeFmEntry <- function(variant_ids = paste0("v", 1:5),
                              pip = seq(0.9, by = -0.15, length.out = 5L),
                              alpha = NULL) {
  if (is.null(alpha)) alpha <- matrix(1/length(variant_ids),
                                       nrow = 1, ncol = length(variant_ids))
  tl <- data.frame(variant_id = variant_ids, pip = pip,
                   stringsAsFactors = FALSE)
  fit <- list(alpha = alpha, pip = setNames(pip, variant_ids),
              V = 0.1)
  FineMappingEntry(variantIds = variant_ids,
                   susieFit = fit,
                   topLoci    = tl)
}

.qep_makeGwasFmr <- function(studies = "G1", n_blocks = 1L,
                              with_sketch = TRUE) {
  entries <- vector("list", n_blocks)
  studyVec <- character(0)
  methodVec <- character(0)
  for (k in seq_len(n_blocks)) {
    # Different variants per block to avoid duplication.
    ids <- paste0("v", (k - 1L) * 3L + (1:3))
    entries[[k]] <- .qep_makeFmEntry(variant_ids = ids,
                                      pip = c(0.5, 0.2, 0.1))
    studyVec <- c(studyVec, studies)
    methodVec <- c(methodVec, "susie")
  }
  GwasFineMappingResult(
    study  = studyVec,
    method = methodVec,
    entry  = entries,
    ldSketch = if (with_sketch) .qep_makeHandle() else NULL)
}

.qep_makeQtlFmr <- function(contexts = "c1", traits = "t1",
                             with_sketch = TRUE) {
  n <- length(contexts) * length(traits)
  studies <- rep("Q1", n)
  ctx <- rep(contexts, length.out = n)
  trs <- rep(traits, each = length(contexts))[seq_len(n)]
  methods <- rep("susie", n)
  entries <- replicate(n,
    .qep_makeFmEntry(variant_ids = paste0("v", 1:5)),
    simplify = FALSE)
  QtlFineMappingResult(
    study   = studies,
    context = ctx,
    trait   = trs,
    method  = methods,
    entry   = entries,
    ldSketch = if (with_sketch) .qep_makeHandle() else NULL)
}

# Mock that returns a plausible enrichment list.
.qep_mockEnrichment <- function(value = 1.5) {
  function(gwasPip, susieQtlRegions, ...) {
    list(enrichment = value,
         enrichmentSe = 0.1,
         enrichmentLogOdds = log(value))
  }
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("qtlEnrichmentPipeline: rejects non-GwasFineMappingResult gwasFmr", {
  qfmr <- .qep_makeQtlFmr()
  expect_error(
    qtlEnrichmentPipeline(gwasFineMappingResult = "no",
                          qtlFineMappingResult  = qfmr),
    "must be a GwasFineMappingResult"
  )
})

test_that("qtlEnrichmentPipeline: rejects non-QtlFineMappingResult qtlFmr", {
  gfmr <- .qep_makeGwasFmr()
  expect_error(
    qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                          qtlFineMappingResult  = "no"),
    "must be a QtlFineMappingResult"
  )
})

test_that("qtlEnrichmentPipeline: NULL ldSketch on the GWAS side errors", {
  gfmr <- .qep_makeGwasFmr(with_sketch = FALSE)
  qfmr <- .qep_makeQtlFmr()
  expect_error(
    qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                          qtlFineMappingResult  = qfmr),
    "must have a non-NULL ldSketch"
  )
})

test_that("qtlEnrichmentPipeline: ldSketch mismatch errors", {
  # Build the QTL with a sketch carrying a different sample set.
  gfmr <- .qep_makeGwasFmr()
  qSketch <- .qep_makeHandle()
  qSketch@sampleIds <- paste0("z", seq_len(qSketch@nSamples))
  qfmr <- QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(.qep_makeFmEntry()),
    ldSketch = qSketch)
  expect_error(
    qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                          qtlFineMappingResult  = qfmr),
    "different sample sets"
  )
})

# ===========================================================================
# Per-study / per-context iteration via mocked qtlEnrichment
# ===========================================================================

test_that("qtlEnrichmentPipeline: returns one row per (gwasStudy, qtlStudy, qtlContext) triple", {
  gfmr <- .qep_makeGwasFmr()
  qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
  local_mocked_bindings(qtlEnrichment = .qep_mockEnrichment(2.0),
                        .package = "pecotmr")
  out <- qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                                qtlFineMappingResult  = qfmr)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)  # 1 GWAS study * 1 QTL study * 2 contexts
  expect_setequal(out$gwasStudy,  "G1")
  expect_setequal(out$qtlStudy,   "Q1")
  expect_setequal(out$qtlContext, c("c1", "c2"))
  expect_equal(out$enrichment, c(2.0, 2.0))
})

test_that("qtlEnrichmentPipeline: distinguishes two QTL studies that share a context label", {
  # Build a QtlFineMappingResult with two studies (Q1, Q2) both
  # tagging the same context "shared_ctx". Per-study filtering must
  # keep them separate (a context-only filter would merge them and
  # produce a single enrichment row instead of two).
  e1 <- .qep_makeFmEntry(variant_ids = paste0("v", 1:5))
  e2 <- .qep_makeFmEntry(variant_ids = paste0("v", 1:5))
  qfmr <- QtlFineMappingResult(
    study    = c("Q1", "Q2"),
    context  = c("shared_ctx", "shared_ctx"),
    trait    = c("t1", "t1"),
    method   = c("susie", "susie"),
    entry    = list(e1, e2),
    ldSketch = .qep_makeHandle())
  gfmr <- .qep_makeGwasFmr()
  capturedRegions <- list()
  local_mocked_bindings(
    qtlEnrichment = function(gwasPip, susieQtlRegions, ...) {
      capturedRegions[[length(capturedRegions) + 1L]] <<- susieQtlRegions
      list(enrichment = 2.0, enrichmentSe = 0.1, enrichmentLogOdds = log(2))
    },
    .package = "pecotmr")
  out <- qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                                qtlFineMappingResult  = qfmr)
  expect_equal(nrow(out), 2L)
  expect_setequal(out$qtlStudy, c("Q1", "Q2"))
  # Each call to qtlEnrichment sees exactly one region (the per-study one),
  # not both regions merged together.
  expect_true(all(lengths(capturedRegions) == 1L))
})

test_that("qtlEnrichmentPipeline: qtlEnrichment failure produces a warning + skip", {
  gfmr <- .qep_makeGwasFmr()
  qfmr <- .qep_makeQtlFmr()
  local_mocked_bindings(
    qtlEnrichment = function(...) stop("synthetic failure"),
    .package = "pecotmr")
  expect_warning(
    out <- qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                                  qtlFineMappingResult  = qfmr),
    "qtlEnrichment failed"
  )
  expect_equal(nrow(out), 0L)
})

test_that("qtlEnrichmentPipeline: empty input collections yield the empty schema", {
  # Build a GwasFineMappingResult whose entries have empty fits so the PIP
  # vector is empty.
  emptyEntry <- FineMappingEntry(
    variantIds = "v1",
    susieFit = list(),  # no pip -> .enrBuildGwasPipVector returns numeric(0)
    topLoci    = data.frame(variant_id = "v1", pip = 0.1,
                            stringsAsFactors = FALSE))
  gfmr <- GwasFineMappingResult(
    study  = "G1", method = "susie",
    entry  = list(emptyEntry),
    ldSketch = .qep_makeHandle())
  qfmr <- .qep_makeQtlFmr()
  expect_warning(
    out <- qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                                  qtlFineMappingResult  = qfmr),
    "no usable PIPs"
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_setequal(colnames(out),
                  c("gwasStudy", "qtlStudy", "qtlContext", "enrichment",
                    "enrichmentSe", "enrichmentLogOdds"))
})

# ===========================================================================
# Internal helpers: .enrBuildGwasPipVector + .enrBuildQtlRegionsList
# ===========================================================================

test_that(".enrBuildGwasPipVector: extracts pip per study", {
  gfmr <- .qep_makeGwasFmr()
  out <- pecotmr:::.enrBuildGwasPipVector(gfmr, "G1")
  expect_equal(length(out), 3L)
  expect_setequal(names(out), paste0("v", 1:3))
})

test_that(".enrBuildGwasPipVector: deduplicates identical PIPs across blocks", {
  # Two rows under the same study, same method, different LD blocks
  # (distinct region_ids auto-supplied by the constructor) — the
  # genome-wide multi-block shape. Both rows share v1 with the same
  # PIP, so dedup must collapse it.
  e1 <- .qep_makeFmEntry(variant_ids = c("v1", "v2"),
                          pip = c(0.5, 0.2))
  e2 <- .qep_makeFmEntry(variant_ids = c("v1", "v3"),
                          pip = c(0.5, 0.4))
  g <- GwasFineMappingResult(
    study = c("G1", "G1"), method = c("susie", "susie"),
    entry = list(e1, e2),
    ldSketch = .qep_makeHandle())
  out <- pecotmr:::.enrBuildGwasPipVector(g, "G1")
  expect_setequal(names(out), c("v1", "v2", "v3"))
})

test_that(".enrBuildGwasPipVector: conflicting PIPs across blocks errors", {
  # Same variant 'v1' appears in two LD blocks with different PIPs —
  # .enrBuildGwasPipVector must refuse to merge them silently.
  e1 <- .qep_makeFmEntry(variant_ids = c("v1"), pip = 0.5,
                          alpha = matrix(0.5, 1, 1))
  e2 <- .qep_makeFmEntry(variant_ids = c("v1"), pip = 0.8,
                          alpha = matrix(0.8, 1, 1))
  g <- GwasFineMappingResult(
    study = c("G1", "G1"), method = c("susie", "susie"),
    entry = list(e1, e2),
    ldSketch = .qep_makeHandle())
  expect_error(
    pecotmr:::.enrBuildGwasPipVector(g, "G1"),
    "conflicting PIPs"
  )
})

test_that(".enrBuildQtlRegionsList: returns per-entry fit shapes for a (study, context) hit", {
  qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
  out <- pecotmr:::.enrBuildQtlRegionsList(qfmr, "Q1", "c1")
  expect_equal(length(out), 1L)
  expect_true(!is.null(out[[1L]]$alpha))
  expect_true(!is.null(out[[1L]]$pip))
})

test_that(".enrBuildQtlRegionsList: returns empty list when the (study, context) tuple is absent", {
  qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
  # Correct context but wrong study -> no hit, even though context exists.
  expect_equal(length(pecotmr:::.enrBuildQtlRegionsList(qfmr, "Q_ghost", "c1")), 0L)
  # Correct study but wrong context.
  expect_equal(length(pecotmr:::.enrBuildQtlRegionsList(qfmr, "Q1", "c_ghost")), 0L)
})

# ===========================================================================
# qtlEnrichment() — kernel wrapper + real C++ integration
# These tests deliberately do NOT mock qtlEnrichmentRcpp so the C++
# kernel in src/qtl_enrichment.cpp gets coverage. The wrapper itself
# (R/qtlEnrichmentPipeline.R::qtlEnrichment) is exercised here directly
# rather than via the deprecated `computeQtlEnrichment` shim (which has
# skip_on_covr()).
# ===========================================================================

# Build a small (gwasPip, susieQtlRegions) fixture with a sparse causal
# signal at known indices so the C++ enrichment routine has something
# meaningful to compute.
.qep_makeRealKernelInputs <- function(seed = 42, nSnps = 50,
                                       causalIdx = c(5, 20, 35),
                                       causalPips = c(0.8, 0.6, 0.9),
                                       L = 2L) {
  set.seed(seed)
  variantNames <- paste0("1:", seq_len(nSnps), ":A:G")
  gwasPip <- rep(0.01, nSnps)
  gwasPip[causalIdx] <- causalPips
  names(gwasPip) <- variantNames

  alpha <- matrix(1 / nSnps, nrow = L, ncol = nSnps)
  alpha[1, ] <- 0.001; alpha[1, causalIdx[1]] <- 0.95
  alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
  alpha[2, ] <- 0.001; alpha[2, causalIdx[2]] <- 0.95
  alpha[2, ] <- alpha[2, ] / sum(alpha[2, ])
  pip <- colSums(alpha)
  names(pip) <- variantNames
  susieFits <- list(
    fit1 = list(pip = pip, alpha = alpha,
                prior_variance = c(0.5, 0.3)))
  list(gwasPip = gwasPip, susieQtlRegions = susieFits,
       variantNames = variantNames)
}

test_that("qtlEnrichment: real C++ kernel returns the expected keys (numGwas + piQtl supplied)", {
  fx <- .qep_makeRealKernelInputs()
  res <- qtlEnrichment(
    gwasPip = fx$gwasPip, susieQtlRegions = fx$susieQtlRegions,
    numGwas = 5000, piQtl = 0.5,
    lambda = 1, impN = 5, numThreads = 1, verbose = FALSE)
  expect_type(res, "list")
  en <- res[[1L]]
  expectedKeys <- c("Intercept", "Enrichment (no shrinkage)",
                    "Enrichment (w/ shrinkage)",
                    "sd (no shrinkage)", "sd (w/ shrinkage)",
                    "Alternative (coloc) p1", "Alternative (coloc) p2",
                    "Alternative (coloc) p12")
  expect_setequal(intersect(expectedKeys, names(en)), expectedKeys)
  expect_true(all(is.finite(unlist(en[expectedKeys]))))
})

test_that("qtlEnrichment: numGwas omitted -> estimates piGwas from data + warns", {
  fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
  expect_warning(
    res <- qtlEnrichment(
      gwasPip = fx$gwasPip, susieQtlRegions = fx$susieQtlRegions,
      piQtl = 0.5, impN = 5, numThreads = 1, verbose = FALSE),
    "numGwas is not provided")
  expect_type(res, "list")
})

test_that("qtlEnrichment: piQtl omitted -> estimates from susieQtlRegions + warns", {
  fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
  expect_warning(
    res <- qtlEnrichment(
      gwasPip = fx$gwasPip, susieQtlRegions = fx$susieQtlRegions,
      numGwas = 3000, impN = 5, numThreads = 1, verbose = FALSE),
    "piQtl is not provided")
  expect_type(res, "list")
})

test_that("qtlEnrichment: errors when piGwas resolves to zero", {
  fx <- .qep_makeRealKernelInputs()
  zeroGwas <- rep(0, length(fx$gwasPip))
  names(zeroGwas) <- names(fx$gwasPip)
  expect_error(
    qtlEnrichment(gwasPip = zeroGwas,
                  susieQtlRegions = fx$susieQtlRegions,
                  piQtl = 0.5, numThreads = 1, verbose = FALSE),
    "No association signal found")
})

test_that("qtlEnrichment: errors when piQtl resolves to zero", {
  fx <- .qep_makeRealKernelInputs()
  expect_error(
    qtlEnrichment(gwasPip = fx$gwasPip,
                  susieQtlRegions = fx$susieQtlRegions,
                  numGwas = 5000, piQtl = 0,
                  numThreads = 1, verbose = FALSE),
    "No QTL associated")
})

test_that("qtlEnrichment: errors when gwasPip has no names", {
  fx <- .qep_makeRealKernelInputs()
  unnamed <- unname(fx$gwasPip)
  expect_error(
    qtlEnrichment(gwasPip = unnamed,
                  susieQtlRegions = fx$susieQtlRegions,
                  numGwas = 5000, piQtl = 0.5,
                  numThreads = 1, verbose = FALSE),
    "Variant names are missing in gwasPip")
})

test_that("qtlEnrichment: errors when susieQtlRegions$pip lacks names", {
  fx <- .qep_makeRealKernelInputs()
  fx$susieQtlRegions$fit1$pip <- unname(fx$susieQtlRegions$fit1$pip)
  expect_error(
    qtlEnrichment(gwasPip = fx$gwasPip,
                  susieQtlRegions = fx$susieQtlRegions,
                  numGwas = 5000, piQtl = 0.5,
                  numThreads = 1, verbose = FALSE),
    "Variant names are missing in susieQtlRegions")
})

test_that("qtlEnrichment: tracks unmatched QTL variants in the output", {
  fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
  # Inject a couple of variant IDs into the QTL fit that don't exist
  # in the GWAS PIP vector.
  newNames <- names(fx$susieQtlRegions$fit1$pip)
  newNames[1:2] <- c("1:9999:A:G", "1:9998:A:G")
  names(fx$susieQtlRegions$fit1$pip) <- newNames
  colnames(fx$susieQtlRegions$fit1$alpha) <- newNames
  res <- qtlEnrichment(
    gwasPip = fx$gwasPip, susieQtlRegions = fx$susieQtlRegions,
    numGwas = 3000, piQtl = 0.5,
    impN = 5, numThreads = 1, verbose = FALSE)
  expect_true("unused_xqtl_variants" %in% names(res))
})


