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
                   trimmedFit = fit,
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

test_that("qtlEnrichmentPipeline: returns one row per (gwasStudy, qtlContext) pair", {
  gfmr <- .qep_makeGwasFmr()
  qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
  local_mocked_bindings(qtlEnrichment = .qep_mockEnrichment(2.0),
                        .package = "pecotmr")
  out <- qtlEnrichmentPipeline(gwasFineMappingResult = gfmr,
                                qtlFineMappingResult  = qfmr)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)  # 1 GWAS study * 2 contexts
  expect_setequal(out$gwasStudy, "G1")
  expect_setequal(out$qtlContext, c("c1", "c2"))
  expect_equal(out$enrichment, c(2.0, 2.0))
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
    trimmedFit = list(),  # no pip -> .enrBuildGwasPipVector returns numeric(0)
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
                  c("gwasStudy", "qtlContext", "enrichment",
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
  # Two rows under the same study but different methods, sharing v1 with
  # an identical PIP value. The (study, method) validity constraint rules
  # out same-method repeats, so we use susie + susieInf for the second.
  e1 <- .qep_makeFmEntry(variant_ids = c("v1", "v2"),
                          pip = c(0.5, 0.2))
  e2 <- .qep_makeFmEntry(variant_ids = c("v1", "v3"),
                          pip = c(0.5, 0.4))
  g <- GwasFineMappingResult(
    study = c("G1", "G1"), method = c("susie", "susieInf"),
    entry = list(e1, e2),
    ldSketch = .qep_makeHandle())
  out <- pecotmr:::.enrBuildGwasPipVector(g, "G1")
  expect_setequal(names(out), c("v1", "v2", "v3"))
})

test_that(".enrBuildGwasPipVector: conflicting PIPs across blocks errors", {
  e1 <- .qep_makeFmEntry(variant_ids = c("v1"), pip = 0.5,
                          alpha = matrix(0.5, 1, 1))
  e2 <- .qep_makeFmEntry(variant_ids = c("v1"), pip = 0.8,
                          alpha = matrix(0.8, 1, 1))
  g <- GwasFineMappingResult(
    study = c("G1", "G1"), method = c("susie", "susieInf"),
    entry = list(e1, e2),
    ldSketch = .qep_makeHandle())
  expect_error(
    pecotmr:::.enrBuildGwasPipVector(g, "G1"),
    "conflicting PIPs"
  )
})

test_that(".enrBuildQtlRegionsList: returns per-entry fit shapes", {
  qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
  out <- pecotmr:::.enrBuildQtlRegionsList(qfmr, "c1")
  expect_equal(length(out), 1L)
  expect_true(!is.null(out[[1L]]$alpha))
  expect_true(!is.null(out[[1L]]$pip))
})


