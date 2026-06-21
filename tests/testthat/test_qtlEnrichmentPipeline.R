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


context("computeQtlEnrichment")

generate_mock_data <- function(seed=1, num_pips = 1000, num_susie_fits = 2) {
  # Simulate fake data for gwas_pip
  n_gwas_pip <- num_pips
  gwas_pip <- runif(n_gwas_pip)
  names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)
  gwas_fit <- list(pip=gwas_pip)

  # Simulate fake data for a single SuSiEFit object
  simulate_susiefit <- function(n, p) {
    pip <- runif(n)
    names(pip) <- paste0("snp", 1:n)
    alpha <- t(matrix(runif(n * p), nrow = n))
    alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
    list(
      pip = pip,
      alpha = alpha,
      prior_variance = runif(p)
    )
  }

  # Simulate multiple SuSiEFit objects
  n_susie_fits <- num_susie_fits
  susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 10), simplify = FALSE)
  # Add these fits to a list, providing names to each element
  names(susie_fits) <- paste0("fit", 1:length(susie_fits))
  return(list(gwas_fit=gwas_fit, susie_fits=susie_fits))
}

test_that("computeQtlEnrichment dummy data single-threaded works",{
  local_mocked_bindings(
      qtlEnrichmentRcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  expect_warning(
    computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, impN = 10, numThreads = 1),
    "numGwas is not provided. Estimating piGwas from the data. Note that this estimate may be biased if the input gwasPip does not contain genome-wide variants.")
  expect_warning(
    computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, impN = 10, numThreads = 1),
    "piQtl is not provided. Estimating piQtl from the data. Note that this estimate may be biased if either 1) the input susieQtlRegions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
  res <- expect_warning(computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, numGwas=5000, piQtl=0.49819, lambda = 1, impN = 10, numThreads = 1))
  expect_true(length(res) > 0)
})

test_that("computeQtlEnrichment dummy data single thread and multi-threaded are equivalent",{
  local_mocked_bindings(
      qtlEnrichmentRcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  res_single <- expect_warning(computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, numGwas=5000, piQtl=0.49819, lambda = 1, impN = 10, numThreads = 1))
  res_multi <- expect_warning(computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, numGwas=5000, piQtl=0.49819, lambda = 1, impN = 10, numThreads = 2))
  expect_equal(res_single, res_multi)
})

# ---- error paths (computeQtlEnrichment.R lines 86, 87, 91) ----
test_that("computeQtlEnrichment errors when pi_gwas is zero", {
  gwas_pip <- rep(0, 10)
  names(gwas_pip) <- paste0("snp", 1:10)
  susie_fits <- list(fit1 = list(pip = setNames(runif(10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    computeQtlEnrichment(gwas_pip, susie_fits, piQtl = 0.5),
    "No association signal found in GWAS data"
  )
})

test_that("computeQtlEnrichment errors when pi_qtl is zero", {
  gwas_pip <- runif(10)
  names(gwas_pip) <- paste0("snp", 1:10)
  susie_fits <- list(fit1 = list(pip = setNames(rep(0, 10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(computeQtlEnrichment(gwas_pip, susie_fits, numGwas = 1000, piQtl = 0)),
    "No QTL associated"
  )
})

test_that("computeQtlEnrichment errors when gwas_pip has no names", {
  gwas_pip <- runif(10)  # no names
  susie_fits <- list(fit1 = list(pip = setNames(runif(10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(computeQtlEnrichment(gwas_pip, susie_fits, numGwas = 1000, piQtl = 0.5)),
    "Variant names are missing in gwasPip"
  )
})

# ---- real C++ qtlEnrichmentRcpp integration test ----
test_that("computeQtlEnrichment calls real C++ enrichment code and returns expected keys", {
  set.seed(42)
  n_snps <- 50
  variantNames <- paste0("1:", 1:n_snps, ":A:G")

  # GWAS PIPs: sparse signal
  gwas_pip <- rep(0.01, n_snps)
  gwas_pip[c(5, 20, 35)] <- c(0.8, 0.6, 0.9)
  names(gwas_pip) <- variantNames

  # SuSiE fit with 2 single effects over same variants
  L <- 2
  alpha <- matrix(1 / n_snps, nrow = L, ncol = n_snps)
  # Concentrate probability on causal variants
  alpha[1, ] <- 0.001; alpha[1, 5] <- 0.95; alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
  alpha[2, ] <- 0.001; alpha[2, 20] <- 0.95; alpha[2, ] <- alpha[2, ] / sum(alpha[2, ])
  pip <- colSums(alpha)
  names(pip) <- variantNames

  susie_fits <- list(
    fit1 = list(pip = pip, alpha = alpha, prior_variance = c(0.5, 0.3))
  )

  # Call without mocking - exercises the real C++ code
  res <- suppressWarnings(
    computeQtlEnrichment(gwas_pip, susie_fits,
                           numGwas = 5000, piQtl = 0.5,
                           lambda = 1, impN = 5, numThreads = 1)
  )
  expect_type(res, "list")
  # The enrichment results are in res[[1]] (the C++ output list)
  en <- res[[1]]
  expected_keys <- c("Intercept", "Enrichment (no shrinkage)", "Enrichment (w/ shrinkage)",
                     "sd (no shrinkage)", "sd (w/ shrinkage)",
                     "Alternative (coloc) p1", "Alternative (coloc) p2", "Alternative (coloc) p12")
  for (key in expected_keys) {
    expect_true(key %in% names(en), info = paste("Missing key:", key))
  }
  # All numeric and finite
  numeric_vals <- unlist(en[expected_keys])
  expect_true(all(is.finite(numeric_vals)))
})

# ---- unmatched variants tracking (computeQtlEnrichment.R line 102) ----
test_that("computeQtlEnrichment tracks unmatched QTL variants", {
  local_mocked_bindings(
    qtlEnrichmentRcpp = function(...) TRUE
  )
  gwas_pip <- runif(10)
  names(gwas_pip) <- paste0("1:", 1:10, ":A:G")
  # QTL has some variants not in GWAS
  qtl_pip <- runif(5)
  names(qtl_pip) <- c(paste0("1:", 1:3, ":A:G"), "1:999:A:G", "1:998:A:G")
  susie_fits <- list(fit1 = list(pip = qtl_pip,
                                  alpha = matrix(runif(5), 1, 5),
                                  prior_variance = 1))
  res <- suppressWarnings(
    computeQtlEnrichment(gwas_pip, susie_fits, numGwas = 1000, piQtl = 0.5)
  )
  expect_true("unused_xqtl_variants" %in% names(res))
})