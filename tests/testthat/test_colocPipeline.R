context("colocPipeline")

# ===========================================================================
# Strategy: mock coloc::coloc.bf_bf to return a tiny fake summary, then
# drive the QTL/GWAS pairing loop on a small fixture. Mock
# fineMappingPipeline so the GwasSumStats input path also runs without the
# heavy susie fits.
# ===========================================================================

.cp_makeHandle <- function(snp_n = 6L, n_samples = 30L,
                           sample_prefix = "s") {
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
    sampleIds = paste0(sample_prefix, seq_len(n_samples)),
    pgenPtr = NULL)
}

.cp_makeFmEntry <- function(variant_ids = paste0("chr1:", 100 * (1:5), ":A:G"),
                             withLbf = TRUE, n_eff = 2L) {
  pip <- seq(0.9, by = -0.15, length.out = length(variant_ids))
  tl <- data.frame(variant_id = variant_ids, pip = pip,
                   stringsAsFactors = FALSE)
  fit <- list(
    alpha = matrix(1/length(variant_ids),
                   nrow = n_eff, ncol = length(variant_ids),
                   dimnames = list(NULL, variant_ids)),
    pip   = setNames(pip, variant_ids),
    V     = rep(0.05, n_eff))
  if (withLbf)
    fit$lbf_variable <- matrix(rnorm(n_eff * length(variant_ids)),
                               nrow = n_eff, ncol = length(variant_ids),
                               dimnames = list(NULL, variant_ids))
  FineMappingEntry(variantIds = variant_ids,
                   trimmedFit = fit,
                   topLoci    = tl)
}

.cp_makeQtlFmr <- function(tuples = list(c("Q1", "c1", "t1", "susie")),
                            entries = NULL, with_sketch = TRUE) {
  if (is.null(entries))
    entries <- replicate(length(tuples), .cp_makeFmEntry(), simplify = FALSE)
  QtlFineMappingResult(
    study   = vapply(tuples, `[[`, character(1), 1),
    context = vapply(tuples, `[[`, character(1), 2),
    trait   = vapply(tuples, `[[`, character(1), 3),
    method  = vapply(tuples, `[[`, character(1), 4),
    entry   = entries,
    ldSketch = if (with_sketch) .cp_makeHandle() else NULL)
}

.cp_makeGwasFmr <- function(tuples = list(c("G1", "susie")),
                             entries = NULL, with_sketch = TRUE) {
  if (is.null(entries))
    entries <- replicate(length(tuples), .cp_makeFmEntry(), simplify = FALSE)
  GwasFineMappingResult(
    study  = vapply(tuples, `[[`, character(1), 1),
    method = vapply(tuples, `[[`, character(1), 2),
    entry  = entries,
    ldSketch = if (with_sketch) .cp_makeHandle() else NULL)
}

.cp_makeGwasSumstats <- function(study = "G1", qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = rnorm(5), N = rep(1000L, 5))
  GwasSumStats(
    study    = study,
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .cp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.cp_mockColocBfBf <- function() {
  function(qLbf, gLbf, p1, p2, p12, ...) {
    list(summary = data.frame(
      idx1 = 1L, idx2 = 1L, nSnps = ncol(qLbf),
      PP.H0.abf = 0.1, PP.H1.abf = 0.2, PP.H2.abf = 0.2,
      PP.H3.abf = 0.2, PP.H4.abf = 0.3,
      stringsAsFactors = FALSE))
  }
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("colocPipeline: rejects non-QtlFineMappingResult qtlFmr", {
  expect_error(
    colocPipeline(qtlFineMappingResult = "no",
                  gwasInput            = .cp_makeGwasFmr()),
    "must be a QtlFineMappingResult"
  )
})

test_that("colocPipeline: rejects gwasInput that is neither GwasSumStats nor GwasFineMappingResult", {
  expect_error(
    colocPipeline(qtlFineMappingResult = .cp_makeQtlFmr(),
                  gwasInput            = 42L),
    "must be a GwasSumStats or a GwasFineMappingResult"
  )
})

test_that("colocPipeline: rejects un-QCd GwasSumStats input", {
  qfmr <- .cp_makeQtlFmr()
  gss <- .cp_makeGwasSumstats(qc = FALSE)
  expect_error(
    colocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = gss),
    "has no QC record"
  )
})

# ===========================================================================
# .colocRequireMatchingLdSketches
# ===========================================================================

test_that(".colocRequireMatchingLdSketches: NULL qtl-side ldSketch is allowed", {
  qfmr <- .cp_makeQtlFmr(with_sketch = FALSE)
  gfmr <- .cp_makeGwasFmr()
  local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(colocPipeline(qtlFineMappingResult = qfmr,
                                         gwasInput            = gfmr))
  expect_s3_class(out, "data.frame")
})

test_that(".colocRequireMatchingLdSketches: non-NULL qtl + NULL gwas errors", {
  qfmr <- .cp_makeQtlFmr()
  gfmr <- .cp_makeGwasFmr(with_sketch = FALSE)
  expect_error(
    colocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = gfmr),
    "ldSketch is NULL"
  )
})

test_that(".colocRequireMatchingLdSketches: panel size mismatch errors", {
  qfmr <- .cp_makeQtlFmr()
  bigSketch <- .cp_makeHandle(snp_n = 7L)
  gfmr <- GwasFineMappingResult(
    study  = "G1", method = "susie",
    entry  = list(.cp_makeFmEntry()),
    ldSketch = bigSketch)
  expect_error(
    colocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = gfmr),
    "ldSketch panels differ in size"
  )
})

# ===========================================================================
# End-to-end with mocked coloc.bf_bf
# ===========================================================================

test_that("colocPipeline: returns one row per (QTL tuple, GWAS tuple) pair", {
  qfmr <- .cp_makeQtlFmr(tuples = list(c("Q1", "c1", "t1", "susie"),
                                       c("Q1", "c2", "t1", "susie")))
  gfmr <- .cp_makeGwasFmr(tuples = list(c("G1", "susie"),
                                        c("G2", "susie")))
  local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(colocPipeline(qtlFineMappingResult = qfmr,
                                         gwasInput            = gfmr))
  expect_equal(nrow(out), 4L)  # 2 QTL tuples * 2 GWAS tuples
  expect_setequal(out$study,      "Q1")
  expect_setequal(out$context,    c("c1", "c2"))
  expect_setequal(out$gwasStudy,  c("G1", "G2"))
})

test_that("colocPipeline: resolves GwasSumStats via fineMappingPipeline (mocked)", {
  qfmr <- .cp_makeQtlFmr()
  gss <- .cp_makeGwasSumstats()
  local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
  local_mocked_bindings(
    fineMappingPipeline = function(data, methods, ...) .cp_makeGwasFmr(),
    .package = "pecotmr")
  out <- suppressWarnings(colocPipeline(qtlFineMappingResult = qfmr,
                                         gwasInput            = gss))
  expect_s3_class(out, "data.frame")
  expect_gte(nrow(out), 1L)
})

test_that("colocPipeline: returnGwasFineMapping attaches the resolved FMR", {
  qfmr <- .cp_makeQtlFmr()
  gss <- .cp_makeGwasSumstats()
  resolvedGfmr <- .cp_makeGwasFmr()
  local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
  local_mocked_bindings(
    fineMappingPipeline = function(data, methods, ...) resolvedGfmr,
    .package = "pecotmr")
  out <- suppressWarnings(colocPipeline(qtlFineMappingResult = qfmr,
                                         gwasInput            = gss,
                                         returnGwasFineMapping = TRUE,
                                         adjustPips           = FALSE))
  expect_identical(attr(out, "gwasFineMapping"), resolvedGfmr)
})

test_that("colocPipeline: coloc.bf_bf failure surfaces as warning and skip", {
  qfmr <- .cp_makeQtlFmr()
  gfmr <- .cp_makeGwasFmr()
  local_mocked_bindings(
    coloc.bf_bf = function(...) stop("synthetic test failure"),
    .package = "coloc")
  expect_warning(
    out <- colocPipeline(qtlFineMappingResult = qfmr,
                         gwasInput            = gfmr),
    "coloc.bf_bf failed"
  )
  expect_equal(nrow(out), 0L)
})

test_that("colocPipeline: empty result has the documented schema", {
  qfmr <- .cp_makeQtlFmr()
  # Build a GWAS FMR whose entry has no usable LBF -> pre-extract returns empty.
  emptyFit <- list(alpha = matrix(0, 1, 1), pip = c(v1 = 0),
                   V = 0, lbf_variable = matrix(NA_real_, 1, 1))
  e <- FineMappingEntry(variantIds = "v1",
                        trimmedFit = emptyFit,
                        topLoci = data.frame(variant_id = "v1", pip = 0,
                                              stringsAsFactors = FALSE))
  gfmr <- GwasFineMappingResult(
    study = "G1", method = "susie",
    entry = list(e),
    ldSketch = .cp_makeHandle())
  out <- suppressWarnings(
    colocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = gfmr))
  expect_equal(nrow(out), 0L)
  expect_setequal(colnames(out),
                  c("study", "context", "trait", "method",
                    "gwasStudy", "gwasMethod", "idx1", "idx2", "nSnps",
                    "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                    "PP.H3.abf", "PP.H4.abf"))
})

# ===========================================================================
# Internal helpers
# ===========================================================================

test_that(".colocExtractLbfFromEntry: entry without trimmedFit returns NULL with warning", {
  e <- FineMappingEntry(
    variantIds = "v1",
    trimmedFit = NULL,
    topLoci    = data.frame(variant_id = "v1", pip = 0.1,
                            stringsAsFactors = FALSE))
  expect_warning(
    out <- pecotmr:::.colocExtractLbfFromEntry(e, FALSE, NULL, 1e-9),
    "has no trimmedFit"
  )
  expect_null(out)
})

test_that(".colocExtractLbfFromEntry: filterLbfCs subsets by cs_index", {
  fit <- list(
    alpha = matrix(0, 3, 4, dimnames = list(NULL, paste0("v", 1:4))),
    pip   = setNames(c(0.9, 0.1, 0.5, 0.2), paste0("v", 1:4)),
    V     = c(0.1, 0.1, 0.1),
    lbf_variable = matrix(1:12, 3, 4, dimnames = list(NULL, paste0("v", 1:4))),
    sets = list(cs_index = c(1L, 3L)))   # keep effects 1 and 3
  e <- FineMappingEntry(variantIds = paste0("v", 1:4),
                        trimmedFit = fit,
                        topLoci = data.frame(variant_id = paste0("v", 1:4),
                                              pip = c(0.9, 0.1, 0.5, 0.2),
                                              stringsAsFactors = FALSE))
  out <- pecotmr:::.colocExtractLbfFromEntry(e, filterLbfCs = TRUE,
                                              filterLbfCsSecondary = NULL,
                                              priorTol = 1e-9)
  expect_equal(nrow(out$lbf), 2L)
})

test_that(".colocAlignLbf: aligned matrices share the common variant set", {
  q <- matrix(0, 2, 4, dimnames = list(NULL,
              c("chr1:10:A:G", "chr1:20:A:G", "chr1:30:A:G", "chr1:40:A:G")))
  g <- matrix(0, 2, 3, dimnames = list(NULL,
              c("chr1:20:A:G", "chr1:30:A:G", "chr1:50:A:G")))
  aligned <- pecotmr:::.colocAlignLbf(q, g)
  expect_setequal(colnames(aligned$qtl),  c("chr1:20:A:G", "chr1:30:A:G"))
  expect_setequal(colnames(aligned$gwas), c("chr1:20:A:G", "chr1:30:A:G"))
})

test_that(".colocRbindLbf: rbinds matrices with union of columns, NAs filled with 0", {
  a <- matrix(1, 2, 2, dimnames = list(NULL, c("v1", "v2")))
  b <- matrix(2, 2, 2, dimnames = list(NULL, c("v2", "v3")))
  out <- pecotmr:::.colocRbindLbf(list(a, b))
  expect_equal(dim(out), c(4L, 3L))
  expect_setequal(colnames(out), c("v1", "v2", "v3"))
  # Column v3 padded with 0 for rows of `a`.
  expect_true(all(out[1:2, "v3"] == 0))
})

test_that(".colocStandardiseRow: fills in missing PP columns with NA", {
  row <- data.frame(idx1 = 1L, stringsAsFactors = FALSE)
  out <- pecotmr:::.colocStandardiseRow(row)
  expect_true(all(c("idx2", "nSnps", paste0("PP.H", 0:4, ".abf")) %in% colnames(out)))
})


context("encoloc")

# The file-path colocalization / enrichment wrappers (colocWrapper,
# xqtlEnrichmentWrapper, colocPostProcessor) and their helpers
# (filterAndOrderColocResults, calculateCumsum, calculate_purity,
# processColocResults, extract_ld_for_variants) have been removed in
# favor of the S4 colocPipeline / qtlEnrichmentPipeline / enlocPipeline
# entry points. The previous tests against the file-path wrappers no
# longer apply; they relied on mocking rssAnalysisPipeline (also
# removed). New tests for colocPipeline / qtlEnrichmentPipeline /
# enlocPipeline live alongside their pipeline implementations.

test_that("xqtlEnrichmentWrapper is a deprecated no-op", {
  expect_warning(
    res <- xqtlEnrichmentWrapper(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})

test_that("colocWrapper is a deprecated no-op", {
  expect_warning(
    res <- colocWrapper(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})

test_that("colocPostProcessor is a deprecated no-op", {
  expect_warning(
    res <- colocPostProcessor(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})


context("enlocPipeline")

# ===========================================================================
# Strategy: enlocPipeline now reuses the colocPipeline LBF + bf_bf path.
# Mock coloc::coloc.bf_bf so the pair loop runs end-to-end on a small
# fixture, with the enrichment-adjusted p12 visible through the captured
# call arguments.
# ===========================================================================

# Reuse fixture builders from test_colocPipeline.R style.
.ep_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
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

.ep_makeFmEntry <- function(variant_ids = paste0("chr1:", 100*(1:5), ":A:G"),
                             n_eff = 2L) {
  pip <- seq(0.9, by = -0.15, length.out = length(variant_ids))
  tl <- data.frame(variant_id = variant_ids, pip = pip,
                   stringsAsFactors = FALSE)
  set.seed(1)
  fit <- list(
    alpha = matrix(1/length(variant_ids),
                   nrow = n_eff, ncol = length(variant_ids),
                   dimnames = list(NULL, variant_ids)),
    pip   = setNames(pip, variant_ids),
    V     = rep(0.05, n_eff),
    lbf_variable = matrix(rnorm(n_eff * length(variant_ids)),
                          nrow = n_eff, ncol = length(variant_ids),
                          dimnames = list(NULL, variant_ids)))
  FineMappingEntry(variantIds = variant_ids,
                   trimmedFit = fit,
                   topLoci    = tl)
}

.ep_mockColocBfBf <- function() {
  function(qLbf, gLbf, p1, p2, p12, ...) {
    list(summary = data.frame(
      idx1 = 1L, idx2 = 1L, nSnps = ncol(qLbf),
      PP.H0.abf = 0.1, PP.H1.abf = 0.2, PP.H2.abf = 0.2,
      PP.H3.abf = 0.2, PP.H4.abf = 0.3,
      p12_actual = p12,
      stringsAsFactors = FALSE))
  }
}

.ep_makeQtlFmr <- function(with_sketch = TRUE) {
  QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(.ep_makeFmEntry()),
    ldSketch = if (with_sketch) .ep_makeHandle() else NULL)
}

.ep_makeGwasFmr <- function(with_sketch = TRUE) {
  GwasFineMappingResult(
    study  = "G1", method = "susie",
    entry  = list(.ep_makeFmEntry()),
    ldSketch = if (with_sketch) .ep_makeHandle() else NULL)
}

.ep_makeGwasSumstats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = rnorm(5), N = rep(1000L, 5))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .ep_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("enlocPipeline: rejects non-QtlFineMappingResult qtlFmr", {
  expect_error(
    enlocPipeline(qtlFineMappingResult = "no",
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment = data.frame(gwasStudy = "G1", qtlContext = "c1",
                                           enrichment = 2.0,
                                           stringsAsFactors = FALSE)),
    "must be a QtlFineMappingResult"
  )
})

test_that("enlocPipeline: rejects gwasInput that is neither GwasSumStats nor GwasFineMappingResult", {
  expect_error(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = 42L,
                  enrichment = data.frame(gwasStudy = "G1", qtlContext = "c1",
                                           enrichment = 2.0,
                                           stringsAsFactors = FALSE)),
    "must be a GwasSumStats or a GwasFineMappingResult"
  )
})

test_that("enlocPipeline: enrichment must be a data.frame", {
  expect_error(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = "not a df"),
    "must be a data.frame"
  )
})

test_that("enlocPipeline: enrichment missing required columns errors", {
  expect_error(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = data.frame(gwasStudy = "G1")),
    "missing column"
  )
})

test_that("enlocPipeline: un-QCd GwasSumStats input is rejected", {
  expect_error(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasSumstats(qc = FALSE),
                  enrichment = data.frame(gwasStudy = "G1", qtlContext = "c1",
                                           enrichment = 2.0,
                                           stringsAsFactors = FALSE)),
    "has no QC record"
  )
})

# ===========================================================================
# .colocLookupEnrichment (formerly .enlocLookupEnrichment, now shared)
# ===========================================================================

test_that(".colocLookupEnrichment: returns the value for a (gwasStudy, qtlContext) hit", {
  enr <- data.frame(gwasStudy = c("G1", "G2"),
                    qtlContext = c("c1", "c1"),
                    enrichment = c(2.0, 3.5),
                    stringsAsFactors = FALSE)
  expect_equal(pecotmr:::.colocLookupEnrichment(enr, "G2", "c1"), 3.5)
})

test_that(".colocLookupEnrichment: returns NA when no row matches", {
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  expect_true(is.na(pecotmr:::.colocLookupEnrichment(enr, "ghost", "c1")))
})

# ===========================================================================
# .colocEmptyResult(enriched = TRUE) — formerly .enlocEmptyResult
# ===========================================================================

test_that(".colocEmptyResult(enriched=TRUE): includes enrichment + p12Used schema", {
  out <- pecotmr:::.colocEmptyResult(enriched = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_true(all(c("enrichment", "p12Used") %in% colnames(out)))
})

# ===========================================================================
# Pair loop — runs end-to-end via the LBF + coloc.bf_bf path
# ===========================================================================

test_that("enlocPipeline: pair loop produces one row per (QTL tuple, GWAS tuple) with adjusted p12", {
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = enr,
                  p12                  = 5e-6,
                  p12Max               = 1e-3))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_equal(out$enrichment, 2.0)
  # 5e-6 * (1 + 2.0) = 1.5e-5 < 1e-3, so capped value is the raw product.
  expect_equal(out$p12Used, 1.5e-5)
  expect_equal(out$p12_actual, 1.5e-5)
})

test_that("enlocPipeline: missing-enrichment pair falls back to baseline p12 with a warning", {
  # An enrichment frame that has no row for (G1, c1).
  enr <- data.frame(gwasStudy = "G_other", qtlContext = "c_other",
                    enrichment = 10.0, stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  expect_warning(
    out <- enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                          gwasInput            = .ep_makeGwasFmr(),
                          enrichment           = enr,
                          p12                  = 5e-6,
                          p12Max               = 1e-3),
    "no enrichment entry"
  )
  # Baseline p12 unchanged because the pair fell back (enRow = 0).
  expect_equal(out$p12Used, 5e-6)
})

test_that("enlocPipeline: p12Max caps the adjusted prior", {
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 1e6,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = enr,
                  p12                  = 5e-6,
                  p12Max               = 1e-4))
  expect_equal(out$p12Used, 1e-4)
})

test_that("enlocPipeline: coloc.bf_bf failures are caught and warned, pair skipped", {
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(
    coloc.bf_bf = function(q, g, ...) stop("boom"),
    .package = "coloc")
  expect_warning(
    out <- enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                          gwasInput            = .ep_makeGwasFmr(),
                          enrichment           = enr),
    "coloc.bf_bf failed"
  )
  expect_equal(nrow(out), 0L)
})

test_that("enlocPipeline: returnGwasFineMapping=TRUE attaches gwasFineMapping attr (non-empty result)", {
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  # Use GwasFineMappingResult input: returnGwasFineMapping has no effect
  # for this branch (only GwasSumStats triggers attachment). To trigger
  # attachment we mock fineMappingPipeline so the GwasSumStats path
  # produces a usable FMR.
  fakeFmr <- .ep_makeGwasFmr()
  local_mocked_bindings(
    fineMappingPipeline = function(data, ...) fakeFmr,
    .package = "pecotmr")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult  = .ep_makeQtlFmr(),
                  gwasInput             = .ep_makeGwasSumstats(),
                  enrichment            = enr,
                  returnGwasFineMapping = TRUE))
  expect_true("gwasFineMapping" %in% names(attributes(out)))
  expect_s4_class(attr(out, "gwasFineMapping"), "GwasFineMappingResult")
})

test_that("enlocPipeline: qLbf NULL (QTL entry's LBF rows drop after priorTol) skips that QTL row", {
  # Build a QTL FMR whose lbf_variable is empty after the V > priorTol filter
  # (V = 0 < default priorTol 1e-9 -> drop all rows -> return NULL).
  emptyFit <- list(
    lbf_variable = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "v1")),
    V = 0.0)
  e <- FineMappingEntry(variantIds = "v1",
                        trimmedFit = emptyFit,
                        topLoci = data.frame(variant_id = "v1", pip = 0,
                                              stringsAsFactors = FALSE))
  qfmr <- QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(e),
    ldSketch = .ep_makeHandle())
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = data.frame(
                    gwasStudy = "G1", qtlContext = "c1",
                    enrichment = 1.0, stringsAsFactors = FALSE)))
  expect_equal(nrow(out), 0L)
})

test_that("enlocPipeline: aligned NULL (disjoint variant sets) skips that pair", {
  # QTL fmr with variant ids that don't overlap the GWAS variant ids.
  qVids <- paste0("chr1:", 100*(1:5), ":A:G")
  gVids <- paste0("chr2:", 200*(1:5), ":A:G")
  qfmr <- QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(.ep_makeFmEntry(variant_ids = qVids)),
    ldSketch = .ep_makeHandle())
  gfmr <- GwasFineMappingResult(
    study   = "G1", method = "susie",
    entry   = list(.ep_makeFmEntry(variant_ids = gVids)),
    ldSketch = .ep_makeHandle())
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = gfmr,
                  enrichment           = data.frame(
                    gwasStudy = "G1", qtlContext = "c1",
                    enrichment = 1.0, stringsAsFactors = FALSE)))
  expect_equal(nrow(out), 0L)
})

test_that("enlocPipeline: empty result schema includes enrichment + p12Used", {
  # Build a GWAS FMR whose entry has no usable LBF -> pre-extract returns empty.
  emptyFit <- list(alpha = matrix(0, 1, 1), pip = c(v1 = 0),
                   V = 0, lbf_variable = matrix(NA_real_, 1, 1))
  e <- FineMappingEntry(variantIds = "v1",
                        trimmedFit = emptyFit,
                        topLoci = data.frame(variant_id = "v1", pip = 0,
                                              stringsAsFactors = FALSE))
  gfmr <- GwasFineMappingResult(
    study = "G1", method = "susie",
    entry = list(e),
    ldSketch = .ep_makeHandle())
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = gfmr,
                  enrichment           = data.frame(
                    gwasStudy = "G1", qtlContext = "c1",
                    enrichment = 1.0, stringsAsFactors = FALSE)))
  expect_equal(nrow(out), 0L)
  expect_true(all(c("enrichment", "p12Used") %in% colnames(out)))
})
