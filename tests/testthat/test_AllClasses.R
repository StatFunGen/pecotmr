context("AllClasses (virtual base classes)")

# Most slots / accessors on the concrete subclasses are exercised in their
# own test files; these tests target the *base-class* behaviors that the
# concrete subclasses inherit without overriding (getStudy on SumStatsBase,
# getQcDiagnostics body branches, and the zero-row adjustPips short-circuit
# on FineMappingResultBase).

# ===========================================================================
# Helpers
# ===========================================================================

.alc_makeHandle <- function(snp_n = 3L) {
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
    nSamples = 10L,
    sampleIds = paste0("s", seq_len(10L)),
    pgenPtr = NULL)
}

.alc_makeGr <- function(n = 3) {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = seq(100L, by = 100L, length.out = n), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("rs", seq_len(n)),
    A1 = rep("A", n), A2 = rep("G", n),
    Z = rnorm(n), N = rep(1000L, n))
  gr
}

.alc_makeGwasSumStats <- function(qcInfo = list()) {
  GwasSumStats(
    study    = "g1",
    entry    = list(.alc_makeGr()),
    genome   = "hg19",
    ldSketch = .alc_makeHandle(),
    qcInfo   = qcInfo)
}

.alc_makeFmEntry <- function(n = 3) {
  tl <- data.frame(
    variant_id     = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    chrom          = rep("1", n),
    pos            = as.integer(100 * seq_len(n)),
    A1             = rep("G", n),
    A2             = rep("A", n),
    N              = rep(1000, n),
    MAF            = rep(0.1, n),
    marginal_beta  = rep(0.1, n),
    marginal_se    = rep(0.05, n),
    marginal_z     = rep(2.0, n),
    marginal_p     = rep(0.05, n),
    pip            = seq(0.9, by = -0.1, length.out = n),
    posterior_mean = rep(0.05, n),
    posterior_sd   = rep(0.02, n),
    stringsAsFactors = FALSE)
  FineMappingEntry(
    variantIds = tl$variant_id,
    susieFit   = list(),
    topLoci    = tl)
}

# ===========================================================================
# SumStatsBase: getStudy (inherited by QtlSumStats / GwasSumStats)
# ===========================================================================

test_that("SumStatsBase: getStudy on a GwasSumStats returns unique study names", {
  ss <- .alc_makeGwasSumStats()
  expect_equal(getStudy(ss), "g1")
})

# ===========================================================================
# SumStatsBase: getQcDiagnostics — every branch
# ===========================================================================

test_that("SumStatsBase: getQcDiagnostics returns NULL on empty qcInfo", {
  ss <- .alc_makeGwasSumStats()  # qcInfo = list() by default
  expect_null(getQcDiagnostics(ss))
})

test_that("SumStatsBase: getQcDiagnostics returns NULL when entryAudit slot is absent", {
  # qcInfo has steps but no entryAudit -> nothing to return.
  ss <- .alc_makeGwasSumStats(qcInfo = list(step1 = "ok"))
  expect_null(getQcDiagnostics(ss))
})

test_that("SumStatsBase: getQcDiagnostics returns the per-entry diagnostics by index", {
  diag1 <- data.frame(SNP = "rs1", outlier = FALSE, stringsAsFactors = FALSE)
  diag2 <- data.frame(SNP = "rs2", outlier = TRUE,  stringsAsFactors = FALSE)
  qc <- list(entryAudit = list(
    list(ldMismatchDiagnostics = diag1),
    list(ldMismatchDiagnostics = diag2)))
  ss <- .alc_makeGwasSumStats(qcInfo = qc)
  expect_identical(getQcDiagnostics(ss, entry = 1L), diag1)
  expect_identical(getQcDiagnostics(ss, entry = 2L), diag2)
})

test_that("SumStatsBase: getQcDiagnostics(entry = NULL) returns the populated entries only", {
  diag1 <- data.frame(SNP = "rs1", outlier = FALSE, stringsAsFactors = FALSE)
  # Entry 2's audit has no ldMismatchDiagnostics field; should be filtered.
  qc <- list(entryAudit = list(
    list(ldMismatchDiagnostics = diag1),
    list(other = "no diagnostics here")))
  ss <- .alc_makeGwasSumStats(qcInfo = qc)
  out <- getQcDiagnostics(ss, entry = NULL)
  expect_type(out, "list")
  expect_equal(length(out), 1L)
  expect_named(out, "1")
  expect_identical(out[["1"]], diag1)
})

test_that("SumStatsBase: getQcDiagnostics(entry = NULL) returns NULL when no entry has diagnostics", {
  qc <- list(entryAudit = list(list(other = 1), list(other = 2)))
  ss <- .alc_makeGwasSumStats(qcInfo = qc)
  expect_null(getQcDiagnostics(ss, entry = NULL))
})

test_that("SumStatsBase: getQcDiagnostics errors on out-of-range entry", {
  qc <- list(entryAudit = list(list(ldMismatchDiagnostics = data.frame(z = 1))))
  ss <- .alc_makeGwasSumStats(qcInfo = qc)
  expect_error(getQcDiagnostics(ss, entry = 0L),  "must be a single integer")
  expect_error(getQcDiagnostics(ss, entry = 99L), "must be a single integer")
  expect_error(getQcDiagnostics(ss, entry = c(1L, 2L)),
               "must be a single integer")
  expect_error(getQcDiagnostics(ss, entry = "first"),
               "must be a single integer")
})

# ===========================================================================
# FineMappingResultBase: adjustPips zero-row short-circuit
# ===========================================================================

test_that("FineMappingResultBase: adjustPips on a zero-row collection returns the input unchanged", {
  e <- .alc_makeFmEntry(3)
  res <- GwasFineMappingResult(study = "g1", method = "susie",
                                entry = list(e))
  empty <- res[integer(0), ]
  expect_s4_class(empty, "GwasFineMappingResult")
  expect_equal(nrow(empty), 0L)
  # Should hit the `if (nrow(x) == 0L) return(x)` early-return.
  out <- adjustPips(empty, character(0))
  expect_identical(out, empty)
})
