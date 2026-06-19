context("qtlSumStats")

# ===========================================================================
# Test helpers
# ===========================================================================

.qtlMakeGenotypeHandle <- function() {
  new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snpInfo = data.frame(),
    nSamples = 0L,
    sampleIds = character(),
    pgenPtr = NULL)
}

.qtlMakeEntryGr <- function(n = 5, chr = "chr1",
                            start_at = 100L, step = 100L,
                            with_maf = FALSE) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep(chr, n),
    ranges = IRanges::IRanges(
      start = seq(start_at, by = step, length.out = n),
      width = 1L)
  )
  mcols_list <- list(
    SNP = paste0("rs", seq_len(n)),
    A1  = rep("A", n),
    A2  = rep("G", n),
    Z   = seq(1.0, by = 0.5, length.out = n),
    N   = rep(1000L, n)
  )
  if (with_maf)
    mcols_list$MAF <- seq(0.05, by = 0.01, length.out = n)
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcols_list)
  gr
}

.qtlMakeOne <- function(study = "study1", context = "Cortex",
                        trait = "ENSG001", n = 5, varY = NA_real_,
                        with_maf = FALSE, qcInfo = list(),
                        genome = "hg19") {
  QtlSumStats(
    study    = study,
    context  = context,
    trait    = trait,
    entry    = list(.qtlMakeEntryGr(n, with_maf = with_maf)),
    genome   = genome,
    ldSketch = .qtlMakeGenotypeHandle(),
    varY     = varY,
    qcInfo   = qcInfo)
}

# ===========================================================================
# Constructor / validity
# ===========================================================================

test_that("QtlSumStats: minimal single-tuple object builds and validates", {
  obj <- .qtlMakeOne()
  expect_s4_class(obj, "QtlSumStats")
  expect_true(methods::validObject(obj))
  expect_equal(nrow(obj), 1L)
})

test_that("QtlSumStats: multi-tuple object is keyed by (study, context, trait)", {
  obj <- QtlSumStats(
    study   = c("s1", "s1", "s2"),
    context = c("c1", "c2", "c1"),
    trait   = c("t1", "t1", "t1"),
    entry   = list(.qtlMakeEntryGr(3),
                   .qtlMakeEntryGr(4),
                   .qtlMakeEntryGr(5)),
    genome   = "hg38",
    ldSketch = .qtlMakeGenotypeHandle())
  expect_equal(nrow(obj), 3L)
  expect_setequal(getContexts(obj), c("c1", "c2"))
  expect_equal(getTraits(obj), "t1")
})

test_that("QtlSumStats: errors when required args are missing", {
  expect_error(QtlSumStats(study = "s1", context = "c1"),
               "all required")
})

test_that("QtlSumStats: errors on length mismatch among study/context/trait/entry", {
  expect_error(
    QtlSumStats(
      study = c("s1", "s2"),
      context = c("c1"),
      trait = c("t1"),
      entry = list(.qtlMakeEntryGr(1)),
      genome = "hg19",
      ldSketch = .qtlMakeGenotypeHandle()),
    "same length"
  )
})

test_that("QtlSumStats: errors when entry is not a list", {
  expect_error(
    QtlSumStats(
      study = "s1", context = "c1", trait = "t1",
      entry = "not_a_list",
      genome = "hg19", ldSketch = .qtlMakeGenotypeHandle()),
    "must be a list"
  )
})

test_that("QtlSumStats: errors on non-singleton genome", {
  expect_error(
    QtlSumStats(
      study = "s1", context = "c1", trait = "t1",
      entry = list(.qtlMakeEntryGr(1)),
      genome = c("hg19", "hg38"),
      ldSketch = .qtlMakeGenotypeHandle()),
    "single character string"
  )
})

test_that("QtlSumStats: errors on duplicate (study, context, trait) tuple", {
  expect_error(
    QtlSumStats(
      study   = c("s1", "s1"),
      context = c("c1", "c1"),
      trait   = c("t1", "t1"),
      entry   = list(.qtlMakeEntryGr(1), .qtlMakeEntryGr(1)),
      genome  = "hg19",
      ldSketch = .qtlMakeGenotypeHandle()),
    "uniqueness violated"
  )
})

test_that("QtlSumStats: errors when entry list contains non-GRanges", {
  expect_error(
    QtlSumStats(
      study = "s1", context = "c1", trait = "t1",
      entry = list("not_a_granges"),
      genome = "hg19",
      ldSketch = .qtlMakeGenotypeHandle()),
    "must be a GRanges"
  )
})

test_that("QtlSumStats: scalar varY is recycled across tuples", {
  obj <- QtlSumStats(
    study   = c("s1", "s2"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    entry   = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
    genome  = "hg19",
    ldSketch = .qtlMakeGenotypeHandle(),
    varY    = 1.5)
  expect_equal(obj$varY, c(1.5, 1.5))
})

test_that("QtlSumStats: errors when varY length is neither 1 nor n", {
  expect_error(
    QtlSumStats(
      study   = c("s1", "s2"),
      context = c("c1", "c1"),
      trait   = c("t1", "t1"),
      entry   = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
      genome  = "hg19",
      ldSketch = .qtlMakeGenotypeHandle(),
      varY    = c(1, 2, 3)),
    "length 1 or length\\(study\\)"
  )
})

test_that("QtlSumStats: accepts and stores extra per-tuple columns via ...", {
  obj <- QtlSumStats(
    study   = c("s1", "s2"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    entry   = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
    genome  = "hg19",
    ldSketch = .qtlMakeGenotypeHandle(),
    cohort  = c("UKB", "FinnGen"))
  expect_equal(as.character(obj$cohort), c("UKB", "FinnGen"))
})

# ===========================================================================
# Accessors
# ===========================================================================

test_that("getSumStats / getZ / getN return values from the single-tuple entry", {
  obj <- .qtlMakeOne(n = 4)
  expect_s4_class(getSumStats(obj), "GRanges")
  expect_equal(length(getSumStats(obj)), 4L)
  expect_equal(getZ(obj), seq(1.0, by = 0.5, length.out = 4))
  expect_equal(getN(obj), rep(1000L, 4))
})

test_that("getMaf returns NULL when MAF column is absent, vector when present", {
  obj_no_maf <- .qtlMakeOne(with_maf = FALSE)
  expect_null(getMaf(obj_no_maf))

  obj_maf <- .qtlMakeOne(n = 4, with_maf = TRUE)
  expect_equal(length(getMaf(obj_maf)), 4L)
})

test_that("nSnps reports the number of variants in the selected entry", {
  obj <- .qtlMakeOne(n = 7)
  expect_equal(nSnps(obj), 7L)
})

test_that("getVarY returns numeric value or NULL when NA", {
  obj_na <- .qtlMakeOne(varY = NA_real_)
  expect_null(getVarY(obj_na))

  obj_v  <- .qtlMakeOne(varY = 2.5)
  expect_equal(getVarY(obj_v), 2.5)
})

test_that("getContexts / getTraits return unique values across tuples", {
  obj <- QtlSumStats(
    study   = c("s1", "s1", "s2"),
    context = c("c1", "c2", "c1"),
    trait   = c("t1", "t1", "t2"),
    entry   = list(.qtlMakeEntryGr(1), .qtlMakeEntryGr(1), .qtlMakeEntryGr(1)),
    genome  = "hg19",
    ldSketch = .qtlMakeGenotypeHandle())
  expect_setequal(getContexts(obj), c("c1", "c2"))
  expect_setequal(getTraits(obj), c("t1", "t2"))
})

test_that("accessors require (study, context, trait) when collection has >1 entry", {
  obj <- QtlSumStats(
    study   = c("s1", "s2"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    entry   = list(.qtlMakeEntryGr(1), .qtlMakeEntryGr(1)),
    genome  = "hg19",
    ldSketch = .qtlMakeGenotypeHandle())
  expect_error(getSumStats(obj),
               "Pass `study`, `context`, and `trait` to select one")

  expect_s4_class(
    getSumStats(obj, study = "s1", context = "c1", trait = "t1"),
    "GRanges")
})

test_that("accessors error on unknown tuple", {
  obj <- .qtlMakeOne(study = "s1", context = "c1", trait = "t1")
  expect_error(
    getSumStats(obj, study = "ghost", context = "c1", trait = "t1"),
    "No entry"
  )
})

test_that("accessors require length-1 selection args", {
  obj <- .qtlMakeOne()
  expect_error(
    getSumStats(obj, study = c("s1", "s2"), context = "c1", trait = "t1"),
    "must each be length 1"
  )
})

# ===========================================================================
# subsetChr
# ===========================================================================

test_that("subsetChr keeps only variants on the requested chromosome", {
  # Multi-chromosome entry: 3 on chr1, 2 on chr2. Built in one shot so
  # the resulting GRanges already shares seqlevels across rows.
  gr <- GenomicRanges::GRanges(
    seqnames = c(rep("chr1", 3), rep("chr2", 2)),
    ranges = IRanges::IRanges(
      start = c(100L, 200L, 300L, 1000L, 1100L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("rs", seq_along(gr)),
    A1  = rep("A", length(gr)),
    A2  = rep("G", length(gr)),
    Z   = seq(1.0, by = 0.5, length.out = length(gr)),
    N   = rep(1000L, length(gr)))
  obj <- QtlSumStats(
    study = "s1", context = "c1", trait = "t1",
    entry = list(gr),
    genome = "hg19",
    ldSketch = .qtlMakeGenotypeHandle())

  chr1_only <- subsetChr(obj, "1")
  expect_equal(nSnps(chr1_only), 3L)

  # The "chr" prefix is accepted equivalently.
  chr2_only <- subsetChr(obj, "chr2")
  expect_equal(nSnps(chr2_only), 2L)
})

test_that("subsetChr returns empty entry when no variants on chromosome", {
  obj <- .qtlMakeOne(n = 4)
  empty <- subsetChr(obj, "22")
  expect_equal(nSnps(empty), 0L)
})

test_that("subsetChr preserves class-level slots (genome, ldSketch, qcInfo)", {
  obj <- .qtlMakeOne(qcInfo = list(step1 = "ok"))
  res <- subsetChr(obj, "1")
  expect_equal(getGenome(res), getGenome(obj))
  expect_equal(getQcInfo(res), getQcInfo(obj))
  expect_identical(getLdSketch(res), getLdSketch(obj))
})

# ===========================================================================
# Show method
# ===========================================================================

test_that("show prints entry count and genome build", {
  obj <- .qtlMakeOne()
  expect_output(show(obj),
                "QtlSumStats: 1 entries, genome build hg19")
})
