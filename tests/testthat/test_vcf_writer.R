# Tests for writeSumstatsVcf (vcf_writer.R)

# =============================================================================
# Test data helpers
# =============================================================================

make_test_gwas_sumstats <- function(n = 5) {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = seq(100, by = 100, length.out = n), width = 1)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("rs", seq_len(n)),
    A1 = rep("A", n),
    A2 = rep("T", n),
    Z = seq(1.5, by = -0.5, length.out = n),
    N = rep(1000L, n)
  )
  new("GwasSumStats",
    sumstats = gr, genome = "hg38",
    traitName = "test_trait", varY = NULL)
}

make_test_finemapping_result <- function(n = 5) {
  tl <- data.frame(
    variant_id = paste0("chr1:", seq(100, by = 100, length.out = n), ":T:A"),
    method = rep("susie", n),
    pip = seq(0.9, by = -0.1, length.out = n),
    cs_index_95 = c(1L, 1L, 0L, 2L, 0L)[seq_len(n)],
    beta = seq(0.5, by = -0.1, length.out = n),
    se = rep(0.1, n),
    z = seq(5.0, by = -1.0, length.out = n),
    stringsAsFactors = FALSE
  )
  new("FineMappingResult",
    variantNames = tl$variant_id,
    trimmedFit = list(),
    topLoci = tl,
    method = "susie",
    sumstats = NULL)
}

# =============================================================================
# GwasSumStats to VCF
# =============================================================================

test_that("writeSumstatsVcf writes GwasSumStats to uncompressed VCF", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")

  ss <- make_test_gwas_sumstats(5)
  out <- tempfile(fileext = ".vcf")
  on.exit(unlink(out), add = TRUE)

  result <- writeSumstatsVcf(ss, out)
  expect_equal(result, out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})

# =============================================================================
# GwasSumStats to bgzipped VCF
# =============================================================================

test_that("writeSumstatsVcf writes GwasSumStats to bgzipped VCF", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")

  ss <- make_test_gwas_sumstats(5)
  out <- tempfile(fileext = ".vcf.bgz")
  on.exit(unlink(c(out, paste0(out, ".tbi")), force = TRUE), add = TRUE)

  result <- writeSumstatsVcf(ss, out)
  expect_equal(result, out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})

# =============================================================================
# FineMappingResult to VCF
# =============================================================================

test_that("writeSumstatsVcf writes FineMappingResult to uncompressed VCF", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")

  fm <- make_test_finemapping_result(5)
  out <- tempfile(fileext = ".vcf")
  on.exit(unlink(out), add = TRUE)

  result <- writeSumstatsVcf(fm, out)
  expect_equal(result, out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})

# =============================================================================
# FineMappingResult to BCF
# =============================================================================

test_that("writeSumstatsVcf writes FineMappingResult to BCF", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("Rsamtools")

  # asBcf may be temporarily disabled in some Rsamtools versions
  asbcf_works <- tryCatch({
    tmp_stem <- tempfile(fileext = ".vcf")
    tmp_bgz <- paste0(tmp_stem, ".bgz")
    gr <- GenomicRanges::GRanges("chr1",
      IRanges::IRanges(start = 100, width = 1, names = "v1"))
    hdr <- VariantAnnotation::VCFHeader(
      header = IRanges::DataFrameList(
        fileformat = S4Vectors::DataFrame(Value = "VCFv4.2",
                                          row.names = "fileformat")),
      sample = "probe")
    cd <- S4Vectors::DataFrame(Samples = "probe", row.names = "probe")
    v <- VariantAnnotation::VCF(rowRanges = gr, colData = cd,
                                exptData = list(header = hdr))
    VariantAnnotation::ref(v) <- Biostrings::DNAStringSet("A")
    VariantAnnotation::alt(v) <- Biostrings::DNAStringSetList(list("T"))
    VariantAnnotation::fixed(v)$FILTER <- "PASS"
    VariantAnnotation::writeVcf(v, tmp_stem, index = TRUE)
    bcf_stem <- tempfile()
    Rsamtools::asBcf(tmp_bgz, dictionary = "chr1", destination = bcf_stem)
    TRUE
  }, error = function(e) FALSE)
  skip_if(!asbcf_works, "Rsamtools::asBcf is not functional in this build")

  fm <- make_test_finemapping_result(5)
  out <- tempfile(fileext = ".bcf")
  on.exit(unlink(out), add = TRUE)

  result <- writeSumstatsVcf(fm, out)
  expect_equal(result, out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})

# =============================================================================
# Empty FineMappingResult errors
# =============================================================================

test_that("writeSumstatsVcf errors on empty FineMappingResult", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")

  empty_tl <- data.frame(
    variant_id = character(0),
    method = character(0),
    stringsAsFactors = FALSE
  )
  fm_empty <- new("FineMappingResult",
    variantNames = character(0),
    trimmedFit = list(),
    topLoci = empty_tl,
    method = "susie",
    sumstats = NULL)

  out <- tempfile(fileext = ".vcf")
  on.exit(unlink(out), add = TRUE)

  expect_error(writeSumstatsVcf(fm_empty, out), "no topLoci")
})
