# Tests for writeSumstatsVcf (vcfWriter.R)

# =============================================================================
# Test data helpers
# =============================================================================

make_test_genotype_handle <- function() {
  new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snpInfo = data.frame(),
    nSamples = 0L,
    sampleIds = character(),
    pgenPtr = NULL)
}

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
  GwasSumStats(
    study = "test_trait",
    entry = list(gr),
    genome = "hg38",
    ldSketch = make_test_genotype_handle())
}

make_test_finemapping_result <- function(n = 5) {
  beta <- seq(0.5, by = -0.1, length.out = n)
  se   <- rep(0.1, n)
  zv   <- seq(5.0, by = -1.0, length.out = n)
  tl <- data.frame(
    variant_id     = paste0("chr1:", seq(100, by = 100, length.out = n), ":T:A"),
    chrom          = rep("1", n),
    pos            = as.integer(seq(100, by = 100, length.out = n)),
    A1             = rep("A", n),
    A2             = rep("T", n),
    N              = rep(1000, n),
    MAF            = rep(0.1, n),
    marginal_beta  = beta,
    marginal_se    = se,
    marginal_z     = zv,
    marginal_p     = 2 * pnorm(-abs(zv)),
    pip            = seq(0.9, by = -0.1, length.out = n),
    posterior_mean = beta * 0.5,
    posterior_sd   = se * 0.5,
    cs_95          = paste0("susie_", c(1L, 1L, 0L, 2L, 0L)[seq_len(n)]),
    method         = rep("susie", n),
    stringsAsFactors = FALSE
  )
  entry <- FineMappingEntry(
    variantIds = tl$variant_id,
    susieFit   = list(),
    topLoci    = tl)
  GwasFineMappingResult(
    study  = "test_study",
    method = "susie",
    entry  = list(entry))
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
    variant_id     = character(0),
    chrom          = character(0),
    pos            = integer(0),
    A1             = character(0),
    A2             = character(0),
    N              = numeric(0),
    MAF            = numeric(0),
    marginal_beta  = numeric(0),
    marginal_se    = numeric(0),
    marginal_z     = numeric(0),
    marginal_p     = numeric(0),
    pip            = numeric(0),
    posterior_mean = numeric(0),
    posterior_sd   = numeric(0),
    stringsAsFactors = FALSE
  )
  entry <- FineMappingEntry(
    variantIds = character(0),
    susieFit   = list(),
    topLoci    = empty_tl)
  fm_empty <- GwasFineMappingResult(
    study  = "test_study",
    method = "susie",
    entry  = list(entry))

  out <- tempfile(fileext = ".vcf")
  on.exit(unlink(out), add = TRUE)

  expect_error(writeSumstatsVcf(fm_empty, out), "no variants to write")
})

# =============================================================================
# splitByContext / splitByTrait — one file per tuple
# =============================================================================

.make_multi_tuple_qtl_fmr <- function() {
  contexts <- c("brain", "blood")
  traits   <- c("ENSG_A", "ENSG_B")
  entries  <- lapply(seq_along(contexts), function(i) {
    ids <- paste0("chr1:", 100 * (1:3), ":T:A")
    tl <- data.frame(
      variant_id     = ids,
      chrom          = rep("1", 3),
      pos            = c(100L, 200L, 300L),
      A1             = rep("A", 3),
      A2             = rep("T", 3),
      N              = rep(1000, 3),
      MAF            = rep(0.1, 3),
      marginal_beta  = c(0.3, 0.1, -0.2) + i / 100,
      marginal_se    = rep(0.05, 3),
      marginal_z     = c(6.0, 2.0, -4.0),
      marginal_p     = c(1e-9, 0.045, 6e-5),
      pip            = c(0.9, 0.5, 0.7),
      posterior_mean = rep(0.05, 3),
      posterior_sd   = rep(0.02, 3),
      cs_95          = paste0("susie_", c(1L, 1L, 0L)),
      stringsAsFactors = FALSE)
    FineMappingEntry(variantIds = ids, susieFit = list(), topLoci = tl)
  })
  QtlFineMappingResult(
    study   = rep("study1", 2),
    context = contexts,
    trait   = traits,
    method  = rep("susie", 2),
    entry   = entries)
}

test_that("writeSumstatsVcf(FineMappingResult): splitByContext emits one VCF per context", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")
  fmr <- .make_multi_tuple_qtl_fmr()
  baseOut <- tempfile(fileext = ".vcf")
  on.exit(unlink(list.files(dirname(baseOut),
                            pattern = basename(tools::file_path_sans_ext(baseOut)),
                            full.names = TRUE)), add = TRUE)
  paths <- writeSumstatsVcf(fmr, baseOut, splitByContext = TRUE)
  expect_length(paths, 2L)
  # Each path is decorated with the context tag.
  expect_true(any(grepl("\\.brain\\.vcf$", paths)))
  expect_true(any(grepl("\\.blood\\.vcf$", paths)))
  for (p in paths) expect_true(file.exists(p))
})

test_that("writeSumstatsVcf(FineMappingResult): splitByTrait emits one VCF per trait", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")
  fmr <- .make_multi_tuple_qtl_fmr()
  baseOut <- tempfile(fileext = ".vcf")
  on.exit(unlink(list.files(dirname(baseOut),
                            pattern = basename(tools::file_path_sans_ext(baseOut)),
                            full.names = TRUE)), add = TRUE)
  paths <- writeSumstatsVcf(fmr, baseOut, splitByTrait = TRUE)
  expect_length(paths, 2L)
  expect_true(any(grepl("\\.ENSG_A\\.vcf$", paths)))
  expect_true(any(grepl("\\.ENSG_B\\.vcf$", paths)))
})

test_that("writeSumstatsVcf(FineMappingResult): splitByContext + splitByTrait combines tags", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")
  fmr <- .make_multi_tuple_qtl_fmr()
  baseOut <- tempfile(fileext = ".vcf")
  on.exit(unlink(list.files(dirname(baseOut),
                            pattern = basename(tools::file_path_sans_ext(baseOut)),
                            full.names = TRUE)), add = TRUE)
  paths <- writeSumstatsVcf(fmr, baseOut,
                            splitByContext = TRUE, splitByTrait = TRUE)
  expect_length(paths, 2L)
  expect_true(any(grepl("\\.brain\\.ENSG_A\\.vcf$", paths)))
  expect_true(any(grepl("\\.blood\\.ENSG_B\\.vcf$", paths)))
})

test_that("writeSumstatsVcf(FineMappingResult): multi-row without split flags requires selectors", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Biostrings")
  fmr <- .make_multi_tuple_qtl_fmr()
  out <- tempfile(fileext = ".vcf")
  expect_error(writeSumstatsVcf(fmr, out),
               "2 matching rows")
})
