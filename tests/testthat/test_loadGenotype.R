# Tests for genotype loading via readGenotypes + extractBlockGenotypes,
# and the loadGenotypeRegion dispatcher.

# Fixtures: 100 samples x 349 variants on chr21:17513228-17592874
test_data_dir <- test_path("test_data")
plink_prefix <- file.path(test_data_dir, "test_variants")
vcf_path <- file.path(test_data_dir, "test_variants.vcf.gz")
gds_path <- file.path(test_data_dir, "test_variants.gds")
region_all <- "chr21:17513228-17592874"
region_sub <- "chr21:17513228-17550000"

# Expected dimensions
n_samples <- 100L
n_variants <- 349L

test_that("format detection supports dotted PLINK2 prefixes", {
  tmp <- tempfile("plink2_dotted_prefix_")
  prefix <- file.path(dirname(tmp), "ADSP.R4.EUR.chr21")
  file.create(paste0(prefix, ".pgen"), paste0(prefix, ".pvar"), paste0(prefix, ".psam"))
  on.exit(unlink(paste0(prefix, c(".pgen", ".pvar", ".psam"))), add = TRUE)

  expect_equal(pecotmr:::.h2DetectFormat(prefix), "plink2")
})

# Shared helper: validate the output structure from loadGenotypeRegion
# (with returnVariantInfo=TRUE)
check_genotype_result <- function(result, expected_nrow = n_samples, expected_ncol = n_variants,
                                  label = "") {
  expect_true(is.list(result), label = paste(label, "is list"))
  expect_named(result, c("X", "variant_info"), ignore.order = TRUE)
  expect_true(is.matrix(result$X))
  expect_true(is.numeric(result$X))
  expect_equal(nrow(result$X), expected_nrow)
  expect_equal(ncol(result$X), expected_ncol)
  expect_true(is.data.frame(result$variant_info))
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result$variant_info)))
  expect_equal(nrow(result$variant_info), expected_ncol)
  # Column names of X match variant IDs
  expect_equal(colnames(result$X), result$variant_info$id)
  # Dosage values should be non-negative integers (0, 1, 2 for biallelic;
  # multiallelic VCF sites can have higher values)
  vals <- result$X[!is.na(result$X)]
  expect_true(all(vals >= 0), label = paste(label, "dosage non-negative"))
  expect_true(all(vals == round(vals)), label = paste(label, "dosage integer-valued"))
}

# --- readGenotypes: PLINK1 (snpStats) ----------------------------------------

test_that("readGenotypes creates plink1 handle", {
  skip_if_not_installed("snpStats")
  handle <- readGenotypes(plink_prefix, format = "plink1")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@format, "plink1")
  expect_equal(handle@nSamples, n_samples)
  expect_equal(nrow(handle@snpInfo), n_variants)
})

test_that("loadGenotypeRegion loads plink1 via dispatch", {
  skip_if_not_installed("snpStats")
  # Use .genotype suffix plink1 files tested elsewhere in test_file_utils
  plink1_path <- file.path(test_data_dir, "protocol_example.genotype")
  skip_if(!file.exists(paste0(plink1_path, ".bed")), "plink1 test fixture missing")
  result <- loadGenotypeRegion(plink1_path, returnVariantInfo = TRUE)
  expect_true(is.list(result))
  expect_true(is.matrix(result$X))
})

# --- readGenotypes: PLINK2 (pgenlibr) ----------------------------------------

test_that("readGenotypes creates plink2 handle", {
  skip_if_not_installed("pgenlibr")
  handle <- readGenotypes(plink_prefix, format = "plink2")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@format, "plink2")
  expect_equal(handle@nSamples, n_samples)
  expect_equal(nrow(handle@snpInfo), n_variants)
})

test_that("extractBlockGenotypes works for plink2", {
  skip_if_not_installed("pgenlibr")
  handle <- readGenotypes(plink_prefix, format = "plink2")
  rse <- extractBlockGenotypes(handle, seq_len(nrow(handle@snpInfo)))
  expect_s4_class(rse, "SummarizedExperiment")
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  expect_equal(nrow(dosage), n_variants)
  expect_equal(ncol(dosage), n_samples)
})

test_that("loadGenotypeRegion filters plink2 by region", {
  skip_if_not_installed("pgenlibr")
  result <- loadGenotypeRegion(plink_prefix, region = region_sub,
                                  returnVariantInfo = TRUE)
  check_genotype_result(result, expected_ncol = 134L, label = "plink2 region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("loadGenotypeRegion filters plink2 indels", {
  skip_if_not_installed("pgenlibr")
  result <- loadGenotypeRegion(plink_prefix, keepIndel = FALSE,
                                  returnVariantInfo = TRUE)
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("loadGenotypeRegion errors on empty region for plink2", {
  skip_if_not_installed("pgenlibr")
  expect_error(loadGenotypeRegion(plink_prefix, region = "chr1:1-2"))
})

# --- readGenotypes: VCF (VariantAnnotation) -----------------------------------

test_that("readGenotypes creates vcf handle", {
  skip_if_not_installed("VariantAnnotation")
  handle <- readGenotypes(vcf_path, format = "vcf")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@format, "vcf")
  expect_equal(handle@nSamples, n_samples)
  expect_equal(nrow(handle@snpInfo), n_variants)
})

test_that("loadGenotypeRegion loads VCF via dispatch", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(loadGenotypeRegion(vcf_path, returnVariantInfo = TRUE))
  check_genotype_result(result, label = "dispatch vcf")
})

test_that("loadGenotypeRegion filters VCF by region", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Rsamtools")
  result <- suppressWarnings(loadGenotypeRegion(vcf_path, region = region_sub,
                                                   returnVariantInfo = TRUE))
  check_genotype_result(result, expected_ncol = 134L, label = "vcf region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("loadGenotypeRegion filters VCF indels", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(loadGenotypeRegion(vcf_path, keepIndel = FALSE,
                                                   returnVariantInfo = TRUE))
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

# --- readGenotypes: GDS (SNPRelate) -------------------------------------------

test_that("readGenotypes creates gds handle", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  handle <- readGenotypes(gds_path, format = "gds")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@format, "gds")
  expect_equal(handle@nSamples, n_samples)
  expect_equal(nrow(handle@snpInfo), n_variants)
})

test_that("loadGenotypeRegion loads GDS via dispatch", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- loadGenotypeRegion(gds_path, returnVariantInfo = TRUE)
  check_genotype_result(result, label = "dispatch gds")
})

test_that("loadGenotypeRegion filters GDS by region", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- loadGenotypeRegion(gds_path, region = region_sub,
                                  returnVariantInfo = TRUE)
  check_genotype_result(result, expected_ncol = 134L, label = "gds region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("loadGenotypeRegion filters GDS indels", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- loadGenotypeRegion(gds_path, keepIndel = FALSE,
                                  returnVariantInfo = TRUE)
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("loadGenotypeRegion errors on empty region for GDS", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  expect_error(loadGenotypeRegion(gds_path, region = "chr1:1-2"))
})

# --- Cross-format consistency -------------------------------------------------

test_that("all formats return same dimensions and positions via loadGenotypeRegion", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SNPRelate")

  p2 <- loadGenotypeRegion(plink_prefix, returnVariantInfo = TRUE)
  vcf <- suppressWarnings(loadGenotypeRegion(vcf_path, returnVariantInfo = TRUE))
  gds <- loadGenotypeRegion(gds_path, returnVariantInfo = TRUE)

  # Same dimensions
  expect_equal(dim(p2$X), dim(vcf$X))
  expect_equal(dim(p2$X), dim(gds$X))

  # Same positions
  expect_equal(p2$variant_info$pos, vcf$variant_info$pos)
  expect_equal(p2$variant_info$pos, gds$variant_info$pos)
})

test_that("PLINK1 and PLINK2 readGenotypes return consistent alleles", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("pgenlibr")
  h1 <- readGenotypes(plink_prefix, format = "plink1")
  h2 <- readGenotypes(plink_prefix, format = "plink2")

  expect_equal(h1@snpInfo$A1, h2@snpInfo$A1)
  expect_equal(h1@snpInfo$A2, h2@snpInfo$A2)
})

# --- loadGenotypeRegion (dispatch) -----------------------------------------

test_that("loadGenotypeRegion dispatches to VCF by extension", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(loadGenotypeRegion(vcf_path, returnVariantInfo = TRUE))
  check_genotype_result(result, label = "dispatch vcf")
})

test_that("loadGenotypeRegion dispatches to GDS by extension", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- loadGenotypeRegion(gds_path, returnVariantInfo = TRUE)
  check_genotype_result(result, label = "dispatch gds")
})

test_that("loadGenotypeRegion dispatches to PLINK2 by prefix", {
  skip_if_not_installed("pgenlibr")
  result <- loadGenotypeRegion(plink_prefix, returnVariantInfo = TRUE)
  check_genotype_result(result, label = "dispatch plink2")
})

test_that("loadGenotypeRegion returns matrix when returnVariantInfo=FALSE", {
  skip_if_not_installed("pgenlibr")
  result <- loadGenotypeRegion(plink_prefix)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_samples)
  expect_equal(ncol(result), n_variants)
})

test_that("loadGenotypeRegion applies region filter", {
  skip_if_not_installed("pgenlibr")
  result <- loadGenotypeRegion(plink_prefix, region = region_sub, returnVariantInfo = TRUE)
  expect_equal(ncol(result$X), 134L)
})

test_that("loadGenotypeRegion errors on unrecognized format", {
  expect_error(loadGenotypeRegion("/nonexistent/file.xyz"), "not found")
})
