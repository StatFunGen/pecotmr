# Tests for genotype loading via readGenotypes + extractBlockGenotypes,
# and the load_genotype_region dispatcher.

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

# Shared helper: validate the output structure from load_genotype_region
# (with return_variant_info=TRUE)
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
  expect_equal(handle@n_samples, n_samples)
  expect_equal(nrow(handle@snp_info), n_variants)
})

test_that("load_genotype_region loads plink1 via dispatch", {
  skip_if_not_installed("snpStats")
  # Use .genotype suffix plink1 files tested elsewhere in test_file_utils
  plink1_path <- file.path(test_data_dir, "protocol_example.genotype")
  skip_if(!file.exists(paste0(plink1_path, ".bed")), "plink1 test fixture missing")
  result <- load_genotype_region(plink1_path, return_variant_info = TRUE)
  expect_true(is.list(result))
  expect_true(is.matrix(result$X))
})

# --- readGenotypes: PLINK2 (pgenlibr) ----------------------------------------

test_that("readGenotypes creates plink2 handle", {
  skip_if_not_installed("pgenlibr")
  handle <- readGenotypes(plink_prefix, format = "plink2")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@format, "plink2")
  expect_equal(handle@n_samples, n_samples)
  expect_equal(nrow(handle@snp_info), n_variants)
})

test_that("extractBlockGenotypes works for plink2", {
  skip_if_not_installed("pgenlibr")
  handle <- readGenotypes(plink_prefix, format = "plink2")
  rse <- extractBlockGenotypes(handle, seq_len(nrow(handle@snp_info)))
  expect_s4_class(rse, "SummarizedExperiment")
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  expect_equal(nrow(dosage), n_variants)
  expect_equal(ncol(dosage), n_samples)
})

test_that("load_genotype_region filters plink2 by region", {
  skip_if_not_installed("pgenlibr")
  result <- load_genotype_region(plink_prefix, region = region_sub,
                                  return_variant_info = TRUE)
  check_genotype_result(result, expected_ncol = 134L, label = "plink2 region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_genotype_region filters plink2 indels", {
  skip_if_not_installed("pgenlibr")
  result <- load_genotype_region(plink_prefix, keep_indel = FALSE,
                                  return_variant_info = TRUE)
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("load_genotype_region errors on empty region for plink2", {
  skip_if_not_installed("pgenlibr")
  expect_error(load_genotype_region(plink_prefix, region = "chr1:1-2"))
})

# --- readGenotypes: VCF (VariantAnnotation) -----------------------------------

test_that("readGenotypes creates vcf handle", {
  skip_if_not_installed("VariantAnnotation")
  handle <- readGenotypes(vcf_path, format = "vcf")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@format, "vcf")
  expect_equal(handle@n_samples, n_samples)
  expect_equal(nrow(handle@snp_info), n_variants)
})

test_that("load_genotype_region loads VCF via dispatch", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(load_genotype_region(vcf_path, return_variant_info = TRUE))
  check_genotype_result(result, label = "dispatch vcf")
})

test_that("load_genotype_region filters VCF by region", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Rsamtools")
  result <- suppressWarnings(load_genotype_region(vcf_path, region = region_sub,
                                                   return_variant_info = TRUE))
  check_genotype_result(result, expected_ncol = 134L, label = "vcf region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_genotype_region filters VCF indels", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(load_genotype_region(vcf_path, keep_indel = FALSE,
                                                   return_variant_info = TRUE))
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
  expect_equal(handle@n_samples, n_samples)
  expect_equal(nrow(handle@snp_info), n_variants)
})

test_that("load_genotype_region loads GDS via dispatch", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_genotype_region(gds_path, return_variant_info = TRUE)
  check_genotype_result(result, label = "dispatch gds")
})

test_that("load_genotype_region filters GDS by region", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_genotype_region(gds_path, region = region_sub,
                                  return_variant_info = TRUE)
  check_genotype_result(result, expected_ncol = 134L, label = "gds region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_genotype_region filters GDS indels", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_genotype_region(gds_path, keep_indel = FALSE,
                                  return_variant_info = TRUE)
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("load_genotype_region errors on empty region for GDS", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  expect_error(load_genotype_region(gds_path, region = "chr1:1-2"))
})

# --- Cross-format consistency -------------------------------------------------

test_that("all formats return same dimensions and positions via load_genotype_region", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SNPRelate")

  p2 <- load_genotype_region(plink_prefix, return_variant_info = TRUE)
  vcf <- suppressWarnings(load_genotype_region(vcf_path, return_variant_info = TRUE))
  gds <- load_genotype_region(gds_path, return_variant_info = TRUE)

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

  expect_equal(h1@snp_info$A1, h2@snp_info$A1)
  expect_equal(h1@snp_info$A2, h2@snp_info$A2)
})

# --- load_genotype_region (dispatch) -----------------------------------------

test_that("load_genotype_region dispatches to VCF by extension", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(load_genotype_region(vcf_path, return_variant_info = TRUE))
  check_genotype_result(result, label = "dispatch vcf")
})

test_that("load_genotype_region dispatches to GDS by extension", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_genotype_region(gds_path, return_variant_info = TRUE)
  check_genotype_result(result, label = "dispatch gds")
})

test_that("load_genotype_region dispatches to PLINK2 by prefix", {
  skip_if_not_installed("pgenlibr")
  result <- load_genotype_region(plink_prefix, return_variant_info = TRUE)
  check_genotype_result(result, label = "dispatch plink2")
})

test_that("load_genotype_region returns matrix when return_variant_info=FALSE", {
  skip_if_not_installed("pgenlibr")
  result <- load_genotype_region(plink_prefix)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_samples)
  expect_equal(ncol(result), n_variants)
})

test_that("load_genotype_region applies region filter", {
  skip_if_not_installed("pgenlibr")
  result <- load_genotype_region(plink_prefix, region = region_sub, return_variant_info = TRUE)
  expect_equal(ncol(result$X), 134L)
})

test_that("load_genotype_region errors on unrecognized format", {
  expect_error(load_genotype_region("/nonexistent/file.xyz"), "not found")
})
