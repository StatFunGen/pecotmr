# Tests for genotype loading functions:
#   load_plink1_data, load_plink2_data, load_vcf_data, load_gds_data,
#   load_genotype_region

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

# Shared helper: validate the output structure from any loader
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

# ─── PLINK1 (snpStats) ───────────────────────────────────────────────────────

test_that("load_plink1_data loads all variants", {
  skip_if_not_installed("snpStats")
  result <- load_plink1_data(plink_prefix)
  check_genotype_result(result, label = "plink1 all")
})

test_that("load_plink1_data filters by region", {
  skip_if_not_installed("snpStats")
  result <- load_plink1_data(plink_prefix, region = region_sub)
  check_genotype_result(result, expected_ncol = 134L, label = "plink1 region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_plink1_data filters indels", {
  skip_if_not_installed("snpStats")
  result <- load_plink1_data(plink_prefix, keep_indel = FALSE)
  # Should have fewer variants than 100 (fixture has indels)
  expect_lt(ncol(result$X), n_variants)
  # All remaining alleles should be single characters
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("load_plink1_data errors on missing files", {
  skip_if_not_installed("snpStats")
  expect_error(load_plink1_data("/nonexistent/path"), "not found")
})

test_that("load_plink1_data errors on empty region", {
  skip_if_not_installed("snpStats")
  expect_error(load_plink1_data(plink_prefix, region = "chr1:1-2"), class = "NoSNPsError")
})

# ─── PLINK2 (pgenlibr) ───────────────────────────────────────────────────────

test_that("load_plink2_data loads all variants", {
  skip_if_not_installed("pgenlibr")
  result <- load_plink2_data(plink_prefix)
  check_genotype_result(result, label = "plink2 all")
})

test_that("load_plink2_data filters by region", {
  skip_if_not_installed("pgenlibr")
  result <- load_plink2_data(plink_prefix, region = region_sub)
  check_genotype_result(result, expected_ncol = 134L, label = "plink2 region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_plink2_data filters indels", {
  skip_if_not_installed("pgenlibr")
  result <- load_plink2_data(plink_prefix, keep_indel = FALSE)
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("load_plink2_data errors on missing files", {
  skip_if_not_installed("pgenlibr")
  expect_error(load_plink2_data("/nonexistent/path"), "not found")
})

test_that("load_plink2_data errors on empty region", {
  skip_if_not_installed("pgenlibr")
  expect_error(load_plink2_data(plink_prefix, region = "chr1:1-2"), class = "NoSNPsError")
})

# ─── VCF (VariantAnnotation) ─────────────────────────────────────────────────

test_that("load_vcf_data loads all variants", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(load_vcf_data(vcf_path))
  check_genotype_result(result, label = "vcf all")
})

test_that("load_vcf_data filters by region", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Rsamtools")
  result <- suppressWarnings(load_vcf_data(vcf_path, region = region_sub))
  check_genotype_result(result, expected_ncol = 134L, label = "vcf region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_vcf_data filters indels", {
  skip_if_not_installed("VariantAnnotation")
  result <- suppressWarnings(load_vcf_data(vcf_path, keep_indel = FALSE))
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("load_vcf_data errors on empty region", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Rsamtools")
  expect_error(suppressWarnings(load_vcf_data(vcf_path, region = "chr1:1-2")))
})

# ─── GDS (SNPRelate) ─────────────────────────────────────────────────────────

test_that("load_gds_data loads all variants", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_gds_data(gds_path)
  check_genotype_result(result, label = "gds all")
})

test_that("load_gds_data filters by region", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_gds_data(gds_path, region = region_sub)
  check_genotype_result(result, expected_ncol = 134L, label = "gds region")
  expect_true(all(result$variant_info$pos >= 17513228 & result$variant_info$pos <= 17550000))
})

test_that("load_gds_data filters indels", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  result <- load_gds_data(gds_path, keep_indel = FALSE)
  expect_lt(ncol(result$X), n_variants)
  expect_true(all(nchar(result$variant_info$A1) == 1))
  expect_true(all(nchar(result$variant_info$A2) == 1))
})

test_that("load_gds_data errors on empty region", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  expect_error(load_gds_data(gds_path, region = "chr1:1-2"), class = "NoSNPsError")
})

# ─── Cross-format consistency ─────────────────────────────────────────────────

test_that("all formats return same dimensions and positions", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SNPRelate")
  p1 <- load_plink1_data(plink_prefix)
  p2 <- load_plink2_data(plink_prefix)
  vcf <- suppressWarnings(load_vcf_data(vcf_path))
  gds <- load_gds_data(gds_path)

  # Same dimensions
  expect_equal(dim(p1$X), dim(p2$X))
  expect_equal(dim(p1$X), dim(vcf$X))
  expect_equal(dim(p1$X), dim(gds$X))

  # Same positions
  expect_equal(p1$variant_info$pos, p2$variant_info$pos)
  expect_equal(p1$variant_info$pos, vcf$variant_info$pos)
  expect_equal(p1$variant_info$pos, gds$variant_info$pos)

  # Same variant IDs between PLINK formats
  expect_equal(p1$variant_info$id, p2$variant_info$id)
})

test_that("PLINK1 and PLINK2 return consistent alleles", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("pgenlibr")
  p1 <- load_plink1_data(plink_prefix)
  p2 <- load_plink2_data(plink_prefix)

  expect_equal(p1$variant_info$A1, p2$variant_info$A1)
  expect_equal(p1$variant_info$A2, p2$variant_info$A2)
})

# ─── load_genotype_region (dispatch) ─────────────────────────────────────────

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
