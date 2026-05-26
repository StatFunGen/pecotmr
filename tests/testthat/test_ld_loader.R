context("ld_loader")

# ===========================================================================
# ld_loader: input validation
# ===========================================================================

test_that("ld_loader errors when no source is provided", {
  expect_error(ld_loader(), "Provide exactly one")
})

test_that("ld_loader errors when multiple sources are provided", {
  R <- list(matrix(1, 2, 2))
  X <- list(matrix(1, 3, 2))
  expect_error(ld_loader(R_list = R, X_list = X), "Provide exactly one")
})

# ===========================================================================
# ld_loader: R_list branch
# ===========================================================================

test_that("ld_loader with R_list returns a function", {
  R <- list(matrix(c(1, 0.5, 0.5, 1), 2, 2))
  loader <- ld_loader(R_list = R)
  expect_type(loader, "closure")
})

test_that("ld_loader R_list returns correct matrix", {
  R1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  R2 <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  loader <- ld_loader(R_list = list(R1, R2))
  expect_equal(loader(1), R1)
  expect_equal(loader(2), R2)
})

test_that("ld_loader R_list with max_variants downsamples", {
  set.seed(42)
  R <- matrix(0.1, 10, 10)
  diag(R) <- 1
  loader <- ld_loader(R_list = list(R), max_variants = 5)
  result <- loader(1)
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
})

test_that("ld_loader R_list without max_variants returns full matrix", {
  R <- matrix(0.1, 10, 10)
  diag(R) <- 1
  loader <- ld_loader(R_list = list(R))
  result <- loader(1)
  expect_equal(nrow(result), 10)
})

test_that("ld_loader R_list max_variants larger than matrix returns full matrix", {
  R <- matrix(0.1, 3, 3)
  diag(R) <- 1
  loader <- ld_loader(R_list = list(R), max_variants = 100)
  result <- loader(1)
  expect_equal(nrow(result), 3)
})

# ===========================================================================
# ld_loader: X_list branch
# ===========================================================================

test_that("ld_loader with X_list returns a function", {
  X <- list(matrix(rnorm(30), 10, 3))
  loader <- ld_loader(X_list = X)
  expect_type(loader, "closure")
})

test_that("ld_loader X_list returns correct matrix", {
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(1:8, 4, 2)
  loader <- ld_loader(X_list = list(X1, X2))
  expect_equal(loader(1), X1)
  expect_equal(loader(2), X2)
})

test_that("ld_loader X_list with max_variants downsamples columns", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  loader <- ld_loader(X_list = list(X), max_variants = 3)
  result <- loader(1)
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
})

test_that("ld_loader X_list max_variants larger than ncol returns full matrix", {
  X <- matrix(rnorm(12), 4, 3)
  loader <- ld_loader(X_list = list(X), max_variants = 100)
  result <- loader(1)
  expect_equal(ncol(result), 3)
})

# ===========================================================================
# ld_loader: ld_meta_path branch validation
# ===========================================================================

test_that("ld_loader with ld_meta_path but no regions errors", {
  expect_error(
    ld_loader(ld_meta_path = "/some/path"),
    "regions.*required"
  )
})

# ===========================================================================
# ld_loader: LD_info branch validation
# ===========================================================================

test_that("ld_loader with LD_info errors when not a data.frame", {
  expect_error(
    ld_loader(LD_info = "not_a_df"),
    "LD_info must be a data.frame"
  )
})

test_that("ld_loader with LD_info errors when missing LD_file column", {
  expect_error(
    ld_loader(LD_info = data.frame(col1 = "a")),
    "LD_info must be a data.frame with column 'LD_file'"
  )
})

# ===========================================================================
# ld_loader: LD_info branch with real genotype fixtures
# ===========================================================================

test_data_dir <- test_path("test_data")

test_that("ld_loader LD_info loads LD from PLINK2 files", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(test_data_dir, "test_variants")
  loader <- ld_loader(LD_info = data.frame(LD_file = plink_prefix))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_equal(ncol(mat), 349L)
  expect_true(isSymmetric(mat))
  expect_true(all(abs(diag(mat) - 1) < 1e-10))
})

test_that("ld_loader LD_info loads LD from VCF file", {
  skip_if_not_installed("VariantAnnotation")
  vcf_path <- file.path(test_data_dir, "test_variants.vcf.gz")
  loader <- ld_loader(LD_info = data.frame(LD_file = vcf_path))
  mat <- suppressWarnings(loader(1))
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_true(isSymmetric(mat))
})

test_that("ld_loader LD_info loads LD from GDS file", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  gds_path <- file.path(test_data_dir, "test_variants.gds")
  loader <- ld_loader(LD_info = data.frame(LD_file = gds_path))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_true(isSymmetric(mat))
})

test_that("ld_loader LD_info loads LD from PLINK1 files", {
  skip_if_not_installed("snpStats")
  plink1_prefix <- file.path(test_data_dir, "protocol_example.genotype")
  loader <- ld_loader(LD_info = data.frame(LD_file = plink1_prefix))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_true(isSymmetric(mat))
})

test_that("ld_loader LD_info loads pre-computed .cor.xz blocks", {
  ld_file <- file.path(test_data_dir, "LD_block_1.chr1_1000_1200.float16.txt.xz")
  bim_file <- file.path(test_data_dir, "LD_block_1.chr1_1000_1200.float16.bim")

  # Mock process_LD_matrix to wrap its result in an LDData S4 object,
  # since extract_ld_matrix now requires an LDData.
  real_process <- pecotmr:::process_LD_matrix
  local_mocked_bindings(
    process_LD_matrix = function(LD_file_path, snp_file_path = NULL) {
      result <- real_process(LD_file_path, snp_file_path)
      mat <- result$LD_matrix
      variant_ids <- result$LD_variants$variants
      ref_panel <- pecotmr:::parse_variant_id(variant_ids)
      ref_panel$variant_id <- variant_ids
      variants_gr <- pecotmr:::.ref_panel_to_granges(ref_panel)
      bm <- data.frame(
        block_id = 1L, chrom = as.character(ref_panel$chrom[1]),
        block_start = min(ref_panel$pos), block_end = max(ref_panel$pos),
        size = length(variant_ids), start_idx = 1L,
        end_idx = length(variant_ids), stringsAsFactors = FALSE
      )
      LDData(correlation = mat, variants = variants_gr, block_metadata = bm)
    },
    .package = "pecotmr"
  )

  loader <- ld_loader(LD_info = data.frame(LD_file = ld_file, SNP_file = bim_file))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_true(isSymmetric(mat))
  expect_true(nrow(mat) > 0)
})

test_that("ld_loader LD_info with max_variants subsamples", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(test_data_dir, "test_variants")
  set.seed(42)
  loader <- ld_loader(LD_info = data.frame(LD_file = plink_prefix), max_variants = 20)
  mat <- loader(1)
  expect_equal(nrow(mat), 20L)
  expect_equal(ncol(mat), 20L)
})

test_that("ld_loader LD_info returns consistent LD across formats", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  plink_prefix <- file.path(test_data_dir, "test_variants")
  gds_path <- file.path(test_data_dir, "test_variants.gds")
  loader_plink <- ld_loader(LD_info = data.frame(LD_file = plink_prefix))
  loader_gds <- ld_loader(LD_info = data.frame(LD_file = gds_path))
  mat_plink <- loader_plink(1)
  mat_gds <- loader_gds(1)
  expect_equal(dim(mat_plink), dim(mat_gds))
})

# ===========================================================================
# ld_loader: ld_meta_path branch with real genotype fixtures
# ===========================================================================

test_that("ld_loader ld_meta_path loads LD from PLINK2 metadata", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(test_data_dir, "ld_meta_plink2_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  region <- "chr21:17513228-17592874"
  loader <- ld_loader(ld_meta_path = meta_file, regions = region)
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_equal(ncol(mat), 349L)
})

test_that("ld_loader ld_meta_path loads LD from VCF metadata", {
  skip_if_not_installed("VariantAnnotation")
  meta_file <- file.path(test_data_dir, "ld_meta_vcf_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  region <- "chr21:17513228-17592874"
  loader <- ld_loader(ld_meta_path = meta_file, regions = region)
  mat <- suppressWarnings(loader(1))
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
})
