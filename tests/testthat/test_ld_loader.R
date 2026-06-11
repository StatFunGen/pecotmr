context("ldLoader")

# ===========================================================================
# ldLoader: input validation
# ===========================================================================

test_that("ldLoader errors when no source is provided", {
  expect_error(ldLoader(), "Provide exactly one")
})

test_that("ldLoader errors when multiple sources are provided", {
  R <- list(matrix(1, 2, 2))
  X <- list(matrix(1, 3, 2))
  expect_error(ldLoader(rList = R, xList = X), "Provide exactly one")
})

# ===========================================================================
# ldLoader: R_list branch
# ===========================================================================

test_that("ldLoader with R_list returns a function", {
  R <- list(matrix(c(1, 0.5, 0.5, 1), 2, 2))
  loader <- ldLoader(rList = R)
  expect_type(loader, "closure")
})

test_that("ldLoader R_list returns correct matrix", {
  R1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  R2 <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  loader <- ldLoader(rList = list(R1, R2))
  expect_equal(loader(1), R1)
  expect_equal(loader(2), R2)
})

test_that("ldLoader R_list with max_variants downsamples", {
  set.seed(42)
  R <- matrix(0.1, 10, 10)
  diag(R) <- 1
  loader <- ldLoader(rList = list(R), maxVariants = 5)
  result <- loader(1)
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
})

test_that("ldLoader R_list without max_variants returns full matrix", {
  R <- matrix(0.1, 10, 10)
  diag(R) <- 1
  loader <- ldLoader(rList = list(R))
  result <- loader(1)
  expect_equal(nrow(result), 10)
})

test_that("ldLoader R_list max_variants larger than matrix returns full matrix", {
  R <- matrix(0.1, 3, 3)
  diag(R) <- 1
  loader <- ldLoader(rList = list(R), maxVariants = 100)
  result <- loader(1)
  expect_equal(nrow(result), 3)
})

# ===========================================================================
# ldLoader: X_list branch
# ===========================================================================

test_that("ldLoader with X_list returns a function", {
  X <- list(matrix(rnorm(30), 10, 3))
  loader <- ldLoader(xList = X)
  expect_type(loader, "closure")
})

test_that("ldLoader X_list returns correct matrix", {
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(1:8, 4, 2)
  loader <- ldLoader(xList = list(X1, X2))
  expect_equal(loader(1), X1)
  expect_equal(loader(2), X2)
})

test_that("ldLoader X_list with max_variants downsamples columns", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  loader <- ldLoader(xList = list(X), maxVariants = 3)
  result <- loader(1)
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
})

test_that("ldLoader X_list max_variants larger than ncol returns full matrix", {
  X <- matrix(rnorm(12), 4, 3)
  loader <- ldLoader(xList = list(X), maxVariants = 100)
  result <- loader(1)
  expect_equal(ncol(result), 3)
})

# ===========================================================================
# ldLoader: ld_meta_path branch validation
# ===========================================================================

test_that("ldLoader with ld_meta_path but no regions errors", {
  expect_error(
    ldLoader(ldMetaPath = "/some/path"),
    "regions.*required"
  )
})

# ===========================================================================
# ldLoader: LD_info branch validation
# ===========================================================================

test_that("ldLoader with LD_info errors when not a data.frame", {
  expect_error(
    ldLoader(ldInfo = "not_a_df"),
    "ldInfo must be a data.frame"
  )
})

test_that("ldLoader with LD_info errors when missing LD_file column", {
  expect_error(
    ldLoader(ldInfo = data.frame(col1 = "a")),
    "ldInfo must be a data.frame with column 'LD_file'"
  )
})

# ===========================================================================
# ldLoader: LD_info branch with real genotype fixtures
# ===========================================================================

test_data_dir <- test_path("test_data")

test_that("ldLoader LD_info loads LD from PLINK2 files", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(test_data_dir, "test_variants")
  loader <- ldLoader(ldInfo = data.frame(LD_file = plink_prefix))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_equal(ncol(mat), 349L)
  expect_true(isSymmetric(mat))
  expect_true(all(abs(diag(mat) - 1) < 1e-10))
})

test_that("ldLoader LD_info loads LD from VCF file", {
  skip_if_not_installed("VariantAnnotation")
  vcf_path <- file.path(test_data_dir, "test_variants.vcf.gz")
  loader <- ldLoader(ldInfo = data.frame(LD_file = vcf_path))
  mat <- suppressWarnings(loader(1))
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_true(isSymmetric(mat))
})

test_that("ldLoader LD_info loads LD from GDS file", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  gds_path <- file.path(test_data_dir, "test_variants.gds")
  loader <- ldLoader(ldInfo = data.frame(LD_file = gds_path))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_true(isSymmetric(mat))
})

test_that("ldLoader LD_info loads LD from PLINK1 files", {
  skip_if_not_installed("snpStats")
  plink1_prefix <- file.path(test_data_dir, "protocol_example.genotype")
  loader <- ldLoader(ldInfo = data.frame(LD_file = plink1_prefix))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_true(isSymmetric(mat))
})

test_that("ldLoader LD_info loads pre-computed .cor.xz blocks", {
  ld_file <- file.path(test_data_dir, "LD_block_1.chr1_1000_1200.float16.txt.xz")
  bim_file <- file.path(test_data_dir, "LD_block_1.chr1_1000_1200.float16.bim")

  # Mock processLdMatrix to wrap its result in an LdData S4 object,
  # since extract_ld_matrix now requires an LdData.
  real_process <- pecotmr:::processLdMatrix
  local_mocked_bindings(
    processLdMatrix = function(LD_file_path, snp_file_path = NULL) {
      result <- real_process(LD_file_path, snp_file_path)
      mat <- result$LD_matrix
      variant_ids <- result$LD_variants$variants
      ref_panel <- pecotmr:::parseVariantId(variant_ids)
      ref_panel$variant_id <- variant_ids
      variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
      bm <- data.frame(
        block_id = 1L, chrom = as.character(ref_panel$chrom[1]),
        block_start = min(ref_panel$pos), block_end = max(ref_panel$pos),
        size = length(variant_ids), start_idx = 1L,
        end_idx = length(variant_ids), stringsAsFactors = FALSE
      )
      LdData(correlation = mat, variants = variants_gr, blockMetadata = bm)
    },
    .package = "pecotmr"
  )

  loader <- ldLoader(ldInfo = data.frame(LD_file = ld_file, SNP_file = bim_file))
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_true(isSymmetric(mat))
  expect_true(nrow(mat) > 0)
})

test_that("ldLoader LD_info with max_variants subsamples", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(test_data_dir, "test_variants")
  set.seed(42)
  loader <- ldLoader(ldInfo = data.frame(LD_file = plink_prefix), maxVariants = 20)
  mat <- loader(1)
  expect_equal(nrow(mat), 20L)
  expect_equal(ncol(mat), 20L)
})

test_that("ldLoader LD_info returns consistent LD across formats", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  plink_prefix <- file.path(test_data_dir, "test_variants")
  gds_path <- file.path(test_data_dir, "test_variants.gds")
  loader_plink <- ldLoader(ldInfo = data.frame(LD_file = plink_prefix))
  loader_gds <- ldLoader(ldInfo = data.frame(LD_file = gds_path))
  mat_plink <- loader_plink(1)
  mat_gds <- loader_gds(1)
  expect_equal(dim(mat_plink), dim(mat_gds))
})

# ===========================================================================
# ldLoader: ld_meta_path branch with real genotype fixtures
# ===========================================================================

test_that("ldLoader ld_meta_path loads LD from PLINK2 metadata", {
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
  loader <- ldLoader(ldMetaPath = meta_file, regions = region)
  mat <- loader(1)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
  expect_equal(ncol(mat), 349L)
})

test_that("ldLoader ld_meta_path loads LD from VCF metadata", {
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
  loader <- ldLoader(ldMetaPath = meta_file, regions = region)
  mat <- suppressWarnings(loader(1))
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 349L)
})
