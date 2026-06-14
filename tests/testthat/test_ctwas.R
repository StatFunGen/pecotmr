context("ctwas")

# ===========================================================================
#  ctwas wrapper tests
# ===========================================================================


# ---------- trimCtwasVariants --------------------------------------------

# Helper: build a minimal region_data structure that trimCtwasVariants expects
make_mock_region_data <- function() {
  # Variant IDs in canonical format (chr:pos:A2:A1)
  variant_ids <- c("chr1:1000:A:G", "chr1:2000:C:T", "chr1:3000:G:A", "chr1:4000:T:C")

  # Weight matrix (4 variants x 1 weight column)
  wgt <- matrix(c(0.5, 0.0001, 0.3, -0.2), nrow = 4, ncol = 1)
  rownames(wgt) <- variant_ids

  gene_id <- "GENE1|ctx1"
  context <- "ctx1"
  study <- "study1"

  weights <- list()
  weights[[gene_id]] <- list()
  weights[[gene_id]][[study]] <- list(
    wgt = wgt,
    context = context,
    p0 = 1000,
    p1 = 4000
  )

  # SuSiE intermediate info
  pip_vals <- c(0.8, 0.05, 0.6, 0.02)
  names(pip_vals) <- variant_ids

  susieWeightsIntermediate <- list()
  susieWeightsIntermediate[["GENE1"]] <- list()
  susieWeightsIntermediate[["GENE1"]][[context]] <- list(
    pip = pip_vals,
    csVariants = list(variant_ids[c(1, 3)]),
    csPurity = list(minAbsCorr = 0.9)
  )

  list(
    weights = weights,
    susieWeightsIntermediate = susieWeightsIntermediate
  )
}

test_that("trimCtwasVariants removes variants below weight cutoff", {
  rd <- make_mock_region_data()
  # Default cutoff 1e-5, variant 2 has weight 0.0001 (above), so all 4 should pass default cutoff
  result <- trimCtwasVariants(rd, twasWeightCutoff = 1e-5)
  expect_true(is.list(result))
  # With a higher cutoff, the near-zero variant should be removed
  result_strict <- trimCtwasVariants(rd, twasWeightCutoff = 0.001)
  # study1 should exist in result
  expect_true("study1" %in% names(result_strict))
  # Get the gene-level result
  gene_weights <- result_strict[["study1"]][["GENE1|ctx1"]]
  # Variant 2 has abs(weight) = 0.0001 < 0.001, so should be removed
  expect_false("chr1:2000:C:T" %in% rownames(gene_weights$wgt))
})

test_that("trimCtwasVariants removes gene when all weights below cutoff", {
  rd <- make_mock_region_data()
  # Set cutoff so high that all variants are dropped
  result <- trimCtwasVariants(rd, twasWeightCutoff = 10)
  # Gene should be removed entirely since no weights pass the cutoff
  # Result should be an empty list
  expect_equal(length(result), 0)
})

test_that("trimCtwasVariants returns result keyed by study", {
  rd <- make_mock_region_data()
  result <- trimCtwasVariants(rd, twasWeightCutoff = 1e-5)
  # merge_by_study reorganizes: weights[[study]][[group]]
  expect_true("study1" %in% names(result))
  expect_true("GENE1|ctx1" %in% names(result[["study1"]]))
})

test_that("trimCtwasVariants updates p0 and p1 positions", {
  rd <- make_mock_region_data()
  # Use a weight cutoff that removes the variant at position 2000
  result <- trimCtwasVariants(rd, twasWeightCutoff = 0.001)
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  # p0 and p1 should reflect the range of remaining variant positions
  remaining_positions <- as.integer(sapply(
    rownames(gene_weights$wgt),
    function(v) strsplit(v, ":")[[1]][2]
  ))
  expect_equal(gene_weights$p0, min(remaining_positions))
  expect_equal(gene_weights$p1, max(remaining_positions))
})

test_that("trimCtwasVariants respects max_num_variants", {
  rd <- make_mock_region_data()
  # Request max 2 variants; since nrow(wgt) == 4 >= max_num_variants == 2,
  # it triggers select_variants which picks by PIP priority
  result <- trimCtwasVariants(rd, twasWeightCutoff = 1e-5, maxNumVariants = 2)
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  expect_true(nrow(gene_weights$wgt) <= 2)
})

test_that("trimCtwasVariants handles NA weights by removing group", {
  rd <- make_mock_region_data()
  # Replace all weights with NA
  rd$weights[["GENE1|ctx1"]][["study1"]]$wgt[, 1] <- NA
  result <- trimCtwasVariants(rd, twasWeightCutoff = 0)
  # The group should be removed because all weights are NA
  expect_equal(length(result), 0)
})

test_that("trimCtwasVariants handles multiple genes", {
  rd <- make_mock_region_data()

  # Add a second gene
  variant_ids2 <- c("chr1:5000:A:G", "chr1:6000:C:T")
  wgt2 <- matrix(c(0.4, -0.3), nrow = 2, ncol = 1)
  rownames(wgt2) <- variant_ids2

  rd$weights[["GENE2|ctx1"]] <- list()
  rd$weights[["GENE2|ctx1"]][["study1"]] <- list(
    wgt = wgt2,
    context = "ctx1",
    p0 = 5000,
    p1 = 6000
  )

  pip_vals2 <- c(0.7, 0.4)
  names(pip_vals2) <- variant_ids2
  rd$susieWeightsIntermediate[["GENE2"]] <- list()
  rd$susieWeightsIntermediate[["GENE2"]][["ctx1"]] <- list(
    pip = pip_vals2,
    csVariants = list(variant_ids2[1]),
    csPurity = list(minAbsCorr = 0.95)
  )

  result <- trimCtwasVariants(rd, twasWeightCutoff = 1e-5)
  expect_true("GENE1|ctx1" %in% names(result[["study1"]]))
  expect_true("GENE2|ctx1" %in% names(result[["study1"]]))
})

test_that("trimCtwasVariants select_variants uses csMinCor to include CS variants", {
  rd <- make_mock_region_data()
  # csPurity minAbsCorr = 0.9, so with csMinCor = 0.8 the CS variants
  # (variant 1 and 3) should be included. Max 2 variants.
  result <- trimCtwasVariants(rd,
    twasWeightCutoff = 1e-5,
    csMinCor = 0.8,
    minPipCutoff = 0.0,
    maxNumVariants = 2
  )
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  included <- rownames(gene_weights$wgt)
  # CS variants chr1:1000:A:G and chr1:3000:G:A have highest PIPs (0.8 and 0.6)
  # and are in the CS, so they should be prioritized
  expect_true("chr1:1000:A:G" %in% included)
  expect_true("chr1:3000:G:A" %in% included)
})

# ===========================================================================
#  Deprecated wrapper: ctwasBimfileLoader
# ===========================================================================

test_that("ctwasBimfileLoader reads .bim and returns legacy column names", {
  bim_path <- tempfile(fileext = ".bim")
  on.exit(unlink(bim_path), add = TRUE)
  cat("1\tchr1:1000:A:G\t0\t1000\tA\tG\n", file = bim_path)
  cat("1\tchr1:2000:C:T\t0\t2000\tC\tT\n", file = bim_path, append = TRUE)

  expect_warning(
    res <- ctwasBimfileLoader(bim_path),
    "deprecated"
  )
  expect_equal(colnames(res), c("chrom", "id", "GD", "pos", "A1", "A2"))
  expect_equal(nrow(res), 2)
  expect_equal(res$pos, c(1000, 2000))
})

test_that("ctwasBimfileLoader accepts .bed path and resolves .bim", {
  bim_path <- tempfile(fileext = ".bim")
  bed_path <- sub("\\.bim$", ".bed", bim_path)
  on.exit(unlink(c(bim_path, bed_path)), add = TRUE)
  cat("22\trs100\t0\t50000\tA\tG\n", file = bim_path)

  expect_warning(
    res <- ctwasBimfileLoader(bed_path),
    "deprecated"
  )
  expect_equal(nrow(res), 1)
  expect_equal(res$pos, 50000)
})

test_that("ctwasBimfileLoader normalizes variant IDs", {
  bim_path <- tempfile(fileext = ".bim")
  on.exit(unlink(bim_path), add = TRUE)
  cat("1\tchr1:1000:A:G\t0\t1000\tA\tG\n", file = bim_path)

  expect_warning(
    res <- ctwasBimfileLoader(bim_path),
    "deprecated"
  )
  # normalizeVariantId should have been applied
  expect_equal(res$id, normalizeVariantId("chr1:1000:A:G"))
})

test_that("ctwasBimfileLoader works with real test fixture", {
  bim_path <- test_path("test_data", "protocol_example.genotype.bim")
  skip_if_not(file.exists(bim_path), "Test fixture not available")

  expect_warning(
    res <- ctwasBimfileLoader(bim_path),
    "deprecated"
  )
  expect_equal(colnames(res), c("chrom", "id", "GD", "pos", "A1", "A2"))
  expect_equal(nrow(res), 100)
})

# ===========================================================================
#  Deprecated wrapper: getCtwasMetaData
# ===========================================================================

test_that("getCtwasMetaData reads LD metadata and returns ldInfo + regionInfo", {
  meta_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("chr1", "1000", "2000", "block1.cor.xz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  cat(paste("chr1", "2000", "3000", "block2.cor.xz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  expect_warning(
    res <- getCtwasMetaData(meta_file),
    "deprecated"
  )
  expect_true(is.list(res))
  expect_true("ldInfo" %in% names(res))
  expect_true("regionInfo" %in% names(res))

  expect_equal(nrow(res$ldInfo), 2)
  expect_equal(colnames(res$ldInfo), c("region_id", "LD_file", "SNP_file"))
  expect_equal(res$ldInfo$region_id, c("1_1000_2000", "1_2000_3000"))

  expect_equal(nrow(res$regionInfo), 2)
  expect_equal(colnames(res$regionInfo), c("chrom", "start", "stop", "region_id"))
  expect_equal(res$regionInfo$chrom, c(1L, 1L))
  expect_equal(res$regionInfo$start, c(1000L, 2000L))
  expect_equal(res$regionInfo$stop, c(2000L, 3000L))
})

test_that("getCtwasMetaData subset_region_ids filters correctly", {
  meta_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("chr1", "1000", "2000", "block1.cor.xz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  cat(paste("chr1", "2000", "3000", "block2.cor.xz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  cat(paste("chr2", "5000", "6000", "block3.cor.xz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  expect_warning(
    res <- getCtwasMetaData(meta_file, subsetRegionIds = "1_1000_2000"),
    "deprecated"
  )
  expect_equal(nrow(res$regionInfo), 1)
  expect_equal(res$regionInfo$region_id, "1_1000_2000")
  # ldInfo is not subsetted (matches original behavior)
  expect_equal(nrow(res$ldInfo), 3)
})

test_that("getCtwasMetaData LD_file paths are relative to metadata directory", {
  tmpdir <- tempdir()
  meta_file <- file.path(tmpdir, "ld_meta.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("chr1", "100", "200", "subdir/block.cor.xz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  expect_warning(
    res <- getCtwasMetaData(meta_file),
    "deprecated"
  )
  expect_equal(res$ldInfo$LD_file, file.path(tmpdir, "subdir/block.cor.xz"))
  expect_equal(res$ldInfo$SNP_file, paste0(file.path(tmpdir, "subdir/block.cor.xz"), ".bim"))
})
