context("LD")
library(tidyverse)

generate_dummy_data <- function() {
  region <- data.frame(
    chrom = "chr1",
    start = c(1000),
    end = c(1190)
  )
  meta_df <- data.frame(
    chrom = "chr1",
    start = c(1000, 1200, 1400, 1600, 1800),
    end = c(1200, 1400, 1600, 1800, 2000),
    path = c(
      "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
      "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
      "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
      "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
      "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
    ))
  return(list(region = region, meta = meta_df))
}

# Generate a wider region that spans multiple blocks for partition testing
generate_multi_block_data <- function() {
  region <- data.frame(
    chrom = "chr1",
    start = c(1000),
    end = c(1500)
  )
  meta_df <- data.frame(
    chrom = "chr1",
    start = c(1000, 1200, 1400, 1600, 1800),
    end = c(1200, 1400, 1600, 1800, 2000),
    path = c(
      "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
      "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
      "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
      "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
      "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
    ))
  return(list(region = region, meta = meta_df))
}

test_that("Check that we correctly retrieve the names from the matrix",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  variants <- unlist(
    c("chr1:1000:A:G", "chr1:1040:A:G", "chr1:1080:A:G", "chr1:1120:A:G", "chr1:1160:A:G"))
  expect_equal(
    unlist(getVariantIds(res)),
    variants)
  file.remove(LD_meta_file_path)
})

test_that("Check that the LD block contains the correct information",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  # Variant names
  variants <- unlist(
    c("chr1:1000:A:G", "chr1:1040:A:G", "chr1:1080:A:G", "chr1:1120:A:G", "chr1:1160:A:G"))
  # Check LD Block 1
  ld_block_one <- getCorrelation(res)
  ld_block_one_original <- as.matrix(
    read_delim(
      "test_data/LD_block_1.chr1_1000_1200.float16.txt.xz",
      delim = " ", col_names = F))
  rownames(ld_block_one_original) <- colnames(ld_block_one_original) <- variants
  expect_equal(ld_block_one, ld_block_one_original)
  file.remove(LD_meta_file_path)
})

# ---- partition_LD_matrix ----

test_that("partition_LD_matrix correctly partitions a single block", {
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix first
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Now partition the matrix
  partitioned <- partition_LD_matrix(ld_data)

  # Expectations for single block case
  expect_equal(length(partitioned$ld_matrices), 1)
  expect_equal(nrow(partitioned$variant_indices), length(getVariantIds(ld_data)))
  expect_equal(unique(partitioned$variant_indices$block_id), 1)
  expect_identical(rownames(partitioned$ld_matrices[[1]]), getVariantIds(ld_data))
  expect_identical(colnames(partitioned$ld_matrices[[1]]), getVariantIds(ld_data))

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix correctly partitions multiple blocks", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Now partition the matrix without merging blocks
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = FALSE)

  # Check if we have the correct number of blocks
  # Should have block 1 (1000-1200), block 2 (1200-1400), and block 3 (1400-1600)
  expected_block_count <- 3
  expect_equal(length(partitioned$ld_matrices), expected_block_count)

  # Check if all variants are assigned to blocks
  expect_equal(nrow(partitioned$variant_indices), length(getVariantIds(ld_data)))

  # Check if block IDs are correct
  expect_setequal(unique(partitioned$variant_indices$block_id), 1:expected_block_count)

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix properly merges small blocks", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Set min_merged_block_size high enough to force merging
  # Each test block likely has 5 variants (based on the existing test)
  min_block_size <- 10

  # Now partition the matrix with block merging
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = TRUE,
                                    min_merged_block_size = min_block_size)

  # We expect fewer blocks after merging
  expect_lt(length(partitioned$ld_matrices), 3)

  # Check if all variants are still assigned to blocks
  expect_equal(nrow(partitioned$variant_indices), length(getVariantIds(ld_data)))

  # Check if merged blocks are larger than min_block_size
  block_sizes <- sapply(partitioned$ld_matrices, nrow)
  expect_true(all(block_sizes >= min_block_size | block_sizes == length(getVariantIds(ld_data))))

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix respects max_merged_block_size", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Set max_merged_block_size to a small value to prevent merging all blocks
  # Each test block likely has 5 variants (based on the existing test)
  max_block_size <- 8

  # Now partition the matrix with restricted block size
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = TRUE,
                                    min_merged_block_size = 2,
                                    max_merged_block_size = max_block_size)

  # Check if no block exceeds max_block_size
  block_sizes <- sapply(partitioned$ld_matrices, nrow)
  expect_true(all(block_sizes <= max_block_size))

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix handles empty matrix gracefully", {
  # Create an empty LD data structure
  empty_ld_data <- list(
    LD_matrix = matrix(0, 0, 0),
    LD_variants = character(0),
    block_metadata = data.frame(
      block_id = integer(0),
      chrom = character(0),
      size = integer(0),
      start_idx = integer(0),
      end_idx = integer(0)
    )
  )

  # Expect an error for empty matrix
  expect_error(partition_LD_matrix(empty_ld_data), "Empty or NULL LD matrix provided")
})

test_that("partition_LD_matrix validates block structure properly", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Create an invalid block structure by converting to a plain list and modifying
  bm <- getBlockMetadata(ld_data)
  vids <- getVariantIds(ld_data)
  ldmat <- getCorrelation(ld_data)
  invalid_ld_data <- list(
    LD_matrix = ldmat,
    LD_variants = vids,
    block_metadata = bm
  )

  # Assuming we have at least 2 blocks:
  if(nrow(invalid_ld_data$block_metadata) >= 2) {
    # Create overlapping blocks with invalid start/end indices
    invalid_ld_data$block_metadata$start_idx[2] <- invalid_ld_data$block_metadata$start_idx[1]
    invalid_ld_data$block_metadata$end_idx[1] <- invalid_ld_data$block_metadata$end_idx[2]

    # Introduce non-zero elements between blocks to trigger validation error
    if(length(invalid_ld_data$LD_variants) >= 2) {
      idx1 <- invalid_ld_data$block_metadata$start_idx[1]
      idx2 <- invalid_ld_data$block_metadata$start_idx[2] + 1
      if(idx1 <= length(invalid_ld_data$LD_variants) &&
         idx2 <= length(invalid_ld_data$LD_variants)) {
        var1 <- invalid_ld_data$LD_variants[idx1]
        var2 <- invalid_ld_data$LD_variants[idx2]
        invalid_ld_data$LD_matrix[var1, var2] <- 0.5
      }
    }

    # Expect an error for invalid block structure
    expect_error(partition_LD_matrix(invalid_ld_data), "Matrix lacks expected block structure")
  }

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix properly maps variants to blocks", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Partition without merging
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = FALSE)

  # Check that each variant is mapped to the correct block
  for(i in seq_along(partitioned$ld_matrices)) {
    # Get variants in this block matrix
    block_variants <- rownames(partitioned$ld_matrices[[i]])

    # Find these variants in the variant_indices dataframe
    variant_block_ids <- partitioned$variant_indices$block_id[
      match(block_variants, partitioned$variant_indices$variant_id)]

    # All variants should be mapped to this block
    expect_true(all(variant_block_ids == i))
  }

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix handles row/column name mismatches", {
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Create a plain list version with mismatched rownames and colnames
  ldmat <- getCorrelation(ld_data)
  vids <- getVariantIds(ld_data)
  rownames(ldmat) <- NULL
  colnames(ldmat) <- NULL
  mismatched_ld_data <- list(
    LD_matrix = ldmat,
    LD_variants = vids,
    block_metadata = getBlockMetadata(ld_data)
  )

  # Should not error and should fix the names
  partitioned <- partition_LD_matrix(mismatched_ld_data)

  # Check if names are fixed
  expect_identical(rownames(partitioned$ld_matrices[[1]]), getVariantIds(ld_data))
  expect_identical(colnames(partitioned$ld_matrices[[1]]), getVariantIds(ld_data))

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix correctly extracts blocks based on metadata", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")

  # Load the LD matrix
  ld_data <- load_LD_matrix(LD_meta_file_path, region)

  # Partition without merging
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = FALSE)

  # For each block, check if the extracted matrix matches the expected submatrix
  ld_variants <- getVariantIds(ld_data)
  ld_matrix <- getCorrelation(ld_data)
  for(i in seq_along(partitioned$ld_matrices)) {
    block_info <- partitioned$block_metadata[i, ]
    start_idx <- block_info$start_idx
    end_idx <- block_info$end_idx

    # Skip if indices are invalid
    if(start_idx > length(ld_variants) ||
       end_idx > length(ld_variants) ||
       end_idx < start_idx) next

    # Get variants for this block
    block_variants <- ld_variants[start_idx:end_idx]

    # Extract expected submatrix
    expected_submatrix <- ld_matrix[block_variants, block_variants, drop = FALSE]

    # Compare with actual block matrix
    expect_equal(partitioned$ld_matrices[[i]], expected_submatrix)
  }

  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix partitions correctly with synthetic data", {
  mat <- matrix(0, 6, 6)
  mat[1:3, 1:3] <- 0.5
  mat[4:6, 4:6] <- 0.5
  diag(mat) <- 1
  variant_ids <- paste0("v", 1:6)
  rownames(mat) <- colnames(mat) <- variant_ids

  ld_data <- list(
    LD_matrix = mat,
    LD_variants = variant_ids,
    block_metadata = data.frame(
      block_id = c(1, 2),
      chrom = c("1", "1"),
      size = c(3, 3),
      start_idx = c(1, 4),
      end_idx = c(3, 6)
    )
  )

  result <- pecotmr:::partition_LD_matrix(ld_data, merge_small_blocks = FALSE)

  expect_type(result, "list")
  expect_true("ld_matrices" %in% names(result))
  expect_true("variant_indices" %in% names(result))
  expect_length(result$ld_matrices, 2)
  expect_equal(nrow(result$ld_matrices[[1]]), 3)
  expect_equal(nrow(result$ld_matrices[[2]]), 3)
})

# ---- order_dedup_regions ----

test_that("order_dedup_regions removes duplicate regions", {
  # Create regions with duplicates
  regions_with_dups <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    start = c(100, 100, 200),  # Note: first two rows are duplicates
    end = c(150, 150, 250)
  )

  result <- order_dedup_regions(regions_with_dups)
  # Should have removed duplicate and return only two rows
  expect_equal(nrow(result), 2)
  expect_equal(result$start, c(100, 200))
})

test_that("order_dedup_regions orders and deduplicates across chromosomes", {
  df <- data.frame(
    chrom = c("chr2", "chr1", "chr1", "chr2"),
    start = c(100, 200, 100, 100),
    end = c(200, 300, 200, 200)
  )
  result <- pecotmr:::order_dedup_regions(df)
  expect_equal(nrow(result), 3)  # one duplicate removed
  expect_true(all(diff(result$start[result$chrom == result$chrom[1]]) >= 0))
})

test_that("order_dedup_regions strips chr prefix", {
  df <- data.frame(chrom = c("chr1", "chr2"), start = c(100, 200), end = c(200, 300))
  result <- pecotmr:::order_dedup_regions(df)
  expect_true(all(result$chrom %in% c(1L, 2L)))
})

# ---- find_intersection_rows ----

test_that("find_intersection_rows correctly identifies start and end rows", {
  # Create a simple genomic dataset
  genomic_data <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(100, 200, 300, 400),
    end = c(150, 250, 350, 450)
  )

  # Region entirely within the dataset
  result <- find_intersection_rows(genomic_data, 1, 220, 330)
  expect_equal(result$start_row$start, 200)
  expect_equal(result$end_row$end, 350)
})

test_that("find_intersection_rows adjusts region bounds if needed", {
  # Create a simple genomic dataset
  genomic_data <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(100, 200, 300, 400),
    end = c(150, 250, 350, 450)
  )

  # Region extends beyond the dataset
  result <- find_intersection_rows(genomic_data, 1, 50, 500)
  # Should adjust to the bounds of the dataset
  expect_equal(result$start_row$start, 100)
  expect_equal(result$end_row$end, 450)
})

test_that("find_intersection_rows errors for non-overlapping regions", {
  # Create a simple genomic dataset
  genomic_data <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(100, 200, 300, 400),
    end = c(150, 250, 350, 450)
  )

  # Region entirely outside the dataset
  expect_error(
    find_intersection_rows(genomic_data, 2, 100, 200),
    "No data for chromosome 2"
  )
})

# ---- validate_selected_region ----

test_that("validate_selected_region passes for valid region", {
  start_row <- data.frame(start = 0)
  end_row <- data.frame(end = 300)
  expect_silent(pecotmr:::validate_selected_region(start_row, end_row, 50, 250))
})

test_that("validate_selected_region errors for uncovered region", {
  start_row <- data.frame(start = 100)
  end_row <- data.frame(end = 200)
  expect_error(
    pecotmr:::validate_selected_region(start_row, end_row, 50, 250),
    "not fully covered"
  )
})

# ---- extract_file_paths ----

test_that("extract_file_paths extracts correct paths", {
  gd <- data.frame(
    chrom = c(1, 1, 1),
    start = c(0, 100, 200),
    end = c(100, 200, 300),
    path = c("f1.ld", "f2.ld", "f3.ld")
  )
  intersection <- list(
    start_row = data.frame(chrom = 1, start = 0),
    end_row = data.frame(start = 200)
  )
  result <- pecotmr:::extract_file_paths(gd, intersection, "path")
  expect_equal(length(result), 3)
})

test_that("extract_file_paths errors on missing column", {
  gd <- data.frame(chrom = 1, start = 0, end = 100)
  intersection <- list(
    start_row = data.frame(chrom = 1, start = 0),
    end_row = data.frame(start = 0)
  )
  expect_error(pecotmr:::extract_file_paths(gd, intersection, "nonexistent"),
               "not found")
})

# ---- partition_LD_matrix: different chromosomes ----

test_that("partition_LD_matrix handles blocks with different chromosomes", {
  # Create test data with blocks on different chromosomes
  test_matrix <- matrix(0, 4, 4)
  diag(test_matrix) <- 1  # Set diagonal to 1
  rownames(test_matrix) <- colnames(test_matrix) <- c("1:100:A:G", "1:200:C:T", "2:100:G:A", "2:200:T:C")

  block_metadata <- data.frame(
    block_id = 1:2,
    chrom = c(1, 2),
    size = c(2, 2),
    start_idx = c(1, 3),
    end_idx = c(2, 4)
  )

  test_ld_data <- list(
    LD_matrix = test_matrix,
    LD_variants = c("1:100:A:G", "1:200:C:T", "2:100:G:A", "2:200:T:C"),
    block_metadata = block_metadata
  )

  # Partition the matrix
  partitioned <- partition_LD_matrix(test_ld_data)

  # Should not merge blocks from different chromosomes
  expect_equal(length(partitioned$ld_matrices), 2)

  # Each block should have the correct variants
  expect_equal(rownames(partitioned$ld_matrices[[1]]), c("1:100:A:G", "1:200:C:T"))
  expect_equal(rownames(partitioned$ld_matrices[[2]]), c("2:100:G:A", "2:200:T:C"))
})

test_that("partition_LD_matrix works with edge case block structures", {
  # Test case: One large block and several tiny blocks that need merging
  large_block_size <- 15
  small_block_size <- 2

  # Create a matrix with blocks of varying sizes
  n_variants <- large_block_size + small_block_size * 3
  test_matrix <- matrix(0, n_variants, n_variants)
  # Set diagonal to 1
  diag(test_matrix) <- 1

  # Generate variant names
  variant_names <- paste0("1:", 100:(100+n_variants-1), ":A:G")
  rownames(test_matrix) <- colnames(test_matrix) <- variant_names

  # Create block metadata
  block_metadata <- data.frame(
    block_id = 1:4,
    chrom = rep(1, 4),
    size = c(large_block_size, small_block_size, small_block_size, small_block_size),
    start_idx = c(1, large_block_size+1, large_block_size+small_block_size+1, large_block_size+small_block_size*2+1),
    end_idx = c(large_block_size, large_block_size+small_block_size, large_block_size+small_block_size*2, n_variants)
  )

  test_ld_data <- list(
    LD_matrix = test_matrix,
    LD_variants = variant_names,
    block_metadata = block_metadata
  )

  # Set minimum block size to force merging of small blocks
  min_merged_size <- small_block_size + 1

  # Partition with merging
  partitioned <- partition_LD_matrix(test_ld_data, merge_small_blocks = TRUE,
                                    min_merged_block_size = min_merged_size)

  # Should merge the small blocks but leave the large block alone
  expect_lt(length(partitioned$ld_matrices), 4)
  expect_gt(length(partitioned$ld_matrices), 1)

  # First block should still be large_block_size
  expect_equal(nrow(partitioned$ld_matrices[[1]]), large_block_size)
})

# ---- extract_LD_for_region ----

test_that("extract_LD_for_region extracts correct region", {
  # Create mock LD matrix and variants
  ld_variants <- data.frame(
    chrom = c(1, 1, 1, 1),
    variants = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
    GD = NA,
    pos = c(100, 200, 300, 400),
    A1 = c("A", "C", "G", "T"),
    A2 = c("G", "T", "A", "C")
  )

  ld_matrix <- matrix(0, 4, 4)
  diag(ld_matrix) <- 1
  rownames(ld_matrix) <- colnames(ld_matrix) <- ld_variants$variants

  # Define a region that should include the middle two variants
  region <- data.frame(
    chrom = 1,
    start = 180,
    end = 320
  )

  result <- extract_LD_for_region(ld_matrix, ld_variants, region, NULL)

  # Should have extracted only the relevant variants
  expect_equal(nrow(result$extracted_LD_variants), 2)
  expect_equal(result$extracted_LD_variants$variants, c("1:200:C:T", "1:300:G:A"))

  # Matrix should be 2x2 with the correct row/column names
  expect_equal(dim(result$extracted_LD_matrix), c(2, 2))
  expect_equal(rownames(result$extracted_LD_matrix), c("1:200:C:T", "1:300:G:A"))
})

test_that("extract_LD_for_region works with extract_coordinates", {
  # Create mock LD matrix and variants
  ld_variants <- data.frame(
    chrom = c(1, 1, 1, 1),
    variants = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
    GD = NA,
    pos = c(100, 200, 300, 400),
    A1 = c("A", "C", "G", "T"),
    A2 = c("G", "T", "A", "C")
  )

  ld_matrix <- matrix(0, 4, 4)
  diag(ld_matrix) <- 1
  rownames(ld_matrix) <- colnames(ld_matrix) <- ld_variants$variants

  # Define a region that should include all variants
  region <- data.frame(
    chrom = 1,
    start = 50,
    end = 450
  )

  # Define specific coordinates to extract
  extract_coordinates <- data.frame(
    chrom = c(1, 1),
    pos = c(100, 300)
  )

  result <- extract_LD_for_region(ld_matrix, ld_variants, region, extract_coordinates)

  # Should have extracted only the specified coordinates
  expect_equal(nrow(result$extracted_LD_variants), 2)
  expect_equal(result$extracted_LD_variants$variants, c("1:100:A:G", "1:300:G:A"))

  # Matrix should be 2x2 with the correct row/column names
  expect_equal(dim(result$extracted_LD_matrix), c(2, 2))
  expect_equal(rownames(result$extracted_LD_matrix), c("1:100:A:G", "1:300:G:A"))
})

# ---- create_LD_matrix ----

test_that("create_LD_matrix correctly combines matrices with overlapping variants", {
  # Create two simple LD matrices with some overlapping variants
  matrix1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(matrix1) <- colnames(matrix1) <- c("1:100:A:G", "1:200:C:T")

  matrix2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  rownames(matrix2) <- colnames(matrix2) <- c("1:200:C:T", "1:300:G:A")

  # Create variants lists
  variants1 <- data.frame(variants = c("1:100:A:G", "1:200:C:T"))
  variants2 <- data.frame(variants = c("1:200:C:T", "1:300:G:A"))

  # Combine matrices
  combined <- create_LD_matrix(
    LD_matrices = list(matrix1, matrix2),
    variants = list(variants1, variants2)
  )

  # Should have created a 3x3 matrix with all unique variants
  expect_equal(dim(combined), c(3, 3))
  expect_equal(rownames(combined), c("1:100:A:G", "1:200:C:T", "1:300:G:A"))

  # Check that values from original matrices are preserved
  expect_equal(combined["1:100:A:G", "1:200:C:T"], 0.5)
  expect_equal(combined["1:200:C:T", "1:300:G:A"], 0.3)

  # Check diagonal values are 1
  expect_equal(combined[1,1], 1)
  expect_equal(combined[2,2], 1)
  expect_equal(combined[3,3], 1)
})

test_that("create_LD_matrix merges non-overlapping blocks", {
  m1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("v1", "v2"), c("v1", "v2")))
  m2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2, dimnames = list(c("v3", "v4"), c("v3", "v4")))

  variants <- list(
    data.frame(variants = c("v1", "v2")),
    data.frame(variants = c("v3", "v4"))
  )
  result <- pecotmr:::create_LD_matrix(list(m1, m2), variants)

  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
  expect_equal(result["v1", "v2"], 0.5)
  expect_equal(result["v3", "v4"], 0.3)
  expect_equal(result["v1", "v3"], 0)  # Cross-block should be 0
})

test_that("create_LD_matrix handles overlapping boundary variant", {
  m1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("v1", "v2"), c("v1", "v2")))
  m2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2, dimnames = list(c("v2", "v3"), c("v2", "v3")))

  variants <- list(
    data.frame(variants = c("v1", "v2")),
    data.frame(variants = c("v2", "v3"))
  )
  result <- pecotmr:::create_LD_matrix(list(m1, m2), variants)

  # v2 is shared, so total should be 3 variants
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
})

# ---- validate_block_structure ----

test_that("validate_block_structure passes for proper block structure", {
  mat <- matrix(0, 6, 6)
  mat[1:3, 1:3] <- 0.5
  mat[4:6, 4:6] <- 0.5
  diag(mat) <- 1

  variant_ids <- paste0("v", 1:6)
  rownames(mat) <- colnames(mat) <- variant_ids

  block_meta <- data.frame(
    block_id = c(1, 2),
    chrom = c("1", "1"),
    size = c(3, 3),
    start_idx = c(1, 4),
    end_idx = c(3, 6)
  )

  expect_silent(pecotmr:::validate_block_structure(mat, block_meta, variant_ids))
})

test_that("validate_block_structure errors on non-block structure", {
  mat <- matrix(0.5, 4, 4)
  diag(mat) <- 1

  variant_ids <- paste0("v", 1:4)
  rownames(mat) <- colnames(mat) <- variant_ids

  block_meta <- data.frame(
    block_id = c(1, 2),
    chrom = c("1", "1"),
    size = c(2, 2),
    start_idx = c(1, 3),
    end_idx = c(2, 4)
  )

  expect_error(pecotmr:::validate_block_structure(mat, block_meta, variant_ids),
               "Matrix lacks expected block structure")
})

# ---- merge_blocks ----

test_that("merge_blocks properly handles blocks at chromosome boundaries", {
  # Create test data with small blocks at chromosome boundaries
  test_matrix <- matrix(0, 6, 6)
  diag(test_matrix) <- 1
  variant_names <- c("1:900:A:G", "1:950:C:T", "2:100:G:A", "2:150:T:C", "3:100:A:G", "3:150:C:T")
  rownames(test_matrix) <- colnames(test_matrix) <- variant_names

  # Create block metadata with small blocks at chromosome boundaries
  block_metadata <- data.frame(
    block_id = 1:3,
    chrom = c(1, 2, 3),
    size = c(2, 2, 2),
    start_idx = c(1, 3, 5),
    end_idx = c(2, 4, 6)
  )

  test_ld_data <- list(
    LD_matrix = test_matrix,
    LD_variants = variant_names,
    block_metadata = block_metadata
  )

  # Set min block size to force merging attempts
  min_block_size <- 3

  # Partition with merging
  partitioned <- partition_LD_matrix(test_ld_data, merge_small_blocks = TRUE,
                                     min_merged_block_size = min_block_size)

  # Should not merge blocks across chromosome boundaries
  expect_equal(length(partitioned$ld_matrices), 3)

  # Each block should match its chromosome
  for (i in 1:3) {
    block_variants <- rownames(partitioned$ld_matrices[[i]])
    chrom_from_variants <- unique(as.integer(sub(":.*", "", block_variants)))
    expect_equal(length(chrom_from_variants), 1)  # Should only have one chromosome per block
    expect_equal(chrom_from_variants, i)  # Should match the expected chromosome
  }
})

test_that("merge_blocks merges small adjacent blocks", {
  block_meta <- data.frame(
    block_id = c(1, 2, 3),
    chrom = c("1", "1", "1"),
    size = c(50, 50, 100),
    start_idx = c(1, 51, 101),
    end_idx = c(50, 100, 200)
  )
  result <- pecotmr:::merge_blocks(block_meta, min_size = 100, max_size = 10000)
  expect_true(nrow(result) < 3)
})

test_that("merge_blocks does not merge cross-chromosome", {
  block_meta <- data.frame(
    block_id = c(1, 2),
    chrom = c("1", "2"),
    size = c(10, 10),
    start_idx = c(1, 11),
    end_idx = c(10, 20)
  )
  result <- pecotmr:::merge_blocks(block_meta, min_size = 50, max_size = 10000)
  expect_equal(nrow(result), 2)  # Cannot merge across chromosomes
})

test_that("merge_blocks returns single block unchanged", {
  block_meta <- data.frame(
    block_id = 1, chrom = "1", size = 10,
    start_idx = 1, end_idx = 10
  )
  result <- pecotmr:::merge_blocks(block_meta, min_size = 100, max_size = 10000)
  expect_equal(nrow(result), 1)
})

# ---- can_merge ----

test_that("can_merge checks chromosome and size", {
  b1 <- data.frame(chrom = "1", size = 100)
  b2 <- data.frame(chrom = "1", size = 200)
  expect_true(pecotmr:::can_merge(b1, b2, max_size = 500))
  expect_false(pecotmr:::can_merge(b1, b2, max_size = 200))

  b3 <- data.frame(chrom = "2", size = 100)
  expect_false(pecotmr:::can_merge(b1, b3, max_size = 500))
})

# ===========================================================================
# check_ld (regularize_ld)
# ===========================================================================

test_that("check_ld reports PD for identity matrix", {
  R <- diag(5)
  result <- check_ld(R)
  expect_true(result$is_pd)
  expect_true(result$is_psd)
  expect_equal(result$method_applied, "none")
  expect_equal(result$R, R)
  expect_equal(result$condition_number, 1)
})

test_that("check_ld reports PD for well-conditioned correlation matrix", {
  R <- matrix(0.3, 4, 4)
  diag(R) <- 1
  result <- check_ld(R)
  expect_true(result$is_pd)
  expect_true(result$is_psd)
  expect_equal(result$n_negative, 0)
  expect_equal(result$method_applied, "none")
})

test_that("check_ld detects non-PSD matrix", {
  R <- matrix(0.9, 3, 3)
  diag(R) <- 1
  R[1, 3] <- R[3, 1] <- -0.5
  result <- check_ld(R)
  expect_false(result$is_psd)
  expect_true(result$n_negative > 0)
  expect_true(result$min_eigenvalue < 0)
  expect_equal(result$method_applied, "none")
})

test_that("check_ld shrink method modifies non-PD matrix", {
  R <- matrix(0.9, 3, 3)
  diag(R) <- 1
  R[1, 3] <- R[3, 1] <- -0.5
  result <- check_ld(R, method = "shrink")
  expect_equal(result$method_applied, "shrink")
  expect_false(identical(result$R, R))
  # With strong enough shrinkage, result should be PD
  result2 <- check_ld(R, method = "shrink", shrinkage = 0.5)
  eig <- eigen(result2$R, symmetric = TRUE)
  expect_true(all(eig$values > 0))
})

test_that("check_ld eigenfix method improves non-PD matrix", {
  R <- matrix(0.9, 3, 3)
  diag(R) <- 1
  R[1, 3] <- R[3, 1] <- -0.5
  original_min_eval <- min(eigen(R, symmetric = TRUE)$values)
  result <- check_ld(R, method = "eigenfix")
  expect_equal(result$method_applied, "eigenfix")
  # Eigenfix should improve the minimum eigenvalue
  fixed_min_eval <- min(eigen(result$R, symmetric = TRUE)$values)
  expect_true(fixed_min_eval > original_min_eval)
  # Unit diagonal preserved
  expect_equal(diag(result$R), rep(1, 3))
  # Symmetry preserved
  expect_equal(result$R, t(result$R))
})

test_that("check_ld shrink does nothing when matrix is already PD", {
  R <- diag(3)
  result <- check_ld(R, method = "shrink")
  expect_equal(result$method_applied, "none")
  expect_equal(result$R, R)
})

test_that("check_ld eigenfix does nothing when matrix is already PD", {
  R <- diag(3)
  result <- check_ld(R, method = "eigenfix")
  expect_equal(result$method_applied, "none")
  expect_equal(result$R, R)
})

# ===========================================================================
# extract_block_matrices: out-of-range blocks
# ===========================================================================

test_that("extract_block_matrices warns and skips out-of-range blocks", {
  mat <- diag(4)
  vnames <- paste0("v", 1:4)
  rownames(mat) <- colnames(mat) <- vnames
  block_metadata <- data.frame(
    block_id = c(1, 2),
    start_idx = c(1, 10),
    end_idx = c(2, 12),
    chrom = c("1", "1"),
    block_start = c(100, 500),
    block_end = c(200, 600),
    size = c(2, 3),
    stringsAsFactors = FALSE
  )
  expect_warning(
    result <- pecotmr:::extract_block_matrices(mat, block_metadata, vnames),
    "outside the range"
  )
  valid_blocks <- result$ld_matrices[!sapply(result$ld_matrices, is.null)]
  expect_equal(length(valid_blocks), 1)
  expect_equal(nrow(valid_blocks[[1]]), 2)
})

# ===========================================================================
# resolve_ld_source: type detection with real fixtures
# ===========================================================================

geno_test_data_dir <- test_path("test_data")
geno_region_all <- "chr21:17513228-17592874"

test_that("resolve_ld_source detects PLINK2 from metadata", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_p2_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- pecotmr:::resolve_ld_source(meta_file)
  expect_equal(result$type, "plink2")
  expect_equal(result$meta_path, meta_file)
})

test_that("resolve_ld_source detects VCF from metadata", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_vcf_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- pecotmr:::resolve_ld_source(meta_file)
  expect_equal(result$type, "vcf")
})

test_that("resolve_ld_source detects GDS from metadata", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_gds_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.gds", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- pecotmr:::resolve_ld_source(meta_file)
  expect_equal(result$type, "gds")
})

test_that("resolve_ld_source detects precomputed from metadata", {
  # Existing LD block metadata with non-zero start/end
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_pre_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("chr1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- pecotmr:::resolve_ld_source(meta_file)
  expect_equal(result$type, "precomputed")
})

test_that("resolve_ld_source errors on missing file", {
  expect_error(pecotmr:::resolve_ld_source("/nonexistent/file.tsv"), "not found")
})

test_that("resolve_ld_source errors on 0:0 sentinel with non-genotype path", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_bad_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "nonexistent_prefix", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  expect_error(pecotmr:::resolve_ld_source(meta_file), "0:0 sentinel")
})

# ===========================================================================
# resolve_genotype_path_for_region
# ===========================================================================

test_that("resolve_genotype_path_for_region resolves correct chromosome path", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_path_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- pecotmr:::resolve_genotype_path_for_region(meta_file, geno_region_all)
  expect_equal(result, file.path(geno_test_data_dir, "test_variants"))
})

test_that("resolve_genotype_path_for_region errors on missing chromosome", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_nochr_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  expect_error(
    pecotmr:::resolve_genotype_path_for_region(meta_file, geno_region_all),
    "No entry for chromosome"
  )
})

# ===========================================================================
# load_LD_from_genotype with real fixtures
# ===========================================================================

test_that("load_LD_from_genotype returns LD matrix with .afreq", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(geno_test_data_dir, "test_variants")
  result <- pecotmr:::load_LD_from_genotype(plink_prefix, geno_region_all)
  expect_true(is(result, "LDData"))
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 349L)
  expect_true(isSymmetric(getCorrelation(result)))
  expect_false(hasGenotypes(result))
  # ref_panel should have allele_freq from .afreq file
  expect_true("allele_freq" %in% names(S4Vectors::mcols(getVariantInfo(result))))
  expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq > 0))
  expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq < 1))
  # block_metadata
  expect_true(is.data.frame(getBlockMetadata(result)))
  expect_equal(nrow(getBlockMetadata(result)), 1L)
})

test_that("load_LD_from_genotype returns genotype matrix when requested", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(geno_test_data_dir, "test_variants")
  result <- pecotmr:::load_LD_from_genotype(plink_prefix, geno_region_all,
                                             return_genotype = TRUE)
  expect_true(hasGenotypes(result))
  X <- getGenotypes(result)
  expect_equal(nrow(X), 100L)  # samples
  expect_equal(ncol(X), 349L)  # variants
})

test_that("load_LD_from_genotype computes variance with n_sample", {
  skip_if_not_installed("pgenlibr")
  plink_prefix <- file.path(geno_test_data_dir, "test_variants")
  result <- pecotmr:::load_LD_from_genotype(plink_prefix, geno_region_all,
                                             n_sample = 100L)
  expect_true("variance" %in% names(S4Vectors::mcols(getVariantInfo(result))))
  expect_true("n_nomiss" %in% names(S4Vectors::mcols(getVariantInfo(result))))
  expect_equal(S4Vectors::mcols(getVariantInfo(result))$n_nomiss[1], 100L)
  expect_true(all(S4Vectors::mcols(getVariantInfo(result))$variance > 0))
})

test_that("load_LD_from_genotype falls back to computed AF without .afreq", {
  skip_if_not_installed("VariantAnnotation")
  vcf_path <- file.path(geno_test_data_dir, "test_variants.vcf.gz")
  result <- suppressWarnings(
    pecotmr:::load_LD_from_genotype(vcf_path, geno_region_all)
  )
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 349L)
  expect_true(isSymmetric(getCorrelation(result)))
  # Allele frequencies computed from genotypes
  expect_true("allele_freq" %in% names(S4Vectors::mcols(getVariantInfo(result))))
  expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq > 0))
  expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq < 1))
})

test_that("load_LD_from_genotype works with GDS files", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  gds_path <- file.path(geno_test_data_dir, "test_variants.gds")
  result <- pecotmr:::load_LD_from_genotype(gds_path, geno_region_all)
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 349L)
  expect_true(isSymmetric(getCorrelation(result)))
})

test_that("load_LD_from_genotype .afreq and computed AF are consistent", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  plink_prefix <- file.path(geno_test_data_dir, "test_variants")
  gds_path <- file.path(geno_test_data_dir, "test_variants.gds")
  res_afreq <- pecotmr:::load_LD_from_genotype(plink_prefix, geno_region_all)
  res_computed <- pecotmr:::load_LD_from_genotype(gds_path, geno_region_all)
  # Allele frequencies should be close (same data, different source)
  expect_true(max(abs(S4Vectors::mcols(getVariantInfo(res_afreq))$allele_freq - S4Vectors::mcols(getVariantInfo(res_computed))$allele_freq)) < 0.01)
})

# ===========================================================================
# load_LD_matrix with real genotype fixtures via metadata
# ===========================================================================

test_that("load_LD_matrix dispatches to PLINK2 genotype source", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_p2_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- load_LD_matrix(meta_file, geno_region_all)
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 349L)
  expect_false(hasGenotypes(result))
})

test_that("load_LD_matrix dispatches to VCF genotype source", {
  skip_if_not_installed("VariantAnnotation")
  meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_vcf_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- suppressWarnings(load_LD_matrix(meta_file, geno_region_all))
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 349L)
})

test_that("load_LD_matrix dispatches to GDS genotype source", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_gds_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.gds", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- load_LD_matrix(meta_file, geno_region_all)
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 349L)
})

test_that("load_LD_matrix return_genotype='auto' returns X for genotype source", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_auto_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- load_LD_matrix(meta_file, geno_region_all, return_genotype = "auto")
  expect_true(hasGenotypes(result))
  X <- getGenotypes(result)
  expect_equal(nrow(X), 100L)  # samples
  expect_equal(ncol(X), 349L)  # variants
})

test_that("load_LD_matrix return_genotype=TRUE errors for precomputed", {
  meta_file <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".tsv"))
  on.exit(unlink(meta_file), add = TRUE)
  meta_df <- data.frame(
    chrom = "chr1", start = 1000, end = 1200,
    path = paste0(
      "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,",
      "./test_data/LD_block_1.chr1_1000_1200.float16.bim"
    )
  )
  write_delim(meta_df, meta_file, delim = "\t")
  region <- data.frame(chrom = "chr1", start = 1000, end = 1190)
  expect_error(
    load_LD_matrix(meta_file, region, return_genotype = TRUE),
    "genotype files"
  )
})

# ===========================================================================
# resolve_ld_source: PLINK1 detection
# ===========================================================================

test_that("resolve_ld_source detects PLINK1 from metadata", {
  skip_if_not_installed("snpStats")
  meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_p1_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("22", "0", "0", "protocol_example.genotype", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- pecotmr:::resolve_ld_source(meta_file)
  expect_equal(result$type, "plink1")
})

# ===========================================================================
# load_LD_matrix: precomputed blocks via real .cor.xz fixtures
# ===========================================================================

test_that("load_LD_matrix loads single precomputed block", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_single_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- load_LD_matrix(meta_file, "chr1:1000-1190")
  expect_true(is.matrix(getCorrelation(result)))
  expect_equal(nrow(getCorrelation(result)), 5L)
  expect_true(isSymmetric(getCorrelation(result)))
  expect_equal(length(getVariantIds(result)), 5L)
  expect_true(all(grepl("^chr1:", getVariantIds(result))))
  expect_false(hasGenotypes(result))
  # block_metadata should have one block
  expect_equal(nrow(getBlockMetadata(result)), 1L)
  # ref_panel (now GRanges) should have variant info via mcols
  ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
  expect_true("variant_id" %in% names(ref_mcols))
})

test_that("load_LD_matrix loads multiple precomputed blocks", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_multi_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  lines <- c(
    paste("chrom", "start", "end", "path", sep = "\t"),
    paste("1", "1000", "1200",
          "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
          sep = "\t"),
    paste("1", "1200", "1400",
          "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
          sep = "\t"),
    paste("1", "1400", "1600",
          "LD_block_3.chr1_1400_1600.float16.txt.xz,LD_block_3.chr1_1400_1600.float16.bim",
          sep = "\t")
  )
  writeLines(lines, meta_file)
  result <- load_LD_matrix(meta_file, "chr1:1000-1500")
  expect_true(is.matrix(getCorrelation(result)))
  # Should span blocks 1-3: 5 + 5 + 5 = 15 unique variants (no overlap in variant IDs)
  expect_true(nrow(getCorrelation(result)) >= 10)
  expect_true(isSymmetric(getCorrelation(result)))
  expect_true(nrow(getBlockMetadata(result)) >= 2)
})

test_that("load_LD_matrix with n_sample for precomputed blocks with freq data", {
  # The 9-column bim format includes allele_freq, variance, n_nomiss;
  # the 6-column bim does not. With 6-col bim and no allele_freq,
  # n_sample cannot compute variance - ref_panel has base columns only.
  meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_nsamp_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- load_LD_matrix(meta_file, "chr1:1000-1190", n_sample = 500L)
  # ref_panel (GRanges) should always have basic variant info in mcols
  ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
  expect_true(all(c("A2", "A1", "variant_id") %in% names(ref_mcols)))
  # 6-col bim lacks allele_freq so variance computation is skipped
  expect_false("variance" %in% names(ref_mcols))
})

# ===========================================================================
# process_LD_matrix: real .cor.xz fixtures
# ===========================================================================

test_that("process_LD_matrix reads .cor.xz with explicit bim path", {
  ld_file <- file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.txt.xz")
  bim_file <- file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.bim")
  result <- pecotmr:::process_LD_matrix(ld_file, bim_file)
  expect_true(is.list(result))
  expect_true(is.matrix(result$LD_matrix))
  expect_equal(nrow(result$LD_matrix), 5L)
  expect_equal(ncol(result$LD_matrix), 5L)
  expect_true(isSymmetric(result$LD_matrix))
  # Diagonal should be 1
  expect_true(all(abs(diag(result$LD_matrix) - 1) < 1e-4))
  # Variant names should be chr:pos:A2:A1 format
  expect_true(all(grepl("^chr1:", rownames(result$LD_matrix))))
  # LD_variants data frame
  expect_true(is.data.frame(result$LD_variants))
  expect_true("variants" %in% names(result$LD_variants))
  expect_equal(nrow(result$LD_variants), 5L)
})

test_that("process_LD_matrix reads different blocks consistently", {
  bim1 <- file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.bim")
  bim2 <- file.path(geno_test_data_dir, "LD_block_2.chr1_1200_1400.float16.bim")
  ld1 <- file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.txt.xz")
  ld2 <- file.path(geno_test_data_dir, "LD_block_2.chr1_1200_1400.float16.txt.xz")
  r1 <- pecotmr:::process_LD_matrix(ld1, bim1)
  r2 <- pecotmr:::process_LD_matrix(ld2, bim2)
  # Different blocks should have different variant positions
  expect_false(any(rownames(r1$LD_matrix) %in% rownames(r2$LD_matrix)))
})

test_that("process_LD_matrix reads 9-column bim with allele_freq/variance/n_nomiss", {
  ld_file <- file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.txt.xz")
  bim_file <- file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.9col.bim")
  result <- pecotmr:::process_LD_matrix(ld_file, bim_file)
  expect_equal(nrow(result$LD_matrix), 5L)
  expect_true(isSymmetric(result$LD_matrix))
  # 9-column bim should include extra columns
  expect_true("allele_freq" %in% names(result$LD_variants))
  expect_true("variance" %in% names(result$LD_variants))
  expect_true("n_nomiss" %in% names(result$LD_variants))
  expect_equal(result$LD_variants$allele_freq, c(0.3, 0.4, 0.2, 0.5, 0.15))
  expect_equal(result$LD_variants$n_nomiss, rep(500, 5))
})

test_that("load_LD_matrix propagates allele_freq/variance/n_nomiss from 9-col bim", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_9col_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.9col.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- load_LD_matrix(meta_file, "chr1:1000-1190")
  # ref_panel (GRanges) should carry the extra columns from the 9-col bim
  ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
  expect_true("allele_freq" %in% names(ref_mcols))
  expect_true("variance" %in% names(ref_mcols))
  expect_true("n_nomiss" %in% names(ref_mcols))
  expect_equal(ref_mcols$allele_freq, c(0.3, 0.4, 0.2, 0.5, 0.15))
  expect_equal(ref_mcols$n_nomiss, rep(500, 5))
})

# ===========================================================================
# get_regional_ld_meta: real .cor.xz fixtures
# ===========================================================================

test_that("get_regional_ld_meta returns correct file paths for single block", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_regional_single_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  lines <- c(
    paste("chrom", "start", "end", "path", sep = "\t"),
    paste("1", "1000", "1200",
          "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
          sep = "\t"),
    paste("1", "1200", "1400",
          "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
          sep = "\t")
  )
  writeLines(lines, meta_file)
  result <- pecotmr:::get_regional_ld_meta(meta_file, "chr1:1050-1150")
  expect_true(is.list(result))
  expect_true(length(result$intersections$LD_file_paths) >= 1)
  # All returned paths should exist
  expect_true(all(file.exists(result$intersections$LD_file_paths)))
  expect_true(all(file.exists(result$intersections$bim_file_paths)))
})

test_that("get_regional_ld_meta spans multiple blocks for wide region", {
  meta_file <- file.path(geno_test_data_dir, "ld_meta_regional_multi_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  lines <- c(
    paste("chrom", "start", "end", "path", sep = "\t"),
    paste("1", "1000", "1200",
          "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
          sep = "\t"),
    paste("1", "1200", "1400",
          "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
          sep = "\t"),
    paste("1", "1400", "1600",
          "LD_block_3.chr1_1400_1600.float16.txt.xz,LD_block_3.chr1_1400_1600.float16.bim",
          sep = "\t")
  )
  writeLines(lines, meta_file)
  result <- pecotmr:::get_regional_ld_meta(meta_file, "chr1:1000-1500")
  expect_true(length(result$intersections$LD_file_paths) >= 2)
})

# ===========================================================================
# drop_collinear_columns: strategy variants
# ===========================================================================

test_that("drop_collinear_columns variance strategy removes lowest-variance column", {
  set.seed(42)
  X <- matrix(rnorm(100 * 4), 100, 4)
  colnames(X) <- c("a", "b", "c", "d")
  # Make column "c" have near-zero variance
  X[, "c"] <- X[1, "c"]
  result <- pecotmr:::drop_collinear_columns(
    X, c("b", "c", "d"), strategy = "variance"
  )
  expect_false("c" %in% colnames(result))
  expect_equal(ncol(result), 3L)
})

test_that("drop_collinear_columns response_correlation strategy works", {
  set.seed(42)
  X <- matrix(rnorm(100 * 3), 100, 3)
  colnames(X) <- c("a", "b", "c")
  y <- X[, "a"] + rnorm(100, sd = 0.1)  # y correlates strongly with "a"
  result <- pecotmr:::drop_collinear_columns(
    X, c("a", "b", "c"), strategy = "response_correlation", response = y
  )
  # Should keep "a" (highest |cor| with response) and remove one of b/c
  expect_true("a" %in% colnames(result))
  expect_equal(ncol(result), 2L)
})

test_that("drop_collinear_columns response_correlation errors without response", {
  X <- matrix(1:12, 4, 3)
  colnames(X) <- c("a", "b", "c")
  expect_error(
    pecotmr:::drop_collinear_columns(X, c("a", "b"), strategy = "response_correlation"),
    "response must be supplied"
  )
})

test_that("drop_collinear_columns with single problematic column removes it", {
  X <- matrix(rnorm(40), 10, 4)
  colnames(X) <- c("a", "b", "c", "d")
  result <- pecotmr:::drop_collinear_columns(X, "b", strategy = "correlation")
  expect_false("b" %in% colnames(result))
  expect_equal(ncol(result), 3L)
})

# ===========================================================================
# enforce_design_full_rank: additional strategies and fallback paths
# ===========================================================================

test_that("enforce_design_full_rank variance strategy produces full rank", {
  set.seed(42)
  X <- matrix(rnorm(100 * 4), 100, 4)
  X[, 4] <- X[, 1] + X[, 2]  # rank deficient
  colnames(X) <- c("a", "b", "c", "d")
  C <- matrix(rnorm(100), 100, 1)
  result <- enforce_design_full_rank(X, C, strategy = "variance")
  full_design <- cbind(1, result, C)
  expect_equal(qr(full_design)$rank, ncol(full_design))
  expect_true(ncol(result) < ncol(X))
})

test_that("enforce_design_full_rank response_correlation strategy works", {
  set.seed(42)
  X <- matrix(rnorm(100 * 4), 100, 4)
  X[, 4] <- X[, 1] + X[, 2]
  colnames(X) <- c("a", "b", "c", "d")
  C <- matrix(rnorm(100), 100, 1)
  y <- X[, "a"] + rnorm(100, sd = 0.1)
  result <- enforce_design_full_rank(X, C, strategy = "response_correlation", response = y)
  full_design <- cbind(1, result, C)
  expect_equal(qr(full_design)$rank, ncol(full_design))
})

test_that("enforce_design_full_rank returns unchanged X when already full rank", {
  set.seed(42)
  X <- matrix(rnorm(100 * 3), 100, 3)
  colnames(X) <- c("a", "b", "c")
  C <- matrix(rnorm(100), 100, 1)
  result <- enforce_design_full_rank(X, C, strategy = "correlation")
  expect_equal(ncol(result), ncol(X))
})

test_that("enforce_design_full_rank fallback to correlation pruning works", {
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * 3), n, 3)
  # Create highly collinear columns that are hard for iterative removal
  X <- cbind(X, X[, 1] + rnorm(n, sd = 1e-10),
                 X[, 2] + rnorm(n, sd = 1e-10),
                 X[, 3] + rnorm(n, sd = 1e-10),
                 X[, 1] + X[, 2] + rnorm(n, sd = 1e-10))
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  C <- matrix(rnorm(n), n, 1)
  result <- enforce_design_full_rank(X, C, strategy = "correlation",
                                      max_iterations = 2L)
  full_design <- cbind(1, result, C)
  expect_equal(qr(full_design)$rank, ncol(full_design))
})

# ===========================================================================
# ld_clump_by_score: edge cases
# ===========================================================================

test_that("ld_clump_by_score errors on empty matrix", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  X <- matrix(numeric(0), nrow = 10, ncol = 0)
  expect_error(ld_clump_by_score(X, score = numeric(0), chr = integer(0), pos = integer(0)),
               "at least one column")
})

test_that("ld_clump_by_score returns 1L for single variant", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  X <- matrix(c(0, 1, 2, 1, 0), ncol = 1)
  result <- ld_clump_by_score(X, score = 1.0, chr = 1L, pos = 100L)
  expect_equal(result, 1L)
})

test_that("ld_clump_by_score errors on mismatched score length", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  X <- matrix(rnorm(20), 5, 4)
  expect_error(ld_clump_by_score(X, score = c(1, 2), chr = rep(1L, 4), pos = 1:4),
               "length\\(score\\)")
})

test_that("ld_clump_by_score errors on mismatched chr/pos length", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  X <- matrix(rnorm(20), 5, 4)
  expect_error(ld_clump_by_score(X, score = runif(4), chr = rep(1L, 2), pos = 1:4),
               "chr and pos")
})

# ===========================================================================
# ld_prune_by_correlation: verbose paths
# ===========================================================================

test_that("ld_prune_by_correlation verbose reports pruning", {
  # Create matrix with correlated columns
  set.seed(42)
  base <- rnorm(100)
  X <- cbind(base, base + rnorm(100, sd = 0.1), rnorm(100), rnorm(100), rnorm(100))
  colnames(X) <- paste0("v", 1:5)
  expect_message(
    ld_prune_by_correlation(X, cor_thres = 0.5, verbose = TRUE),
    "pruned"
  )
})

test_that("ld_prune_by_correlation verbose reports no pruning", {
  # Create a small matrix with no correlated columns
  set.seed(42)
  X <- matrix(rnorm(500), 100, 5)
  colnames(X) <- paste0("v", 1:5)
  expect_message(
    ld_prune_by_correlation(X, cor_thres = 0.999, verbose = TRUE),
    "no columns pruned"
  )
})
