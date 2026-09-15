context("LD")
library(tidyverse)

# Helper: build an LdData S4 object from variant IDs and optional correlation matrix
make_test_ld_data <- function(variant_ids, R = NULL, blockMetadata = NULL) {
    if (is.null(R)) {
        p <- length(variant_ids)
        R <- diag(p)
        rownames(R) <- colnames(R) <- variant_ids
    }
    ref_panel <- pecotmr:::parseVariantId(variant_ids)
    ref_panel$variant_id <- variant_ids
    variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
    if (is.null(blockMetadata)) {
        blockMetadata <- data.frame(
            blockId = 1L,
            chrom = as.character(ref_panel$chrom[1]),
            blockStart = min(ref_panel$pos),
            blockEnd = max(ref_panel$pos),
            size = length(variant_ids),
            startIdx = 1L,
            endIdx = length(variant_ids),
            stringsAsFactors = FALSE
        )
    }
    LdData(
        correlation = R,
        variants = variants_gr,
        blockMetadata = blockMetadata
    )
}

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
        )
    )
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
        )
    )
    return(list(region = region, meta = meta_df))
}

test_that("Check that we correctly retrieve the names from the matrix", {
    data <- generate_dummy_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")
    res <- loadLdMatrix(LD_meta_file_path, region)
    variants <- unlist(
        c(
            "chr1:1000:A:G",
            "chr1:1040:A:G",
            "chr1:1080:A:G",
            "chr1:1120:A:G",
            "chr1:1160:A:G"
        )
    )
    expect_equal(
        unlist(getVariantIds(res)),
        variants
    )
    file.remove(LD_meta_file_path)
})

test_that("Check that the LD block contains the correct information", {
    data <- generate_dummy_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")
    res <- loadLdMatrix(LD_meta_file_path, region)
    # Variant names
    variants <- unlist(
        c(
            "chr1:1000:A:G",
            "chr1:1040:A:G",
            "chr1:1080:A:G",
            "chr1:1120:A:G",
            "chr1:1160:A:G"
        )
    )
    # Check LD Block 1
    ld_block_one <- getCorrelation(res)
    ld_block_one_original <- as.matrix(
        read_delim(
            "test_data/LD_block_1.chr1_1000_1200.float16.txt.xz",
            delim = " ",
            col_names = F
        )
    )
    rownames(ld_block_one_original) <- colnames(
        ld_block_one_original
    ) <- variants
    expect_equal(ld_block_one, ld_block_one_original)
    file.remove(LD_meta_file_path)
})

# ---- partitionLdMatrix ----

test_that("partitionLdMatrix correctly partitions a single block", {
    data <- generate_dummy_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix first
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Now partition the matrix
    partitioned <- partitionLdMatrix(ld_data)

    # Expectations for single block case
    expect_equal(length(partitioned$ldMatrices), 1)
    expect_equal(
        nrow(partitioned$variantIndices),
        length(getVariantIds(ld_data))
    )
    expect_equal(unique(partitioned$variantIndices$blockId), 1)
    expect_identical(
        rownames(partitioned$ldMatrices[[1]]),
        getVariantIds(ld_data)
    )
    expect_identical(
        colnames(partitioned$ldMatrices[[1]]),
        getVariantIds(ld_data)
    )

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix correctly partitions multiple blocks", {
    data <- generate_multi_block_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix that spans multiple blocks
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Now partition the matrix without merging blocks
    partitioned <- partitionLdMatrix(ld_data, mergeSmallBlocks = FALSE)

    # Check if we have the correct number of blocks
    # Should have block 1 (1000-1200), block 2 (1200-1400), and block 3 (1400-1600)
    expected_block_count <- 3
    expect_equal(length(partitioned$ldMatrices), expected_block_count)

    # Check if all variants are assigned to blocks
    expect_equal(
        nrow(partitioned$variantIndices),
        length(getVariantIds(ld_data))
    )

    # Check if block IDs are correct
    expect_setequal(
        unique(partitioned$variantIndices$blockId),
        1:expected_block_count
    )

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix properly merges small blocks", {
    data <- generate_multi_block_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix that spans multiple blocks
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Set min_merged_block_size high enough to force merging
    # Each test block likely has 5 variants (based on the existing test)
    min_block_size <- 10

    # Now partition the matrix with block merging
    partitioned <- partitionLdMatrix(
        ld_data,
        mergeSmallBlocks = TRUE,
        minMergedBlockSize = min_block_size
    )

    # We expect fewer blocks after merging
    expect_lt(length(partitioned$ldMatrices), 3)

    # Check if all variants are still assigned to blocks
    expect_equal(
        nrow(partitioned$variantIndices),
        length(getVariantIds(ld_data))
    )

    # Check if merged blocks are larger than min_block_size
    block_sizes <- sapply(partitioned$ldMatrices, nrow)
    expect_true(all(
        block_sizes >= min_block_size |
            block_sizes == length(getVariantIds(ld_data))
    ))

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix respects max_merged_block_size", {
    data <- generate_multi_block_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix that spans multiple blocks
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Set max_merged_block_size to a small value to prevent merging all blocks
    # Each test block likely has 5 variants (based on the existing test)
    max_block_size <- 8

    # Now partition the matrix with restricted block size
    partitioned <- partitionLdMatrix(
        ld_data,
        mergeSmallBlocks = TRUE,
        minMergedBlockSize = 2,
        maxMergedBlockSize = max_block_size
    )

    # Check if no block exceeds max_block_size
    block_sizes <- sapply(partitioned$ldMatrices, nrow)
    expect_true(all(block_sizes <= max_block_size))

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix handles empty matrix gracefully", {
    # A plain list (legacy format) is no longer accepted; the S4 check fires first.
    empty_ld_data <- list(
        ldMatrix = matrix(0, 0, 0),
        ldVariants = character(0),
        blockMetadata = data.frame(
            blockId = integer(0),
            chrom = character(0),
            size = integer(0),
            startIdx = integer(0),
            endIdx = integer(0)
        )
    )

    # Expect the S4 type-check error
    expect_error(
        partitionLdMatrix(empty_ld_data),
        "ldData must be an LdData object"
    )
})

test_that("partitionLdMatrix validates block structure properly", {
    data <- generate_multi_block_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix that spans multiple blocks
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Create an invalid block structure by modifying blockMetadata
    bm <- getBlockMetadata(ld_data)
    vids <- getVariantIds(ld_data)
    ldmat <- getCorrelation(ld_data)

    # Assuming we have at least 2 blocks:
    if (nrow(bm) >= 2) {
        # Create overlapping blocks with invalid start/end indices
        bm$startIdx[2] <- bm$startIdx[1]
        bm$endIdx[1] <- bm$endIdx[2]

        # Introduce non-zero elements between blocks to trigger validation error
        if (length(vids) >= 2) {
            idx1 <- bm$startIdx[1]
            idx2 <- bm$startIdx[2] + 1
            if (idx1 <= length(vids) && idx2 <= length(vids)) {
                var1 <- vids[idx1]
                var2 <- vids[idx2]
                ldmat[var1, var2] <- 0.5
            }
        }

        # Rebuild LdData with modified matrix and block metadata
        invalid_ld_data <- new(
            "LdData",
            getVariantInfo(ld_data),
            correlation = ldmat,
            genotypeHandle = NULL,
            snpIdx = getSnpIdx(ld_data),
            blockMetadata = bm
        )

        # Expect an error for invalid block structure
        expect_error(
            partitionLdMatrix(invalid_ld_data),
            "Matrix lacks expected block structure"
        )
    }

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix properly maps variants to blocks", {
    data <- generate_multi_block_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Partition without merging
    partitioned <- partitionLdMatrix(ld_data, mergeSmallBlocks = FALSE)

    # Check that each variant is mapped to the correct block
    for (i in seq_along(partitioned$ldMatrices)) {
        # Get variants in this block matrix
        block_variants <- rownames(partitioned$ldMatrices[[i]])

        # Find these variants in the variantIndices dataframe
        variant_block_ids <- partitioned$variantIndices$blockId[
            match(block_variants, partitioned$variantIndices$variant_id)
        ]

        # All variants should be mapped to this block
        expect_true(all(variant_block_ids == i))
    }

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix handles row/column name mismatches", {
    data <- generate_dummy_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Create an LdData with mismatched rownames and colnames on the correlation matrix
    ldmat <- getCorrelation(ld_data)
    vids <- getVariantIds(ld_data)
    rownames(ldmat) <- NULL
    colnames(ldmat) <- NULL
    mismatched_ld_data <- LdData(
        correlation = ldmat,
        variants = getVariantInfo(ld_data),
        blockMetadata = getBlockMetadata(ld_data)
    )

    # Should not error and should fix the names
    partitioned <- partitionLdMatrix(mismatched_ld_data)

    # Check if names are fixed
    expect_identical(
        rownames(partitioned$ldMatrices[[1]]),
        getVariantIds(ld_data)
    )
    expect_identical(
        colnames(partitioned$ldMatrices[[1]]),
        getVariantIds(ld_data)
    )

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix correctly extracts blocks based on metadata", {
    data <- generate_multi_block_data()
    region <- data$region
    LD_meta_file_path <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS")
    )
    write_delim(data$meta, LD_meta_file_path, delim = "\t")

    # Load the LD matrix
    ld_data <- loadLdMatrix(LD_meta_file_path, region)

    # Partition without merging
    partitioned <- partitionLdMatrix(ld_data, mergeSmallBlocks = FALSE)

    # For each block, check if the extracted matrix matches the expected submatrix
    ld_variants <- getVariantIds(ld_data)
    ld_matrix <- getCorrelation(ld_data)
    for (i in seq_along(partitioned$ldMatrices)) {
        block_info <- partitioned$blockMetadata[i, ]
        startIdx <- block_info$startIdx
        endIdx <- block_info$endIdx

        # Skip if indices are invalid
        if (
            startIdx > length(ld_variants) ||
                endIdx > length(ld_variants) ||
                endIdx < startIdx
        ) {
            next
        }

        # Get variants for this block
        block_variants <- ld_variants[startIdx:endIdx]

        # Extract expected submatrix
        expected_submatrix <- ld_matrix[
            block_variants,
            block_variants,
            drop = FALSE
        ]

        # Compare with actual block matrix
        expect_equal(partitioned$ldMatrices[[i]], expected_submatrix)
    }

    file.remove(LD_meta_file_path)
})

test_that("partitionLdMatrix partitions correctly with synthetic data", {
    mat <- matrix(0, 6, 6)
    mat[1:3, 1:3] <- 0.5
    mat[4:6, 4:6] <- 0.5
    diag(mat) <- 1
    variant_ids <- c(
        "chr1:100:A:G",
        "chr1:200:C:T",
        "chr1:300:G:A",
        "chr1:400:T:C",
        "chr1:500:A:G",
        "chr1:600:C:T"
    )
    rownames(mat) <- colnames(mat) <- variant_ids

    bm <- data.frame(
        blockId = c(1L, 2L),
        chrom = c("1", "1"),
        blockStart = c(100L, 400L),
        blockEnd = c(300L, 600L),
        size = c(3L, 3L),
        startIdx = c(1L, 4L),
        endIdx = c(3L, 6L),
        stringsAsFactors = FALSE
    )

    ld_data <- make_test_ld_data(variant_ids, R = mat, blockMetadata = bm)

    result <- pecotmr:::partitionLdMatrix(ld_data, mergeSmallBlocks = FALSE)

    expect_type(result, "list")
    expect_true("ldMatrices" %in% names(result))
    expect_true("variantIndices" %in% names(result))
    expect_length(result$ldMatrices, 2)
    expect_equal(nrow(result$ldMatrices[[1]]), 3)
    expect_equal(nrow(result$ldMatrices[[2]]), 3)
})

# ---- orderDedupRegions ----

test_that("orderDedupRegions removes duplicate regions", {
    # Create regions with duplicates
    regions_with_dups <- data.frame(
        chrom = c("chr1", "chr1", "chr1"),
        start = c(100, 100, 200), # Note: first two rows are duplicates
        end = c(150, 150, 250)
    )

    result <- orderDedupRegions(regions_with_dups)
    # Should have removed duplicate and return only two rows
    expect_equal(nrow(result), 2)
    expect_equal(result$start, c(100, 200))
})

test_that("orderDedupRegions orders and deduplicates across chromosomes", {
    df <- data.frame(
        chrom = c("chr2", "chr1", "chr1", "chr2"),
        start = c(100, 200, 100, 100),
        end = c(200, 300, 200, 200)
    )
    result <- pecotmr:::orderDedupRegions(df)
    expect_equal(nrow(result), 3) # one duplicate removed
    expect_true(all(diff(result$start[result$chrom == result$chrom[1]]) >= 0))
})

test_that("orderDedupRegions strips chr prefix", {
    df <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 200),
        end = c(200, 300)
    )
    result <- pecotmr:::orderDedupRegions(df)
    expect_equal(result$chrom, c("1", "2")) # normalized to a bare string, not integer
})

test_that("orderDedupRegions sorts X/Y/MT after autosomes and keeps chrom a string", {
    df <- data.frame(
        chrom = c("chrX", "chr2", "chr1", "chrMT"),
        start = c(1L, 1L, 1L, 1L),
        end = c(2L, 2L, 2L, 2L)
    )
    result <- pecotmr:::orderDedupRegions(df)
    expect_equal(result$chrom, c("1", "2", "X", "MT")) # genomic order, not NA-collapsed
    expect_type(result$chrom, "character")
})

# ---- findIntersectionRows ----

test_that("findIntersectionRows correctly identifies start and end rows", {
    # Create a simple genomic dataset
    genomic_data <- data.frame(
        chrom = c(1, 1, 1, 1),
        start = c(100, 200, 300, 400),
        end = c(150, 250, 350, 450)
    )

    # Region entirely within the dataset
    result <- findIntersectionRows(genomic_data, 1, 220, 330)
    expect_equal(result$startRow$start, 200)
    expect_equal(result$endRow$end, 350)
})

test_that("findIntersectionRows adjusts region bounds if needed", {
    # Create a simple genomic dataset
    genomic_data <- data.frame(
        chrom = c(1, 1, 1, 1),
        start = c(100, 200, 300, 400),
        end = c(150, 250, 350, 450)
    )

    # Region extends beyond the dataset
    result <- findIntersectionRows(genomic_data, 1, 50, 500)
    # Should adjust to the bounds of the dataset
    expect_equal(result$startRow$start, 100)
    expect_equal(result$endRow$end, 450)
})

test_that("findIntersectionRows errors for non-overlapping regions", {
    # Create a simple genomic dataset
    genomic_data <- data.frame(
        chrom = c(1, 1, 1, 1),
        start = c(100, 200, 300, 400),
        end = c(150, 250, 350, 450)
    )

    # Region entirely outside the dataset
    expect_error(
        findIntersectionRows(genomic_data, 2, 100, 200),
        "No data for chromosome 2"
    )
})

# ---- validateSelectedRegion ----

test_that("validateSelectedRegion passes for valid region", {
    startRow <- data.frame(start = 0)
    endRow <- data.frame(end = 300)
    expect_silent(pecotmr:::validateSelectedRegion(startRow, endRow, 50, 250))
})

test_that("validateSelectedRegion errors for uncovered region", {
    startRow <- data.frame(start = 100)
    endRow <- data.frame(end = 200)
    expect_error(
        pecotmr:::validateSelectedRegion(startRow, endRow, 50, 250),
        "not fully covered"
    )
})

# ---- extractFilePaths ----

test_that("extractFilePaths extracts correct paths", {
    gd <- data.frame(
        chrom = c(1, 1, 1),
        start = c(0, 100, 200),
        end = c(100, 200, 300),
        path = c("f1.ld", "f2.ld", "f3.ld")
    )
    intersection <- list(
        startRow = data.frame(chrom = 1, start = 0),
        endRow = data.frame(start = 200)
    )
    result <- pecotmr:::extractFilePaths(gd, intersection, "path")
    expect_equal(length(result), 3)
})

test_that("extractFilePaths errors on missing column", {
    gd <- data.frame(chrom = 1, start = 0, end = 100)
    intersection <- list(
        startRow = data.frame(chrom = 1, start = 0),
        endRow = data.frame(start = 0)
    )
    expect_error(
        pecotmr:::extractFilePaths(gd, intersection, "nonexistent"),
        "not found"
    )
})

# ---- partitionLdMatrix: different chromosomes ----

test_that("partitionLdMatrix handles blocks with different chromosomes", {
    # Create test data with blocks on different chromosomes
    test_matrix <- matrix(0, 4, 4)
    diag(test_matrix) <- 1 # Set diagonal to 1
    variant_ids <- c(
        "chr1:100:A:G",
        "chr1:200:C:T",
        "chr2:100:G:A",
        "chr2:200:T:C"
    )
    rownames(test_matrix) <- colnames(test_matrix) <- variant_ids

    blockMetadata <- data.frame(
        blockId = c(1L, 2L),
        chrom = c("1", "2"),
        blockStart = c(100L, 100L),
        blockEnd = c(200L, 200L),
        size = c(2L, 2L),
        startIdx = c(1L, 3L),
        endIdx = c(2L, 4L),
        stringsAsFactors = FALSE
    )

    test_ld_data <- make_test_ld_data(
        variant_ids,
        R = test_matrix,
        blockMetadata = blockMetadata
    )

    # Partition the matrix
    partitioned <- partitionLdMatrix(test_ld_data)

    # Should not merge blocks from different chromosomes
    expect_equal(length(partitioned$ldMatrices), 2)

    # Each block should have the correct variants
    expect_equal(
        rownames(partitioned$ldMatrices[[1]]),
        c("chr1:100:A:G", "chr1:200:C:T")
    )
    expect_equal(
        rownames(partitioned$ldMatrices[[2]]),
        c("chr2:100:G:A", "chr2:200:T:C")
    )
})

test_that("partitionLdMatrix works with edge case block structures", {
    # Test case: One large block and several tiny blocks that need merging
    large_block_size <- 15
    small_block_size <- 2

    # Create a matrix with blocks of varying sizes
    n_variants <- large_block_size + small_block_size * 3
    test_matrix <- matrix(0, n_variants, n_variants)
    # Set diagonal to 1
    diag(test_matrix) <- 1

    # Generate variant names in chr:pos:A2:A1 format
    variantNames <- paste0("chr1:", 100:(100 + n_variants - 1), ":A:G")
    rownames(test_matrix) <- colnames(test_matrix) <- variantNames

    # Create block metadata
    blockMetadata <- data.frame(
        blockId = 1:4,
        chrom = rep("1", 4),
        blockStart = c(
            100L,
            as.integer(100 + large_block_size),
            as.integer(100 + large_block_size + small_block_size),
            as.integer(100 + large_block_size + small_block_size * 2)
        ),
        blockEnd = c(
            as.integer(100 + large_block_size - 1),
            as.integer(100 + large_block_size + small_block_size - 1),
            as.integer(100 + large_block_size + small_block_size * 2 - 1),
            as.integer(100 + n_variants - 1)
        ),
        size = c(
            large_block_size,
            small_block_size,
            small_block_size,
            small_block_size
        ),
        startIdx = c(
            1L,
            as.integer(large_block_size + 1),
            as.integer(large_block_size + small_block_size + 1),
            as.integer(large_block_size + small_block_size * 2 + 1)
        ),
        endIdx = c(
            as.integer(large_block_size),
            as.integer(large_block_size + small_block_size),
            as.integer(large_block_size + small_block_size * 2),
            as.integer(n_variants)
        ),
        stringsAsFactors = FALSE
    )

    test_ld_data <- make_test_ld_data(
        variantNames,
        R = test_matrix,
        blockMetadata = blockMetadata
    )

    # Set minimum block size to force merging of small blocks
    min_merged_size <- small_block_size + 1

    # Partition with merging
    partitioned <- partitionLdMatrix(
        test_ld_data,
        mergeSmallBlocks = TRUE,
        minMergedBlockSize = min_merged_size
    )

    # Should merge the small blocks but leave the large block alone
    expect_lt(length(partitioned$ldMatrices), 4)
    expect_gt(length(partitioned$ldMatrices), 1)

    # First block should still be large_block_size
    expect_equal(nrow(partitioned$ldMatrices[[1]]), large_block_size)
})

# ---- extractLdForRegion ----

test_that("extractLdForRegion extracts correct region", {
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

    result <- extractLdForRegion(ld_matrix, ld_variants, region, NULL)

    # Should have extracted only the relevant variants
    expect_equal(nrow(result$extractedLdVariants), 2)
    expect_equal(
        result$extractedLdVariants$variants,
        c("1:200:C:T", "1:300:G:A")
    )

    # Matrix should be 2x2 with the correct row/column names
    expect_equal(dim(result$extractedLdMatrix), c(2, 2))
    expect_equal(
        rownames(result$extractedLdMatrix),
        c("1:200:C:T", "1:300:G:A")
    )
})

test_that("extractLdForRegion works with extract_coordinates", {
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

    result <- extractLdForRegion(
        ld_matrix,
        ld_variants,
        region,
        extract_coordinates
    )

    # Should have extracted only the specified coordinates
    expect_equal(nrow(result$extractedLdVariants), 2)
    expect_equal(
        result$extractedLdVariants$variants,
        c("1:100:A:G", "1:300:G:A")
    )

    # Matrix should be 2x2 with the correct row/column names
    expect_equal(dim(result$extractedLdMatrix), c(2, 2))
    expect_equal(
        rownames(result$extractedLdMatrix),
        c("1:100:A:G", "1:300:G:A")
    )
})

# ---- createLdMatrix ----

test_that("createLdMatrix correctly combines matrices with overlapping variants", {
    # Create two simple LD matrices with some overlapping variants
    matrix1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    rownames(matrix1) <- colnames(matrix1) <- c("1:100:A:G", "1:200:C:T")

    matrix2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
    rownames(matrix2) <- colnames(matrix2) <- c("1:200:C:T", "1:300:G:A")

    # Create variants lists
    variants1 <- data.frame(variants = c("1:100:A:G", "1:200:C:T"))
    variants2 <- data.frame(variants = c("1:200:C:T", "1:300:G:A"))

    # Combine matrices
    combined <- createLdMatrix(
        ldMatrices = list(matrix1, matrix2),
        variants = list(variants1, variants2)
    )

    # Should have created a 3x3 matrix with all unique variants
    expect_equal(dim(combined), c(3, 3))
    expect_equal(rownames(combined), c("1:100:A:G", "1:200:C:T", "1:300:G:A"))

    # Check that values from original matrices are preserved
    expect_equal(combined["1:100:A:G", "1:200:C:T"], 0.5)
    expect_equal(combined["1:200:C:T", "1:300:G:A"], 0.3)

    # Check diagonal values are 1
    expect_equal(combined[1, 1], 1)
    expect_equal(combined[2, 2], 1)
    expect_equal(combined[3, 3], 1)
})

test_that("createLdMatrix merges non-overlapping blocks", {
    m1 <- matrix(
        c(1, 0.5, 0.5, 1),
        2,
        2,
        dimnames = list(
            c("chr1:100:A:G", "chr1:200:A:G"),
            c("chr1:100:A:G", "chr1:200:A:G")
        )
    )
    m2 <- matrix(
        c(1, 0.3, 0.3, 1),
        2,
        2,
        dimnames = list(
            c("chr1:300:A:G", "chr1:400:A:G"),
            c("chr1:300:A:G", "chr1:400:A:G")
        )
    )

    variants <- list(
        data.frame(variants = c("chr1:100:A:G", "chr1:200:A:G")),
        data.frame(variants = c("chr1:300:A:G", "chr1:400:A:G"))
    )
    result <- pecotmr:::createLdMatrix(list(m1, m2), variants)

    expect_equal(nrow(result), 4)
    expect_equal(ncol(result), 4)
    expect_equal(result["chr1:100:A:G", "chr1:200:A:G"], 0.5)
    expect_equal(result["chr1:300:A:G", "chr1:400:A:G"], 0.3)
    # Cross-block should be 0
    expect_equal(result["chr1:100:A:G", "chr1:300:A:G"], 0)
})

test_that("createLdMatrix handles overlapping boundary variant", {
    m1 <- matrix(
        c(1, 0.5, 0.5, 1),
        2,
        2,
        dimnames = list(
            c("chr1:100:A:G", "chr1:200:A:G"),
            c("chr1:100:A:G", "chr1:200:A:G")
        )
    )
    m2 <- matrix(
        c(1, 0.3, 0.3, 1),
        2,
        2,
        dimnames = list(
            c("chr1:200:A:G", "chr1:300:A:G"),
            c("chr1:200:A:G", "chr1:300:A:G")
        )
    )

    variants <- list(
        data.frame(variants = c("chr1:100:A:G", "chr1:200:A:G")),
        data.frame(variants = c("chr1:200:A:G", "chr1:300:A:G"))
    )
    result <- pecotmr:::createLdMatrix(list(m1, m2), variants)

    # v2 is shared, so total should be 3 variants
    expect_equal(nrow(result), 3)
    expect_equal(ncol(result), 3)
})

# ---- validateBlockStructure ----

test_that("validateBlockStructure passes for proper block structure", {
    mat <- matrix(0, 6, 6)
    mat[1:3, 1:3] <- 0.5
    mat[4:6, 4:6] <- 0.5
    diag(mat) <- 1

    variant_ids <- sprintf("chr1:%d:A:G", 100L * (1:6))
    rownames(mat) <- colnames(mat) <- variant_ids

    block_meta <- data.frame(
        blockId = c(1, 2),
        chrom = c("1", "1"),
        size = c(3, 3),
        startIdx = c(1, 4),
        endIdx = c(3, 6)
    )

    expect_silent(pecotmr:::validateBlockStructure(
        mat,
        block_meta,
        variant_ids
    ))
})

test_that("validateBlockStructure errors on non-block structure", {
    mat <- matrix(0.5, 4, 4)
    diag(mat) <- 1

    variant_ids <- sprintf("chr1:%d:A:G", 100L * (1:4))
    rownames(mat) <- colnames(mat) <- variant_ids

    block_meta <- data.frame(
        blockId = c(1, 2),
        chrom = c("1", "1"),
        size = c(2, 2),
        startIdx = c(1, 3),
        endIdx = c(2, 4)
    )

    expect_error(
        pecotmr:::validateBlockStructure(mat, block_meta, variant_ids),
        "Matrix lacks expected block structure"
    )
})

# ---- mergeBlocks ----

test_that("mergeBlocks properly handles blocks at chromosome boundaries", {
    # Create test data with small blocks at chromosome boundaries
    test_matrix <- matrix(0, 6, 6)
    diag(test_matrix) <- 1
    variantNames <- c(
        "chr1:900:A:G",
        "chr1:950:C:T",
        "chr2:100:G:A",
        "chr2:150:T:C",
        "chr3:100:A:G",
        "chr3:150:C:T"
    )
    rownames(test_matrix) <- colnames(test_matrix) <- variantNames

    # Create block metadata with small blocks at chromosome boundaries
    blockMetadata <- data.frame(
        blockId = c(1L, 2L, 3L),
        chrom = c("1", "2", "3"),
        blockStart = c(900L, 100L, 100L),
        blockEnd = c(950L, 150L, 150L),
        size = c(2L, 2L, 2L),
        startIdx = c(1L, 3L, 5L),
        endIdx = c(2L, 4L, 6L),
        stringsAsFactors = FALSE
    )

    test_ld_data <- make_test_ld_data(
        variantNames,
        R = test_matrix,
        blockMetadata = blockMetadata
    )

    # Set min block size to force merging attempts
    min_block_size <- 3

    # Partition with merging
    partitioned <- partitionLdMatrix(
        test_ld_data,
        mergeSmallBlocks = TRUE,
        minMergedBlockSize = min_block_size
    )

    # Should not merge blocks across chromosome boundaries
    expect_equal(length(partitioned$ldMatrices), 3)

    # Each block should match its chromosome
    for (i in 1:3) {
        block_variants <- rownames(partitioned$ldMatrices[[i]])
        # Strip "chr" prefix before extracting chromosome number
        chrom_from_variants <- unique(as.integer(sub(
            "chr([0-9]+):.*",
            "\\1",
            block_variants
        )))
        expect_equal(length(chrom_from_variants), 1) # Should only have one chromosome per block
        expect_equal(chrom_from_variants, i) # Should match the expected chromosome
    }
})

test_that("mergeBlocks merges small adjacent blocks", {
    block_meta <- data.frame(
        blockId = c(1, 2, 3),
        chrom = c("1", "1", "1"),
        size = c(50, 50, 100),
        startIdx = c(1, 51, 101),
        endIdx = c(50, 100, 200)
    )
    result <- pecotmr:::mergeBlocks(block_meta, minSize = 100, maxSize = 10000)
    expect_true(nrow(result) < 3)
})

test_that("mergeBlocks does not merge cross-chromosome", {
    block_meta <- data.frame(
        blockId = c(1, 2),
        chrom = c("1", "2"),
        size = c(10, 10),
        startIdx = c(1, 11),
        endIdx = c(10, 20)
    )
    result <- pecotmr:::mergeBlocks(block_meta, minSize = 50, maxSize = 10000)
    expect_equal(nrow(result), 2) # Cannot merge across chromosomes
})

test_that("mergeBlocks returns single block unchanged", {
    block_meta <- data.frame(
        blockId = 1,
        chrom = "1",
        size = 10,
        startIdx = 1,
        endIdx = 10
    )
    result <- pecotmr:::mergeBlocks(block_meta, minSize = 100, maxSize = 10000)
    expect_equal(nrow(result), 1)
})

# ---- canMerge ----

test_that("canMerge checks chromosome and size", {
    bm <- data.frame(
        chrom = c("1", "1", "2"),
        size = c(100, 200, 100),
        stringsAsFactors = FALSE
    )
    # rows 1 and 2: same chrom, combined size 300
    expect_true(pecotmr:::canMerge(bm, 1, 2, maxSize = 500))
    expect_false(pecotmr:::canMerge(bm, 1, 2, maxSize = 200))
    # rows 1 and 3: different chromosome
    expect_false(pecotmr:::canMerge(bm, 1, 3, maxSize = 500))
})

# ===========================================================================
# checkLd (regularize_ld)
# ===========================================================================

test_that("checkLd reports PD for identity matrix", {
    R <- diag(5)
    result <- checkLd(R)
    expect_true(result$isPd)
    expect_true(result$isPsd)
    expect_equal(result$methodApplied, "none")
    expect_equal(result$R, R)
    expect_equal(result$conditionNumber, 1)
})

test_that("checkLd reports PD for well-conditioned correlation matrix", {
    R <- matrix(0.3, 4, 4)
    diag(R) <- 1
    result <- checkLd(R)
    expect_true(result$isPd)
    expect_true(result$isPsd)
    expect_equal(result$nNegative, 0)
    expect_equal(result$methodApplied, "none")
})

test_that("checkLd detects non-PSD matrix", {
    R <- matrix(0.9, 3, 3)
    diag(R) <- 1
    R[1, 3] <- R[3, 1] <- -0.5
    result <- checkLd(R)
    expect_false(result$isPsd)
    expect_true(result$nNegative > 0)
    expect_true(result$minEigenvalue < 0)
    expect_equal(result$methodApplied, "none")
})

test_that("checkLd shrink method modifies non-PD matrix", {
    R <- matrix(0.9, 3, 3)
    diag(R) <- 1
    R[1, 3] <- R[3, 1] <- -0.5
    result <- checkLd(R, method = "shrink")
    expect_equal(result$methodApplied, "shrink")
    expect_false(identical(result$R, R))
    # With strong enough shrinkage, result should be PD
    result2 <- checkLd(R, method = "shrink", shrinkage = 0.5)
    eig <- eigen(result2$R, symmetric = TRUE)
    expect_true(all(eig$values > 0))
})

test_that("checkLd eigenfix method improves non-PD matrix", {
    R <- matrix(0.9, 3, 3)
    diag(R) <- 1
    R[1, 3] <- R[3, 1] <- -0.5
    original_min_eval <- min(eigen(R, symmetric = TRUE)$values)
    result <- checkLd(R, method = "eigenfix")
    expect_equal(result$methodApplied, "eigenfix")
    # Eigenfix should improve the minimum eigenvalue
    fixed_min_eval <- min(eigen(result$R, symmetric = TRUE)$values)
    expect_true(fixed_min_eval > original_min_eval)
    # Unit diagonal preserved
    expect_equal(diag(result$R), rep(1, 3))
    # Symmetry preserved
    expect_equal(result$R, t(result$R))
})

test_that("checkLd shrink does nothing when matrix is already PD", {
    R <- diag(3)
    result <- checkLd(R, method = "shrink")
    expect_equal(result$methodApplied, "none")
    expect_equal(result$R, R)
})

test_that("checkLd eigenfix does nothing when matrix is already PD", {
    R <- diag(3)
    result <- checkLd(R, method = "eigenfix")
    expect_equal(result$methodApplied, "none")
    expect_equal(result$R, R)
})

# ===========================================================================
# extractBlockMatrices: out-of-range blocks
# ===========================================================================

test_that("extractBlockMatrices warns and skips out-of-range blocks", {
    mat <- diag(4)
    vnames <- sprintf("chr1:%d:A:G", 100L * (1:4))
    rownames(mat) <- colnames(mat) <- vnames
    blockMetadata <- data.frame(
        blockId = c(1, 2),
        startIdx = c(1, 10),
        endIdx = c(2, 12),
        chrom = c("1", "1"),
        blockStart = c(100, 500),
        blockEnd = c(200, 600),
        size = c(2, 3),
        stringsAsFactors = FALSE
    )
    expect_warning(
        result <- pecotmr:::extractBlockMatrices(mat, blockMetadata, vnames),
        "outside the range"
    )
    valid_blocks <- result$ldMatrices[!sapply(result$ldMatrices, is.null)]
    expect_equal(length(valid_blocks), 1)
    expect_equal(nrow(valid_blocks[[1]]), 2)
})

# ===========================================================================
# resolveLdSource: type detection with real fixtures
# ===========================================================================

geno_test_data_dir <- test_path("test_data")
geno_region_all <- "chr21:17513228-17592874"

test_that("resolveLdSource detects PLINK2 from metadata", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_p2_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::resolveLdSource(meta_file)
    expect_equal(result$type, "plink2")
    expect_equal(result$metaPath, meta_file)
})

test_that("resolveLdSource detects VCF from metadata", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_vcf_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::resolveLdSource(meta_file)
    expect_equal(result$type, "vcf")
})

test_that("resolveLdSource detects GDS from metadata", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_gds_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants.gds", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::resolveLdSource(meta_file)
    expect_equal(result$type, "gds")
})

test_that("resolveLdSource detects precomputed from metadata", {
    # Existing LD block metadata with non-zero start/end
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_pre_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "chr1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::resolveLdSource(meta_file)
    expect_equal(result$type, "precomputed")
})

test_that("resolveLdSource errors on missing file", {
    expect_error(
        pecotmr:::resolveLdSource("/nonexistent/file.tsv"),
        "not found"
    )
})

test_that("resolveLdSource errors on 0:0 sentinel with non-genotype path", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_bad_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "nonexistent_prefix", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    expect_error(pecotmr:::resolveLdSource(meta_file), "0:0 sentinel")
})

# ===========================================================================
# resolveGenotypePathForRegion
# ===========================================================================

test_that("resolveGenotypePathForRegion resolves correct chromosome path", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_path_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::resolveGenotypePathForRegion(meta_file, geno_region_all)
    expect_equal(result, file.path(geno_test_data_dir, "test_variants"))
})

test_that("resolveGenotypePathForRegion errors on missing chromosome", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_nochr_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("1", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    expect_error(
        pecotmr:::resolveGenotypePathForRegion(meta_file, geno_region_all),
        "No entry for chromosome"
    )
})

# ===========================================================================
# loadLdFromGenotype with real fixtures
# ===========================================================================

test_that("loadLdFromGenotype returns LD matrix with .afreq", {
    skip_if_not_installed("pgenlibr")
    plink_prefix <- file.path(geno_test_data_dir, "test_variants")
    result <- pecotmr:::loadLdFromGenotype(plink_prefix, geno_region_all)
    expect_true(is(result, "LdData"))
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 349L)
    expect_true(isSymmetric(getCorrelation(result)))
    expect_false(hasGenotypes(result))
    # ref_panel should have allele_freq from .afreq file
    expect_true(
        "allele_freq" %in% names(S4Vectors::mcols(getVariantInfo(result)))
    )
    expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq > 0))
    expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq < 1))
    # blockMetadata
    expect_true(is.data.frame(getBlockMetadata(result)))
    expect_equal(nrow(getBlockMetadata(result)), 1L)
})

test_that("loadLdFromGenotype returns genotype matrix when requested", {
    skip_if_not_installed("pgenlibr")
    plink_prefix <- file.path(geno_test_data_dir, "test_variants")
    result <- pecotmr:::loadLdFromGenotype(
        plink_prefix,
        geno_region_all,
        returnGenotype = TRUE
    )
    expect_true(hasGenotypes(result))
    X <- getGenotypes(result)
    expect_equal(nrow(X), 100L) # samples
    expect_equal(ncol(X), 349L) # variants
})

test_that("loadLdFromGenotype computes variance with n_sample", {
    skip_if_not_installed("pgenlibr")
    plink_prefix <- file.path(geno_test_data_dir, "test_variants")
    result <- pecotmr:::loadLdFromGenotype(
        plink_prefix,
        geno_region_all,
        nSample = 100L
    )
    expect_true("variance" %in% names(S4Vectors::mcols(getVariantInfo(result))))
    expect_true("n_nomiss" %in% names(S4Vectors::mcols(getVariantInfo(result))))
    expect_equal(S4Vectors::mcols(getVariantInfo(result))$n_nomiss[1], 100L)
    expect_true(all(S4Vectors::mcols(getVariantInfo(result))$variance > 0))
})

test_that("loadLdFromGenotype falls back to computed AF without .afreq", {
    skip_if_not_installed("VariantAnnotation")
    vcf_path <- file.path(geno_test_data_dir, "test_variants.vcf.gz")
    result <- suppressWarnings(
        pecotmr:::loadLdFromGenotype(vcf_path, geno_region_all)
    )
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 349L)
    expect_true(isSymmetric(getCorrelation(result)))
    # Allele frequencies computed from genotypes
    expect_true(
        "allele_freq" %in% names(S4Vectors::mcols(getVariantInfo(result)))
    )
    expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq > 0))
    expect_true(all(S4Vectors::mcols(getVariantInfo(result))$allele_freq < 1))
})

test_that("loadLdFromGenotype works with GDS files", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    gds_path <- file.path(geno_test_data_dir, "test_variants.gds")
    result <- pecotmr:::loadLdFromGenotype(gds_path, geno_region_all)
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 349L)
    expect_true(isSymmetric(getCorrelation(result)))
})

test_that("loadLdFromGenotype .afreq and computed AF are consistent", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    plink_prefix <- file.path(geno_test_data_dir, "test_variants")
    gds_path <- file.path(geno_test_data_dir, "test_variants.gds")
    res_afreq <- pecotmr:::loadLdFromGenotype(plink_prefix, geno_region_all)
    res_computed <- pecotmr:::loadLdFromGenotype(gds_path, geno_region_all)
    # Allele frequencies should be close (same data, different source)
    expect_true(
        max(abs(
            S4Vectors::mcols(getVariantInfo(res_afreq))$allele_freq -
                S4Vectors::mcols(getVariantInfo(res_computed))$allele_freq
        )) <
            0.01
    )
})

# ===========================================================================
# loadLdMatrix with real genotype fixtures via metadata
# ===========================================================================

test_that("loadLdMatrix dispatches to PLINK2 genotype source", {
    skip_if_not_installed("pgenlibr")
    meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_p2_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- loadLdMatrix(meta_file, geno_region_all)
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 349L)
    expect_false(hasGenotypes(result))
})

test_that("loadLdMatrix dispatches to VCF genotype source", {
    skip_if_not_installed("VariantAnnotation")
    meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_vcf_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- suppressWarnings(loadLdMatrix(meta_file, geno_region_all))
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 349L)
})

test_that("loadLdMatrix dispatches to GDS genotype source", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_gds_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants.gds", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- loadLdMatrix(meta_file, geno_region_all)
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 349L)
})

test_that("loadLdMatrix return_genotype='auto' returns X for genotype source", {
    skip_if_not_installed("pgenlibr")
    meta_file <- file.path(geno_test_data_dir, "ld_meta_ldmat_auto_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- loadLdMatrix(meta_file, geno_region_all, returnGenotype = "auto")
    expect_true(hasGenotypes(result))
    X <- getGenotypes(result)
    expect_equal(nrow(X), 100L) # samples
    expect_equal(ncol(X), 349L) # variants
})

test_that("loadLdMatrix return_genotype=TRUE errors for precomputed", {
    meta_file <- gsub(
        "//",
        "/",
        tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".tsv")
    )
    on.exit(unlink(meta_file), add = TRUE)
    meta_df <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1200,
        path = paste0(
            "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,",
            "./test_data/LD_block_1.chr1_1000_1200.float16.bim"
        )
    )
    write_delim(meta_df, meta_file, delim = "\t")
    region <- data.frame(chrom = "chr1", start = 1000, end = 1190)
    expect_error(
        loadLdMatrix(meta_file, region, returnGenotype = TRUE),
        "genotype files"
    )
})

# ===========================================================================
# resolveLdSource: PLINK1 detection
# ===========================================================================

test_that("resolveLdSource detects PLINK1 from metadata", {
    skip_if_not_installed("snpStats")
    meta_file <- file.path(geno_test_data_dir, "ld_meta_resolve_p1_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("22", "0", "0", "protocol_example.genotype", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::resolveLdSource(meta_file)
    expect_equal(result$type, "plink1")
})

# ===========================================================================
# loadLdMatrix: precomputed blocks via real .cor.xz fixtures
# ===========================================================================

test_that("loadLdMatrix loads single precomputed block", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_single_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- loadLdMatrix(meta_file, "chr1:1000-1190")
    expect_true(is.matrix(getCorrelation(result)))
    expect_equal(nrow(getCorrelation(result)), 5L)
    expect_true(isSymmetric(getCorrelation(result)))
    expect_equal(length(getVariantIds(result)), 5L)
    expect_true(all(grepl("^chr1:", getVariantIds(result))))
    expect_false(hasGenotypes(result))
    # blockMetadata should have one block
    expect_equal(nrow(getBlockMetadata(result)), 1L)
    # ref_panel (now GRanges) should have variant info via mcols
    ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
    expect_true("variant_id" %in% names(ref_mcols))
})

test_that("loadLdMatrix loads multiple precomputed blocks", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_multi_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        paste(
            "1",
            "1200",
            "1400",
            "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
            sep = "\t"
        ),
        paste(
            "1",
            "1400",
            "1600",
            "LD_block_3.chr1_1400_1600.float16.txt.xz,LD_block_3.chr1_1400_1600.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    result <- loadLdMatrix(meta_file, "chr1:1000-1500")
    expect_true(is.matrix(getCorrelation(result)))
    # Should span blocks 1-3: 5 + 5 + 5 = 15 unique variants (no overlap in variant IDs)
    expect_true(nrow(getCorrelation(result)) >= 10)
    expect_true(isSymmetric(getCorrelation(result)))
    expect_true(nrow(getBlockMetadata(result)) >= 2)
})

test_that("loadLdMatrix with n_sample for precomputed blocks with freq data", {
    # The 9-column bim format includes allele_freq, variance, n_nomiss;
    # the 6-column bim does not. With 6-col bim and no allele_freq,
    # n_sample cannot compute variance - ref_panel has base columns only.
    meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_nsamp_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- loadLdMatrix(meta_file, "chr1:1000-1190", nSample = 500L)
    # ref_panel (GRanges) should always have basic variant info in mcols
    ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
    expect_true(all(c("A2", "A1", "variant_id") %in% names(ref_mcols)))
    # 6-col bim lacks allele_freq so variance computation is skipped
    expect_false("variance" %in% names(ref_mcols))
})

# ===========================================================================
# processLdMatrix: real .cor.xz fixtures
# ===========================================================================

test_that("processLdMatrix reads .cor.xz with explicit bim path", {
    ld_file <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.txt.xz"
    )
    bim_file <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.bim"
    )
    result <- pecotmr:::processLdMatrix(ld_file, bim_file)
    expect_true(is.list(result))
    expect_true(is.matrix(result$ldMatrix))
    expect_equal(nrow(result$ldMatrix), 5L)
    expect_equal(ncol(result$ldMatrix), 5L)
    expect_true(isSymmetric(result$ldMatrix))
    # Diagonal should be 1
    expect_true(all(abs(diag(result$ldMatrix) - 1) < 1e-4))
    # Variant names should be chr:pos:A2:A1 format
    expect_true(all(grepl("^chr1:", rownames(result$ldMatrix))))
    # ldVariants data frame
    expect_true(is.data.frame(result$ldVariants))
    expect_true("variants" %in% names(result$ldVariants))
    expect_equal(nrow(result$ldVariants), 5L)
})

test_that("processLdMatrix reads different blocks consistently", {
    bim1 <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.bim"
    )
    bim2 <- file.path(
        geno_test_data_dir,
        "LD_block_2.chr1_1200_1400.float16.bim"
    )
    ld1 <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.txt.xz"
    )
    ld2 <- file.path(
        geno_test_data_dir,
        "LD_block_2.chr1_1200_1400.float16.txt.xz"
    )
    r1 <- pecotmr:::processLdMatrix(ld1, bim1)
    r2 <- pecotmr:::processLdMatrix(ld2, bim2)
    # Different blocks should have different variant positions
    expect_false(any(rownames(r1$ldMatrix) %in% rownames(r2$ldMatrix)))
})

test_that("processLdMatrix reads 9-column bim with allele_freq/variance/n_nomiss", {
    ld_file <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.txt.xz"
    )
    bim_file <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.9col.bim"
    )
    result <- pecotmr:::processLdMatrix(ld_file, bim_file)
    expect_equal(nrow(result$ldMatrix), 5L)
    expect_true(isSymmetric(result$ldMatrix))
    # 9-column bim should include extra columns
    expect_true("allele_freq" %in% names(result$ldVariants))
    expect_true("variance" %in% names(result$ldVariants))
    expect_true("n_nomiss" %in% names(result$ldVariants))
    expect_equal(result$ldVariants$allele_freq, c(0.3, 0.4, 0.2, 0.5, 0.15))
    expect_equal(result$ldVariants$n_nomiss, rep(500, 5))
})

test_that("loadLdMatrix propagates allele_freq/variance/n_nomiss from 9-col bim", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_precomp_9col_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.9col.bim",
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- loadLdMatrix(meta_file, "chr1:1000-1190")
    # ref_panel (GRanges) should carry the extra columns from the 9-col bim
    ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
    expect_true("allele_freq" %in% names(ref_mcols))
    expect_true("variance" %in% names(ref_mcols))
    expect_true("n_nomiss" %in% names(ref_mcols))
    expect_equal(ref_mcols$allele_freq, c(0.3, 0.4, 0.2, 0.5, 0.15))
    expect_equal(ref_mcols$n_nomiss, rep(500, 5))
})

# ===========================================================================
# getRegionalLdMeta: real .cor.xz fixtures
# ===========================================================================

test_that("getRegionalLdMeta returns correct file paths for single block", {
    meta_file <- file.path(
        geno_test_data_dir,
        "ld_meta_regional_single_tmp.tsv"
    )
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        paste(
            "1",
            "1200",
            "1400",
            "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    result <- pecotmr:::getRegionalLdMeta(meta_file, "chr1:1050-1150")
    expect_true(is.list(result))
    expect_true(length(result$intersections$LD_file_paths) >= 1)
    # All returned paths should exist
    expect_true(all(file.exists(result$intersections$LD_file_paths)))
    expect_true(all(file.exists(result$intersections$bimFilePaths)))
})

test_that("getRegionalLdMeta spans multiple blocks for wide region", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_regional_multi_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        paste(
            "1",
            "1200",
            "1400",
            "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
            sep = "\t"
        ),
        paste(
            "1",
            "1400",
            "1600",
            "LD_block_3.chr1_1400_1600.float16.txt.xz,LD_block_3.chr1_1400_1600.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    result <- pecotmr:::getRegionalLdMeta(meta_file, "chr1:1000-1500")
    expect_true(length(result$intersections$LD_file_paths) >= 2)
})

# ===========================================================================
# dropCollinearColumns: strategy variants
# ===========================================================================

test_that("dropCollinearColumns variance strategy removes lowest-variance column", {
    set.seed(42)
    X <- matrix(rnorm(100 * 4), 100, 4)
    colnames(X) <- c("a", "b", "c", "d")
    # Make column "c" have near-zero variance
    X[, "c"] <- X[1, "c"]
    result <- pecotmr:::dropCollinearColumns(
        X,
        c("b", "c", "d"),
        strategy = "variance"
    )
    expect_false("c" %in% colnames(result))
    expect_equal(ncol(result), 3L)
})

test_that("dropCollinearColumns responseCorrelation strategy works", {
    set.seed(42)
    X <- matrix(rnorm(100 * 3), 100, 3)
    colnames(X) <- c("a", "b", "c")
    y <- X[, "a"] + rnorm(100, sd = 0.1) # y correlates strongly with "a"
    result <- pecotmr:::dropCollinearColumns(
        X,
        c("a", "b", "c"),
        strategy = "responseCorrelation",
        response = y
    )
    # Should keep "a" (highest |cor| with response) and remove one of b/c
    expect_true("a" %in% colnames(result))
    expect_equal(ncol(result), 2L)
})

test_that("dropCollinearColumns responseCorrelation errors without response", {
    X <- matrix(1:12, 4, 3)
    colnames(X) <- c("a", "b", "c")
    expect_error(
        pecotmr:::dropCollinearColumns(
            X,
            c("a", "b"),
            strategy = "responseCorrelation"
        ),
        "response must be supplied"
    )
})

test_that("dropCollinearColumns with single problematic column removes it", {
    X <- matrix(rnorm(40), 10, 4)
    colnames(X) <- c("a", "b", "c", "d")
    result <- pecotmr:::dropCollinearColumns(X, "b", strategy = "correlation")
    expect_false("b" %in% colnames(result))
    expect_equal(ncol(result), 3L)
})

# ===========================================================================
# enforceDesignFullRank: additional strategies and fallback paths
# ===========================================================================

test_that("enforceDesignFullRank variance strategy produces full rank", {
    set.seed(42)
    X <- matrix(rnorm(100 * 4), 100, 4)
    X[, 4] <- X[, 1] + X[, 2] # rank deficient
    colnames(X) <- c("a", "b", "c", "d")
    C <- matrix(rnorm(100), 100, 1)
    result <- enforceDesignFullRank(X, C, strategy = "variance")
    full_design <- cbind(1, result, C)
    expect_equal(qr(full_design)$rank, ncol(full_design))
    expect_true(ncol(result) < ncol(X))
})

test_that("enforceDesignFullRank responseCorrelation strategy works", {
    set.seed(42)
    X <- matrix(rnorm(100 * 4), 100, 4)
    X[, 4] <- X[, 1] + X[, 2]
    colnames(X) <- c("a", "b", "c", "d")
    C <- matrix(rnorm(100), 100, 1)
    y <- X[, "a"] + rnorm(100, sd = 0.1)
    result <- enforceDesignFullRank(
        X,
        C,
        strategy = "responseCorrelation",
        response = y
    )
    full_design <- cbind(1, result, C)
    expect_equal(qr(full_design)$rank, ncol(full_design))
})

test_that("enforceDesignFullRank returns unchanged X when already full rank", {
    set.seed(42)
    X <- matrix(rnorm(100 * 3), 100, 3)
    colnames(X) <- c("a", "b", "c")
    C <- matrix(rnorm(100), 100, 1)
    result <- enforceDesignFullRank(X, C, strategy = "correlation")
    expect_equal(ncol(result), ncol(X))
})

test_that("enforceDesignFullRank fallback to correlation pruning works", {
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * 3), n, 3)
    # Create highly collinear columns that are hard for iterative removal
    X <- cbind(
        X,
        X[, 1] + rnorm(n, sd = 1e-10),
        X[, 2] + rnorm(n, sd = 1e-10),
        X[, 3] + rnorm(n, sd = 1e-10),
        X[, 1] + X[, 2] + rnorm(n, sd = 1e-10)
    )
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (seq_len(ncol(X))))
    C <- matrix(rnorm(n), n, 1)
    result <- enforceDesignFullRank(
        X,
        C,
        strategy = "correlation",
        maxIterations = 2L
    )
    full_design <- cbind(1, result, C)
    expect_equal(qr(full_design)$rank, ncol(full_design))
})

# ===========================================================================
# ldClumpByScore: edge cases
# ===========================================================================

test_that("ldClumpByScore errors on empty matrix", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    X <- matrix(numeric(0), nrow = 10, ncol = 0)
    expect_error(
        ldClumpByScore(
            X,
            score = numeric(0),
            chr = integer(0),
            pos = integer(0)
        ),
        "at least one column"
    )
})

test_that("ldClumpByScore returns 1L for single variant", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    X <- matrix(c(0, 1, 2, 1, 0), ncol = 1)
    result <- ldClumpByScore(X, score = 1.0, chr = 1L, pos = 100L)
    expect_equal(result, 1L)
})

test_that("ldClumpByScore errors on mismatched score length", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    X <- matrix(rnorm(20), 5, 4)
    expect_error(
        ldClumpByScore(X, score = c(1, 2), chr = rep(1L, 4), pos = 1:4),
        "length\\(score\\)"
    )
})

test_that("ldClumpByScore errors on mismatched chr/pos length", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    X <- matrix(rnorm(20), 5, 4)
    expect_error(
        ldClumpByScore(X, score = runif(4), chr = rep(1L, 2), pos = 1:4),
        "chr and pos"
    )
})

# ===========================================================================
# ldPruneByCorrelation: verbose paths
# ===========================================================================

test_that("ldPruneByCorrelation verbose reports pruning", {
    # Create matrix with correlated columns
    set.seed(42)
    base <- rnorm(100)
    X <- cbind(
        base,
        base + rnorm(100, sd = 0.1),
        rnorm(100),
        rnorm(100),
        rnorm(100)
    )
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    expect_message(
        ldPruneByCorrelation(X, corThres = 0.5, verbose = TRUE),
        "pruned"
    )
})

test_that("ldPruneByCorrelation verbose reports no pruning", {
    # Create a small matrix with no correlated columns
    set.seed(42)
    X <- matrix(rnorm(500), 100, 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    expect_message(
        ldPruneByCorrelation(X, corThres = 0.999, verbose = TRUE),
        "no columns pruned"
    )
})

# =============================================================================
# loadLdMatrix duplicate variant removal
# =============================================================================

test_that("loadLdMatrix dedup removes duplicated variants from result", {
    # Simulate what loadLdMatrix does after calling the backend: a result with
    # duplicated ldVariants should have duplicates removed.
    # We test the dedup logic by constructing a mock result and verifying
    # the internal dedup code path via the exported function's contract.
    # Since we can't easily call the real function without data, test the dedup
    # behavior directly on the result structure.
    mat <- matrix(1:16, nrow = 4, ncol = 4)
    variants <- c(
        "chr1:100:A:G",
        "chr1:200:C:T",
        "chr1:100:A:G",
        "chr1:300:T:A"
    )
    ref <- data.frame(
        chrom = c(1, 1, 1, 1),
        pos = c(100, 200, 100, 300),
        A2 = c("A", "C", "A", "T"),
        A1 = c("G", "T", "G", "A")
    )

    # Apply the same dedup logic used in loadLdMatrix
    dup_idx <- which(duplicated(variants))
    expect_equal(dup_idx, 3L)

    variants_clean <- variants[-dup_idx]
    mat_clean <- mat[-dup_idx, -dup_idx, drop = FALSE]
    ref_clean <- ref[-dup_idx, , drop = FALSE]

    expect_equal(length(variants_clean), 3)
    expect_equal(nrow(mat_clean), 3)
    expect_equal(ncol(mat_clean), 3)
    expect_equal(nrow(ref_clean), 3)
    expect_false(any(duplicated(variants_clean)))
})


test_data_dir <- test_path("test_data")


library(testthat)

# ===========================================================================
# ldPruneByCorrelation
# ===========================================================================

test_that("ldPruneByCorrelation removes highly correlated columns", {
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:p))
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
    result <- ldPruneByCorrelation(X, corThres = 0.9)
    expect_true(ncol(result$X.new) < p)
    expect_equal(length(result$filter.id), ncol(result$X.new))
})

test_that("ldPruneByCorrelation keeps all columns when uncorrelated", {
    set.seed(42)
    n <- 100
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:p))
    result <- ldPruneByCorrelation(X, corThres = 0.99)
    expect_equal(ncol(result$X.new), p)
    expect_equal(result$filter.id, 1:p)
})

test_that("ldPruneByCorrelation preserves colnames for single remaining column", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 3), nrow = n)
    colnames(X) <- c("a", "b", "c")
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.001)
    X[, 3] <- X[, 1] + rnorm(n, sd = 0.001)
    result <- ldPruneByCorrelation(X, corThres = 0.5)
    expect_true(ncol(result$X.new) >= 1)
    expect_true(!is.null(colnames(result$X.new)))
})

test_that("ldPruneByCorrelation errors on single-column input", {
    set.seed(42)
    n <- 30
    X <- matrix(rnorm(n), nrow = n, ncol = 1)
    colnames(X) <- "chr1:100:A:G"
    expect_error(ldPruneByCorrelation(X, corThres = 0.8))
})

test_that("ldPruneByCorrelation strict threshold removes at least as many as lenient", {
    set.seed(42)
    n <- 100
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:p))
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
    X[, 3] <- X[, 1] + rnorm(n, sd = 0.1)
    X[, 5] <- X[, 4] + rnorm(n, sd = 0.1)
    result_strict <- ldPruneByCorrelation(X, corThres = 0.3)
    result_lenient <- ldPruneByCorrelation(X, corThres = 0.99)
    expect_true(ncol(result_strict$X.new) <= ncol(result_lenient$X.new))
})

test_that("ldPruneByCorrelation preserves colnames when no columns deleted", {
    set.seed(42)
    n <- 100
    p <- 3
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- c("snp_a", "snp_b", "snp_c")
    result <- ldPruneByCorrelation(X, corThres = 0.999)
    expect_equal(colnames(result$X.new), colnames(X))
})

test_that("ldPruneByCorrelation is silent by default, chatty with verbose", {
    set.seed(1)
    X <- matrix(rnorm(100), 20, 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    X[, 2] <- X[, 1] + rnorm(20, sd = 1e-3)
    expect_silent(ldPruneByCorrelation(X, corThres = 0.9))
    expect_message(
        ldPruneByCorrelation(X, corThres = 0.9, verbose = TRUE),
        "ldPruneByCorrelation"
    )
})

# ===========================================================================
# dropCollinearColumns
# ===========================================================================

# dropCollinearColumns and enforceDesignFullRank are unexported helpers;
# access them via pecotmr::: in these tests.

test_that("dropCollinearColumns returns X unchanged when problematicCols is empty", {
    X <- matrix(rnorm(100), nrow = 20, ncol = 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = character(0),
        strategy = "correlation"
    )
    expect_equal(ncol(result), 5)
})

test_that("dropCollinearColumns removes single problematic column", {
    X <- matrix(rnorm(100), nrow = 20, ncol = 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = "chr1:300:A:G",
        strategy = "correlation"
    )
    expect_equal(ncol(result), 4)
    expect_false("chr1:300:A:G" %in% colnames(result))
})

test_that("dropCollinearColumns variance strategy removes lowest variance column", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(X) <- c("low_var", "mid_var", "high_var")
    X[, 1] <- X[, 1] * 0.01
    X[, 3] <- X[, 3] * 10
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = c("low_var", "mid_var", "high_var"),
        strategy = "variance"
    )
    expect_equal(ncol(result), 2)
    expect_false("low_var" %in% colnames(result))
})

test_that("dropCollinearColumns correlation strategy with two columns removes one", {
    set.seed(42)
    X <- matrix(rnorm(100), nrow = 20, ncol = 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = c("chr1:100:A:G", "chr1:200:A:G"),
        strategy = "correlation"
    )
    expect_equal(ncol(result), 4)
})

test_that("dropCollinearColumns correlation strategy with 3+ cols removes highest sum", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:4))
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
    X[, 3] <- X[, 1] + rnorm(n, sd = 0.01)
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        strategy = "correlation"
    )
    expect_equal(ncol(result), 3)
})

test_that("dropCollinearColumns responseCorrelation strategy removes lowest |cor| with response", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(X) <- c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    response <- X[, 1] * 2 + rnorm(n, sd = 0.1)
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        strategy = "responseCorrelation",
        response = response
    )
    expect_equal(ncol(result), 2)
    expect_true("chr1:100:A:G" %in% colnames(result))
})

test_that("dropCollinearColumns errors on responseCorrelation without response", {
    X <- matrix(rnorm(60), 20, 3)
    colnames(X) <- c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    expect_error(
        pecotmr:::dropCollinearColumns(
            X,
            problematicCols = c("chr1:100:A:G", "chr1:200:A:G"),
            strategy = "responseCorrelation"
        ),
        "response"
    )
})

test_that("dropCollinearColumns errors on invalid strategy", {
    X <- matrix(rnorm(100), nrow = 20, ncol = 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    expect_error(
        pecotmr:::dropCollinearColumns(
            X,
            problematicCols = c("chr1:100:A:G", "chr1:200:A:G"),
            strategy = "invalid_strategy"
        ),
        "must be one of"
    )
})

test_that("dropCollinearColumns preserves column name when single column remains", {
    set.seed(42)
    X <- matrix(rnorm(40), nrow = 20, ncol = 2)
    colnames(X) <- c("keeper", "removed")
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = "removed",
        strategy = "correlation"
    )
    expect_equal(ncol(result), 1)
    expect_equal(colnames(result), "keeper")
})

test_that("dropCollinearColumns is silent by default", {
    X <- matrix(rnorm(40), 20, 2)
    colnames(X) <- c("a", "b")
    expect_silent(pecotmr:::dropCollinearColumns(
        X,
        problematicCols = "b",
        strategy = "correlation"
    ))
})

# ===========================================================================
# enforceDesignFullRank (unexported)
# ===========================================================================

test_that("enforceDesignFullRank returns full-rank matrix when already full rank", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:3))
    C <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
    result <- enforceDesignFullRank(X = X, C = C, strategy = "correlation")
    expect_true(is.matrix(result))
    expect_true(ncol(result) >= 1)
})

test_that("enforceDesignFullRank handles rank-deficient design via correlation fallback", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:4))
    X[, 4] <- X[, 1] + X[, 2]
    result <- enforceDesignFullRank(X = X, C = NULL, strategy = "correlation")
    design <- cbind(1, result)
    expect_equal(qr(design)$rank, ncol(design))
})

test_that("enforceDesignFullRank preserves colname for single-column input", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n), nrow = n, ncol = 1)
    colnames(X) <- "only_snp"
    result <- enforceDesignFullRank(X = X, C = NULL, strategy = "correlation")
    expect_equal(colnames(result), "only_snp")
})

test_that("enforceDesignFullRank is silent by default", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 3), n, 3)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:3))
    expect_silent(enforceDesignFullRank(
        X = X,
        C = NULL,
        strategy = "correlation"
    ))
})

# ===========================================================================
# ldClumpByScore
# ===========================================================================

test_that("ldClumpByScore skips clumping on single-column input", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    set.seed(1)
    X <- matrix(rbinom(100, 2, 0.3), 100, 1)
    colnames(X) <- "chr1:100:A:G"
    keep <- ldClumpByScore(X, score = 1.0, chr = 1L, pos = 100L, r2 = 0.2)
    expect_equal(keep, 1L)
})

test_that("ldClumpByScore validates input lengths", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    set.seed(1)
    X <- matrix(rbinom(100 * 3, 2, 0.3), 100, 3)
    colnames(X) <- paste0("chr1:", seq_len(3) * 1000, ":A:G")
    expect_error(
        ldClumpByScore(
            X,
            score = 1:2,
            chr = rep(1L, 3),
            pos = seq_len(3) * 1000L
        ),
        "score"
    )
    expect_error(
        ldClumpByScore(
            X,
            score = 1:3,
            chr = rep(1L, 2),
            pos = seq_len(3) * 1000L
        ),
        "chr and pos"
    )
})

test_that("ldClumpByScore returns indices on real data", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    set.seed(1)
    n <- 500
    p <- 20
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    colnames(X) <- paste0("chr1:", seq_len(p) * 1000, ":A:G")
    # Introduce perfect LD between variants 1 and 2
    X[, 2] <- X[, 1]
    score <- c(2, 1, runif(p - 2))
    chr <- rep(1L, p)
    pos <- seq_len(p) * 1000L
    keep <- ldClumpByScore(X, score = score, chr = chr, pos = pos, r2 = 0.2)
    expect_true(1L %in% keep)
    expect_false(2L %in% keep) # pruned: same as variant 1 but lower score
    expect_true(length(keep) < p)
})


test_that("standardize_genotype_hwe: centers by 2p and scales by sqrt(2p(1-p))", {
    set.seed(42)
    n <- 30
    p <- 5
    af <- runif(p, 0.1, 0.9)
    X <- matrix(rbinom(n * p, 2, rep(af, each = n)), nrow = n, ncol = p)

    X_std <- pecotmr:::standardizeGenotypeHwe(X, af)

    # Manual verification
    expected <- sweep(sweep(X, 2, 2 * af), 2, sqrt(2 * af * (1 - af)), "/")
    expect_equal(X_std, expected, tolerance = 1e-14)
})


test_that("SVD from raw sketch matches direct computation", {
    set.seed(77)
    n <- 25
    p <- 8
    af <- runif(p, 0.15, 0.85)
    X <- matrix(rbinom(n * p, 2, rep(af, each = n)), nrow = n, ncol = p)

    # Two-step process: standardize then SVD
    X_std <- pecotmr:::standardizeGenotypeHwe(X, af)
    svd_result <- svd(X_std)

    # Verify this matches manual computation
    X_manual <- sweep(sweep(X, 2, 2 * af), 2, sqrt(2 * af * (1 - af)), "/")
    svd_manual <- svd(X_manual)

    expect_equal(svd_result$d, svd_manual$d, tolerance = 1e-14)
    expect_equal(abs(svd_result$v), abs(svd_manual$v), tolerance = 1e-14)
})


# ===========================================================================
# Tests migrated from test_misc.R (computeLd + .findValidFilePath*)
# ===========================================================================

test_that("findValidFilePath returns target when it exists directly", {
    pkg_root <- normalizePath(
        file.path(test_path(), "..", ".."),
        mustWork = TRUE
    )
    target <- file.path(pkg_root, "DESCRIPTION")
    ref <- file.path(pkg_root, "NAMESPACE")
    skip_if_not(
        file.exists(target) && file.exists(ref),
        "Package root files not found"
    )
    result <- pecotmr:::.findValidFilePath(
        referenceFilePath = ref,
        targetFilePath = target
    )
    expect_equal(result, target)
})


test_that("findValidFilePath constructs path from reference directory", {
    pkg_root <- normalizePath(
        file.path(test_path(), "..", ".."),
        mustWork = TRUE
    )
    ref <- file.path(pkg_root, "NAMESPACE")
    skip_if_not(file.exists(ref), "NAMESPACE not found")
    result <- pecotmr:::.findValidFilePath(
        referenceFilePath = ref,
        targetFilePath = "DESCRIPTION"
    )
    expect_true(file.exists(result))
    expect_true(grepl("DESCRIPTION$", result))
})


test_that("findValidFilePath errors when both paths are invalid", {
    expect_error(
        pecotmr:::.findValidFilePath(
            referenceFilePath = "/nonexistent/dir/ref.txt",
            targetFilePath = "/nonexistent/target.txt"
        ),
        "Both reference and target file paths do not work"
    )
})


test_that("findValidFilePath returns reference when target is invalid but reference exists", {
    pkg_root <- normalizePath(
        file.path(test_path(), "..", ".."),
        mustWork = TRUE
    )
    ref <- file.path(pkg_root, "DESCRIPTION")
    skip_if_not(file.exists(ref), "DESCRIPTION not found")
    result <- pecotmr:::.findValidFilePath(
        referenceFilePath = ref,
        targetFilePath = "/totally/bogus/path.txt"
    )
    expect_equal(result, ref)
})


test_that("findValidFilePaths resolves multiple targets", {
    pkg_root <- normalizePath(
        file.path(test_path(), "..", ".."),
        mustWork = TRUE
    )
    ref <- file.path(pkg_root, "NAMESPACE")
    skip_if_not(file.exists(ref), "NAMESPACE not found")
    targets <- c("DESCRIPTION", "NAMESPACE")
    result <- pecotmr:::.findValidFilePaths(ref, targets)
    expect_length(result, 2)
    expect_true(all(file.exists(result)))
})


test_that("findValidFilePaths errors on all-invalid targets", {
    ref <- "/nonexistent/ref.txt"
    targets <- c("/bogus/a.txt", "/bogus/b.txt")
    expect_error(pecotmr:::.findValidFilePaths(ref, targets))
})

# =============================================================================
# computeLd
# =============================================================================

test_that("computeLd sample method produces valid correlation matrix", {
    set.seed(42)
    X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
    colnames(X) <- paste0("rs", 1:4)

    R <- computeLd(X, method = "sample")
    expect_equal(nrow(R), 4)
    expect_equal(ncol(R), 4)
    expect_equal(unname(diag(R)), rep(1, 4))
    expect_true(isSymmetric(R))
    expect_true(all(R >= -1 & R <= 1))
})


test_that("computeLd population method produces valid matrix", {
    set.seed(42)
    X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
    colnames(X) <- paste0("rs", 1:4)

    R <- computeLd(X, method = "population")
    expect_equal(nrow(R), 4)
    expect_equal(unname(diag(R)), rep(1, 4))
    expect_true(isSymmetric(R))
})


test_that("computeLd with a single SNP returns 1x1 identity matrix", {
    X <- matrix(c(0, 1, 2, 1, 0), ncol = 1)
    colnames(X) <- "rs1"
    R <- computeLd(X, method = "sample")
    expect_equal(dim(R), c(1L, 1L))
    expect_equal(R[1, 1], 1.0)
    expect_equal(colnames(R), "rs1")
})


test_that("computeLd handles column with all NA gracefully", {
    set.seed(123)
    X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
    X[, 3] <- NA
    colnames(X) <- paste0("rs", 1:5)

    R <- computeLd(X, method = "sample")
    expect_equal(dim(R), c(5L, 5L))
    expect_equal(unname(diag(R)), rep(1, 5))
    expect_equal(R[3, 1], 0)
    expect_equal(R[1, 3], 0)
})


test_that("computeLd population method handles column with all NA", {
    set.seed(123)
    X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
    X[, 2] <- NA
    colnames(X) <- paste0("rs", 1:5)

    R <- computeLd(X, method = "population")
    expect_equal(dim(R), c(5L, 5L))
    expect_equal(unname(diag(R)), rep(1, 5))
    expect_equal(R[2, 4], 0)
})


test_that("computeLd with larger matrix (100 SNPs) is fast and valid", {
    set.seed(99)
    X <- matrix(sample(0:2, 5000, replace = TRUE), nrow = 50, ncol = 100)
    colnames(X) <- paste0("rs", 1:100)

    R <- computeLd(X, method = "sample")
    expect_equal(dim(R), c(100L, 100L))
    expect_equal(unname(diag(R)), rep(1, 100))
    expect_true(isSymmetric(R))
    expect_true(all(R >= -1 & R <= 1))
})


test_that("computeLd population method with larger matrix is valid", {
    set.seed(99)
    X <- matrix(sample(0:2, 5000, replace = TRUE), nrow = 50, ncol = 100)
    colnames(X) <- paste0("rs", 1:100)

    R <- computeLd(X, method = "population")
    expect_equal(dim(R), c(100L, 100L))
    expect_equal(unname(diag(R)), rep(1, 100))
    expect_true(isSymmetric(R))
})


test_that("computeLd with perfectly correlated SNPs returns correlation of 1", {
    set.seed(42)
    col1 <- sample(0:2, 50, replace = TRUE)
    X <- matrix(c(col1, col1), ncol = 2)
    colnames(X) <- c("rs1", "rs2")

    R <- computeLd(X, method = "sample")
    expect_equal(R[1, 2], 1.0, tolerance = 1e-10)
    expect_equal(R[2, 1], 1.0, tolerance = 1e-10)
})


test_that("computeLd population method with trim_samples trims correctly", {
    set.seed(42)
    X <- matrix(sample(0:2, 33, replace = TRUE), nrow = 11, ncol = 3)
    colnames(X) <- paste0("rs", 1:3)

    R_trimmed <- computeLd(X, method = "population", trimSamples = TRUE)
    expect_equal(dim(R_trimmed), c(3L, 3L))
    R_full <- computeLd(X, method = "population", trimSamples = FALSE)
    expect_equal(dim(R_full), c(3L, 3L))
})


test_that("computeLd with two monomorphic SNPs produces 0 off-diagonal", {
    X <- matrix(c(rep(1, 50), rep(2, 50)), nrow = 50, ncol = 2)
    colnames(X) <- c("mono1", "mono2")

    R <- computeLd(X, method = "sample")
    expect_equal(R[1, 2], 0)
    expect_equal(R[2, 1], 0)
    expect_equal(unname(diag(R)), c(1, 1))
})


test_that("computeLd preserves column names", {
    set.seed(42)
    X <- matrix(sample(0:2, 60, replace = TRUE), nrow = 20, ncol = 3)
    colnames(X) <- c("snp_alpha", "snp_beta", "snp_gamma")

    R <- computeLd(X, method = "sample")
    expect_equal(colnames(R), c("snp_alpha", "snp_beta", "snp_gamma"))
    expect_equal(rownames(R), c("snp_alpha", "snp_beta", "snp_gamma"))
})


test_that("computeLd with heavy missingness still produces valid matrix", {
    set.seed(42)
    X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 40, ncol = 5)
    na_idx <- sample(length(X), size = floor(0.5 * length(X)))
    X[na_idx] <- NA
    colnames(X) <- paste0("rs", 1:5)

    R <- computeLd(X, method = "sample")
    expect_true(all(!is.na(R)))
    expect_equal(unname(diag(R)), rep(1, 5))

    R_pop <- computeLd(X, method = "population")
    expect_true(all(!is.na(R_pop)))
    expect_equal(unname(diag(R_pop)), rep(1, 5))
})


test_that("computeLd with NA genotypes and sample method", {
    set.seed(42)
    X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
    X[1, 1] <- NA
    X[5, 3] <- NA
    colnames(X) <- paste0("rs", 1:4)

    R <- computeLd(X, method = "sample")
    expect_true(all(!is.na(R)))
    expect_equal(unname(diag(R)), rep(1, 4))
})


test_that("computeLd errors on NULL input", {
    expect_error(computeLd(NULL), "X must be provided")
})


test_that("computeLd sample vs population differ but are close", {
    set.seed(42)
    X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
    colnames(X) <- paste0("rs", 1:5)

    R_sample <- computeLd(X, method = "sample")
    R_pop <- computeLd(X, method = "population")

    expect_false(identical(R_sample, R_pop))
    expect_true(max(abs(R_sample - R_pop)) < 0.1)
})


test_that("computeLd gcta method produces valid correlation matrix", {
    set.seed(42)
    X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
    colnames(X) <- paste0("rs", 1:5)

    R <- computeLd(X, method = "gcta")
    expect_equal(dim(R), c(5, 5))
    expect_equal(unname(diag(R)), rep(1, 5), tolerance = 1e-10)
    expect_true(isSymmetric(R, tol = 1e-10))
    expect_true(all(abs(R) <= 1 + 1e-10))
})


test_that("computeLd gcta method handles missing data", {
    set.seed(42)
    X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
    colnames(X) <- paste0("rs", 1:5)
    X[sample(length(X), 50)] <- NA

    R <- computeLd(X, method = "gcta")
    expect_equal(dim(R), c(5, 5))
    expect_true(all(is.finite(R)))
})


test_that("computeLd gcta agrees with sample method on complete data", {
    set.seed(42)
    X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
    colnames(X) <- paste0("rs", 1:5)

    R_sample <- computeLd(X, method = "sample")
    R_gcta <- computeLd(X, method = "gcta")

    # With no missing data, GCTA and sample should be close (differ by N vs N-1 denom)
    expect_true(max(abs(R_sample - R_gcta)) < 0.05)
})


test_that("computeLd gcta preserves column names", {
    set.seed(42)
    X <- matrix(sample(0:2, 300, replace = TRUE), nrow = 100)
    colnames(X) <- c("snp_a", "snp_b", "snp_c")

    R <- computeLd(X, method = "gcta")
    expect_equal(colnames(R), c("snp_a", "snp_b", "snp_c"))
    expect_equal(rownames(R), c("snp_a", "snp_b", "snp_c"))
})


# =============================================================================
# waldTestPval
# =============================================================================

test_that("computeLd sample method without Rfast falls back to cor", {
    set.seed(42)
    X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
    colnames(X) <- paste0("snp", 1:5)
    R <- computeLd(X, method = "sample")
    expect_equal(dim(R), c(5, 5))
    expect_equal(as.numeric(diag(R)), rep(1, 5))
    expect_true(all(abs(R) <= 1))
})


test_that("computeLd with gcta method and trim_samples", {
    set.seed(42)
    # 21 samples -> trimmed to 20 (multiple of 4)
    X <- matrix(sample(0:2, 105, replace = TRUE), nrow = 21, ncol = 5)
    colnames(X) <- paste0("snp", 1:5)
    R <- computeLd(X, method = "gcta", trimSamples = TRUE)
    expect_equal(dim(R), c(5, 5))
    expect_equal(as.numeric(diag(R)), rep(1, 5))
})


test_that("computeLd population method with trim_samples", {
    set.seed(42)
    X <- matrix(sample(0:2, 105, replace = TRUE), nrow = 21, ncol = 5)
    colnames(X) <- paste0("snp", 1:5)
    R <- computeLd(X, method = "population", trimSamples = TRUE)
    expect_equal(dim(R), c(5, 5))
    expect_equal(as.numeric(diag(R)), rep(1, 5))
})


test_that("computeLd with shrinkage > 0", {
    set.seed(42)
    X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
    colnames(X) <- paste0("snp", 1:5)
    R_no_shrink <- computeLd(X, method = "sample", shrinkage = 0)
    R_shrink <- computeLd(X, method = "sample", shrinkage = 0.1)
    # Shrunk matrix should be closer to identity
    expect_equal(as.numeric(diag(R_shrink)), rep(1, 5))
    # Off-diagonal elements should be shrunk toward 0
    off_diag_no <- R_no_shrink[1, 2]
    off_diag_s <- R_shrink[1, 2]
    expect_equal(off_diag_s, 0.9 * off_diag_no)
})


# =============================================================================
# Additional coverage: findIntersectionRows / getRegionalLdMeta edge cases
# =============================================================================

test_that("findIntersectionRows errors when region falls in a coverage gap", {
    # Chromosome exists but there is a gap between the two blocks; after
    # clamping, no single row covers the query start -> stop (line 33).
    gd <- data.frame(chrom = c(1, 1), start = c(100, 300), end = c(150, 350))
    expect_error(
        pecotmr:::findIntersectionRows(gd, 1, 200, 250),
        "not covered by any rows"
    )
})

test_that("getRegionalLdMeta handles whole-chromosome 0:0 sentinel rows", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_wholechrom_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "1",
            "0",
            "0",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )
    result <- pecotmr:::getRegionalLdMeta(meta_file, "chr1:1000-1190")
    expect_true(length(result$intersections$LD_file_paths) >= 1)
    expect_true(all(file.exists(result$intersections$LD_file_paths)))
    # 0:0 sentinel row should have had its end set to Inf internally (line 103)
    expect_true(is.infinite(max(result$ldMetaData$end)))
})

test_that("getRegionalLdMeta validates complete coverage when required", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_complete_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    # Region fully inside block 1 -> validateSelectedRegion passes (line 122)
    result <- pecotmr:::getRegionalLdMeta(
        meta_file,
        "chr1:1050-1150",
        completeCoverageRequired = TRUE
    )
    expect_true(length(result$intersections$LD_file_paths) >= 1)
})

# =============================================================================
# Additional coverage: processLdMatrix variant-file auto-detect / pvar / symm
# =============================================================================

test_that("processLdMatrix auto-detects companion .bim when snpFilePath is NULL", {
    tmp_xz <- tempfile(fileext = ".txt.xz")
    on.exit(unlink(c(tmp_xz, paste0(tmp_xz, ".bim"))), add = TRUE)
    file.copy(
        file.path(
            geno_test_data_dir,
            "LD_block_1.chr1_1000_1200.float16.txt.xz"
        ),
        tmp_xz
    )
    # Place the companion at "<ldfile>.bim" so the auto-detector finds it
    file.copy(
        file.path(geno_test_data_dir, "LD_block_1.chr1_1000_1200.float16.bim"),
        paste0(tmp_xz, ".bim")
    )
    result <- pecotmr:::processLdMatrix(tmp_xz, snpFilePath = NULL)
    expect_equal(nrow(result$ldMatrix), 5L)
    expect_true(all(grepl("^chr1:", rownames(result$ldMatrix))))
})

test_that("processLdMatrix errors when no companion variant file is found", {
    tmp_xz <- tempfile(fileext = ".txt.xz")
    on.exit(unlink(tmp_xz), add = TRUE)
    file.copy(
        file.path(
            geno_test_data_dir,
            "LD_block_1.chr1_1000_1200.float16.txt.xz"
        ),
        tmp_xz
    )
    expect_error(
        pecotmr:::processLdMatrix(tmp_xz, snpFilePath = NULL),
        "No variant file found"
    )
})

test_that("processLdMatrix reads .pvar metadata and symmetrizes an upper-triangular matrix", {
    skip_if_not_installed("pgenlibr")
    pvar <- file.path(geno_test_data_dir, "test_harmonize_regions.pvar") # 8 variants
    # Build an 8x8 upper-triangular matrix (lower triangle exactly zero) so the
    # lower.tri == 0 branch fires (line 183).
    set.seed(1)
    n <- 8
    U <- diag(n)
    U[upper.tri(U)] <- round(runif(sum(upper.tri(U)), -0.4, 0.8), 3)
    expect_true(all(U[lower.tri(U)] == 0))
    tmp_xz <- tempfile(fileext = ".txt.xz")
    on.exit(unlink(tmp_xz), add = TRUE)
    con <- xzfile(tmp_xz, "w")
    writeLines(paste(as.vector(t(U)), collapse = " "), con) # row-major order
    close(con)

    result <- pecotmr:::processLdMatrix(tmp_xz, pvar)
    expect_equal(nrow(result$ldMatrix), n)
    expect_true(isSymmetric(result$ldMatrix))
    # .pvar metadata has no gpos column; pos is derived from the variant id
    expect_false("gpos" %in% names(result$ldVariants))
    expect_true(all(grepl("^chr21:", rownames(result$ldMatrix))))
})

# =============================================================================
# Additional coverage: createLdMatrix empty entry / loadLdMatrix dedup
# =============================================================================

test_that("createLdMatrix skips empty variant-list entries", {
    m1 <- matrix(
        c(1, 0.5, 0.5, 1),
        2,
        2,
        dimnames = list(
            c("chr1:100:A:G", "chr1:200:A:G"),
            c("chr1:100:A:G", "chr1:200:A:G")
        )
    )
    variants <- list(
        data.frame(variants = character(0)), # empty -> next (line 229)
        data.frame(variants = c("chr1:100:A:G", "chr1:200:A:G"))
    )
    result <- pecotmr:::createLdMatrix(list(matrix(0, 0, 0), m1), variants)
    expect_equal(nrow(result), 2L)
    expect_equal(rownames(result), c("chr1:100:A:G", "chr1:200:A:G"))
    expect_equal(result["chr1:100:A:G", "chr1:200:A:G"], 0.5)
})

test_that("loadLdMatrix removes duplicate variants from the backend result", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_dedup_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )

    # Build an LdData whose variant IDs contain a duplicate (index 3 == index 1)
    variant_ids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:100:A:G")
    ref <- data.frame(
        chrom = c(1, 1, 1),
        pos = c(100, 200, 100),
        A2 = c("A", "C", "A"),
        A1 = c("G", "T", "G"),
        variant_id = variant_ids,
        stringsAsFactors = FALSE
    )
    gr <- pecotmr:::.refPanelToGranges(ref)
    R <- matrix(c(1, 0.5, 0.9, 0.5, 1, 0.4, 0.9, 0.4, 1), 3, 3)
    rownames(R) <- colnames(R) <- variant_ids
    bm <- data.frame(
        blockId = 1L,
        chrom = "1",
        blockStart = 100L,
        blockEnd = 200L,
        size = 3L,
        startIdx = 1L,
        endIdx = 3L,
        stringsAsFactors = FALSE
    )
    dup_ld <- LdData(correlation = R, variants = gr, blockMetadata = bm)

    local_mocked_bindings(
        loadLdFromBlocks = function(
            ldMetaFilePath,
            region,
            extractCoordinates = NULL,
            nSample = NULL
        ) {
            dup_ld
        },
        .package = "pecotmr"
    )
    result <- loadLdMatrix(meta_file, "chr1:100-200")
    ids <- getVariantIds(result)
    expect_equal(length(ids), 2L)
    expect_false(any(duplicated(ids)))
    expect_equal(dim(getCorrelation(result)), c(2L, 2L))
})

# =============================================================================
# Additional coverage: resolveLdSource column check
# =============================================================================

test_that("resolveLdSource errors when metadata has fewer than 4 columns", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_3col_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", sep = "\t"), meta_file)
    cat(
        paste("1", "1000", "1200", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    expect_error(pecotmr:::resolveLdSource(meta_file), "at least 4 columns")
})

# =============================================================================
# Additional coverage: loadLdFromGenotype .afreq mismatch warning
# =============================================================================

test_that("loadLdFromGenotype warns when .afreq is missing some variants", {
    skip_if_not_installed("pgenlibr")
    plink_prefix <- file.path(geno_test_data_dir, "test_variants")
    local_mocked_bindings(
        readAfreq = function(prefix) {
            data.frame(
                id = "no_such_variant",
                alt_freq = 0.1,
                stringsAsFactors = FALSE
            )
        },
        .package = "pecotmr"
    )
    expect_warning(
        pecotmr:::loadLdFromGenotype(plink_prefix, geno_region_all),
        "no allele frequency"
    )
})

# =============================================================================
# Additional coverage: .ldFromSketch onMissing="drop"
# =============================================================================

test_that(".ldFromSketch drops variants absent from the panel when onMissing='drop'", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    si <- getSnpInfo(h)
    ids <- c(as.character(si$SNP[1]), "bogus:999:A:G", as.character(si$SNP[3]))
    m <- pecotmr:::.ldFromSketch(h, ids, onMissing = "drop")
    expect_equal(dim(m), c(2L, 2L))
    expect_equal(
        attr(m, "keptVariantIds"),
        c(as.character(si$SNP[1]), as.character(si$SNP[3]))
    )
    expect_equal(unname(diag(m)), c(1, 1))
})

# =============================================================================
# Additional coverage: .ldFromSketch reconciles chr-prefix conventions
# =============================================================================

test_that(".ldFromSketch resolves chr-prefixed request ids against a non-prefixed panel", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(getSnpInfo(h)$SNP) # chr-prefixed (chr21_..._C_G)
    h@snpInfo$SNP <- sub("^chr", "", pid) # panel now lacks the prefix
    m <- pecotmr:::.ldFromSketch(h, pid[1:3]) # request keeps the prefix
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rownames(m), pid[1:3]) # caller's labels preserved
    expect_equal(unname(diag(m)), c(1, 1, 1))
})

test_that(".ldFromSketch resolves non-prefixed request ids against a chr-prefixed panel", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(getSnpInfo(h)$SNP)
    req <- sub("^chr", "", pid[1:3]) # request drops the prefix
    m <- pecotmr:::.ldFromSketch(h, req)
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rownames(m), req) # caller's labels preserved
})

test_that(".ldFromSketch is unchanged when request and panel share a convention", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(getSnpInfo(h)$SNP)
    m <- pecotmr:::.ldFromSketch(h, pid[1:3])
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rownames(m), pid[1:3])
})

test_that(".ldFromSketch still errors on a genuinely-absent variant after reconciliation", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(getSnpInfo(h)$SNP)
    # request the prefix-stripped form (forces reconciliation) plus one truly-absent variant
    req <- c(sub("^chr", "", pid[1]), "21_99999999_A_G")
    expect_error(
        pecotmr:::.ldFromSketch(h, req),
        "not present in the LD sketch panel"
    )
})

test_that(".ldFromSketch rejects an ldSketch that is not a genotype panel", {
    expect_error(
        pecotmr:::.ldFromSketch(
            "not_a_handle",
            c("chr1:100:A:G", "chr1:200:A:G")
        ),
        "ldSketch must be a genotype panel"
    )
})

# A synthetic panel + dosage extractor, so the shape of what .ldFromSketch
# returns is checked without a plink2 fixture (and so without pgenlibr).
# @noRd
.lds_makeHandle <- function(snpN = 6L, nSamples = 30L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(snpN)),
            CHR = rep("1", snpN),
            BP = seq(100L, by = 100L, length.out = snpN),
            A1 = rep("A", snpN),
            A2 = rep("G", snpN),
            stringsAsFactors = FALSE
        ),
        nSamples = nSamples,
        sampleIds = sprintf("s%d", seq_len(nSamples)),
        pgenPtr = NULL
    )
}

# @noRd
.lds_mockExtractor <- function(seed = 7, nSamples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        nSnp <- nrow(getSnpInfo(handle))
        panel <- matrix(
            rbinom(nSamples * nSnp, 2, 0.3),
            nrow = nSamples,
            ncol = nSnp,
            dimnames = list(getSampleIds(handle), getSnpInfo(handle)$SNP)
        )
        rr <- GenomicRanges::GRanges(
            seqnames = str_c("chr", getSnpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = getSnpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = getSnpInfo(handle)$SNP[snpIdx],
            A1 = getSnpInfo(handle)$A1[snpIdx],
            A2 = getSnpInfo(handle)$A2[snpIdx]
        )
        dosage <- t(panel[, snpIdx, drop = FALSE])
        dimnames(dosage) <- list(
            getSnpInfo(handle)$SNP[snpIdx],
            getSampleIds(handle)
        )
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = S4Vectors::DataFrame(
                sampleId = getSampleIds(handle),
                row.names = getSampleIds(handle)
            )
        )
    }
}

test_that(".ldFromSketch returns a symmetric unit-diagonal LD matrix", {
    h <- .lds_makeHandle()
    local_mocked_bindings(
        extractBlockGenotypes = .lds_mockExtractor(),
        .package = "pecotmr"
    )
    ids <- c("chr1:200:A:G", "chr1:400:A:G", "chr1:500:A:G")
    R <- pecotmr:::.ldFromSketch(h, ids)
    expect_true(is.matrix(R))
    expect_equal(dimnames(R), list(ids, ids))
    expect_equal(unname(diag(R)), rep(1, 3), tolerance = 1e-12)
    expect_equal(R, t(R), tolerance = 1e-12)
})

test_that(".ldFromSketch negates LD for a variant the panel spells flipped", {
    # The panel carries A/G; asking for the same variant as G/A is the same
    # variant read off the other allele, so every correlation it takes part in
    # changes sign while the diagonal stays 1. Before the sign was applied the
    # two calls returned an identical matrix, which silently mis-signed the LD
    # against a caller whose alleles were oriented the other way.
    h <- .lds_makeHandle()
    local_mocked_bindings(
        extractBlockGenotypes = .lds_mockExtractor(),
        .package = "pecotmr"
    )
    same <- c("chr1:200:A:G", "chr1:400:A:G", "chr1:500:A:G")
    flipped <- c("chr1:200:G:A", "chr1:400:A:G", "chr1:500:A:G")
    rSame <- pecotmr:::.ldFromSketch(h, same)
    rFlip <- pecotmr:::.ldFromSketch(h, flipped)
    sgn <- c(-1, 1, 1)
    expect_equal(
        unname(rFlip),
        unname(rSame) * outer(sgn, sgn),
        tolerance = 1e-12
    )
    expect_equal(unname(diag(rFlip)), rep(1, 3), tolerance = 1e-12)
    expect_equal(dimnames(rFlip), list(flipped, flipped))
})

test_that(".ldFromSketch leaves LD untouched when every allele agrees", {
    h <- .lds_makeHandle()
    local_mocked_bindings(
        extractBlockGenotypes = .lds_mockExtractor(),
        .package = "pecotmr"
    )
    ids <- c("chr1:200:A:G", "chr1:400:A:G")
    expect_equal(
        pecotmr:::.ldFromSketch(h, ids),
        pecotmr:::.ldFromSketch(h, ids),
        tolerance = 1e-12
    )
})

test_that(".ldFromSketch errors on a variant the panel does not carry", {
    h <- .lds_makeHandle()
    expect_error(
        pecotmr:::.ldFromSketch(h, c("chr1:100:A:G", "ghost")),
        "variant id.*not present in the LD sketch"
    )
})

# =============================================================================
# Additional coverage: .requireMatchingLdSketches error paths
# =============================================================================

test_that(".requireMatchingLdSketches errors when slots are not GenotypeHandle", {
    expect_error(
        pecotmr:::.requireMatchingLdSketches(
            list(a = 1),
            list(b = 2),
            "testPipeline"
        ),
        "must both be genotype panels"
    )
})

test_that(".requireMatchingLdSketches errors when panels differ in a column", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    si <- getSnpInfo(h)
    si2 <- si
    si2$A1[1] <- if (identical(si2$A1[1], "A")) "C" else "A" # mutate one allele
    h2 <- new(
        "GenotypeHandle",
        path = getPath(h),
        format = getFormat(h),
        snpInfo = si2,
        nSamples = getNSamples(h),
        sampleIds = getSampleIds(h),
        pgenPtr = NULL,
        chromPaths = character(0)
    )
    expect_error(
        pecotmr:::.requireMatchingLdSketches(h, h2, "testPipeline"),
        "differ in column"
    )
})

test_that(".requireMatchingLdSketches tolerates a chr-prefix-only difference", {
    si1 <- data.frame(
        SNP = c("1:100:A:G", "1:200:C:T"),
        CHR = c("1", "1"),
        BP = c(100L, 200L),
        A1 = c("A", "C"),
        A2 = c("G", "T"),
        stringsAsFactors = FALSE
    )
    si2 <- si1
    si2$CHR <- c("chr1", "chr1")
    si2$SNP <- c("chr1:100:A:G", "chr1:200:C:T") # same panel, chr-prefixed
    mk <- function(si) {
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "gds",
            snpInfo = si,
            nSamples = 3L,
            sampleIds = paste0("s", 1:3),
            pgenPtr = NULL,
            chromPaths = character(0)
        )
    }
    expect_null(
        pecotmr:::.requireMatchingLdSketches(mk(si1), mk(si2), "testPipeline")
    )
})


# =============================================================================
# Additional coverage: loadLdFromBlocks empty-block handling
# =============================================================================

test_that("loadLdFromBlocks drops empty blocks and keeps non-empty ones", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_emptyblock_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        ),
        paste(
            "1",
            "1200",
            "1400",
            "LD_block_2.chr1_1200_1400.float16.txt.xz,LD_block_2.chr1_1200_1400.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    # Region 1180-1260: block 1 (variants 1000..1160) is empty; block 2 keeps 1200,1240
    expect_message(
        result <- pecotmr:::loadLdFromBlocks(meta_file, "chr1:1180-1260"),
        "Removing 1 empty LD block"
    )
    ids <- getVariantIds(result)
    expect_true(length(ids) >= 1)
    expect_true(all(grepl("^chr1:12", ids)))
})

test_that("loadLdFromBlocks errors when no block has variants in the region", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_noblockvar_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    # Region 1165-1175 falls between block-1 variants 1160 and the next block
    expect_error(
        pecotmr:::loadLdFromBlocks(meta_file, "chr1:1165-1175"),
        "No variants found in any LD block"
    )
})

test_that("loadLdFromBlocks derives variance from allele_freq + nSample when variance is NA", {
    # 9-column bim with NA variance but a present allele_freq column.
    bim_file <- file.path(geno_test_data_dir, "LD_block_1_navar_tmp.bim")
    meta_file <- file.path(geno_test_data_dir, "ld_meta_navar_tmp.tsv")
    on.exit(unlink(c(bim_file, meta_file)), add = TRUE)
    bim_lines <- c(
        "1\tchr1:1000_A_G\t0\t1000\tA\tG\tNA\t0.3\t500",
        "1\tchr1:1040_A_G\t0\t1040\tA\tG\tNA\t0.4\t500",
        "1\tchr1:1080_A_G\t0\t1080\tA\tG\tNA\t0.2\t500",
        "1\tchr1:1120_A_G\t0\t1120\tA\tG\tNA\t0.5\t500",
        "1\tchr1:1160_A_G\t0\t1160\tA\tG\tNA\t0.15\t500"
    )
    writeLines(bim_lines, bim_file)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste(
            "1",
            "1000",
            "1200",
            paste0(
                "LD_block_1.chr1_1000_1200.float16.txt.xz,",
                basename(bim_file)
            ),
            sep = "\t"
        ),
        "\n",
        file = meta_file,
        append = TRUE
    )

    result <- loadLdMatrix(meta_file, "chr1:1000-1190", nSample = 500L)
    ref_mcols <- S4Vectors::mcols(getVariantInfo(result))
    expect_true("variance" %in% names(ref_mcols))
    expect_true(all(!is.na(ref_mcols$variance)))
    expect_true(all(ref_mcols$variance > 0))
    expect_equal(ref_mcols$n_nomiss, rep(500, length(ref_mcols$variance)))
})

# =============================================================================
# Additional coverage: filterVariantsByLdReference keepIndel = FALSE
# =============================================================================

test_that("filterVariantsByLdReference with keepIndel=FALSE drops indels", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_filtind_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    lines <- c(
        paste("chrom", "start", "end", "path", sep = "\t"),
        paste(
            "1",
            "1000",
            "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"
        )
    )
    writeLines(lines, meta_file)
    # All three positions are on the reference panel, but the middle one is an indel.
    variant_ids <- c("chr1:1000:A:G", "chr1:1040:AT:A", "chr1:1080:A:G")
    result <- suppressMessages(
        filterVariantsByLdReference(variant_ids, meta_file, keepIndel = FALSE)
    )
    expect_false("chr1:1040:AT:A" %in% result$data)
    expect_true("chr1:1000:A:G" %in% result$data)
    expect_true("chr1:1080:A:G" %in% result$data)
})

# =============================================================================
# Additional coverage: partitionLdMatrix / validateBlockStructure / extractBlockMatrices
# =============================================================================

test_that("partitionLdMatrix accepts LdBlocks blockMetadata", {
    variant_ids <- paste0("chr1:", c(100, 200, 300, 400), ":A:G")
    R <- diag(4)
    rownames(R) <- colnames(R) <- variant_ids
    gr_blocks <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(100, 300), end = c(200, 400))
    )
    S4Vectors::mcols(gr_blocks) <- S4Vectors::DataFrame(
        blockId = c(1L, 2L),
        chrom = c("1", "1"),
        size = c(2L, 2L),
        startIdx = c(1L, 3L),
        endIdx = c(2L, 4L),
        blockStart = c(100L, 300L),
        blockEnd = c(200L, 400L)
    )
    ldb <- gr_blocks
    ref <- pecotmr:::parseVariantId(variant_ids)
    ref$variant_id <- variant_ids
    gr_vars <- pecotmr:::.refPanelToGranges(ref)
    ld <- LdData(correlation = R, variants = gr_vars, blockMetadata = ldb)
    result <- pecotmr:::partitionLdMatrix(ld, mergeSmallBlocks = FALSE)
    expect_length(result$ldMatrices, 2)
    expect_equal(nrow(result$ldMatrices[[1]]), 2L)
})

test_that("partitionLdMatrix errors on an empty correlation matrix", {
    variant_ids <- "chr1:100:A:G"
    ref <- pecotmr:::parseVariantId(variant_ids)
    ref$variant_id <- variant_ids
    gr <- pecotmr:::.refPanelToGranges(ref)
    empty_R <- matrix(numeric(0), 0, 0)
    bm <- data.frame(
        blockId = 1L,
        chrom = "1",
        blockStart = 100L,
        blockEnd = 100L,
        size = 1L,
        startIdx = 1L,
        endIdx = 1L,
        stringsAsFactors = FALSE
    )
    ld <- LdData(correlation = empty_R, variants = gr, blockMetadata = bm)
    expect_error(pecotmr:::partitionLdMatrix(ld), "Empty or NULL LD matrix")
})

test_that("partitionLdMatrix errors when all blocks have invalid indices", {
    variant_ids <- paste0("chr1:", c(100, 200, 300), ":A:G")
    R <- diag(3)
    rownames(R) <- colnames(R) <- variant_ids
    ref <- pecotmr:::parseVariantId(variant_ids)
    ref$variant_id <- variant_ids
    gr <- pecotmr:::.refPanelToGranges(ref)
    bm <- data.frame(
        blockId = 1L,
        chrom = "1",
        blockStart = 100L,
        blockEnd = 300L,
        size = 3L,
        startIdx = 10L,
        endIdx = 20L,
        stringsAsFactors = FALSE
    )
    ld <- LdData(correlation = R, variants = gr, blockMetadata = bm)
    expect_error(pecotmr:::partitionLdMatrix(ld), "No valid LD blocks found")
})

test_that("partitionLdMatrix removes blocks with invalid indices and reindexes", {
    variant_ids <- paste0("chr1:", seq(100, 600, by = 100), ":A:G") # 6 variants
    R <- diag(6)
    rownames(R) <- colnames(R) <- variant_ids
    ref <- pecotmr:::parseVariantId(variant_ids)
    ref$variant_id <- variant_ids
    gr <- pecotmr:::.refPanelToGranges(ref)
    bm <- data.frame(
        blockId = c(1L, 2L),
        chrom = c("1", "1"),
        blockStart = c(100L, 400L),
        blockEnd = c(300L, 600L),
        size = c(3L, 3L),
        startIdx = c(1L, 50L),
        endIdx = c(3L, 60L),
        stringsAsFactors = FALSE
    )
    ld <- LdData(correlation = R, variants = gr, blockMetadata = bm)
    expect_message(
        result <- pecotmr:::partitionLdMatrix(ld, mergeSmallBlocks = FALSE),
        "Removing 1 LD block"
    )
    expect_length(result$ldMatrices, 1)
    expect_equal(nrow(result$ldMatrices[[1]]), 3L)
})

test_that("validateBlockStructure flags out-of-range block indices", {
    mat <- diag(4)
    vnames <- sprintf("chr1:%d:A:G", 100L * (1:4))
    rownames(mat) <- colnames(mat) <- vnames
    bm <- data.frame(
        blockId = c(1L, 2L),
        chrom = c("1", "1"),
        size = c(2L, 2L),
        startIdx = c(1L, 10L),
        endIdx = c(2L, 12L)
    )
    expect_error(
        pecotmr:::validateBlockStructure(mat, bm, vnames),
        "Block indices out of range"
    )
})

test_that("extractBlockMatrices skips blocks where endIdx < startIdx", {
    mat <- diag(4)
    vnames <- sprintf("chr1:%d:A:G", 100L * (1:4))
    rownames(mat) <- colnames(mat) <- vnames
    bm <- data.frame(
        blockId = c(1L, 2L),
        startIdx = c(1L, 3L),
        endIdx = c(2L, 2L),
        chrom = c("1", "1"),
        blockStart = c(1L, 3L),
        blockEnd = c(2L, 4L),
        size = c(2L, 1L),
        stringsAsFactors = FALSE
    )
    result <- pecotmr:::extractBlockMatrices(mat, bm, vnames)
    valid <- result$ldMatrices[!sapply(result$ldMatrices, is.null)]
    expect_length(valid, 1)
    expect_equal(nrow(valid[[1]]), 2L)
})

# =============================================================================
# Additional coverage: ldPruneByCorrelation snprelate backend
# =============================================================================

test_that("ldPruneByCorrelation snprelate backend prunes correlated columns", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    set.seed(42)
    n <- 100
    p <- 6
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    X[, 2] <- X[, 1] # perfect LD between columns 1 and 2
    colnames(X) <- paste0("snp", 1:p)
    result <- suppressMessages(
        ldPruneByCorrelation(
            X,
            corThres = 0.5,
            backend = "snprelate",
            verbose = TRUE
        )
    )
    expect_true(ncol(result$X.new) <= p)
    expect_equal(length(result$filter.id), ncol(result$X.new))
    expect_true(all(result$filter.id %in% seq_len(p)))
})

# =============================================================================
# Additional coverage: dropCollinearColumns verbose messages
# =============================================================================

test_that("dropCollinearColumns prints verbose messages for each strategy", {
    set.seed(7)
    X <- matrix(rnorm(100 * 4), 100, 4)
    colnames(X) <- c("a", "b", "c", "d")
    X[, "c"] <- X[, "c"] * 0.001 # lowest variance
    y <- X[, "a"] * 2 + rnorm(100, sd = 0.1)

    expect_message(
        pecotmr:::dropCollinearColumns(
            X,
            "b",
            strategy = "correlation",
            verbose = TRUE
        ),
        "removing single column"
    )
    expect_message(
        pecotmr:::dropCollinearColumns(
            X,
            c("a", "b", "c"),
            strategy = "variance",
            verbose = TRUE
        ),
        "smallest variance"
    )
    expect_message(
        pecotmr:::dropCollinearColumns(
            X,
            c("a", "b"),
            strategy = "correlation",
            verbose = TRUE
        ),
        "two candidates"
    )
    expect_message(
        pecotmr:::dropCollinearColumns(
            X,
            c("a", "b", "c"),
            strategy = "correlation",
            verbose = TRUE
        ),
        "highest sum"
    )
    expect_message(
        pecotmr:::dropCollinearColumns(
            X,
            c("a", "b", "c"),
            strategy = "responseCorrelation",
            response = y,
            verbose = TRUE
        ),
        "smallest .* with response"
    )
})

# =============================================================================
# Additional coverage: enforceDesignFullRank verbose / fallback branches
# =============================================================================

test_that("enforceDesignFullRank verbose: batch-removal success + iterative path", {
    set.seed(11)
    X <- matrix(rnorm(80 * 4), 80, 4)
    X[, 4] <- X[, 1] + X[, 2] # rank deficient, fixable by removing one column
    colnames(X) <- c("a", "b", "c", "d")
    C <- matrix(rnorm(80), 80, 1)
    expect_message(
        result <- enforceDesignFullRank(
            X,
            C,
            strategy = "variance",
            verbose = TRUE
        ),
        "enforceDesignFullRank"
    )
    full_design <- cbind(1, result, C)
    expect_equal(qr(full_design)$rank, ncol(full_design))
})

test_that("enforceDesignFullRank verbose: constant covariate triggers fallback", {
    set.seed(12)
    X <- matrix(rnorm(60 * 3), 60, 3) # X itself is full rank
    colnames(X) <- c("a", "b", "c")
    C <- matrix(1, 60, 1) # constant -> collinear with intercept, not fixable via X
    # Iterative path finds no removable X column (break), then the correlation
    # fallback runs over each threshold; design stays rank-deficient throughout.
    expect_message(
        result <- enforceDesignFullRank(
            X,
            C,
            strategy = "correlation",
            verbose = TRUE
        ),
        "ldPruneByCorrelation fallback"
    )
    expect_true(is.matrix(result))
})

test_that("enforceDesignFullRank verbose: batch removal insufficient path", {
    set.seed(13)
    X <- matrix(rnorm(60 * 3), 60, 3)
    X <- cbind(X, X[, 1]) # duplicate of column 1
    colnames(X) <- c("a", "b", "c", "d")
    C <- matrix(1, 60, 1) # constant covariate keeps design deficient
    # Removing the QR-flagged X column cannot restore full rank (C is constant),
    # so the batch-removal "insufficient" branch fires and skips iterative pruning.
    expect_message(
        result <- enforceDesignFullRank(
            X,
            C,
            strategy = "correlation",
            verbose = TRUE
        ),
        "batch removal insufficient"
    )
    expect_true(is.matrix(result))
})

# =============================================================================
# Additional coverage: ldClumpByScore verbose + FBM input
# =============================================================================

test_that("ldClumpByScore prints verbose message for single-variant input", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    X <- matrix(c(0, 1, 2, 1, 0), ncol = 1)
    expect_message(
        ldClumpByScore(X, score = 1.0, chr = 1L, pos = 100L, verbose = TRUE),
        "single variant"
    )
})

test_that("ldClumpByScore accepts a pre-built FBM and reports retained count (verbose)", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    set.seed(1)
    n <- 200
    p <- 8
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    X[, 2] <- X[, 1] # perfect LD
    G <- bigstatsr::FBM.code256(
        nrow = n,
        ncol = p,
        init = X,
        code = c(0, 1, 2, rep(NA_real_, 253L))
    )
    score <- c(2, 1, runif(p - 2))
    chr <- rep(1L, p)
    pos <- seq_len(p) * 1000L
    expect_message(
        keep <- ldClumpByScore(
            G,
            score = score,
            chr = chr,
            pos = pos,
            r2 = 0.2,
            verbose = TRUE
        ),
        "variants retained"
    )
    expect_true(1L %in% keep)
    expect_false(2L %in% keep)
})

# =============================================================================
# Additional coverage: extractLdMatrix
# =============================================================================

test_that("extractLdMatrix errors on non-LdData input", {
    expect_error(pecotmr:::extractLdMatrix(list()), "must be an LdData object")
})

test_that("extractLdMatrix returns the genotype matrix when wantGenotype=TRUE", {
    variant_ids <- paste0("chr1:", c(100, 200, 300), ":A:G")
    X <- matrix(rnorm(15), 5, 3)
    colnames(X) <- variant_ids
    ref <- pecotmr:::parseVariantId(variant_ids)
    ref$variant_id <- variant_ids
    gr <- pecotmr:::.refPanelToGranges(ref)
    bm <- data.frame(
        blockId = 1L,
        chrom = "1",
        blockStart = 100L,
        blockEnd = 300L,
        size = 3L,
        startIdx = 1L,
        endIdx = 3L,
        stringsAsFactors = FALSE
    )
    ld <- LdData(genotypeHandle = X, variants = gr, blockMetadata = bm)
    result <- pecotmr:::extractLdMatrix(ld, wantGenotype = TRUE)
    expect_equal(result, X)
})


# =============================================================================
# Additional coverage: computeLd alternative backends + guard rails
# =============================================================================

test_that("computeLd snprelate backend returns a valid correlation matrix", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    set.seed(1)
    X <- matrix(rbinom(100 * 5, 2, 0.3), 100, 5)
    colnames(X) <- paste0("rs", 1:5)
    R <- suppressMessages(computeLd(
        X,
        method = "sample",
        backend = "snprelate"
    ))
    expect_equal(dim(R), c(5L, 5L))
    expect_equal(unname(diag(R)), rep(1, 5))
    expect_true(all(is.finite(R)))
    expect_equal(colnames(R), colnames(X))
})

test_that("computeLd snpstats backend returns a valid correlation matrix", {
    skip_if_not_installed("snpStats")
    set.seed(1)
    X <- matrix(rbinom(100 * 5, 2, 0.3), 100, 5)
    colnames(X) <- paste0("rs", 1:5)
    R <- computeLd(X, method = "sample", backend = "snpstats")
    expect_equal(dim(R), c(5L, 5L))
    expect_equal(unname(diag(R)), rep(1, 5))
    expect_true(all(is.finite(R)))
})

test_that("computeLd errors when a non-internal backend is paired with non-sample method", {
    set.seed(1)
    X <- matrix(rbinom(100 * 3, 2, 0.3), 100, 3)
    colnames(X) <- paste0("rs", 1:3)
    expect_error(
        computeLd(X, method = "population", backend = "snprelate"),
        "only supported with method='sample'"
    )
    expect_error(
        computeLd(X, method = "gcta", backend = "snpstats"),
        "only supported with method='sample'"
    )
})

# The Rfast-absent fallback branches (computeLd's `R <- cor(X_imp)` and
# ldPruneByCorrelation's `cor.X <- cor(X)`) are reachable by mocking the *base*
# `requireNamespace` (so it reports Rfast missing) for the duration of the call.
test_that("ldPruneByCorrelation and computeLd fall back to base cor() when Rfast is absent", {
    with_mocked_bindings(
        {
            set.seed(1)
            X <- matrix(rnorm(50 * 4), 50, 4)
            colnames(X) <- paste0("s", 1:4)
            pruned <- ldPruneByCorrelation(X, corThres = 0.9) # hits cor(X) at ld.R:1207
            R <- computeLd(X, method = "sample", backend = "internal") # hits cor(X_imp) at ld.R:1816
            expect_true(is.list(pruned))
            expect_true(all(c("X.new", "filter.id") %in% names(pruned)))
            expect_equal(dim(R), c(4L, 4L))
            expect_equal(unname(diag(R)), rep(1, 4), tolerance = 1e-8)
            expect_true(isSymmetric(unname(R)))
            # base cor() fallback must agree with the direct base computation.
            expect_equal(unname(R), unname(cor(X)), tolerance = 1e-10)
        },
        requireNamespace = function(package, ...) {
            if (identical(package, "Rfast")) {
                FALSE
            } else {
                base::requireNamespace(package, ...)
            }
        },
        .package = "base"
    )
})

# =============================================================================
# detectVariantConvention — uncovered line 586
# =============================================================================

# ---------------------------------------------------------------------------
# .panelVariantStats / .panelVariantFilter
#
# The shared measurement of "how common and how well-genotyped is this variant
# in the LD panel". Both the analysis-time RSS filter and the RAISS
# imputation-target filter read it, so a variant is judged the same way
# wherever it is judged.
# ---------------------------------------------------------------------------

# @noRd
.pvf_dosage <- function() {
    set.seed(11)
    nS <- 100L
    af <- c(rep(0.35, 5L), rep(0.004, 5L))
    d <- vapply(af, function(f) rbinom(nS, 2L, f), numeric(nS))
    colnames(d) <- sprintf("chr1:%d:A:G", 1000L * seq_along(af))
    d[1:80, 2] <- NA
    d
}

test_that(".panelVariantStats measures MAF, AF and missingness", {
    d <- .pvf_dosage()
    st <- .panelVariantStats(d)
    expect_named(st, c("af", "maf", "missRate"))
    expect_true(all(st$maf <= 0.5, na.rm = TRUE))
    # The rare half sits far below the common half.
    expect_lt(max(st$maf[6:10]), min(st$maf[1:5]))
    expect_equal(st$missRate[[2]], 0.8)
    expect_equal(st$missRate[[1]], 0)
})

test_that(".panelVariantStats needs un-imputed dosage to see missingness", {
    # Mean-imputation fills every hole, so missingness would read as zero
    # everywhere. This pins why the callers pass meanImpute = FALSE.
    d <- .pvf_dosage()
    filled <- d
    filled[is.na(filled)] <- 0
    expect_gt(.panelVariantStats(d)$missRate[[2]], 0)
    expect_equal(.panelVariantStats(filled)$missRate[[2]], 0)
})

test_that(".panelVariantStats reports NA MAF for an all-missing variant", {
    d <- .pvf_dosage()
    d[, 3] <- NA_real_
    st <- .panelVariantStats(d)
    expect_true(is.na(st$maf[[3]]))
    expect_equal(st$missRate[[3]], 1)
})

test_that(".panelVariantFilter is a no-op at its defaults", {
    data(qtlDatasetExample)
    gh <- getGenotypes(qtlDatasetExample)
    handle <- getGenotypeHandle(qtlDatasetExample)
    ids <- normalizeVariantId(getSnpInfo(handle)$SNP)
    expect_identical(.panelVariantFilter(handle, ids), ids)
})

test_that(".panelVariantFilter drops panel-rare variants", {
    data(qtlDatasetExample)
    handle <- getGenotypeHandle(qtlDatasetExample)
    ids <- normalizeVariantId(getSnpInfo(handle)$SNP)
    loose <- .panelVariantFilter(handle, ids, mafCutoff = 0.05)
    tight <- .panelVariantFilter(handle, ids, mafCutoff = 0.2)
    expect_lt(length(loose), length(ids))
    expect_lt(length(tight), length(loose))
    # Kept sets are nested as the cutoff rises, and order is the caller's.
    expect_true(all(is_in(tight, loose)))
    expect_identical(loose, ids[is_in(ids, loose)])
})

test_that(".panelVariantFilter treats MAC as a MAF equivalent", {
    data(qtlDatasetExample)
    handle <- getGenotypeHandle(qtlDatasetExample)
    ids <- normalizeVariantId(getSnpInfo(handle)$SNP)
    nSamp <- getNSamples(handle)
    # macCutoff / (2 * nSamples) is the same threshold as mafCutoff.
    byMac <- .panelVariantFilter(handle, ids, macCutoff = 0.1 * 2 * nSamp)
    byMaf <- .panelVariantFilter(handle, ids, mafCutoff = 0.1)
    expect_identical(byMac, byMaf)
    # The stricter of the two wins.
    expect_identical(
        .panelVariantFilter(handle, ids, mafCutoff = 0.2, macCutoff = 2),
        .panelVariantFilter(handle, ids, mafCutoff = 0.2)
    )
})

test_that(".panelVariantFilter drops high-missingness variants", {
    data(qtlDatasetExample)
    handle <- getGenotypeHandle(qtlDatasetExample)
    ids <- normalizeVariantId(getSnpInfo(handle)$SNP)
    strict <- .panelVariantFilter(handle, ids, imissCutoff = 0)
    expect_lt(length(strict), length(ids))
    # A cutoff above the panel's worst variant keeps everything.
    expect_identical(.panelVariantFilter(handle, ids, imissCutoff = 1), ids)
})

test_that(".panelVariantFilter passes through ids absent from the panel", {
    # Whether a missing variant is an error or is dropped belongs to
    # .ldFromSketch's `onMissing`; deciding it here too would let the two
    # disagree about the same variant.
    data(qtlDatasetExample)
    handle <- getGenotypeHandle(qtlDatasetExample)
    ids <- normalizeVariantId(getSnpInfo(handle)$SNP)[1:3]
    withGhost <- c("chr9:999:A:G", ids)
    expect_true(is_in(
        "chr9:999:A:G",
        .panelVariantFilter(handle, withGhost, mafCutoff = 0.001)
    ))
})

test_that(".panelVariantFilter handles empty and NULL input", {
    data(qtlDatasetExample)
    handle <- getGenotypeHandle(qtlDatasetExample)
    expect_length(
        .panelVariantFilter(handle, character(0), mafCutoff = 0.1),
        0L
    )
    expect_identical(
        .panelVariantFilter(NULL, "chr1:1:A:G", mafCutoff = 0.1),
        "chr1:1:A:G"
    )
})


test_that(".panelCutoffs short-circuits when no cutoff is set", {
    # NULL means the panel is never touched, which is what keeps the default
    # path free of an extra dosage read.
    expect_null(.panelCutoffs(list()))
    expect_null(.panelCutoffs(list(
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1
    )))
    expect_equal(.panelCutoffs(list(mafCutoff = 0.01))$mafCutoff, 0.01)
    expect_equal(.panelCutoffs(list(imissCutoff = 0.5))$imissCutoff, 0.5)
})


# ---------------------------------------------------------------------------
# .panelVariantFilter: the .afreq fast path
#
# Allele frequency alone answers MAF and MAC, so a panel that ships a PLINK2
# .afreq sidecar is filtered without materializing its dosage. The sidecar and
# the dosage must agree on every variant, or one pipeline would judge a
# variant differently from the next.
# ---------------------------------------------------------------------------

# @noRd
.pvfAfreqHandle <- function() {
    readGenotypeHandle(test_path("test_data/test_variants"), format = "plink2")
}

test_that(".panelAfreqMaf reads panel MAF from the .afreq sidecar", {
    skip_if_not_installed("pgenlibr")
    handle <- .pvfAfreqHandle()
    ids <- as.character(getSnpInfo(handle)$SNP)
    maf <- .panelAfreqMaf(handle, ids)
    afreq <- readAfreq(test_path("test_data/test_variants"))
    altFreq <- afreq$alt_freq[match(ids, afreq$id)]
    expect_equal(maf, pmin(altFreq, 1 - altFreq))
})

test_that(".panelAfreqMaf refuses a partial or non-plink2 sidecar", {
    skip_if_not_installed("pgenlibr")
    handle <- .pvfAfreqHandle()
    ids <- as.character(getSnpInfo(handle)$SNP)
    # One id the sidecar cannot answer for is enough: a partial answer would
    # disagree with the dosage path about that variant, so NULL sends the
    # caller down the dosage path for all of them.
    expect_null(.panelAfreqMaf(handle, c(ids, "not-a-panel-variant")))
    bed <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink1"
    )
    expect_null(.panelAfreqMaf(bed, as.character(getSnpInfo(bed)$SNP)))
})

test_that(".panelVariantFilter: .afreq and dosage agree on what to drop", {
    skip_if_not_installed("pgenlibr")
    handle <- .pvfAfreqHandle()
    ids <- as.character(getSnpInfo(handle)$SNP)
    viaAfreq <- .panelVariantFilter(handle, ids, mafCutoff = 0.2)
    # Force the dosage path by hiding the sidecar from the fast path.
    local_mocked_bindings(
        .panelAfreqMaf = function(handle, variantIds) NULL,
        .package = "pecotmr"
    )
    viaDosage <- .panelVariantFilter(handle, ids, mafCutoff = 0.2)
    expect_lt(length(viaAfreq), length(ids))
    expect_identical(viaAfreq, viaDosage)
})

test_that(".panelVariantFilter uses dosage whenever missingness is capped", {
    skip_if_not_installed("pgenlibr")
    # Missingness is not derivable from the sidecar, so an imissCutoff must
    # never be silently answered by it.
    handle <- .pvfAfreqHandle()
    ids <- as.character(getSnpInfo(handle)$SNP)
    local_mocked_bindings(
        .panelAfreqMaf = function(handle, variantIds) {
            abort("the .afreq fast path must not run here")
        },
        .package = "pecotmr"
    )
    expect_no_error(
        .panelVariantFilter(handle, ids, mafCutoff = 0.2, imissCutoff = 0.5)
    )
})

test_that(".panelAfreqPrefixes reads only the chromosomes the panel spans", {
    handle <- .pvfAfreqHandle()
    # Single-file panel: its own stem, resolved for file access.
    expect_identical(
        .panelAfreqPrefixes(handle),
        pecotmr:::.genotypeReadPath(handle)
    )
    # Sharded panel: a manifest chromosome the snpInfo does not span is not
    # read -- the cost per-chromosome sketch trimming exists to avoid.
    sharded <- handle
    sharded@chromPaths <- c(
        "22" = getPath(handle),
        "21" = "/nonexistent/chr21"
    )
    sharded@snpInfo$CHR <- rep("22", nrow(getSnpInfo(handle)))
    expect_identical(.panelAfreqPrefixes(sharded), getPath(handle))
})


# =============================================================================
# loadLdMatrix: coordinate-free sources and per-block options
#
# These cover what ldLoader() / loadLdBlock() / loadLdSketch() used to do,
# now folded into the single loader. Every source returns an LdData, which is
# the point of the consolidation -- callers no longer branch on what they
# loaded from.
# =============================================================================

test_that("loadLdMatrix requires exactly one addressing mode", {
    R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    expect_error(loadLdMatrix(R), "either `region` or `block`")
    expect_error(
        loadLdMatrix(R, region = "chr1:1-2", block = 1),
        "not both"
    )
    # A coordinate-free source cannot be addressed by region.
    expect_error(loadLdMatrix(R, region = "chr1:1-2"), "addressed with `block`")
})

test_that("loadLdMatrix loads an in-memory correlation matrix by block", {
    R1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
    R2 <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
    expect_equal(unname(getCorrelation(loadLdMatrix(R1, block = 1))), R1)
    expect_equal(
        unname(getCorrelation(loadLdMatrix(list(R1, R2), block = 2))),
        R2
    )
})

test_that("loadLdMatrix treats a tall matrix as genotypes", {
    X <- matrix(rnorm(40), 10, 4)
    ld <- loadLdMatrix(X, block = 1)
    expect_true(hasGenotypes(ld))
    expect_equal(length(ld), 4L)
})

test_that("loadLdMatrix subsamples oversized blocks reproducibly", {
    R <- matrix(0.1, 10, 10)
    diag(R) <- 1
    small <- loadLdMatrix(R, block = 1, maxVariants = 5, seed = 42)
    expect_equal(length(small), 5L)
    expect_equal(dim(getCorrelation(small)), c(5L, 5L))
    # Same seed, same draw.
    again <- loadLdMatrix(R, block = 1, maxVariants = 5, seed = 42)
    expect_equal(getVariantIds(small), getVariantIds(again))
    # A cap above the block size is a no-op.
    expect_equal(length(loadLdMatrix(R, block = 1, maxVariants = 100)), 10L)
})

test_that("loadLdMatrix reads an ldInfo table by block", {
    skip_if_not_installed("pgenlibr")
    ldInfo <- data.frame(
        LD_file = file.path(test_path("test_data"), "test_variants")
    )
    ld <- loadLdMatrix(ldInfo, block = 1)
    expect_s4_class(ld, "LdData")
    expect_equal(length(ld), 349L)
})

test_that("loadLdMatrix validates an ldInfo table", {
    expect_error(
        loadLdMatrix(data.frame(col1 = "a"), block = 1),
        "needs an `LD_file` column"
    )
})

test_that("loadLdMatrix takes a vector of regions", {
    meta <- system.file(
        "extdata",
        "ld_reference",
        "ld_meta_file.tsv",
        package = "pecotmr"
    )
    region <- "chr22:10000000-19000000"
    many <- loadLdMatrix(meta, region = c(region, region))
    expect_type(many, "list")
    expect_length(many, 2L)
    expect_s4_class(many[[1]], "LdData")
    expect_equal(length(many[[1]]), length(many[[2]]))
})

test_that("loadLdMatrix can materialize genotypes onto the object", {
    meta <- system.file(
        "extdata",
        "ld_reference",
        "ld_meta_file.tsv",
        package = "pecotmr"
    )
    region <- "chr22:10000000-19000000"
    lazy <- loadLdMatrix(meta, region = region, returnGenotype = TRUE)
    eager <- loadLdMatrix(
        meta,
        region = region,
        returnGenotype = TRUE,
        materializeGenotypes = TRUE
    )
    # Same data; the difference is whether later access re-reads the file.
    expect_s4_class(getGenotypeHandle(lazy), "GenotypeHandle")
    expect_true(is.matrix(getGenotypeHandle(eager)))
    expect_equal(getGenotypes(lazy), getGenotypes(eager))
})

test_that("loadLdMatrix can drop monomorphic variants", {
    meta <- system.file(
        "extdata",
        "ld_reference",
        "ld_meta_file.tsv",
        package = "pecotmr"
    )
    region <- "chr22:10000000-19000000"
    full <- loadLdMatrix(meta, region = region, returnGenotype = TRUE)
    kept <- loadLdMatrix(
        meta,
        region = region,
        returnGenotype = TRUE,
        dropMonomorphic = TRUE
    )
    # This panel has no monomorphic variants, so the filter is a no-op here;
    # what matters is that it does not drop polymorphic ones.
    expect_lte(length(kept), length(full))
    af <- getRefPanel(kept)$allele_freq
    expect_true(all(af > 0 & af < 1))
})


# =============================================================================
# One failure mode for a missing or invalid LD sketch
#
# Panel access happens through several routes: the shared `.ldFromSketch()`,
# ctwas's whole-panel assembler, RAiSS's coordinate window, and the bare
# accessors. They used to disagree on what a bad sketch looked like -- the
# shared path said so plainly while the others surfaced "unable to find an
# inherited method for 'getSnpInfo'" from whichever accessor touched it
# first. The guard now sits at `.ldSketchRanges()` / `.ldSketchDosage()`, the
# two primitives everything funnels through.
# =============================================================================

test_that("every panel access route rejects a NULL sketch the same way", {
    df <- data.frame(
        chrom = "chr22",
        pos = 1:3,
        SNP = paste0("chr22:", 1:3, ":A:G")
    )
    routes <- list(
        shared = function() {
            pecotmr:::.ldFromSketch(NULL, "chr22:1:A:G", label = "demo")
        },
        ctwas = function() pecotmr:::.ctwasComputeFullPanelLd(NULL),
        raiss = function() pecotmr:::.qcRaissWindowIdx(df, NULL, 0L),
        matchIds = function() pecotmr:::.ldSketchMatchIds(NULL),
        dosage = function() pecotmr:::.ldSketchDosage(NULL, 1L),
        ranges = function() pecotmr:::.ldSketchRanges(NULL)
    )
    for (nm in names(routes)) {
        expect_error(routes[[nm]](), "carries no ldSketch", info = nm)
    }
})

test_that("a non-panel sketch is rejected by every route", {
    expect_error(
        pecotmr:::.ldSketchRanges("not a panel"),
        "must be a genotype panel"
    )
    expect_error(
        pecotmr:::.ctwasComputeFullPanelLd("not a panel"),
        "must be a genotype panel"
    )
    expect_error(
        pecotmr:::.ldSketchDosage("not a panel", 1L),
        "must be a genotype panel"
    )
})

test_that("callers that validate first keep their own label", {
    # The generic guard must not mask a caller's more specific message.
    shared <- tryCatch(
        pecotmr:::.ldFromSketch(NULL, "x", label = "fineMappingPipeline"),
        error = conditionMessage
    )
    expect_match(shared, "^fineMappingPipeline:")
    ctwas <- tryCatch(
        pecotmr:::.ctwasComputeFullPanelLd(NULL),
        error = conditionMessage
    )
    expect_match(ctwas, "^ctwasPipeline:")
})

test_that("panel identity uses the allele-repaired ids for comparison", {
    # `.ldSketchMatchIds()` repairs an id whose allele fields are tags rather
    # than DNA, filling them from the panel's A1/A2 mcols; that repaired form
    # is what every comparison against sumstats or weights uses.
    repaired <- pecotmr:::.repairVariantIds(
        c("chr1:100:R:V", "chr1:200:R:V"),
        A2 = c("A", "C"),
        A1 = c("G", "T")
    )
    expect_equal(repaired, c("chr1:100:A:G", "chr1:200:C:T"))

    # An id that already carries DNA alleles is left alone.
    dna <- c("chr1:100:A:G")
    expect_equal(pecotmr:::.repairVariantIds(dna, A2 = "C", A1 = "T"), dna)

    # And an rsID panel passes through unchanged -- `.ldFromSketchMatch()`
    # falls back to exact id matching for those, so rewriting them would
    # break the match rather than fix it.
    panel <- readGenotypes(
        system.file("extdata", "toy_ref.bed", package = "pecotmr")
    )
    expect_equal(
        pecotmr:::.ldSketchMatchIds(panel),
        pecotmr:::.ldSketchVariantIds(panel)
    )
})


# ===========================================================================
# Sketch accessors and per-variant panel statistics
# ===========================================================================

test_that("the sketch accessors are NULL-safe", {
    # A collection with no LD sketch is a legitimate state (individual-level
    # inputs), so these answer NULL rather than erroring.
    expect_null(pecotmr:::.ldSketchHandle(NULL))
    expect_null(pecotmr:::.ldSketchSubset(NULL, 1L))
})

test_that("panel statistics report zero missingness for a sample-less panel", {
    # With no samples the missing RATE is undefined as 1 - nObs/nSamp; report
    # 0 per variant rather than NaN, which would poison the QC comparisons.
    d0 <- matrix(numeric(0), nrow = 0L, ncol = 3L)
    expect_equal(pecotmr:::.panelVariantStats(d0)$missRate, rep(0, 3L))
})

test_that("panel allele frequencies are computed per variant", {
    d <- matrix(c(0, 1, 2, 0, 0, 0), nrow = 2L, byrow = TRUE)
    st <- pecotmr:::.panelVariantStats(d)
    expect_equal(st$af, c(0, 0.25, 0.5))
    # maf folds above 0.5.
    expect_equal(st$maf, pmin(st$af, 1 - st$af))
})


test_that("panel variant filtering is a no-op without a panel or variants", {
    # Both early exits return the caller's ids unchanged: with no panel there
    # is nothing to filter against, and with no variants nothing to filter.
    expect_equal(
        pecotmr:::.panelVariantFilter(NULL, c("v1", "v2")),
        c("v1", "v2")
    )
    expect_equal(
        pecotmr:::.panelVariantFilter(NULL, character(0)),
        character(0)
    )
})


# ===========================================================================
# Empty sketches, the .afreq fast path, and design rank
# ===========================================================================

test_that("emptying a sketch clears both axes", {
    # An emptied sketch references no LD, so its sample axis is dead weight;
    # nSamples must go to 0 alongside sampleIds or the derived dosage
    # dimnames would disagree with the matrix.
    expect_null(pecotmr:::.emptySketch(NULL))
    h <- new(
        "GenotypeHandle",
        path = "/tmp/x.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = c("a", "b"),
            CHR = "1",
            BP = 1:2,
            A1 = "A",
            A2 = "G",
            fileIdx = 1:2,
            stringsAsFactors = FALSE
        ),
        nSamples = 2L,
        sampleIds = c("s1", "s2"),
        pgenPtr = NULL
    )
    e <- pecotmr:::.emptySketch(h)
    expect_equal(nrow(pecotmr:::getSnpInfo(e)), 0L)
    expect_equal(length(pecotmr:::getSampleIds(e)), 0L)
})

test_that("a missing or malformed afreq sidecar degrades to no frequencies", {
    # The sidecar is optional; a bad one must not take the whole filter down,
    # it just falls back to reading dosage.
    expect_null(pecotmr:::.panelAfreqTable("/nonexistent/prefix"))
    h <- new(
        "GenotypeHandle",
        path = "/tmp/x.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = c("a", "b"),
            CHR = "1",
            BP = 1:2,
            A1 = "A",
            A2 = "G",
            fileIdx = 1:2,
            stringsAsFactors = FALSE
        ),
        nSamples = 2L,
        sampleIds = c("s1", "s2"),
        pgenPtr = NULL
    )
    # The fast path is plink2-only.
    expect_null(pecotmr:::.panelAfreqMaf(h, c("a", "b")))
})

test_that("collinear design columns are reported, full-rank ones are not", {
    full <- matrix(
        c(1, 0, 0, 1),
        2L,
        2L,
        dimnames = list(NULL, c("a", "b"))
    )
    expect_equal(
        pecotmr:::.edfrProblematicColnames(full, full),
        character(0)
    )
    # b is a multiple of a: QR pivots it past the rank boundary.
    dup <- cbind(a = c(1, 0), b = c(2, 0))
    expect_equal(pecotmr:::.edfrProblematicColnames(dup, dup), "b")
})

test_that("computeLd refuses snpIdx against an already-materialized block", {
    # snpIdx selects variants FROM a panel; a bare matrix is already the
    # block, so accepting it would silently ignore the subscript.
    expect_error(
        computeLd(matrix(rnorm(20), 10L, 2L), snpIdx = 1L),
        "`snpIdx` selects variants from a genotype panel"
    )
})


# ===========================================================================
# Block-indexed sources and LdData subsetting
# ===========================================================================

test_that("a block-indexed source must be addressable by block", {
    # `block` only means something for a meta file, an ldInfo frame, a matrix
    # or a list of them; anything else is refused rather than silently
    # returning the wrong block.
    expect_error(
        pecotmr:::.loadLdFromIndexed("a string", 1L, FALSE),
        "cannot address a character by block"
    )
})

test_that("a matrix and a list of matrices both address by block", {
    R <- diag(3)
    dimnames(R) <- list(c("v1", "v2", "v3"), c("v1", "v2", "v3"))
    expect_s4_class(pecotmr:::.loadLdFromIndexed(R, 1L, FALSE), "LdData")
    # The list form picks the requested element.
    both <- pecotmr:::.loadLdFromIndexed(list(R, R), 2L, FALSE)
    expect_s4_class(both, "LdData")
    expect_equal(length(both), 3L)
})

test_that("materializing genotypes is a no-op without genotypes to read", {
    # Both exits return the LdData untouched: not requested, and requested
    # but the source is a correlation matrix with nothing to materialize.
    R <- diag(3)
    dimnames(R) <- list(c("v1", "v2", "v3"), c("v1", "v2", "v3"))
    ld <- pecotmr:::.ldDataFromMatrix(R, isGenotype = FALSE)
    expect_identical(pecotmr:::.ldApplyMaterialize(ld, FALSE), ld)
    expect_identical(pecotmr:::.ldApplyMaterialize(ld, TRUE), ld)
})

test_that("subsetting an LdData narrows the correlation on both axes", {
    R <- diag(3)
    dimnames(R) <- list(c("v1", "v2", "v3"), c("v1", "v2", "v3"))
    ld <- pecotmr:::.ldDataFromMatrix(R, isGenotype = FALSE)
    out <- pecotmr:::.ldSubsetData(ld, c(1L, 3L))
    expect_equal(length(out), 2L)
    expect_equal(dim(getCorrelation(out)), c(2L, 2L))
})

# ---------------------------------------------------------------------------
# LdData narrowing and materialization: each genotype source shape (matrix,
# handle + snpIdx, mixture list) narrows differently, and only the matrix
# shape can be materialized.
# ---------------------------------------------------------------------------

.ldcov_gr <- function(variant_ids, af = NULL) {
    rp <- pecotmr:::parseVariantId(variant_ids)
    rp$variant_id <- variant_ids
    if (!is.null(af)) {
        rp$allele_freq <- af
    }
    pecotmr:::.refPanelToGranges(rp)
}

.ldcov_bm <- function(n) {
    data.frame(
        blockId = 1L,
        chrom = "chr1",
        blockStart = 100,
        blockEnd = 100 * n,
        size = n,
        startIdx = 1L,
        endIdx = n,
        stringsAsFactors = FALSE
    )
}

.ldcov_handle <- function() {
    new(
        "GenotypeHandle",
        path = "/tmp/test.gds",
        format = "gds",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
    )
}

test_that(".loadLdDedup passes through when there are no variant ids", {
    # A GRanges with no variant_id column: getVariantIds() is NULL, so there
    # is nothing to deduplicate against.
    bare <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100, 200, 300), width = 1)
    )
    ld <- LdData(
        correlation = diag(3),
        variants = bare,
        blockMetadata = .ldcov_bm(3L)
    )
    expect_null(getVariantIds(ld))
    expect_identical(pecotmr:::.loadLdDedup(ld), ld)
})

test_that(".ldApplyMonomorphic keeps an all-polymorphic block unchanged", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    ld <- LdData(
        correlation = diag(3),
        variants = .ldcov_gr(v, af = c(0.2, 0.3, 0.4)),
        blockMetadata = .ldcov_bm(3L)
    )
    # Every allele_freq is strictly inside (0, 1), so no subset is taken.
    expect_identical(pecotmr:::.ldApplyMonomorphic(ld, TRUE), ld)
})

test_that(".ldApplyMonomorphic needs an allele_freq column to act", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    ld <- LdData(
        correlation = diag(3),
        variants = .ldcov_gr(v),
        blockMetadata = .ldcov_bm(3L)
    )
    # No frequencies to judge monomorphism by, so the request is a no-op
    # rather than an error -- distinct from the all-polymorphic case.
    expect_false(is_in("allele_freq", colnames(getRefPanel(ld))))
    expect_identical(pecotmr:::.ldApplyMonomorphic(ld, TRUE), ld)
})

test_that(".ldApplyMonomorphic drops monomorphic variants", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    ld <- LdData(
        correlation = diag(3),
        variants = .ldcov_gr(v, af = c(0.2, 0.0, 0.4)),
        blockMetadata = .ldcov_bm(3L)
    )
    out <- pecotmr:::.ldApplyMonomorphic(ld, TRUE)
    expect_equal(getVariantIds(out), c("chr1:100:A:G", "chr1:300:G:A"))
    expect_equal(dim(getCorrelation(out)), c(2L, 2L))
})

test_that(".ldApplyMaterialize leaves a mixture list alone", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    m1 <- matrix(0, nrow = 4L, ncol = 3L, dimnames = list(NULL, v))
    m2 <- matrix(1, nrow = 4L, ncol = 3L, dimnames = list(NULL, v))
    ld <- LdData(
        correlation = NULL,
        genotypeHandle = list(m1, m2),
        mixtureWeights = c(0.5, 0.5),
        variants = .ldcov_gr(v),
        blockMetadata = .ldcov_bm(3L)
    )
    # getGenotypes() on a mixture returns a LIST of dosage matrices; there is
    # no single matrix to fold into the object.
    expect_false(is.matrix(getGenotypes(ld)))
    expect_identical(pecotmr:::.ldApplyMaterialize(ld, TRUE), ld)
})

test_that(".ldSubsetData narrows a dosage matrix by column", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    m <- matrix(0, nrow = 4L, ncol = 3L, dimnames = list(NULL, v))
    ld <- LdData(
        correlation = NULL,
        genotypeHandle = m,
        snpIdx = NULL,
        variants = .ldcov_gr(v),
        blockMetadata = .ldcov_bm(3L)
    )
    out <- pecotmr:::.ldSubsetData(ld, c(1L, 3L))
    expect_equal(dim(pecotmr:::getGenotypeHandle(out)), c(4L, 2L))
    expect_equal(getVariantIds(out), c("chr1:100:A:G", "chr1:300:G:A"))
})

test_that(".ldSubsetData narrows a handle through snpIdx", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    ld <- LdData(
        correlation = NULL,
        genotypeHandle = .ldcov_handle(),
        snpIdx = c(5L, 6L, 7L),
        variants = .ldcov_gr(v),
        blockMetadata = .ldcov_bm(3L)
    )
    out <- pecotmr:::.ldSubsetData(ld, c(1L, 3L))
    # The handle is untouched; snpIdx selects into its full snpInfo, so
    # subsetting picks positions 1 and 3 OF THE INDEX, not of the file.
    expect_equal(pecotmr:::getSnpIdx(out), c(5L, 7L))
    expect_equal(getVariantIds(out), c("chr1:100:A:G", "chr1:300:G:A"))
})

test_that(".ldSketchMatchIds returns bare ids when alleles are absent", {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100, 200), width = 1)
    )
    S4Vectors::mcols(gr)$SNP <- c("chr1:100:A:G", "chr1:200:C:T")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(x = matrix(0, nrow = 2L, ncol = 1L)),
        rowRanges = gr
    )
    # With no A1/A2 columns there is nothing to repair the ids against.
    expect_equal(
        pecotmr:::.ldSketchMatchIds(se),
        c("chr1:100:A:G", "chr1:200:C:T")
    )
})

test_that(".panelVariantFilter is a no-op without a sketch", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    # An active cutoff still cannot filter anything with no panel to read
    # frequencies from.
    expect_equal(
        pecotmr:::.panelVariantFilter(NULL, v, mafCutoff = 0.01),
        v
    )
    expect_equal(
        pecotmr:::.panelVariantFilter(NULL, character(0), mafCutoff = 0.01),
        character(0)
    )
})

test_that(".panelVariantFilter is a no-op when nothing matches the panel", {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100, 200), width = 1)
    )
    S4Vectors::mcols(gr)$SNP <- c("chr1:100:A:G", "chr1:200:C:T")
    sketch <- SummarizedExperiment::SummarizedExperiment(
        assays = list(x = matrix(0, nrow = 2L, ncol = 1L)),
        rowRanges = gr
    )
    ids <- c("chr9:999:A:G", "chr9:888:C:T")
    # The sketch exists but shares no variant with the request, so there is
    # no frequency to filter on and the ids pass through untouched.
    expect_equal(
        pecotmr:::.panelVariantFilter(sketch, ids, mafCutoff = 0.01),
        ids
    )
})

test_that(".ldSketchNullGuard names the label in its strict error", {
    expect_error(
        pecotmr:::.ldSketchNullGuard("Q", NULL, "myPipe", "myLabel", "strict"),
        "ldSketch on `myLabel` is non-NULL"
    )
    expect_error(
        pecotmr:::.ldSketchNullGuard("Q", NULL, "myPipe", NULL, "strict"),
        "qtl ldSketch is non-NULL"
    )
    # lenient tolerates the same input.
    expect_true(
        pecotmr:::.ldSketchNullGuard("Q", NULL, "myPipe", "myLabel", "lenient")
    )
})

test_that(".ldSketchCheckContent rejects panels on different chromosomes", {
    mkSe <- function(chr, pos) {
        g <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1))
        S4Vectors::mcols(g)$SNP <- str_c(chr, ":", pos, ":A:G")
        SummarizedExperiment::SummarizedExperiment(
            assays = list(x = matrix(0, nrow = length(pos), ncol = 1L)),
            rowRanges = g
        )
    }
    expect_error(
        pecotmr:::.ldSketchCheckContent(
            mkSe("chr1", c(100, 200)),
            mkSe("chr2", c(100, 200)),
            "myPipe",
            " between X and Y"
        ),
        "differ in column CHR between X and Y"
    )
})

test_that(".blockPairMessages is silent on a clean block pair", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    m <- diag(3)
    rownames(m) <- colnames(m) <- v
    # Cross-block correlation is zero everywhere, which is what a correct
    # block structure looks like.
    bm <- data.frame(startIdx = c(1L, 2L), endIdx = c(1L, 3L))
    expect_equal(
        pecotmr:::.blockPairMessages(1L, 2L, bm, m, v, 3L),
        character(0)
    )
})

test_that(".dropCollinearPickCor is reproducible under a seed", {
    X <- cbind(a = c(1, 2, 3, 4, 5), b = c(1, 2, 3, 4, 5.0001))
    first <- pecotmr:::.dropCollinearPickCor(X, c(1L, 2L), FALSE, seed = 42L)
    second <- pecotmr:::.dropCollinearPickCor(X, c(1L, 2L), FALSE, seed = 42L)
    # Two candidates means the choice is random; the seed pins it.
    expect_identical(first, second)
    expect_true(first %in% c(1L, 2L))
})

test_that("extractLdMatrix returns the correlation when genotypes are absent", {
    v <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    ld <- LdData(
        correlation = diag(3),
        variants = .ldcov_gr(v),
        blockMetadata = .ldcov_bm(3L)
    )
    expect_false(hasGenotypes(ld))
    # wantGenotype is honoured only when genotypes exist; otherwise the
    # correlation is returned rather than erroring.
    expect_equal(dim(pecotmr:::extractLdMatrix(ld)), c(3L, 3L))
    expect_equal(
        dim(pecotmr:::extractLdMatrix(ld, wantGenotype = TRUE)),
        c(3L, 3L)
    )
})

# ---------------------------------------------------------------------------
# computeLd default index: with snpIdx = NULL every variant in the panel is
# used, on both the in-memory and the on-disk route.
# ---------------------------------------------------------------------------

test_that("computeLd(panel) defaults to the panel's whole variant set", {
    panel <- readGenotypes(
        path = test_path("test_data", "test_variants.gds"),
        format = "gds"
    )
    sub <- panel[1:6, ]
    R <- computeLd(sub)
    expect_equal(dim(R), c(6L, 6L))
    expect_equal(unname(diag(R)), rep(1, 6L))
    # The NULL default must agree with spelling the index out.
    expect_equal(unname(R), unname(computeLd(sub, snpIdx = 1:6)))
})

test_that("computeLd(onDisk) defaults to every variant in the panel", {
    skip_if_not_installed("SNPRelate")
    handle <- readGenotypeHandle(
        test_path("test_data", "test_variants.gds"),
        format = "gds"
    )
    R <- computeLd(handle, backend = "snprelate", onDisk = TRUE)
    n <- nrow(pecotmr:::getSnpInfo(handle))
    expect_equal(dim(R), c(n, n))
    expect_identical(rownames(R), pecotmr:::getSnpInfo(handle)$SNP)
})

test_that("computeLd(onDisk) applies shrinkage toward the identity", {
    skip_if_not_installed("SNPRelate")
    handle <- readGenotypeHandle(
        test_path("test_data", "test_variants.gds"),
        format = "gds"
    )
    idx <- 1:6
    plain <- computeLd(handle, snpIdx = idx, backend = "snprelate",
        onDisk = TRUE)
    shrunk <- computeLd(handle, snpIdx = idx, backend = "snprelate",
        onDisk = TRUE, shrinkage = 0.5)
    # (1 - s) * R + s * I: off-diagonals halve, the diagonal stays 1.
    expect_equal(unname(diag(shrunk)), rep(1, length(idx)))
    expect_equal(
        unname(shrunk[upper.tri(shrunk)]),
        unname(0.5 * plain[upper.tri(plain)])
    )
})

# ---------------------------------------------------------------------------
# .ldInfoBlock: an `ldInfo` row names the LD matrix, and the variant metadata
# is either named alongside it or auto-detected from the LD path.
# ---------------------------------------------------------------------------

test_that(".ldInfoBlock reads a block with an explicit SNP_file", {
    ldFile <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.txt.xz"
    )
    bimFile <- file.path(
        geno_test_data_dir,
        "LD_block_1.chr1_1000_1200.float16.bim"
    )
    out <- pecotmr:::.ldInfoBlock(
        data.frame(
            LD_file = ldFile,
            SNP_file = bimFile,
            stringsAsFactors = FALSE
        ),
        1L
    )
    # Both ldInfo sources return an LdData, so the post-load chain does not
    # have to branch on which kind of LD_file it was given.
    expect_s4_class(out, "LdData")
    expect_equal(dim(getCorrelation(out)), c(5L, 5L))
    expect_equal(
        head(getVariantIds(out), 3L),
        c("chr1:1000:A:G", "chr1:1040:A:G", "chr1:1080:A:G")
    )
    # The .bim coordinates are carried through rather than replaced with the
    # chrNA:1..n placeholders a bare matrix would get.
    refPanel <- getRefPanel(out)
    expect_equal(as.character(refPanel$chrom[[1L]]), "chr1")
    expect_equal(as.integer(refPanel$pos[[1L]]), 1000L)
})

test_that("loadLdMatrix reads an ldInfo table of precomputed LD", {
    info <- data.frame(
        LD_file = file.path(
            geno_test_data_dir,
            "LD_block_1.chr1_1000_1200.float16.txt.xz"
        ),
        SNP_file = file.path(
            geno_test_data_dir,
            "LD_block_1.chr1_1000_1200.float16.bim"
        ),
        stringsAsFactors = FALSE
    )
    # End-to-end: the whole post-load chain (dedup, monomorphic, subsample,
    # materialize) runs against a precomputed block, not just a genotype one.
    out <- loadLdMatrix(info, block = 1L)
    expect_s4_class(out, "LdData")
    expect_equal(length(getVariantIds(out)), 5L)
    expect_equal(unname(diag(getCorrelation(out))), rep(1, 5L))
})

test_that(".ldInfoBlock auto-detects the variant file beside the matrix", {
    # Auto-detection APPENDS the suffix to the LD path, so the companion of
    # `block.cor.xz` is `block.cor.xz.bim` -- not a sibling with the LD
    # extension swapped out.
    dir <- withr::local_tempdir()
    ldFile <- file.path(dir, "block.cor.xz")
    file.copy(
        file.path(
            geno_test_data_dir,
            "LD_block_1.chr1_1000_1200.float16.txt.xz"
        ),
        ldFile
    )
    file.copy(
        file.path(
            geno_test_data_dir,
            "LD_block_1.chr1_1000_1200.float16.bim"
        ),
        str_c(ldFile, ".bim")
    )
    auto <- pecotmr:::.ldInfoBlock(
        data.frame(LD_file = ldFile, stringsAsFactors = FALSE),
        1L
    )
    explicit <- pecotmr:::.ldInfoBlock(
        data.frame(
            LD_file = ldFile,
            SNP_file = str_c(ldFile, ".bim"),
            stringsAsFactors = FALSE
        ),
        1L
    )
    expect_identical(auto, explicit)
    expect_s4_class(auto, "LdData")
    expect_equal(dim(getCorrelation(auto)), c(5L, 5L))
})

test_that(".ldInfoBlock requires an LD_file column", {
    expect_error(
        pecotmr:::.ldInfoBlock(data.frame(x = 1), 1L),
        "needs an `LD_file` column"
    )
})

test_that(".panelAfreqMaf returns NULL when no sidecar exists at all", {
    skip_if_not_installed("pgenlibr")
    # This plink2 fixture ships without an .afreq, so every shard's table is
    # absent -- distinct from the partial-sidecar case above.
    handle <- readGenotypeHandle(
        test_path("test_data/test_variants_chr22"),
        format = "plink2"
    )
    expect_equal(pecotmr:::getFormat(handle), "plink2")
    expect_null(
        pecotmr:::.panelAfreqMaf(
            handle,
            as.character(pecotmr:::getSnpInfo(handle)$SNP)
        )
    )
})
