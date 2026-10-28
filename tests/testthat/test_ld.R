context("LD")
library(tidyverse)

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

test_that("orderDedupRegions sorts X/Y/MT after autosomes and keeps chrom a string", {
    df <- data.frame(
        chrom = c("chrX", "chr2", "chr1", "chrM"),
        start = c(1L, 1L, 1L, 1L),
        end = c(2L, 2L, 2L, 2L)
    )
    result <- pecotmr:::orderDedupRegions(df)
    # Genomic order, not alphabetical and not NA-collapsed.
    expect_equal(result$chrom, c("chr1", "chr2", "chrX", "chrM"))
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
# resolveLdSource: type detection with real fixtures
# ===========================================================================

geno_test_data_dir <- test_path("test_data")
geno_region_all <- GenomicRanges::GRanges("chr21:17513228-17592874")

test_that("resolveLdSource tells the format of the genotypes metadata names", {
    # One row per chromosome naming its genotype file or PLINK prefix,
    # resolved beside the metadata file.
    rows <- c(
        plink2 = "21\t0\t0\ttest_variants",
        vcf = "21\t0\t0\ttest_variants.vcf.gz",
        gds = "21\t0\t0\ttest_variants.gds",
        plink1 = "22\t0\t0\tprotocol_example.genotype"
    )
    for (format in names(rows)) {
        meta_file <- withr::local_tempfile(
            tmpdir = geno_test_data_dir,
            fileext = ".tsv"
        )
        writeLines(c("chrom\tstart\tend\tpath", rows[[format]]), meta_file)
        result <- pecotmr:::resolveLdSource(meta_file)
        expect_equal(result$type, format, info = format)
        expect_equal(result$metaPath, meta_file, info = format)
    }
})

test_that("resolveLdSource tells a file's format as readGenotypes does", {
    # A .vcf.bgz is a VCF, and an extension's case does not matter.
    dir <- withr::local_tempdir()
    files <- file.path(dir, c("x.vcf.bgz", "x.VCF.GZ", "x.GDS", "x.Bgen"))
    file.create(files)
    expect_equal(
        purrr::map_chr(files, function(f) pecotmr:::resolveLdSource(f)$type),
        c("vcf", "vcf", "gds", "bgen")
    )
})

test_that("resolveLdSource refuses a row pointing at pre-computed LD", {
    # pecotmr reads LD only from genotypes, so a row naming a .cor.xz matrix
    # is an error rather than a source type.
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
    expect_error(
        pecotmr:::resolveLdSource(meta_file),
        "does not resolve to genotype files"
    )
})

test_that("resolveLdSource takes a genotype source given directly", {
    # A genotype file or PLINK prefix is its own source, with no metadata.
    stem <- str_remove(
        system.file("extdata", "toy_ref.bed", package = "pecotmr"),
        "\\.bed$"
    )
    res <- pecotmr:::resolveLdSource(stem)
    expect_equal(res$type, "plink1")
    expect_equal(res$dataPath, stem)
    expect_null(res$metaPath)
    vcf <- test_path("test_data", "test_variants.vcf.gz")
    expect_equal(pecotmr:::resolveLdSource(vcf)$type, "vcf")
})

test_that("resolveLdSource errors on missing file", {
    expect_error(
        pecotmr:::resolveLdSource("/nonexistent/file.tsv"),
        "LD metadata file: File does not exist"
    )
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

test_that("resolveGenotypePathForRegion takes a window on one chromosome", {
    # A per-chromosome meta names one file per chromosome, so a window over
    # two cannot resolve to one path; and a window is a GRanges, not a string.
    expect_error(
        pecotmr:::resolveGenotypePathForRegion(
            "unused.tsv",
            GenomicRanges::GRanges(c("chr21:1-10", "chr22:1-10"))
        ),
        "`region` spans chromosomes chr21, chr22"
    )
    expect_error(
        pecotmr:::resolveGenotypePathForRegion("unused.tsv", "chr21:1-10"),
        "`region` must be a GRanges"
    )
})

# ===========================================================================
# regionalLdMeta: real .cor.xz fixtures
# ===========================================================================

test_that("regionalLdMeta returns correct file paths for single block", {
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
    result <- pecotmr:::regionalLdMeta(
        meta_file,
        GenomicRanges::GRanges("chr1:1050-1150")
    )
    expect_true(is.list(result))
    expect_true(length(result$intersections$LD_file_paths) >= 1)
    # All returned paths should exist
    expect_true(all(file.exists(result$intersections$LD_file_paths)))
    expect_true(all(file.exists(result$intersections$bimFilePaths)))
})

test_that("regionalLdMeta spans multiple blocks for wide region", {
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
    result <- pecotmr:::regionalLdMeta(
        meta_file,
        GenomicRanges::GRanges("chr1:1000-1500")
    )
    expect_true(length(result$intersections$LD_file_paths) >= 2)
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
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(ncol(X)))
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
        "score.*Must have length 4"
    )
})

test_that("ldClumpByScore errors on mismatched chr/pos length", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    X <- matrix(rnorm(20), 5, 4)
    expect_error(
        ldClumpByScore(X, score = runif(4), chr = rep(1L, 2), pos = 1:4),
        "chr.*Must have length 4"
    )
})

# ===========================================================================
# ldPruneByCorrelation: verbose paths
# ===========================================================================

test_that("ldPruneByCorrelation reports pruning only when verbose", {
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
    expect_silent(ldPruneByCorrelation(X, corThres = 0.5))
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
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
    result <- ldPruneByCorrelation(X, corThres = 0.9)
    expect_true(ncol(result$X.new) < p)
    expect_equal(length(result$filter.id), ncol(result$X.new))
    expect_true(all(result$filter.id %in% seq_len(p)))
})

test_that("ldPruneByCorrelation keeps all columns when uncorrelated", {
    set.seed(42)
    n <- 100
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    result <- ldPruneByCorrelation(X, corThres = 0.99)
    expect_equal(ncol(result$X.new), p)
    expect_equal(result$filter.id, seq_len(p))
    expect_equal(colnames(result$X.new), colnames(X))
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

test_that("ldPruneByCorrelation keeps a single column, with nothing to prune", {
    # hclust() could not cluster one object, so this used to error.
    set.seed(42)
    n <- 30
    X <- matrix(rnorm(n), nrow = n, ncol = 1)
    colnames(X) <- "chr1:100:A:G"
    res <- ldPruneByCorrelation(X, corThres = 0.8)
    expect_identical(res$X.new, X)
    expect_identical(res$filter.id, 1L)
})

test_that("ldPruneByCorrelation strict threshold removes at least as many as lenient", {
    set.seed(42)
    n <- 100
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
    X[, 3] <- X[, 1] + rnorm(n, sd = 0.1)
    X[, 5] <- X[, 4] + rnorm(n, sd = 0.1)
    result_strict <- ldPruneByCorrelation(X, corThres = 0.3)
    result_lenient <- ldPruneByCorrelation(X, corThres = 0.99)
    expect_true(ncol(result_strict$X.new) <= ncol(result_lenient$X.new))
})

# ===========================================================================
# dropCollinearColumns
# ===========================================================================

# dropCollinearColumns and enforceDesignFullRank are unexported helpers;
# access them via pecotmr::: in these tests.

test_that("dropCollinearColumns returns X unchanged with no candidates", {
    X <- matrix(rnorm(100), nrow = 20, ncol = 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = character(0),
        strategy = "correlation"
    )
    expect_identical(result, X)
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

test_that("dropCollinearColumns correlation drops one of two candidates", {
    set.seed(42)
    X <- matrix(rnorm(100), nrow = 20, ncol = 5)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:5))
    candidates <- c("chr1:100:A:G", "chr1:200:A:G")
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = candidates,
        strategy = "correlation"
    )
    # Which of the two goes is a random pick, but it is always one of them,
    # and every other column comes through untouched.
    dropped <- setdiff(colnames(X), colnames(result))
    expect_length(dropped, 1L)
    expect_true(is_in(dropped, candidates))
    expect_identical(result, X[, colnames(X) != dropped])
})

test_that("dropCollinearColumns correlation drops the highest-sum candidate", {
    set.seed(42)
    n <- 50
    X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:4))
    # chr1:200 is built from both of the other candidates, which are
    # independent of each other, so its sum of |cor| with them is the
    # largest by a wide margin. It sits mid-list, so neither the first nor
    # the last candidate would pass for it.
    X[, 2] <- X[, 1] + X[, 3] + rnorm(n, sd = 0.1)
    result <- pecotmr:::dropCollinearColumns(
        X,
        problematicCols = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        strategy = "correlation"
    )
    expect_identical(result, X[, -2L])
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
        "response must be supplied"
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

test_that("ldClumpByScore drops the lower-scoring copy of a perfect-LD pair", {
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

# ===========================================================================
# Tests migrated from test_misc.R (computeLd + .findValidFilePath*)
# ===========================================================================

# The fixtures sit in test_data, which R CMD check copies with the tests,
# unlike the package root's own files.
.ldBed <- test_path("test_data", "protocol_example.genotype.bed")

test_that("findValidFilePath returns target when it exists directly", {
    target <- test_path("test_data", "protocol_example.genotype.bim")
    result <- pecotmr:::.findValidFilePath(
        referenceFilePath = .ldBed,
        targetFilePath = target
    )
    expect_equal(result, target)
})

test_that("findValidFilePath constructs path from reference directory", {
    result <- pecotmr:::.findValidFilePath(
        referenceFilePath = .ldBed,
        targetFilePath = "protocol_example.genotype.bim"
    )
    expect_equal(
        result,
        test_path("test_data", "protocol_example.genotype.bim")
    )
})

test_that("findValidFilePath errors when both paths are invalid", {
    expect_error(
        pecotmr:::.findValidFilePath(
            referenceFilePath = "/nonexistent/dir/ref.txt",
            targetFilePath = "target.txt"
        ),
        str_c(
            "Both reference and target file paths do not work. Tried paths: ",
            "'/nonexistent/dir/ref.txt' and '/nonexistent/dir/target.txt'"
        ),
        fixed = TRUE
    )
})

test_that("findValidFilePath returns reference when target is invalid but reference exists", {
    result <- pecotmr:::.findValidFilePath(
        referenceFilePath = .ldBed,
        targetFilePath = "/totally/bogus/path.txt"
    )
    expect_equal(result, .ldBed)
})

test_that("findValidFilePaths resolves multiple targets", {
    targets <- c(
        "protocol_example.genotype.bim",
        "protocol_example.genotype.fam"
    )
    result <- pecotmr:::.findValidFilePaths(.ldBed, targets)
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

test_that("computeLd with perfectly correlated SNPs returns correlation of 1", {
    set.seed(42)
    col1 <- sample(0:2, 50, replace = TRUE)
    X <- matrix(c(col1, col1), ncol = 2)
    colnames(X) <- c("rs1", "rs2")

    R <- computeLd(X, method = "sample")
    expect_equal(R[1, 2], 1.0, tolerance = 1e-10)
    expect_equal(R[2, 1], 1.0, tolerance = 1e-10)
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
# computeLd: trimSamples and shrinkage
# =============================================================================

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
    # 21 samples trim to the first 20, a multiple of four.
    expect_equal(R, computeLd(X[1:20, ], method = "population"))
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

test_that("computeLd on a GDS panel uses the in-memory path", {
    # The native SNPRelate route is no longer chosen implicitly by format:
    # it applies a different missing-data policy, so it is opt-in via
    # onDisk = TRUE rather than a silent consequence of the panel's format.
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    panel <- readGenotypes(test_path("test_data", "test_variants.gds"))
    R <- computeLd(panel, snpIdx = 1:6, backend = "internal")
    expect_equal(dim(R), c(6L, 6L))
    expect_true(all(is.finite(R)))
    expect_false(anyNA(R))
    expect_equal(unname(diag(R)), rep(1, 6), tolerance = 1e-6)
    expect_identical(rownames(R), rownames(panel)[1:6])
})

test_that("computeLd returns identity for a single-variant block", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    R <- computeLd(handle, snpIdx = 1L, backend = "internal")
    expect_equal(R, diag(1))
})

test_that("computeLd on a panel computes via the general path", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    R <- computeLd(handle, snpIdx = 1:5, backend = "internal")
    expect_equal(dim(R), c(5L, 5L))
    expect_equal(unname(diag(R)), rep(1, 5), tolerance = 1e-6)
})

# =============================================================================
# Additional coverage: findIntersectionRows / regionalLdMeta edge cases
# =============================================================================

test_that("findIntersectionRows errors when region falls in a coverage gap", {
    # Chromosome exists but there is a gap between the two blocks; after
    # clamping, no single row covers the query start -> error.
    gd <- data.frame(chrom = c(1, 1), start = c(100, 300), end = c(150, 350))
    expect_error(
        pecotmr:::findIntersectionRows(gd, 1, 200, 250),
        "not covered by any rows"
    )
})

test_that("regionalLdMeta handles whole-chromosome 0:0 sentinel rows", {
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
    result <- pecotmr:::regionalLdMeta(
        meta_file,
        GenomicRanges::GRanges("chr1:1000-1190")
    )
    expect_true(length(result$intersections$LD_file_paths) >= 1)
    expect_true(all(file.exists(result$intersections$LD_file_paths)))
    # 0:0 sentinel row should have had its end set to Inf internally
    expect_true(is.infinite(max(result$ldMetaData$end)))
})

test_that("regionalLdMeta validates complete coverage when required", {
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
    # Region fully inside block 1 -> validateSelectedRegion passes
    result <- pecotmr:::regionalLdMeta(
        meta_file,
        GenomicRanges::GRanges("chr1:1050-1150"),
        completeCoverageRequired = TRUE
    )
    expect_true(length(result$intersections$LD_file_paths) >= 1)
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
# Additional coverage: .ldFromRef onMissing="drop"
# =============================================================================

test_that(".ldFromRef drops variants the panel lacks under onMissing='drop'", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    si <- snpInfo(h)
    ids <- c(as.character(si$SNP[1]), "bogus:999:A:G", as.character(si$SNP[3]))
    m <- pecotmr:::.ldFromRef(h, ids, onMissing = "drop")
    expect_equal(dim(m), c(2L, 2L))
    expect_equal(
        attr(m, "keptVariantIds"),
        c(as.character(si$SNP[1]), as.character(si$SNP[3]))
    )
    expect_equal(unname(diag(m)), c(1, 1))
})

# =============================================================================
# Additional coverage: .ldFromRef reconciles chr-prefix conventions
# =============================================================================

test_that(".ldFromRef matches chr-prefixed ids to a non-prefixed panel", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(snpInfo(h)$SNP) # chr-prefixed (chr21_..._C_G)
    h@snpInfo$SNP <- sub("^chr", "", pid) # panel now lacks the prefix
    m <- pecotmr:::.ldFromRef(h, pid[1:3]) # request keeps the prefix
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rownames(m), pid[1:3]) # caller's labels preserved
    expect_equal(unname(diag(m)), c(1, 1, 1))
})

test_that(".ldFromRef matches non-prefixed ids to a chr-prefixed panel", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(snpInfo(h)$SNP)
    req <- sub("^chr", "", pid[1:3]) # request drops the prefix
    m <- pecotmr:::.ldFromRef(h, req)
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rownames(m), req) # caller's labels preserved
})

test_that(".ldFromRef still errors on a variant absent after reconciling", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    pid <- as.character(snpInfo(h)$SNP)
    # request the prefix-stripped form (forces reconciliation) plus one truly-absent variant
    req <- c(sub("^chr", "", pid[1]), "21_99999999_A_G")
    expect_error(
        pecotmr:::.ldFromRef(h, req),
        "not present in the LD reference"
    )
})

test_that(".ldFromRef rejects an ldReference that is not a genotype panel", {
    expect_error(
        pecotmr:::.ldFromRef(
            "not_a_handle",
            c("chr1:100:A:G", "chr1:200:A:G")
        ),
        "ldReference must be a genotype panel"
    )
})

# A synthetic panel + dosage extractor, so the shape of what .ldFromRef
# returns is checked without a plink2 fixture (and so without pgenlibr).
# @noRd
.lds_makeHandle <- function(snpN = 6L, nSamples = 30L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch",
        format = "plink2",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(snpN)),
            CHR = rep("chr1", snpN),
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
        nSnp <- nrow(snpInfo(handle))
        panel <- matrix(
            rbinom(nSamples * nSnp, 2, 0.3),
            nrow = nSamples,
            ncol = nSnp,
            dimnames = list(sampleIds(handle), snpInfo(handle)$SNP)
        )
        rr <- GenomicRanges::GRanges(
            seqnames = str_c("chr", snpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = snpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = snpInfo(handle)$SNP[snpIdx],
            A1 = snpInfo(handle)$A1[snpIdx],
            A2 = snpInfo(handle)$A2[snpIdx]
        )
        dosage <- t(panel[, snpIdx, drop = FALSE])
        dimnames(dosage) <- list(
            snpInfo(handle)$SNP[snpIdx],
            sampleIds(handle)
        )
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = S4Vectors::DataFrame(
                sampleId = sampleIds(handle),
                row.names = sampleIds(handle)
            )
        )
    }
}

test_that(".ldFromRef returns a symmetric unit-diagonal LD matrix", {
    h <- .lds_makeHandle()
    local_mocked_bindings(
        extractBlockGenotypes = .lds_mockExtractor(),
        .package = "pecotmr"
    )
    ids <- c("chr1:200:A:G", "chr1:400:A:G", "chr1:500:A:G")
    R <- pecotmr:::.ldFromRef(h, ids)
    expect_true(is.matrix(R))
    expect_equal(dimnames(R), list(ids, ids))
    expect_equal(unname(diag(R)), rep(1, 3), tolerance = 1e-12)
    expect_equal(R, t(R), tolerance = 1e-12)
})

test_that(".ldFromRef negates LD for a variant the panel spells flipped", {
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
    rSame <- pecotmr:::.ldFromRef(h, same)
    rFlip <- pecotmr:::.ldFromRef(h, flipped)
    sgn <- c(-1, 1, 1)
    expect_equal(
        unname(rFlip),
        unname(rSame) * outer(sgn, sgn),
        tolerance = 1e-12
    )
    expect_equal(unname(diag(rFlip)), rep(1, 3), tolerance = 1e-12)
    expect_equal(dimnames(rFlip), list(flipped, flipped))
})

test_that(".ldFromRef leaves LD untouched when every allele agrees", {
    h <- .lds_makeHandle()
    local_mocked_bindings(
        extractBlockGenotypes = .lds_mockExtractor(),
        .package = "pecotmr"
    )
    ids <- c("chr1:200:A:G", "chr1:400:A:G")
    # The sample correlation of the very dosages the panel serves for those
    # variants: with no allele to flip, no sign may change.
    served <- .lds_mockExtractor()(h, c(2L, 4L))
    direct <- cor(t(SummarizedExperiment::assay(served, "dosage")))
    expect_equal(pecotmr:::.ldFromRef(h, ids), direct, tolerance = 1e-12)
})

test_that(".ldFromRef errors on a variant the panel does not carry", {
    h <- .lds_makeHandle()
    expect_error(
        pecotmr:::.ldFromRef(h, c("chr1:100:A:G", "ghost")),
        "variant id.*not present in the LD reference"
    )
})

# =============================================================================
# Additional coverage: .requireMatchingLdRefs error paths
# =============================================================================

test_that(".requireMatchingLdRefs errors when slots are not GenotypeHandle", {
    expect_error(
        pecotmr:::.requireMatchingLdRefs(
            list(a = 1),
            list(b = 2),
            "testPipeline"
        ),
        "must both be genotype panels"
    )
})

test_that(".requireMatchingLdRefs errors on swapped alleles", {
    skip_if_not_installed("pgenlibr")
    h <- readGenotypeHandle(
        file.path(geno_test_data_dir, "test_variants"),
        format = "plink2"
    )
    si <- snpInfo(h)
    si2 <- si
    si2$A1[1] <- si$A2[1] # swap the coding of one shared variant
    si2$A2[1] <- si$A1[1]
    h2 <- new(
        "GenotypeHandle",
        path = path(h),
        format = genotypeFormat(h),
        snpInfo = si2,
        nSamples = nSamples(h),
        sampleIds = sampleIds(h),
        pgenPtr = NULL,
        chromPaths = character(0)
    )
    expect_error(
        pecotmr:::.requireMatchingLdRefs(h, h2, "testPipeline"),
        "swapped A1/A2"
    )
})

test_that(".requireMatchingLdRefs accepts panels trimmed differently", {
    # The reported bug: one LD sketch shared by a QtlSumStats and a
    # GwasSumStats, each trimmed by its own summaryStatsQc to a different
    # surviving variant set. The panels overlap but are not identical.
    si <- data.frame(
        SNP = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
        CHR = rep("chr1", 4L),
        BP = c(100L, 200L, 300L, 400L),
        A1 = c("A", "C", "G", "T"),
        A2 = c("G", "T", "A", "C"),
        stringsAsFactors = FALSE
    )
    mk <- function(rows) {
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "plink2",
            snpInfo = si[rows, , drop = FALSE],
            nSamples = 3L,
            sampleIds = str_c("s", 1:3),
            pgenPtr = NULL,
            chromPaths = character(0)
        )
    }
    expect_warning(
        expect_null(
            pecotmr:::.requireMatchingLdRefs(
                mk(1:3),
                mk(2:4),
                "trimPipelineA"
            )
        ),
        "share 2 variant"
    )
})

test_that(".requireMatchingLdRefs errors when panels share no variant", {
    mkSi <- function(bp) {
        data.frame(
            SNP = str_c("1:", bp, ":A:G"),
            CHR = rep("chr1", length(bp)),
            BP = as.integer(bp),
            A1 = rep("A", length(bp)),
            A2 = rep("G", length(bp)),
            stringsAsFactors = FALSE
        )
    }
    mk <- function(bp) {
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "plink2",
            snpInfo = mkSi(bp),
            nSamples = 3L,
            sampleIds = str_c("s", 1:3),
            pgenPtr = NULL,
            chromPaths = character(0)
        )
    }
    expect_error(
        pecotmr:::.requireMatchingLdRefs(
            mk(c(100, 200)),
            mk(c(900, 950)),
            "disjointPipeline"
        ),
        "share no variant"
    )
})

test_that(".requireMatchingLdRefs errors on a different sample set", {
    si <- data.frame(
        SNP = "1:100:A:G",
        CHR = "chr1",
        BP = 100L,
        A1 = "A",
        A2 = "G",
        stringsAsFactors = FALSE
    )
    mk <- function(ids) {
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "plink2",
            snpInfo = si,
            nSamples = length(ids),
            sampleIds = ids,
            pgenPtr = NULL,
            chromPaths = character(0)
        )
    }
    expect_error(
        pecotmr:::.requireMatchingLdRefs(
            mk(str_c("s", 1:3)),
            mk(str_c("t", 1:3)),
            "samplePipeline"
        ),
        "different sample sets"
    )
})

test_that(".requireMatchingLdRefs tolerates chr-prefixed SNP labels", {
    si1 <- data.frame(
        SNP = c("1:100:A:G", "1:200:C:T"),
        CHR = c("chr1", "chr1"),
        BP = c(100L, 200L),
        A1 = c("A", "C"),
        A2 = c("G", "T"),
        stringsAsFactors = FALSE
    )
    si2 <- si1
    si2$SNP <- c("chr1:100:A:G", "chr1:200:C:T") # same panel, chr-prefixed
    mk <- function(si) {
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "plink2",
            snpInfo = si,
            nSamples = 3L,
            sampleIds = paste0("s", 1:3),
            pgenPtr = NULL,
            chromPaths = character(0)
        )
    }
    expect_null(
        pecotmr:::.requireMatchingLdRefs(mk(si1), mk(si2), "testPipeline")
    )
})

# =============================================================================
# Additional coverage: filterVariantsByLdReference keepIndel = FALSE
# =============================================================================

test_that("filterVariantsByLdReference with keepIndel=FALSE drops indels", {
    meta_file <- file.path(geno_test_data_dir, "ld_meta_filtind_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(
        c(
            paste("chrom", "start", "end", "path", sep = "\t"),
            paste("21", "0", "0", "test_variants", sep = "\t")
        ),
        meta_file
    )
    # All three are on the reference panel, but the middle one is an indel.
    variant_ids <- c(
        "chr21:17513228:C:G",
        "chr21:17592874:G:GACTTT",
        "chr21:17513580:C:T"
    )
    result <- suppressMessages(
        filterVariantsByLdReference(variant_ids, meta_file, keepIndel = FALSE)
    )
    expect_false("chr21:17592874:G:GACTTT" %in% result$data)
    expect_true("chr21:17513228:C:G" %in% result$data)
    expect_true("chr21:17513580:C:T" %in% result$data)
})

test_that("hclust pruning keeps what single-linkage clustering keeps", {
    # The groups are found as threshold-graph components a stripe of the
    # correlation matrix at a time; 600 columns span three stripes, and the
    # AR chains carry LD across their boundaries.
    set.seed(9)
    n <- 60
    p <- 600
    X <- matrix(0, n, p)
    X[, 1] <- rnorm(n)
    for (j in 2:p) {
        X[, j] <- 0.93 * X[, j - 1] + rnorm(n, sd = 0.4)
    }
    X[, 400] <- X[, 3]
    X[, 401] <- 0
    for (thr in c(0.5, 0.8, 0.95)) {
        r <- suppressWarnings(cor(X))
        r[is.na(r)] <- 0
        tree <- hclust(as.dist(1 - abs(r)), method = "single")
        want <- which(!duplicated(cutree(tree, h = 1 - thr)))
        got <- ldPruneByCorrelation(X, corThres = thr)
        expect_equal(got$filter.id, want, label = str_c("corThres ", thr))
        expect_equal(got$X.new, X[, want, drop = FALSE])
    }
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
    colnames(X) <- paste0("snp", seq_len(p))
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

test_that("the genotype-code backends refuse dosages rather than round them", {
    # GDS and SnpMatrix hold 2-bit codes: a fractional dosage would be
    # rounded and its LD changed, and a 3 would read as missing.
    set.seed(1)
    X <- matrix(rbinom(100 * 4, 2, 0.3), 100, 4)
    dosage <- replace(X, 5L, 0.7)
    outOfRange <- replace(X, 7L, 3)
    expect_error(
        pecotmr:::.ldRequireIntegerGenotypes(dosage, "snpstats", "internal"),
        "1 value\\(s\\) here are not one of them"
    )
    expect_error(
        pecotmr:::.ldRequireIntegerGenotypes(outOfRange, "snpstats", "x"),
        "not one of them"
    )
    expect_true(pecotmr:::.ldRequireIntegerGenotypes(
        replace(X, 3L, NA),
        "snpstats",
        "internal"
    ))
    skip_if_not_installed("snpStats")
    expect_error(
        computeLd(dosage, backend = "snpstats"),
        "Use backend = \"internal\""
    )
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    expect_error(
        computeLd(dosage, backend = "snprelate"),
        "backend = \"snprelate\" stores genotypes"
    )
    expect_error(
        ldPruneByCorrelation(dosage, backend = "snprelate"),
        "Use backend = \"hclust\""
    )
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
    d <- exec(cbind, !!!map(af, function(f) rbinom(nS, 2L, f)))
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
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    ids <- rownames(panel)
    expect_identical(.panelVariantFilter(panel, ids), ids)
})

test_that(".panelVariantFilter drops panel-rare variants", {
    data(qtlDatasetExample)
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    ids <- rownames(panel)
    loose <- .panelVariantFilter(
        panel,
        ids,
        PanelFilterParam(mafCutoff = 0.05)
    )
    tight <- .panelVariantFilter(
        panel,
        ids,
        PanelFilterParam(mafCutoff = 0.2)
    )
    expect_lt(length(loose), length(ids))
    expect_lt(length(tight), length(loose))
    # Kept sets are nested as the cutoff rises, and order is the caller's.
    expect_true(all(is_in(tight, loose)))
    expect_identical(loose, ids[is_in(ids, loose)])
})

test_that(".panelVariantFilter treats MAC as a MAF equivalent", {
    data(qtlDatasetExample)
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    ids <- rownames(panel)
    nSamp <- ncol(panel)
    # macCutoff / (2 * nSamples) is the same threshold as mafCutoff.
    byMac <- .panelVariantFilter(
        panel,
        ids,
        PanelFilterParam(macCutoff = 0.1 * 2 * nSamp)
    )
    byMaf <- .panelVariantFilter(
        panel,
        ids,
        PanelFilterParam(mafCutoff = 0.1)
    )
    expect_identical(byMac, byMaf)
    # The stricter of the two wins.
    expect_identical(
        .panelVariantFilter(
            panel,
            ids,
            PanelFilterParam(mafCutoff = 0.2, macCutoff = 2)
        ),
        .panelVariantFilter(panel, ids, PanelFilterParam(mafCutoff = 0.2))
    )
})

test_that(".panelVariantFilter drops high-missingness variants", {
    data(qtlDatasetExample)
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    ids <- rownames(panel)
    strict <- .panelVariantFilter(
        panel,
        ids,
        PanelFilterParam(imissCutoff = 0)
    )
    expect_lt(length(strict), length(ids))
    # A cutoff above the panel's worst variant keeps everything.
    expect_identical(
        .panelVariantFilter(panel, ids, PanelFilterParam(imissCutoff = 1)),
        ids
    )
})

test_that(".panelVariantFilter passes through ids absent from the panel", {
    # Whether a missing variant is an error or is dropped belongs to
    # .ldFromRef's `onMissing`; deciding it here too would let the two
    # disagree about the same variant.
    data(qtlDatasetExample)
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    ids <- rownames(panel)[1:3]
    withGhost <- c("chr9:999:A:G", ids)
    expect_true(is_in(
        "chr9:999:A:G",
        .panelVariantFilter(
            panel,
            withGhost,
            PanelFilterParam(mafCutoff = 0.001)
        )
    ))
})

test_that(".panelVariantFilter handles empty and NULL input", {
    data(qtlDatasetExample)
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    expect_length(
        .panelVariantFilter(
            panel,
            character(0),
            PanelFilterParam(mafCutoff = 0.1)
        ),
        0L
    )
    expect_identical(
        .panelVariantFilter(
            NULL,
            "chr1:1:A:G",
            PanelFilterParam(mafCutoff = 0.1)
        ),
        "chr1:1:A:G"
    )
})

test_that(".panelCutoffs short-circuits when no cutoff is set", {
    # NULL means the panel is never touched, which is what keeps the default
    # path free of an extra dosage read.
    expect_null(.panelCutoffs())
    expect_null(.panelCutoffs(PanelFilterParam(
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1
    )))
    expect_equal(
        .panelCutoffs(PanelFilterParam(mafCutoff = 0.01))$mafCutoff,
        0.01
    )
    expect_equal(
        .panelCutoffs(PanelFilterParam(imissCutoff = 0.5))$imissCutoff,
        0.5
    )
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

# .panelAfreqMaf() for an LD reference -- a panel or a bare handle: its
# genotype source and chromosomes.
# @noRd
.pvfAfreqMaf <- function(ref, ids) {
    .panelAfreqMaf(
        pecotmr:::.ldRefSource(ref),
        pecotmr:::.ldRefChrom(ref),
        ids
    )
}

test_that(".panelAfreqMaf reads panel MAF from the .afreq sidecar", {
    skip_if_not_installed("pgenlibr")
    handle <- .pvfAfreqHandle()
    ids <- as.character(snpInfo(handle)$SNP)
    maf <- .pvfAfreqMaf(handle, ids)
    afreq <- readAfreq(test_path("test_data/test_variants"))
    altFreq <- afreq$alt_freq[match(ids, afreq$id)]
    expect_equal(maf, pmin(altFreq, 1 - altFreq))
})

test_that(".panelAfreqMaf refuses a partial sidecar and reads any format's", {
    skip_if_not_installed("pgenlibr")
    handle <- .pvfAfreqHandle()
    ids <- as.character(snpInfo(handle)$SNP)
    # One id the sidecar cannot answer for is enough: a partial answer would
    # disagree with the dosage path about that variant, so NULL sends the
    # caller down the dosage path for all of them.
    expect_null(.pvfAfreqMaf(handle, c(ids, "not-a-panel-variant")))
    # The sidecar is the file's, not the format's: the PLINK1 fileset beside
    # the same .afreq reads the same frequencies.
    skip_if_not_installed("PlinkMatrix")
    bed <- readGenotypes(plink1Prefix = test_path("test_data/test_variants"))
    bedIds <- as.character(SummarizedExperiment::rowRanges(bed)$SNP)
    expect_equal(.pvfAfreqMaf(bed, bedIds), .pvfAfreqMaf(handle, bedIds))
})

test_that(".panelVariantFilter: .afreq and dosage agree on what to drop", {
    skip_if_not_installed("pgenlibr")
    handle <- .pvfAfreqHandle()
    ids <- as.character(snpInfo(handle)$SNP)
    viaAfreq <- .panelVariantFilter(
        handle,
        ids,
        PanelFilterParam(mafCutoff = 0.2)
    )
    # Force the dosage path by hiding the sidecar from the fast path.
    local_mocked_bindings(
        .panelAfreqMaf = function(source, chroms, variantIds) NULL,
        .package = "pecotmr"
    )
    viaDosage <- .panelVariantFilter(
        handle,
        ids,
        PanelFilterParam(mafCutoff = 0.2)
    )
    expect_lt(length(viaAfreq), length(ids))
    expect_identical(viaAfreq, viaDosage)
})

test_that(".panelVariantFilter uses dosage whenever missingness is capped", {
    skip_if_not_installed("pgenlibr")
    # Missingness is not derivable from the sidecar, so an imissCutoff must
    # never be silently answered by it.
    handle <- .pvfAfreqHandle()
    ids <- as.character(snpInfo(handle)$SNP)
    local_mocked_bindings(
        .panelAfreqMaf = function(handle, variantIds) {
            abort("the .afreq fast path must not run here")
        },
        .package = "pecotmr"
    )
    expect_no_error(
        .panelVariantFilter(
            handle,
            ids,
            PanelFilterParam(mafCutoff = 0.2, imissCutoff = 0.5)
        )
    )
})

test_that(".afreqPrefixes reads only the chromosomes the panel spans", {
    handle <- .pvfAfreqHandle()
    source <- pecotmr:::.genotypeSourceRecord(handle)
    # Single-file panel: its own stem, resolved for file access.
    expect_identical(
        .afreqPrefixes(source, snpInfo(handle)$CHR),
        pecotmr:::.genotypeReadPath(handle)
    )
    # Sharded panel: a manifest chromosome the panel does not span is not
    # read -- the cost per-chromosome sketch trimming exists to avoid.
    sharded <- list_assign(
        source,
        chromPaths = c("22" = path(handle), "21" = "/nonexistent/chr21")
    )
    expect_identical(.afreqPrefixes(sharded, "22"), path(handle))
    # A panel held in memory reads no sidecar.
    expect_identical(.afreqPrefixes(NULL, "22"), character(0))
})

# =============================================================================
# One failure mode for a missing or invalid LD sketch
#
# Panel access happens through several routes: the shared `.ldFromRef()`,
# ctwas's whole-panel assembler, RAiSS's coordinate window, and the bare
# accessors. They used to disagree on what a bad sketch looked like -- the
# shared path said so plainly while the others surfaced "unable to find an
# inherited method for 'snpInfo'" from whichever accessor touched it
# first. The guard now sits at `.ldRefRanges()` / `.ldRefDosage()`, the
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
            pecotmr:::.ldFromRef(NULL, "chr22:1:A:G", label = "demo")
        },
        ctwas = function() pecotmr:::.ctwasPanelFor(NULL),
        raiss = function() pecotmr:::.qcRaissWindowIdx(df, NULL, 0L),
        matchIds = function() pecotmr:::.ldRefMatchIds(NULL),
        dosage = function() pecotmr:::.ldRefDosage(NULL, 1L),
        ranges = function() pecotmr:::.ldRefRanges(NULL)
    )
    for (nm in names(routes)) {
        expect_error(routes[[nm]](), "carries no ldReference", info = nm)
    }
})

test_that("a non-panel sketch is rejected by every route", {
    expect_error(
        pecotmr:::.ldRefRanges("not a panel"),
        "must be a genotype panel"
    )
    expect_error(
        pecotmr:::.ctwasPanelFor("not a panel"),
        "must be a genotype panel"
    )
    expect_error(
        pecotmr:::.ldRefDosage("not a panel", 1L),
        "must be a genotype panel"
    )
})

test_that("callers that validate first keep their own label", {
    # The generic guard must not mask a caller's more specific message.
    shared <- tryCatch(
        pecotmr:::.ldFromRef(NULL, "x", label = "fineMappingPipeline"),
        error = conditionMessage
    )
    expect_match(shared, "^fineMappingPipeline:")
    ctwas <- tryCatch(
        pecotmr:::.ctwasPanelFor(NULL),
        error = conditionMessage
    )
    expect_match(ctwas, "^ctwasPipeline:")
})

test_that("panel identity uses the allele-repaired ids for comparison", {
    # `.ldRefMatchIds()` repairs an id whose allele fields are tags rather
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

    # And an rsID panel passes through unchanged -- `.ldFromRefMatch()`
    # falls back to exact id matching for those, so rewriting them would
    # break the match rather than fix it.
    panel <- readGenotypes(
        system.file("extdata", "toy_ref.bed", package = "pecotmr")
    )
    expect_equal(
        pecotmr:::.ldRefMatchIds(panel),
        pecotmr:::.ldRefVariantIds(panel)
    )
})

# ===========================================================================
# Sketch accessors and per-variant panel statistics
# ===========================================================================

test_that("the sketch accessors are NULL-safe", {
    # A collection with no LD sketch is a legitimate state (individual-level
    # inputs), so these answer NULL rather than erroring.
    expect_null(pecotmr:::.ldRefHandle(NULL))
    expect_null(pecotmr:::.ldRefSubset(NULL, 1L))
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

# ===========================================================================
# Empty sketches, the .afreq fast path, and design rank
# ===========================================================================

test_that("emptying a sketch clears both axes", {
    # An emptied sketch references no LD, so its sample axis is dead weight;
    # nSamples must go to 0 alongside sampleIds or the derived dosage
    # dimnames would disagree with the matrix.
    expect_null(pecotmr:::.emptyLdRef(NULL))
    h <- new(
        "GenotypeHandle",
        path = "/tmp/x",
        format = "plink2",
        snpInfo = data.frame(
            SNP = c("a", "b"),
            CHR = "chr1",
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
    e <- pecotmr:::.emptyLdRef(h)
    expect_equal(nrow(pecotmr:::snpInfo(e)), 0L)
    expect_equal(length(pecotmr:::sampleIds(e)), 0L)
})

test_that(".emptyLdRef drops the seed, not just the rows", {
    # A panel carries its handle twice: as rowRanges and inside the dosage
    # assay's DelayedArray seed. Dropping rows alone leaves the seed holding
    # the whole panel, which is the weight this exists to shed.
    h <- .sk_panel()
    e <- .emptyLdRef(h)
    expect_equal(nrow(e), 0L)
    expect_null(.panelSeed(e))
    expect_lt(
        length(serialize(e, NULL)),
        length(serialize(h, NULL)) / 2
    )
})

test_that(".emptyLdRef is NULL-safe and idempotent", {
    expect_null(.emptyLdRef(NULL))
    e <- .emptyLdRef(.sk_panel())
    expect_equal(nrow(.emptyLdRef(e)), 0L)
})

test_that(".emptyLdRef drops both axes and the seed on an RSE sketch", {
    # The seed carries every variant's metadata; an empty panel keeps none
    # of it, only the record of what it was read from.
    gh <- .gh_pgenHandle()
    empty <- pecotmr:::.emptyLdRef(pecotmr:::.genotypeExperiment(gh))
    expect_equal(dim(empty), c(0L, 0L))
    expect_null(pecotmr:::.panelSeed(empty))
    expect_null(pecotmr:::.ldRefHandle(empty))
    expect_equal(pecotmr:::.ldRefFormat(empty), genotypeFormat(gh))
})

test_that("a missing or malformed afreq sidecar degrades to no frequencies", {
    # The sidecar is optional; a bad one must not take the whole filter down,
    # it just falls back to reading dosage.
    expect_null(pecotmr:::.panelAfreqTable("/nonexistent/prefix"))
    prefix <- file.path(withr::local_tempdir(), "panel")
    sidecar <- str_c(prefix, ".afreq")
    h <- new(
        "GenotypeHandle",
        path = prefix,
        format = "plink2",
        snpInfo = data.frame(
            SNP = c("a", "b"),
            CHR = "chr1",
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
    # No sidecar beside the file: nothing to read, in any format.
    expect_null(.pvfAfreqMaf(h, c("a", "b")))
    # A well-formed sidecar there is read ...
    writeLines(
        c(
            "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT",
            "1\ta\tG\tA\t0.1\t4",
            "1\tb\tG\tA\t0.8\t4"
        ),
        sidecar
    )
    expect_equal(.pvfAfreqMaf(h, c("a", "b")), c(0.1, 0.2))
    # ... and a malformed one degrades to no frequencies: one without
    # ALT_FREQS, which the reader refuses, or one that is not a table.
    writeLines(
        c(
            "#CHROM\tID\tREF\tALT\tOBS_CT",
            "1\ta\tG\tA\t4",
            "1\tb\tG\tA\t4"
        ),
        sidecar
    )
    expect_error(readAfreq(prefix), "ALT_FREQS")
    expect_null(pecotmr:::.panelAfreqTable(prefix))
    expect_null(.pvfAfreqMaf(h, c("a", "b")))
    writeLines("not an afreq file", sidecar)
    expect_null(pecotmr:::.panelAfreqTable(prefix))
    expect_null(.pvfAfreqMaf(h, c("a", "b")))
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

test_that(".ldRefMatchIds returns bare ids when alleles are absent", {
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
        pecotmr:::.ldRefMatchIds(se),
        c("chr1:100:A:G", "chr1:200:C:T")
    )
})

test_that(".ldRefMatchIds repairs tag ids and leaves clean labels alone", {
    h <- .ssQ_makeHandleVid() # ids chr1:100:G:A .. chr1:800:G:A, A1=A A2=G
    h@snpInfo$SNP[h@snpInfo$BP == 200L] <- "chr1:200:INS:A"
    ids <- pecotmr:::.ldRefMatchIds(h)
    # the tag id is re-rendered from the panel's own A2/A1 columns ...
    expect_equal(ids[2L], "chr1:200:G:A")
    # ... and every other label is untouched.
    expect_equal(ids[-2L], pecotmr:::.ldRefVariantIds(h)[-2L])
    # the raw accessor still answers the panel FILE's label (.afreq keys on it)
    expect_equal(pecotmr:::.ldRefVariantIds(h)[2L], "chr1:200:INS:A")
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
        pecotmr:::.panelVariantFilter(
            sketch,
            ids,
            PanelFilterParam(mafCutoff = 0.01)
        ),
        ids
    )
})

test_that(".ldRefNullGuard names the label in its strict error", {
    expect_error(
        pecotmr:::.ldRefNullGuard("Q", NULL, "myPipe", "myLabel", "strict"),
        "ldReference on `myLabel` is non-NULL"
    )
    expect_error(
        pecotmr:::.ldRefNullGuard("Q", NULL, "myPipe", NULL, "strict"),
        "qtl ldReference is non-NULL"
    )
    # lenient tolerates the same input.
    expect_true(
        pecotmr:::.ldRefNullGuard("Q", NULL, "myPipe", "myLabel", "lenient")
    )
})

test_that(".ldRefCheckContent rejects panels on different chromosomes", {
    # Disjoint chromosomes leave the two panels with nothing in common, which
    # is a different LD reference rather than a trimming difference.
    mkSe <- function(chr, pos) {
        g <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1))
        S4Vectors::mcols(g)$SNP <- str_c(chr, ":", pos, ":A:G")
        SummarizedExperiment::SummarizedExperiment(
            assays = list(x = matrix(0, nrow = length(pos), ncol = 1L)),
            rowRanges = g
        )
    }
    expect_error(
        pecotmr:::.ldRefCheckContent(
            mkSe("chr1", c(100, 200)),
            mkSe("chr2", c(100, 200)),
            "myPipe",
            " between X and Y"
        ),
        "share no variant between X and Y"
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
    skip_if_not_installed("GDSArray")
    panel <- readGenotypes(test_path("test_data", "test_variants.gds"))
    R <- computeLd(panel, backend = "snprelate", onDisk = TRUE)
    expect_equal(dim(R), c(nrow(panel), nrow(panel)))
    expect_identical(
        rownames(R),
        SummarizedExperiment::rowRanges(panel)$SNP
    )
})

test_that("computeLd(onDisk) refuses a SeqArray GDS panel", {
    # SNPRelate's on-disk routine reads its own SNP GDS files; a SeqArray
    # file is refused by name rather than inside snpgdsOpen().
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("SeqArray")
    skip_if_not_installed("GDSArray")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), c("s1", "s2", "s3"))
    panel <- suppressWarnings(readGenotypes(.writeTestSeqArray(vcf)))
    expect_error(
        computeLd(panel, backend = "snprelate", onDisk = TRUE),
        "is a SeqArray GDS file"
    )
})

test_that("cutting a VCF panel's rows cuts its seed too", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), c("s1", "s2", "s3"))
    panel <- suppressWarnings(readGenotypes(vcf))
    dosageSeed <- function(x) {
        DelayedArray::seed(SummarizedExperiment::assay(x, "dosage"))
    }
    dose <- function(x) {
        unname(as.matrix(SummarizedExperiment::assay(x, "dosage")))
    }
    cut <- pecotmr:::.panelSubsetRows(panel, c(2L, 6L, 7L))
    expect_s4_class(dosageSeed(cut), "GenotypeDosageSeed")
    expect_equal(dim(dosageSeed(cut)), c(3L, 3L))
    expect_equal(dose(cut), .vcfOverlapGt()[c(2L, 6L, 7L), ])
    expect_equal(S4Vectors::metadata(cut), S4Vectors::metadata(panel))
    # A cut that reorders rows is taken from the assay, over the whole seed.
    reordered <- pecotmr:::.panelSubsetRows(panel, c(6L, 2L))
    expect_equal(dim(dosageSeed(reordered)), c(7L, 3L))
    expect_equal(dose(reordered), .vcfOverlapGt()[c(6L, 2L), ])
})

test_that("computeLd(onDisk) applies shrinkage toward the identity", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("GDSArray")
    handle <- readGenotypes(test_path("test_data", "test_variants.gds"))
    idx <- 1:6
    plain <- computeLd(
        handle,
        snpIdx = idx,
        backend = "snprelate",
        onDisk = TRUE
    )
    shrunk <- computeLd(
        handle,
        snpIdx = idx,
        backend = "snprelate",
        onDisk = TRUE,
        shrinkage = 0.5
    )
    # (1 - s) * R + s * I: off-diagonals halve, the diagonal stays 1.
    expect_equal(unname(diag(shrunk)), rep(1, length(idx)))
    expect_equal(
        unname(shrunk[upper.tri(shrunk)]),
        unname(0.5 * plain[upper.tri(plain)])
    )
})

test_that("onDisk LD matches the in-memory path on complete data", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("GDSArray")
    handle <- readGenotypes(test_path("test_data", "test_variants.gds"))
    idx <- c(4L, 1L, 6L, 2L, 5L, 3L)
    onDisk <- computeLd(
        handle,
        snpIdx = idx,
        backend = "snprelate",
        onDisk = TRUE
    )
    inMem <- computeLd(handle, snpIdx = idx)
    # The fixture has no missing calls, which is the only regime where the
    # two are meant to agree.
    expect_equal(unname(onDisk), unname(inMem), tolerance = 1e-8)
})

test_that("onDisk LD labels and orders by the request", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("GDSArray")
    handle <- readGenotypes(test_path("test_data", "test_variants.gds"))
    idx <- c(4L, 1L, 6L, 2L, 5L, 3L)
    R <- computeLd(handle, snpIdx = idx, backend = "snprelate", onDisk = TRUE)
    # snpgdsLDMat() returns file order and labels nothing; both are corrected.
    expect_identical(
        rownames(R),
        SummarizedExperiment::rowRanges(handle)$SNP[idx]
    )
    expect_identical(colnames(R), rownames(R))
})

test_that("onDisk LD is refused where it cannot be honoured", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("GDSArray")
    skip_if_not_installed("PlinkMatrix")
    gds <- readGenotypes(test_path("test_data", "test_variants.gds"))
    plink <- readGenotypes(
        plink1Prefix = test_path("test_data", "test_variants")
    )
    expect_error(
        computeLd(plink, backend = "snprelate", onDisk = TRUE),
        "GDS-backed panel"
    )
    expect_error(
        computeLd(gds, backend = "internal", onDisk = TRUE),
        "only available for"
    )
    expect_error(
        computeLd(matrix(1, 4, 3), backend = "snprelate", onDisk = TRUE),
        "must be one"
    )
    expect_error(
        computeLd(gds, method = "gcta", backend = "snprelate", onDisk = TRUE),
        "sample correlation"
    )
})

test_that("LD from a sketch does not depend on the request order", {
    skip_if_not_installed("pgenlibr")
    handle <- readGenotypeHandle(
        test_path("test_data", "test_variants"),
        format = "plink2"
    )
    ids <- snpInfo(handle)$SNP[1:6]
    shuffled <- ids[c(4L, 1L, 6L, 2L, 5L, 3L)]
    sorted <- pecotmr:::.ldFromRef(handle, ids, label = "t")
    mixed <- pecotmr:::.ldFromRef(handle, shuffled, label = "t")
    expect_identical(rownames(mixed), shuffled)
    expect_equal(sorted[ids, ids], mixed[ids, ids])
})

test_that("a PLINK1 panel cut to a few variants keeps only those", {
    # Its seed is cut with it, so a saved panel carries the kept variants'
    # metadata only, and reads them back from the file.
    skip_if_not_installed("PlinkMatrix")
    panel <- readGenotypes(
        plink1Prefix = test_path("test_data", "test_variants")
    )
    cut <- pecotmr:::.panelSubsetRows(panel, c(3L, 10L, 200L))
    # 100 samples in the panel
    expect_equal(dim(pecotmr:::.panelSeed(cut)), c(3L, 100L))
    expect_equal(.dose(cut), .dose(panel)[c(3L, 10L, 200L), ])
    saved <- tempfile(fileext = ".rds")
    saveRDS(cut, saved)
    expect_equal(.dose(readRDS(saved)), .dose(cut))
})

test_that("ldReference cross-checks are no-ops when a side has no variants", {
    # A zero-variant panel -- what an emptied object carries -- contradicts
    # nothing, so there is nothing to compare and the check returns quietly.
    local_mocked_bindings(.ldRefA1 = function(x) character(0))
    expect_null(pecotmr:::.ldRefCheckAlleleCoding(
        NULL,
        NULL,
        integer(0),
        integer(0),
        "testPipeline",
        ""
    ))
})

test_that(".ldRefCheckOverlap is a no-op when a side has no variants", {
    local_mocked_bindings(.ldRefVariantKeys = function(x) character(0))
    expect_null(
        pecotmr:::.ldRefCheckOverlap(NULL, NULL, "testPipeline", "")
    )
})

test_that("LdPruningOptions refuses what the snprelate backend owns", {
    expect_error(LdPruningOptions(gdsobj = 1), "the temporary GDS")
    expect_error(LdPruningOptions(method = "r"), "fixed at 'corr'")
    expect_error(
        LdPruningOptions(ld.threshold = 0.5),
        "the caller's `corThres`"
    )
    expect_error(LdPruningOptions(verbose = TRUE), "the caller's `verbose`")
    expect_error(LdPruningOptions(nosuch = 1), "unknown argument")
    expect_equal(LdPruningOptions(slide.max.bp = 1e6)$slide.max.bp, 1e6)
})

test_that("ldPruneByCorrelation refuses options its backend never uses", {
    set.seed(1)
    X <- matrix(rnorm(100 * 5), 100, 5)
    expect_error(
        ldPruneByCorrelation(
            X,
            backend = "hclust",
            methodArgs = LdPruningOptions(slide.max.n = 10L)
        ),
        "backend 'hclust' does not call"
    )
    # No options is the default, so the hclust path still works untouched.
    expect_length(ldPruneByCorrelation(X, corThres = 0.9)$filter.id, 5L)
    expect_error(
        ldPruneByCorrelation(X, methodArgs = list(slide.max.n = 10L)),
        "LdPruningOptions"
    )
})

test_that("ldPruneByCorrelation forwards methodArgs to SNPRelate", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    seen <- NULL
    real <- SNPRelate::snpgdsLDpruning
    set.seed(2)
    X <- matrix(sample(0:2, 100 * 6, replace = TRUE), 100, 6)
    colnames(X) <- str_c("snp", seq_len(6))
    suppressMessages(with_mocked_bindings(
        ldPruneByCorrelation(
            X,
            corThres = 0.9,
            backend = "snprelate",
            methodArgs = LdPruningOptions(slide.max.n = 3L)
        ),
        snpgdsLDpruning = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "SNPRelate"
    ))
    expect_equal(seen$slide.max.n, 3L)
    expect_equal(seen$method, "corr")
    expect_equal(seen$ld.threshold, 0.9)
})

# =============================================================================
# The LD reference dispatch: a genotype panel and an LdMatrix answer alike
# =============================================================================

# A panel of real genotypes, and the LdMatrix holding exactly its LD.
.ldrefPair <- function(nSamples = 100L) {
    panel <- readGenotypes(file.path(geno_test_data_dir, "test_variants.bed"))
    ld <- LdMatrix(
        unname(computeLd(panel)),
        variants = SummarizedExperiment::rowRanges(panel),
        nSamples = nSamples
    )
    list(panel = panel, ld = ld)
}

# The stochastic-sketch fixture: 100 pseudo-samples projected from 50 real
# ones, with an .afreq that carries U_MIN / U_MAX.
.ldrefSketch <- function() {
    readGenotypes(file.path(geno_test_data_dir, "test_harmonize_regions.pgen"))
}

test_that("an LdMatrix gives the same LD as the panel it was built from", {
    pr <- .ldrefPair()
    ids <- names(SummarizedExperiment::rowRanges(pr$panel))[c(9, 2, 40, 5)]
    expect_equal(
        pecotmr:::.ldFromRef(pr$ld, ids, label = "t"),
        pecotmr:::.ldFromRef(pr$panel, ids, label = "t")
    )
    # A request in the opposite allele coding is re-signed the same way.
    p <- parseVariantId(ids[[1L]])
    swapped <- replace(
        ids,
        1L,
        str_c("chr", p$chrom, ":", p$pos, ":", p$A1, ":", p$A2)
    )
    flipped <- pecotmr:::.ldFromRef(pr$ld, swapped, label = "t")
    expect_equal(
        unname(flipped[1, -1]),
        -unname(pecotmr:::.ldFromRef(pr$panel, ids, label = "t")[1, -1])
    )
    expect_error(
        pecotmr:::.ldFromRef(pr$ld, c(ids, "chr21:1:A:G"), label = "t"),
        "not present in the LD reference"
    )
})

test_that("a genotype reference's factor reproduces its LD", {
    # crossprod(F) / (nrow(F) - 1) is the LD itself, and F's columns have
    # zero mean, so an engine that correlates F's rows sees the same LD.
    pr <- .ldrefPair()
    sketch <- .ldrefSketch()
    for (ref in list(pr$panel, sketch)) {
        idx <- c(7L, 1L, 8L, 2L, 3L)
        f <- pecotmr:::.ldRefFactor(ref, idx)
        expect_equal(
            unname(crossprod(f)) / (nrow(f) - 1),
            unname(pecotmr:::.ldRefLdAt(ref, idx))
        )
        expect_equal(unname(colMeans(f)), rep(0, length(idx)))
    }
})

test_that("a monomorphic variant's factor column is padded to unit variance", {
    # computeLd() treats a monomorphic variant as independent of every other
    # with unit variance; its standardized column is zero, so the factor adds
    # two rows of +/- sqrt((m - 1) / 2) for it and rescales the rest.
    set.seed(5)
    X <- matrix(sample(0:2, 60, replace = TRUE), nrow = 15)
    X[, 3] <- 1
    X[c(2, 9), 1] <- NA
    f <- pecotmr:::.ldFactorFromDosage(X)
    expect_equal(nrow(f), 17L)
    expect_equal(crossprod(f) / (nrow(f) - 1), unname(computeLd(X)))
    expect_equal(colMeans(f), rep(0, 4))
    # Nothing to pad: the factor is the standardized dosages themselves.
    plain <- pecotmr:::.ldFactorFromDosage(X[, -3])
    expect_equal(nrow(plain), 15L)
    expect_equal(crossprod(plain) / 14, unname(computeLd(X[, -3])))
})

test_that(".ldFactorFromRef aligns the factor as .ldFromRef aligns the LD", {
    pr <- .ldrefPair()
    ids <- names(SummarizedExperiment::rowRanges(pr$panel))[c(9, 2, 40, 5)]
    p <- parseVariantId(ids[[1L]])
    swapped <- replace(
        ids,
        1L,
        str_c("chr", p$chrom, ":", p$pos, ":", p$A1, ":", p$A2)
    )
    f <- pecotmr:::.ldFactorFromRef(pr$panel, swapped, label = "t")
    expect_identical(colnames(f), swapped)
    expect_equal(
        crossprod(f) / (nrow(f) - 1),
        pecotmr:::.ldFromRef(pr$panel, swapped, label = "t")
    )
    dropped <- pecotmr:::.ldFactorFromRef(
        pr$panel,
        c(ids, "chr21:1:A:G"),
        label = "t",
        onMissing = "drop"
    )
    expect_identical(attr(dropped, "keptVariantIds"), ids)
    expect_null(
        pecotmr:::.ldFactorFromRef(
            pr$panel,
            "chr21:1:A:G",
            label = "t",
            onMissing = "drop"
        )
    )
    expect_error(
        pecotmr:::.ldFactorFromRef(pr$ld, ids, label = "t"),
        "Only a single genotype panel"
    )
})

test_that("an LdMatrix answers the variant-table primitives as a panel does", {
    pr <- .ldrefPair()
    expect_identical(
        pecotmr:::.ldRefMatchIds(pr$ld),
        pecotmr:::.ldRefMatchIds(pr$panel)
    )
    expect_identical(
        pecotmr:::.ldRefChrom(pr$ld),
        pecotmr:::.ldRefChrom(pr$panel)
    )
    expect_identical(
        pecotmr:::.ldRefVariantKeys(pr$ld),
        pecotmr:::.ldRefVariantKeys(pr$panel)
    )
    expect_null(pecotmr:::.ldRefHandle(pr$ld))
    expect_identical(pecotmr:::.asLdRef(pr$ld), pr$ld)
    expect_true(pecotmr:::.ldIsReference(pr$ld))
    expect_false(pecotmr:::.ldIsPanel(pr$ld))
    expect_match(pecotmr:::.ldRefLabel(pr$ld), "LdMatrix of 349 variants")
    sub <- pecotmr:::.ldRefSubset(pr$ld, c(TRUE, FALSE, TRUE, rep(FALSE, 346)))
    expect_equal(length(sub), 2L)
    expect_equal(length(pecotmr:::.emptyLdRef(pr$ld)), 0L)
    expect_error(pecotmr:::.asLdRef(1:3), "or an LdMatrix")
})

test_that("an LdMatrix's sample size is what it records", {
    pr <- .ldrefPair(nSamples = 80L)
    expect_identical(pecotmr:::.ldRefNSamples(pr$ld), 80L)
    expect_identical(pecotmr:::.ldRefRealNSamples(pr$ld), 80L)
    expect_identical(pecotmr:::.ldRefNoiseN(pr$ld), 80L)
    unknown <- LdMatrix(ldMatrix(pr$ld))
    expect_true(is.na(pecotmr:::.ldRefRealNSamples(unknown)))
    expect_error(
        pecotmr:::.ldRefRequireNSamples(unknown, "DENTIST"),
        "DENTIST needs the LD reference's sample size"
    )
    expect_equal(pecotmr:::.ldRefNoiseN(pr$panel), 100L)
})

test_that("a stochastic sketch records its real size and frequencies", {
    skip_if_not_installed("pgenlibr")
    sketch <- .ldrefSketch()
    expect_true(pecotmr:::.ldRefIsProjection(sketch))
    # 100 pseudo-samples, and an .afreq that carries U_MIN / U_MAX (so the
    # panel is a stochastic sketch) with OBS_CT 100: 50 real people. A MAC
    # cutoff and DENTIST's truncation rank count those, not the columns.
    expect_equal(pecotmr:::.ldRefNSamples(sketch), 100L)
    expect_identical(pecotmr:::.ldRefRealNSamples(sketch), 50L)
    # 1 / (1/50 + 1/100): the noise of the real samples and of the
    # projection's 100 pseudo-samples together.
    expect_identical(pecotmr:::.ldRefNoiseN(sketch), 33L)
    af <- readAfreq(
        file.path(geno_test_data_dir, "test_harmonize_regions")
    )
    expect_equal(
        S4Vectors::mcols(SummarizedExperiment::rowRanges(sketch))$AF,
        af$alt_freq[match(pecotmr:::.ldRefVariantIds(sketch), af$id)]
    )
    # Subsetting keeps the record with the variants.
    sub <- sketch[3:5, ]
    expect_true(pecotmr:::.ldRefIsProjection(sub))
    expect_equal(pecotmr:::.ldRefAf(sub), pecotmr:::.ldRefAf(sketch)[3:5])
    # A panel of real genotypes records nothing.
    real <- readGenotypes(file.path(geno_test_data_dir, "test_variants.pgen"))
    expect_false(pecotmr:::.ldRefIsProjection(real))
    expect_null(pecotmr:::.ldRefAf(real))
})

test_that("frequencies come from the record, never from projected dosage", {
    skip_if_not_installed("pgenlibr")
    sketch <- .ldrefSketch()
    st <- pecotmr:::.ldRefVariantStats(sketch, 1:3, FALSE, "t")
    af <- pecotmr:::.ldRefAf(sketch)[1:3]
    expect_equal(st$maf, pmin(af, 1 - af))
    expect_null(st$missRate)
    expect_error(
        pecotmr:::.ldRefVariantStats(sketch, 1:3, TRUE, "The LD-panel filter"),
        "needs genotype missingness.*stochastic sketch"
    )
    expect_equal(
        pecotmr:::.ldRefVariance(sketch, 1:3),
        2 * af * (1 - af)
    )
})

test_that("an LdMatrix answers frequency questions only from its AF", {
    pr <- .ldrefPair()
    m <- pecotmr:::.ldFromRefMatch(
        pr$ld,
        pecotmr:::.ldRefMatchIds(pr$ld)[1:5],
        "t",
        "drop"
    )
    cutoffs <- list(mafCutoff = 0.1, macCutoff = 0, imissCutoff = 1)
    expect_error(
        pecotmr:::.panelDropMask(pr$ld, m, cutoffs),
        "needs allele frequencies, and the LD reference records none"
    )
    expect_true(all(is.na(pecotmr:::.ldRefVariance(pr$ld, 1:3))))
    expect_error(
        pecotmr:::.panelDropMask(
            pr$ld,
            m,
            list(mafCutoff = 0, macCutoff = 0, imissCutoff = 0.1)
        ),
        "needs genotype missingness.*an LdMatrix"
    )

    gr <- SummarizedExperiment::rowRanges(pr$panel)
    S4Vectors::mcols(gr)$AF <- rep(c(0.02, 0.3), length.out = length(gr))
    withAf <- LdMatrix(ldMatrix(pr$ld), variants = gr, nSamples = 100)
    expect_equal(
        pecotmr:::.panelDropMask(withAf, m, cutoffs),
        c(TRUE, FALSE, TRUE, FALSE, TRUE)
    )
    # MAC 20 over 100 samples is MAF 0.1.
    expect_equal(
        pecotmr:::.panelDropMask(
            withAf,
            m,
            list(mafCutoff = 0, macCutoff = 20, imissCutoff = 1)
        ),
        c(TRUE, FALSE, TRUE, FALSE, TRUE)
    )
    mix <- LdMatrix(list(withAf, withAf), mixtureWeights = c(0.5, 0.5))
    expect_error(
        pecotmr:::.panelDropMask(
            mix,
            m,
            list(mafCutoff = 0, macCutoff = 20, imissCutoff = 1)
        ),
        "LdMatrix is a mixture"
    )
    expect_equal(
        pecotmr:::.ldRefVariance(withAf, 1:2),
        2 * c(0.02, 0.3) * (1 - c(0.02, 0.3))
    )
})

test_that("the same-reference check compares LdMatrix objects by their LD", {
    pr <- .ldrefPair()
    a <- pr$ld[1:200]
    b <- pr$ld[150:349]
    expect_no_error(
        pecotmr:::.requireMatchingLdRefs(a, b, "pipe")
    )
    other <- LdMatrix(
        diag(200),
        variants = SummarizedExperiment::rowRanges(pr$panel)[150:349]
    )
    expect_error(
        pecotmr:::.requireMatchingLdRefs(a, other, "pipe"),
        "disagree on the LD among the variants they share"
    )
    expect_error(
        pecotmr:::.requireMatchingLdRefs(a, pr$panel, "pipe"),
        "one side's ldReference is an LdMatrix and the other's a genotype panel"
    )
})

test_that("the LD's eigenvalues come from its factor's thin SVD", {
    set.seed(4)
    X <- matrix(rbinom(30 * 50, 2, 0.3), 30, 50)
    F <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(F) / (nrow(F) - 1)
    dec <- pecotmr:::.ldFactorEigen(F)
    # The SVD lists min(m, p) directions; the rest of R's are exactly zero.
    expect_equal(dec$nNull, 50L - nrow(F))
    expect_equal(length(dec$values) + dec$nNull, 50L)
    expect_equal(dec$values[1:25], eigen(R, symmetric = TRUE)$values[1:25])
    expect_equal(dec$vectors %*% (dec$values * t(dec$vectors)), R)
})

test_that("R o R products read the same from a matrix, factor or eigenpairs", {
    set.seed(5)
    X <- matrix(rbinom(40 * 600, 2, 0.4), 40, 600)
    F <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(F) / (nrow(F) - 1)
    e <- eigen(R, symmetric = TRUE)
    eig <- list(vectors = e$vectors, values = e$values)
    A <- cbind(1, runif(600), rbinom(600, 1, 0.3))
    want <- (R^2) %*% A
    # 600 variants are read in three stripes, the last one partial.
    expect_equal(pecotmr:::.ldSquaredProduct(list(F = F), A), want)
    expect_equal(pecotmr:::.ldSquaredProduct(list(R = R), A), want)
    expect_equal(pecotmr:::.ldSquaredProduct(eig, A), want)
    expect_equal(pecotmr:::.ldSourceRows(list(F = F), 3:5), R[3:5, ])
    expect_equal(pecotmr:::.ldSourceRows(eig, 3:5), R[3:5, ])
    expect_equal(pecotmr:::.ldSourceSize(list(F = F)), 600L)
    expect_equal(pecotmr:::.ldSourceSize(eig), 600L)
})

test_that("a sketch's .afreq sidecar is read beside a file of any format", {
    # The sidecar is named after the file's stem, whatever the format; a
    # stochastic one (U_MIN / U_MAX, OBS_CT) records the real sample count
    # and A1 frequencies on the panel.
    expect_equal(
        pecotmr:::.genotypeStem(c(
            "a/x.vcf.gz",
            "x.VCF.BGZ",
            "x.bcf",
            "x.gds",
            "x.bgen",
            "plink/stem"
        )),
        c("a/x", "x", "x", "x", "x", "plink/stem")
    )
    skip_if_not_installed("RSQLite")
    fx <- .bgw_fixture(n = 6L, nVar = 5L)
    ids <- c("chr1:100:A:G", sprintf("rs%d", 1:4))
    freq <- c(0.1, 0.2, 0.3, 0.4, 0.5)
    afreq <- tibble::tibble(
        `#CHROM` = "1",
        ID = ids,
        REF = "A",
        ALT = "G",
        ALT_FREQS = freq,
        OBS_CT = c(200, 200, 198, 200, 200),
        U_MIN = -1,
        U_MAX = 1
    )
    readr::write_tsv(afreq, str_c(str_remove(fx$path, "\\.bgen$"), ".afreq"))
    panel <- readGenotypes(fx$path)
    expect_equal(pecotmr:::.ldRefProjection(panel)$nSamples, 100L)
    expect_equal(SummarizedExperiment::rowRanges(panel)$AF, freq)
    # The same sidecar also answers the panel filter's MAF.
    expect_equal(
        pecotmr:::.panelAfreqMaf(
            pecotmr:::.ldRefSource(panel),
            pecotmr:::.ldRefChrom(panel),
            ids
        ),
        pmin(freq, 1 - freq)
    )
})

test_that("a panel held in memory serves as an LD reference", {
    # Any RSE with a dosage assay and SNP / A1 / A2 rowRanges is a panel: it
    # reads through no handle and records no file, and every LD route reads
    # it as it reads the same genotypes from disk.
    disk <- readGenotypes(test_path("test_data", "test_variants.bed"))[1:20, ]
    mem <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            dosage = as.matrix(SummarizedExperiment::assay(disk, "dosage"))
        ),
        rowRanges = SummarizedExperiment::rowRanges(disk)
    )
    expect_null(pecotmr:::.ldRefHandle(mem))
    expect_null(pecotmr:::.ldRefSource(mem))
    expect_match(pecotmr:::.ldRefLabel(mem), "held in memory")
    expect_equal(computeLd(mem), computeLd(disk))
    expect_equal(
        pecotmr:::.ldRefFactor(mem, 1:20),
        pecotmr:::.ldRefFactor(disk, 1:20)
    )
    expect_equal(dim(pecotmr:::.ldRefSubset(mem, 3:7)), c(5L, ncol(mem)))
    expect_equal(dim(pecotmr:::.emptyLdRef(mem)), c(0L, 0L))
})

# ===========================================================================
# Per-row LD references
# ===========================================================================

# An LdMatrix over five chr1 variants at 100 bp spacing.
.ldr_ldm <- function(rho = 0.2, chrom = "chr1") {
    ids <- sprintf("%s:%d:A:G", chrom, 100L * (1:5))
    r <- rho^abs(outer(1:5, 1:5, "-"))
    LdMatrix(`dimnames<-`(r, list(ids, ids)))
}

# A QtlSumStats of one tuple per reference in `refs`, every tuple over the same
# five variants on `chrom`.
.ldr_qss <- function(refs, chrom = "chr1") {
    n <- length(refs)
    gr <- GenomicRanges::GRanges(
        chrom,
        IRanges::IRanges(100L * (1:5), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("%s:%d:A:G", chrom, 100L * (1:5)),
        A1 = "A",
        A2 = "G",
        Z = seq_len(5) / 10,
        N = 100L
    )
    QtlSumStats(
        studyName = rep("s1", n),
        context = str_c("c", seq_len(n)),
        trait = rep("t1", n),
        entry = rep(list(gr), n),
        genome = "hg19",
        ldReference = refs
    )
}

test_that(".ldRefPerTuple reads one reference or one per tuple", {
    a <- .ldr_ldm()
    b <- .ldr_ldm(rho = 0.5)
    expect_identical(pecotmr:::.ldRefPerTuple(a, 3L), list(a, a, a))
    # What ldReference() returns reads the same way: a one-entry list is that
    # reference for every tuple, an empty one is LD-free.
    expect_identical(pecotmr:::.ldRefPerTuple(list(ld1 = a), 2L), list(a, a))
    expect_identical(pecotmr:::.ldRefPerTuple(list(), 2L), list(NULL, NULL))
    expect_identical(pecotmr:::.ldRefPerTuple(NULL, 2L), list(NULL, NULL))
    expect_identical(pecotmr:::.ldRefPerTuple(list(a, b), 2L), list(a, b))
    expect_error(pecotmr:::.ldRefPerTuple(list(a, b), 3L), "one per tuple")
})

test_that(".ldRefDict stores each distinct reference once and keys every row", {
    a <- .ldr_ldm()
    b <- .ldr_ldm(rho = 0.5)
    # An equal copy built separately is the same reference, not a second one.
    dict <- pecotmr:::.ldRefDict(list(a, NULL, b, a, .ldr_ldm()))
    expect_identical(dict$refs, list(ld1 = a, ld2 = b))
    expect_identical(dict$keys, c("ld1", NA, "ld2", "ld1", "ld1"))
    expect_identical(pecotmr:::.ldRefDict(list(NULL))$refs, list())
})

test_that("a collection keeps each row on its own reference, stored once", {
    a <- .ldr_ldm()
    b <- .ldr_ldm(rho = 0.5)
    qss <- .ldr_qss(list(a, b, a))
    expect_identical(ldReference(qss), list(ld1 = a, ld2 = b))
    expect_identical(qss$ldReference, c("ld1", "ld2", "ld1"))
    expect_identical(pecotmr:::.rowLdRefs(qss), list(a, b, a))
    # A subset reports only the references its rows still use.
    expect_identical(ldReference(qss[c(1L, 3L)]), list(ld1 = a))
    expect_output(show(qss), "2 references")
})

test_that("validity refuses a row keyed to a reference the slot lacks", {
    qss <- .ldr_qss(list(.ldr_ldm()))
    expect_error(
        {
            qss$ldReference <- "ld9"
            methods::validObject(qss)
        },
        "not in ldReference: ld9"
    )
})

test_that("an object saved before per-row references is refused", {
    # Read as a list, its single reference would pass for LD-free.
    old <- .ldr_qss(list(.ldr_ldm()))
    attr(old, "ldReference") <- .ldr_ldm()
    md <- attr(old, "elementMetadata")
    md[["ldReference"]] <- NULL
    attr(old, "elementMetadata") <- md
    expect_error(ldReference(old), "saved before LD references became")
})

test_that(".withRowLdRefs puts every row on one reference", {
    a <- .ldr_ldm()
    b <- .ldr_ldm(rho = 0.5)
    moved <- pecotmr:::.withRowLdRefs(.ldr_qss(list(a, b)), b)
    expect_identical(ldReference(moved), list(ld1 = b))
    expect_identical(moved$ldReference, c("ld1", "ld1"))
    expect_true(methods::validObject(moved))
})

test_that(".requireMatchingRowLdRefs checks only the pairs it is given", {
    a <- .ldr_ldm()
    qtl <- .ldr_qss(list(a, .ldr_ldm(rho = 0.5)))
    gwas <- .ldr_qss(list(a))
    # QTL row 1 shares the GWAS row's reference; row 2 does not, but is only
    # refused when it is paired with it.
    expect_null(pecotmr:::.requireMatchingRowLdRefs(qtl, gwas, 1L, 1L, "p"))
    expect_error(
        pecotmr:::.requireMatchingRowLdRefs(qtl, gwas, 2L, 1L, "p"),
        "disagree on the LD"
    )
})

test_that(".overlappingRowPairs pairs rows whose regions overlap", {
    a <- .ldr_ldm()
    qtl <- combineQtlSumStats(
        .ldr_qss(list(a)),
        .ldr_qss(list(.ldr_ldm(chrom = "chr2")), chrom = "chr2")
    )
    pairs <- pecotmr:::.overlappingRowPairs(qtl, .ldr_qss(list(a)))
    expect_identical(pairs, list(qtlRows = 1L, gwasRows = 1L))
    # Collections on no common chromosome have no pairs, and say nothing.
    expect_silent(
        none <- pecotmr:::.overlappingRowPairs(
            qtl,
            .ldr_qss(list(.ldr_ldm(chrom = "chr3")), chrom = "chr3")
        )
    )
    expect_identical(none, list(qtlRows = integer(0), gwasRows = integer(0)))
})

test_that(".overlapsOnShared compares ranges on the chromosomes both have", {
    q <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(100L, 100L), width = 10L)
    )
    s <- GenomicRanges::GRanges("chr2", IRanges::IRanges(105L, width = 1L))
    expect_identical(pecotmr:::.overlapsOnShared(q, s), c(FALSE, TRUE))
    # With no chromosome in common nothing overlaps, and nothing is said.
    far <- GenomicRanges::GRanges("chr3", IRanges::IRanges(100L, width = 1L))
    expect_silent(none <- pecotmr:::.overlapsOnShared(q, far))
    expect_identical(none, c(FALSE, FALSE))
    # Strand says nothing about where a variant or a block sits.
    minus <- GenomicRanges::GRanges(
        "chr2",
        IRanges::IRanges(105L, width = 1L),
        strand = "-"
    )
    plus <- GenomicRanges::`strand<-`(q, value = "+")
    expect_identical(pecotmr:::.overlapsOnShared(plus, minus), c(FALSE, TRUE))
})

test_that(".locusRanges keeps where ranges sit and the build, nothing else", {
    gr <- GenomicRanges::GRanges(
        c("chr1", "chrM"),
        IRanges::IRanges(c(100L, 5L), width = 1L),
        strand = c("+", "-"),
        SNP = c("a", "b")
    )
    GenomeInfoDb::genome(gr) <- "hg19"
    out <- pecotmr:::.locusRanges(gr)
    expect_identical(
        as.character(GenomicRanges::seqnames(out)),
        c("chr1", "chrM")
    )
    expect_identical(as.character(GenomicRanges::strand(out)), c("*", "*"))
    expect_identical(ncol(S4Vectors::mcols(out)), 0L)
    expect_identical(unname(GenomeInfoDb::genome(out)), c("hg19", "hg19"))
})

test_that("ranges on two genome builds refuse to meet", {
    q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, width = 10L))
    s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(105L, width = 1L))
    GenomeInfoDb::genome(q) <- "hg19"
    GenomeInfoDb::genome(s) <- "hg38"
    expect_error(pecotmr:::.overlapsOnShared(q, s), "incompatible genomes")
})

test_that(".rowsSharedLdRef returns the rows' one reference or refuses", {
    a <- .ldr_ldm()
    qss <- .ldr_qss(list(a, .ldr_ldm(rho = 0.5), a))
    expect_identical(pecotmr:::.rowsSharedLdRef(qss, c(1L, 3L), "rows"), a)
    expect_error(
        pecotmr:::.rowsSharedLdRef(qss, 1:2, "these rows"),
        "these rows use 2 different LD references"
    )
})

# ===========================================================================
# On a real genotype panel: checkLd, ldClumpByScore, enforceDesignFullRank,
# filterVariantsByLdReference
# ===========================================================================

# Fixtures: 100 samples x 349 variants on chr21
plink_prefix <- test_path("test_data", "test_variants")

# The panel's genotype matrix and variant info
loadTestGenotype <- function() {
    loadGenotypeRegion(plink_prefix, returnVariantInfo = TRUE)
}

# --- checkLd ----------------------------------------------------------------

test_that("checkLd reports a rank-deficient sample LD matrix as PSD", {
    skip_if_not_installed("pgenlibr")
    X <- loadTestGenotype()$X
    # 100 samples cap the rank at 99 of 349 variants.
    R <- computeLd(X, method = "sample")
    result <- checkLd(R)
    expect_true(is.list(result))
    expect_false(result$isPd)
    expect_true(result$isPsd)
    expect_equal(result$methodApplied, "none")
    # min eigenvalue may be near-zero (numerically PSD, not strictly PD)
    expect_true(result$minEigenvalue > -1e-7)
    expect_true(result$nNegative == 0)
    expect_true(is.finite(result$conditionNumber))
})

# --- ldClumpByScore -------------------------------------------------------

test_that("ldClumpByScore takes UCSC chromosome labels on a real panel", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("bigstatsr")
    geno <- loadTestGenotype()
    set.seed(42)
    score <- runif(ncol(geno$X))
    # UCSC labels (chr21) are taken as they are.
    chr <- geno$variant_info$chrom
    pos <- geno$variant_info$pos
    keep <- ldClumpByScore(
        geno$X,
        score = score,
        chr = chr,
        pos = pos,
        r2 = 0.2
    )
    expect_true(is.integer(keep))
    expect_true(length(keep) > 0)
    expect_true(length(keep) <= ncol(geno$X))
    expect_true(all(keep %in% seq_len(ncol(geno$X))))
})

# --- enforceDesignFullRank ------------------------------------------------

test_that("enforceDesignFullRank handles genotype matrix with covariates", {
    skip_if_not_installed("pgenlibr")
    X <- loadTestGenotype()$X
    # Create a simple covariate matrix (e.g., first 2 PCs of X)
    pca <- prcomp(X, rank. = 2)
    C <- pca$x
    result <- enforceDesignFullRank(X, C, strategy = "correlation")
    expect_true(is.matrix(result))
    expect_equal(nrow(result), nrow(X))
    # Should produce full-rank design
    full_design <- cbind(1, result, C)
    expect_equal(qr(full_design)$rank, ncol(full_design))
})

# --- filterVariantsByLdReference -----------------------------------------

test_that("filterVariantsByLdReference filters via a PLINK metadata file", {
    skip_if_not_installed("pgenlibr")
    geno <- loadTestGenotype()
    vi <- geno$variant_info
    variant_ids <- paste0(vi$chrom, ":", vi$pos, ":", vi$A2, ":", vi$A1)
    fake_ids <- c("21:999999:A:G", "21:888888:C:T")
    all_ids <- c(variant_ids, fake_ids)

    # Create a metadata TSV in the same directory as the PLINK files
    # so the relative path resolves correctly
    meta_file <- file.path(test_data_dir, "ld_meta_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(
        paste("chrom", "start", "end", "path", sep = "\t"),
        meta_file
    )
    cat(
        paste("21", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )

    result <- suppressMessages(
        filterVariantsByLdReference(all_ids, meta_file, keepIndel = TRUE)
    )
    expect_true(is.list(result))
    expect_true(length(result$data) <= length(all_ids))
    expect_true(length(result$idx) == length(result$data))
    # Fake variants should be filtered out
    expect_true(length(result$data) <= length(variant_ids))
})
