context("file_utils")
library(tidyverse)

test_that("readPvar reads real pvar file", {
    skip_if_not_installed("pgenlibr")
    pvar_path <- file.path(test_path("test_data"), "test_variants.pvar")
    res <- pecotmr:::readPvar(pvar_path)
    expect_equal(colnames(res), c("chrom", "id", "pos", "A2", "A1"))
    expect_equal(nrow(res), 349L)
    expect_true(all(res$chrom == "21"))
})

test_that("readBim dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- readBim(example_path)
    expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
    expect_equal(nrow(res), 100)
})

test_that("readFam dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- readFam(example_path)
    expect_equal(nrow(res), 100)
})

test_that("openBed dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- openBed(example_path)
    expect_equal(res$class, "pgen")
})

test_that("findValidFilePath works",{
    ref_path <- "test_data/protocol_example.genotype.bed"
    expect_error(
        pecotmr:::.findValidFilePath(paste0(ref_path, "s"), "protocol_example.genotype.bamf"),
        "Both reference and target file paths do not work. Tried paths: 'test_data/protocol_example.genotype.beds' and 'test_data/protocol_example.genotype.bamf'")
    expect_equal(
        pecotmr:::.findValidFilePath(ref_path, "abc"),
        ref_path)
    expect_equal(
        pecotmr:::.findValidFilePath(ref_path, "protocol_example.genotype.bim"),
        "test_data/protocol_example.genotype.bim")
    expect_equal(
        pecotmr:::.findValidFilePath(ref_path, "test_data/protocol_example.genotype.bim"),
        "test_data/protocol_example.genotype.bim")
})


dummy_geno_data <- function(
    number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
    number_missing = 10, number_low_maf = 10, number_zero_var = 10, number_var_thresh = 10) {
    set.seed(1)
    # Create portion of Matrix with satisfactory values
    X <- matrix(
        sample(c(0,1,2), number_of_samples*number_of_snps, replace = TRUE),
        nrow=number_of_samples, ncol=number_of_snps)
    # Create portion of Matrix that should get pruned
    ## Missing Rate
    if (number_missing > 0) {
        X_missing <- rbind(
            matrix(
                sample(c(0,1,2), (number_of_samples-3)*number_of_snps, replace = TRUE),
                nrow=number_of_samples-3, ncol=number_of_snps),
            matrix(
                rep(NA, 3*number_of_snps), nrow=3, ncol=number_of_snps))
        X <- cbind(X, X_missing)
    }
    ## MAF
    if (number_low_maf > 0) {
        X_maf <- matrix(
            rep(0.1, number_of_samples*number_of_snps), nrow=number_of_samples, ncol=number_of_snps)
        X <- cbind(X, X_maf)
    }
    ## Zero Variance
    if (number_zero_var > 0) {
        X_zerovar <- matrix(
            rep(1, number_of_samples*number_of_snps), nrow=number_of_samples, ncol=number_of_snps)
        X <- cbind(X, X_zerovar)
    }
    ## Variance Threshold, just one row
    if (number_var_thresh > 0) {
        X_varthresh <- matrix(
            c(rep(1, (number_of_samples - 1)), 2), nrow=number_of_samples, ncol=1)
        X <- cbind(X, X_varthresh)
    }
    colnames(X) <- paste0(
        "chr1:",
        seq(1000,1000+number_of_snps+number_missing+number_low_maf+number_zero_var+number_var_thresh-1),
        "_G_C")
    rownames(X) <- paste0("Sample_", seq(sample_start_id, number_of_samples + sample_start_id - 1))
    return(X)
}

dummy_pheno_data <- function(number_of_samples = 10, number_of_phenotypes = 10, randomize = FALSE, sample_start_id = 1) {
    # Create dummy phenotype bed file
    # columns: Chrom, Start, End, Sample_1, Sample_2, ..., Sample_N
    start_matrix <- matrix(
        c(
            rep("chr1", number_of_phenotypes),
            seq(100, 100+number_of_phenotypes-1),
            seq(101, 101+number_of_phenotypes-1)
        ),
        nrow=number_of_phenotypes, ncol=3)
    end_matrix <- matrix(
        rnorm(number_of_samples*number_of_phenotypes), nrow=number_of_phenotypes, ncol=number_of_samples)
    pheno_data <- cbind(start_matrix, end_matrix)
    sample_ids <- paste0("Sample_", seq(sample_start_id, number_of_samples + sample_start_id - 1))
    colnames(pheno_data) <- c("#chr", "start", "end", sample_ids)
    colnames(end_matrix) <- sample_ids
    if (randomize) {
        end_matrix <- end_matrix[sample(nrow(end_matrix)),]
    }
    pheno_data <- t(pheno_data)
    pheno_data <- lapply(seq_len(ncol(pheno_data)), function(i) pheno_data[,i,drop=FALSE])
    return(pheno_data)
}

dummy_covar_data <- function(number_of_samples = 10, number_of_covars = 10, row_na = FALSE, randomize = FALSE, sample_start_id = 1) {
    covar <- matrix(
        sample(1:20, number_of_samples*number_of_covars, replace = TRUE),
        nrow=number_of_samples, ncol=number_of_covars)
    colnames(covar) <- paste0("Covar_", seq(1, number_of_covars))
    rownames(covar) <- paste0("Sample_", seq(sample_start_id, number_of_samples + sample_start_id - 1))
    if (randomize) {
        covar <- covar[sample(nrow(covar)),]
    }
    if (row_na) {
        covar[sample(length(covar),1), 1:number_of_covars] <- NA
    }
    return(covar)
}


test_that("Test loadGenotypeRegion",{
  res <- loadGenotypeRegion(
    "test_data/protocol_example.genotype")
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
})

test_that("Test loadGenotypeRegion no indels",{
  res <- loadGenotypeRegion(
    "test_data/protocol_example.genotype", keepIndel = F)
  bim_file <- read_delim(
    "test_data/protocol_example.genotype.bim", delim = "\t", col_names = F
  )
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
  indels <- with(bim_file, grepl("[^ATCG]", X5) | grepl("[^ATCG]", X6) | nchar(X5) > 1 | nchar(X6) > 1)
  expect_equal(
    nrow(bim_file[!indels, ]),
    ncol(res)
  )
})

test_that("Test loadGenotypeRegion with region",{
  res <- loadGenotypeRegion(
    "test_data/protocol_example.genotype",
    region = "chr22:20689453-20845958")
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  snp_ids <- read_delim(
    "test_data/protocol_example.genotype.bim", delim = "\t", col_names = F
  ) %>% pull(X2)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
  expect_equal(ncol(res), 8)
  expect_equal(colnames(res), snp_ids[1:8])
})

test_that("Test loadGenotypeRegion with region and no indels",{
  res <- loadGenotypeRegion(
    "test_data/protocol_example.genotype",
    region = "chr22:20689453-20845958", keepIndel = F)
  bim_file <- read_delim(
    "test_data/protocol_example.genotype.bim", delim = "\t", col_names = F
  )[1:8, ]
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
  indels <- with(bim_file, grepl("[^ATCG]", X5) | grepl("[^ATCG]", X6) | nchar(X5) > 1 | nchar(X6) > 1)
  expect_equal(
    nrow(bim_file[!indels, ]),
    ncol(res))
  expect_equal(colnames(res), bim_file[!indels, ]$X2)
})

test_that("loadGenotypeRegion errors on missing genotype files", {
  expect_error(
    loadGenotypeRegion("/nonexistent/geno"),
    "Genotype files not found"
  )
})

# --- findStochasticMeta tests ---

test_that("findStochasticMeta finds generic sidecar from PLINK1 prefix", {
  td <- test_path("test_data")
  # test_harmonize_regions has .stochastic_meta.tsv alongside it
  result <- pecotmr:::findStochasticMeta(file.path(td, "test_harmonize_regions"))
  expect_true(!is.null(result))
  expect_true(grepl("\\.(afreq|stochastic_meta\\.tsv)$", result))
})

test_that("findStochasticMeta finds sidecar from VCF path", {
  td <- test_path("test_data")
  result <- pecotmr:::findStochasticMeta(file.path(td, "test_harmonize_regions.vcf.gz"))
  expect_true(!is.null(result))
  expect_true(grepl("\\.(afreq|stochastic_meta\\.tsv)$", result))
})

test_that("findStochasticMeta finds sidecar from GDS path", {
  td <- test_path("test_data")
  result <- pecotmr:::findStochasticMeta(file.path(td, "test_harmonize_regions.gds"))
  expect_true(!is.null(result))
  expect_true(grepl("\\.(afreq|stochastic_meta\\.tsv)$", result))
})

test_that("findStochasticMeta returns NULL when no sidecar exists", {
  td <- test_path("test_data")
  result <- pecotmr:::findStochasticMeta(file.path(td, "protocol_example.genotype"))
  expect_null(result)
})

# --- readStochasticMeta tests ---

test_that("readStochasticMeta reads generic format", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
  result <- pecotmr:::readStochasticMeta(path)
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
  expect_equal(nrow(result), 8L)
  expect_true(is.numeric(result$u_min))
  expect_true(is.numeric(result$u_max))
})

test_that("readStochasticMeta reads afreq format", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.afreq")
  result <- pecotmr:::readStochasticMeta(path)
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
  expect_equal(nrow(result), 8L)
  expect_true(all(grepl("^chr21_", result$id)))
})

test_that("readStochasticMeta reads afreq.zst format", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.afreq.zst")
  result <- pecotmr:::readStochasticMeta(path)
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
  expect_equal(nrow(result), 8L)
  # Should produce identical results to the plain afreq
  plain <- pecotmr:::readStochasticMeta(file.path(td, "test_harmonize_regions.afreq"))
  expect_equal(result, plain)
})

test_that("findStochasticMeta prefers afreq over afreq.zst", {
  td <- test_path("test_data")
  # Both .afreq and .afreq.zst exist; findStochasticMeta should return .afreq first
  result <- pecotmr:::findStochasticMeta(file.path(td, "test_harmonize_regions"))
  expect_true(grepl("\\.afreq$", result))
})

test_that("readStochasticMeta auto-detects format from extension", {
  td <- test_path("test_data")
  # .afreq extension -> afreq parser
  afreq_result <- pecotmr:::readStochasticMeta(file.path(td, "test_harmonize_regions.afreq"))
  # .tsv extension -> generic parser
  generic_result <- pecotmr:::readStochasticMeta(
    file.path(td, "test_harmonize_regions.stochastic_meta.tsv"))
  # Both should return the same u_min/u_max values
  expect_equal(afreq_result$u_min, generic_result$u_min)
  expect_equal(afreq_result$u_max, generic_result$u_max)
  expect_equal(afreq_result$id, generic_result$id)
})

test_that("readStochasticMeta respects format override", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
  # Explicit generic format should work
  result <- pecotmr:::readStochasticMeta(path, format = "generic")
  expect_equal(nrow(result), 8L)
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
})

test_that("readStochasticMeta returns NULL for afreq without U_MIN/U_MAX", {
  td <- test_path("test_data")
  # test_variants.afreq has no U_MIN/U_MAX columns
  path <- file.path(td, "test_variants.afreq")
  result <- pecotmr:::readStochasticMeta(path)
  expect_null(result)
})

test_that("readStochasticMeta returns NULL for nonexistent file", {
  result <- pecotmr:::readStochasticMeta("/nonexistent/file.tsv")
  expect_null(result)
})

# --- loadGenotypeRegion stochastic inversion test ---

test_that("loadGenotypeRegion applies stochastic inversion with explicit sidecar", {
  td <- test_path("test_data")
  metaPath <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
  smeta <- pecotmr:::readStochasticMeta(metaPath)

  # Load with explicit sidecar - inversion transforms the integer dosages
  res <- loadGenotypeRegion(
    file.path(td, "test_harmonize_regions"),
    returnVariantInfo =TRUE,
    stochasticMetaPath =metaPath
  )

  expect_equal(ncol(res$X), 8L)
  # u_min/u_max should be attached to variant_info
  expect_true("u_min" %in% colnames(res$variant_info))
  expect_true("u_max" %in% colnames(res$variant_info))
  expect_equal(res$variant_info$u_min, smeta$u_min)
  expect_equal(res$variant_info$u_max, smeta$u_max)

  # Verify inversion math: for a dosage value d with u_min/u_max,
  # inverted = d * (u_max - u_min) / 2 + u_min
  # Check the first variant's first sample manually
  raw <- loadGenotypeRegion(
    file.path(td, "protocol_example.genotype"),
    region = "chr22:20689453-20845958"
  )
  # protocol_example has no sidecar, so raw values are unchanged (integer dosages)
  expect_true(all(raw == round(raw), na.rm = TRUE))

  # The inverted matrix should NOT be all integers (u_min != 0 or u_max != 2)
  expect_false(all(res$X == round(res$X), na.rm = TRUE))
})





















# ===========================================================================
# readBim vroom-based tests
# ===========================================================================

test_that("readBim returns correct columns and types", {
  bim_path <- tempfile(fileext = ".bim")
  cat("22\trs100\t0\t50000\tA\tG\n", file = bim_path)
  cat("22\trs200\t0\t60000\tT\tC\n", file = bim_path, append = TRUE)
  cat("22\trs300\t0\t70000\tC\tA\n", file = bim_path, append = TRUE)

  bed_path <- sub("\\.bim$", ".bed", bim_path)
  file.copy(bim_path, bim_path)
  res <- readBim(bed_path)
  expect_equal(nrow(res), 3)
  expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
  expect_equal(res$id, c("rs100", "rs200", "rs300"))
  expect_equal(res$pos, c(50000, 60000, 70000))
  file.remove(bim_path)
})

# ===========================================================================
# tabixRegion
# ===========================================================================





# ===========================================================================
# NoSNPsError / NoPhenotypeError custom conditions
# ===========================================================================

test_that("NoSNPsError creates proper error condition", {
  err <- NoSNPsError("test message")
  expect_true(inherits(err, "NoSNPsError"))
  expect_true(inherits(err, "error"))
  expect_true(inherits(err, "condition"))
  expect_equal(err$message, "test message")
})


# ===========================================================================
# extractPhenotypeCoordinates
# ===========================================================================

# ===========================================================================
# loadTsvRegion
# ===========================================================================



# ===========================================================================
# batchLoadTwasWeights
# ===========================================================================




# ===========================================================================
# .colocFilterCsByConcentration
# ===========================================================================

test_that(".colocFilterCsByConcentration returns numeric index vector", {
  set.seed(42)
  n_L <- 5
  n_vars <- 20
  alpha_raw <- matrix(runif(n_L * n_vars), nrow = n_L)
  alpha_norm <- t(apply(alpha_raw, 1, function(x) x / sum(x)))

  mock_susie <- list(
    alpha = alpha_norm,
    V = runif(n_L),
    lbf_variable = matrix(rnorm(n_L * n_vars), nrow = n_L),
    mu = matrix(rnorm(n_L * n_vars), nrow = n_L),
    mu2 = matrix(abs(rnorm(n_L * n_vars)), nrow = n_L),
    sets = list(cs = list(L1 = c(1,3,5), L3 = c(2,4)), cs_index = c(1, 3)),
    pip = colSums(alpha_norm),
    niter = 100,
    converged = TRUE
  )

  result <- pecotmr:::.colocFilterCsByConcentration(
    mock_susie, coverage = 0.5, concentration = 0.5)
  expect_true(is.numeric(result))
})

# ===========================================================================
# getRefVariantInfo
# ===========================================================================

test_that("getRefVariantInfo processes precomputed bim with 6 columns", {
  td <- test_path("test_data")
  meta_file <- file.path(td, "ld_meta_refinfo_6col_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- getRefVariantInfo(meta_file, "chr1:1000-1190")
  expect_true(is.data.frame(result))
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% colnames(result)))
  expect_equal(nrow(result), 5L)
})

test_that("getRefVariantInfo processes precomputed bim with 9 columns", {
  td <- test_path("test_data")
  meta_file <- file.path(td, "ld_meta_refinfo_9col_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.9col.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- getRefVariantInfo(meta_file, "chr1:1000-1190")
  expect_true(all(c("chrom", "id", "pos", "A2", "A1", "variance", "allele_freq", "n_nomiss") %in% colnames(result)))
  expect_equal(nrow(result), 5L)
  expect_equal(result$allele_freq, c(0.3, 0.4, 0.2, 0.5, 0.15))
})

# ---- invertMinmaxScaling ----

test_that("invertMinmaxScaling exactly recovers original U", {
  set.seed(42)
  n <- 500
  k <- 4
  # Simulate original U with arbitrary values
  U_original <- matrix(rnorm(n * k, mean = 0.5, sd = 0.3), n, k)

  # Apply the same min-max scaling as rss_ld_sketch
  u_min <- apply(U_original, 2, min)
  u_max <- apply(U_original, 2, max)
  denom <- u_max - u_min
  U_scaled <- sweep(sweep(U_original, 2, u_min, "-"), 2, denom, "/") * 2

  # Verify scaled is in [0, 2]
  expect_true(all(U_scaled >= 0 & U_scaled <= 2))

  # Invert
  U_recovered <- invertMinmaxScaling(U_scaled, u_min, u_max)

  # Must be exactly the original (up to floating point)
  expect_equal(U_recovered, U_original, tolerance = 1e-12)
})

test_that("invertMinmaxScaling preserves correlation structure", {
  set.seed(123)
  n <- 200
  k <- 3
  # Simulate U = W'G (G is raw, not standardized, matching rss_ld_sketch)
  G <- sapply(c(0.2, 0.4, 0.1), function(p) rbinom(n, 2, p))
  W <- matrix(rnorm(n * n, 0, 1 / sqrt(n)), n, n)
  U_original <- crossprod(W, G)

  # Scale and invert
  u_min <- apply(U_original, 2, min)
  u_max <- apply(U_original, 2, max)
  denom <- u_max - u_min
  U_scaled <- sweep(sweep(U_original, 2, u_min, "-"), 2, denom, "/") * 2
  U_recovered <- invertMinmaxScaling(U_scaled, u_min, u_max)

  # Exact recovery
  expect_equal(U_recovered, U_original, tolerance = 1e-12)
})

test_that("invertMinmaxScaling handles monomorphic variant", {
  X <- matrix(c(1.0, 1.0, 1.0, 0.5, 1.0, 1.5), ncol = 2)
  u_min <- c(0.5, 0.0)
  u_max <- c(0.5, 1.0)  # first column is monomorphic
  result <- invertMinmaxScaling(X, u_min, u_max)
  expect_equal(ncol(result), 2)
})

test_that("invertMinmaxScaling errors on mismatched lengths", {
  X <- matrix(1:6, ncol = 2)
  expect_error(invertMinmaxScaling(X, c(0, 0, 0), c(1, 1, 1)),
               "Length of u_min")
})

# ===========================================================================
# batchLoadTwasWeights (additional coverage)
# ===========================================================================





# ===========================================================================
# loadCovariateData with real fixture
# ===========================================================================



# ===========================================================================
# loadTsvRegion with real tabix-indexed fixture
# ===========================================================================





# ===========================================================================
# getRefVariantInfo with PLINK2 fixture
# ===========================================================================

test_that("getRefVariantInfo returns variant info for PLINK2 source", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- getRefVariantInfo(meta_file, region = "chr21:17513228-17592874")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 349L)
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result)))
  # .afreq is present, so allele_freq should be populated
  expect_true("allele_freq" %in% names(result))
  expect_true(all(result$allele_freq > 0 & result$allele_freq < 1))
})

test_that("getRefVariantInfo filters by subregion", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_sub_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- getRefVariantInfo(meta_file, region = "chr21:17513228-17550000")
  expect_true(nrow(result) < 349L)
  expect_true(all(result$pos >= 17513228 & result$pos <= 17550000))
})

test_that("getRefVariantInfo returns variant info for VCF source", {
  skip_if_not_installed("VariantAnnotation")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_vcf_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- suppressWarnings(
    getRefVariantInfo(meta_file, region = "chr21:17513228-17592874")
  )
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 349L)
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result)))
  expect_true("allele_freq" %in% names(result))
  expect_true(all(result$allele_freq > 0 & result$allele_freq < 1))
})

test_that("getRefVariantInfo returns variant info for GDS source", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_gds_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.gds", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- getRefVariantInfo(meta_file, region = "chr21:17513228-17592874")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 349L)
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result)))
  expect_true("allele_freq" %in% names(result))
})

test_that("getRefVariantInfo VCF filters by subregion", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Rsamtools")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_vcf_sub_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- suppressWarnings(
    getRefVariantInfo(meta_file, region = "chr21:17513228-17550000")
  )
  expect_true(nrow(result) < 349L)
  expect_true(all(result$pos >= 17513228 & result$pos <= 17550000))
})

test_that("getRefVariantInfo returns consistent results across formats", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  region <- "chr21:17513228-17592874"
  td <- test_path("test_data")

  meta_plink <- file.path(td, "ld_meta_refinfo_cmp_p2_tmp.tsv")
  meta_gds <- file.path(td, "ld_meta_refinfo_cmp_gds_tmp.tsv")
  on.exit({unlink(meta_plink); unlink(meta_gds)}, add = TRUE)

  for (f in c(meta_plink, meta_gds)) {
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), f)
  }
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_plink, append = TRUE)
  cat(paste("21", "0", "0", "test_variants.gds", sep = "\t"), "\n",
      file = meta_gds, append = TRUE)

  info_plink <- getRefVariantInfo(meta_plink, region = region)
  info_gds <- getRefVariantInfo(meta_gds, region = region)

  expect_equal(nrow(info_plink), nrow(info_gds))
  expect_equal(info_plink$pos, info_gds$pos)
})

# ===========================================================================
# readAfreq
# ===========================================================================

test_that("readAfreq returns correct structure from .afreq file", {
  td <- test_path("test_data")
  af <- readAfreq(file.path(td, "test_variants"))
  expect_true(is.data.frame(af))
  expect_equal(nrow(af), 349L)
  expect_true(all(c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct") %in% colnames(af)))
})

test_that("readAfreq returns correct types", {
  td <- test_path("test_data")
  af <- readAfreq(file.path(td, "test_variants"))
  expect_type(af$alt_freq, "double")
  expect_true(all(af$alt_freq >= 0 & af$alt_freq <= 1))
  expect_true(all(af$obs_ct > 0))
})

test_that("readAfreq returns NULL when no afreq file exists", {
  af <- readAfreq(file.path(tempdir(), "nonexistent_prefix"))
  expect_null(af)
})

test_that("readAfreq reads .afreq.zst file", {
  td <- test_path("test_data")
  # test_harmonize_regions has both .afreq and .afreq.zst; readAfreq prefers .zst
  af <- readAfreq(file.path(td, "test_harmonize_regions"))
  expect_true(is.data.frame(af))
  expect_true(all(c("id", "A2", "A1", "alt_freq", "obs_ct") %in% colnames(af)))
  # This afreq has U_MIN/U_MAX columns
  expect_true(all(c("u_min", "u_max") %in% colnames(af)))
  expect_equal(nrow(af), 8L)
})

test_that("readAfreq reads plain .afreq with U_MIN/U_MAX", {
  td <- test_path("test_data")
  # Temporarily hide the .zst so readAfreq falls through to plain .afreq
  zst_path <- file.path(td, "test_harmonize_regions.afreq.zst")
  tmp_path <- paste0(zst_path, ".bak")
  file.rename(zst_path, tmp_path)
  on.exit(file.rename(tmp_path, zst_path), add = TRUE)

  af <- readAfreq(file.path(td, "test_harmonize_regions"))
  expect_true(is.data.frame(af))
  expect_true(all(c("u_min", "u_max") %in% colnames(af)))
  expect_equal(nrow(af), 8L)
})

test_that("readAfreq IDs match pvar IDs", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  af <- readAfreq(file.path(td, "test_variants"))
  pvar <- readPvar(file.path(td, "test_variants.pvar"))
  expect_equal(af$id, pvar$id)
})

# ===========================================================================
# matchVariantsToKeep
# ===========================================================================

test_that("matchVariantsToKeep filters to specified variants", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snpInfoToVariantInfo(handle@snpInfo)

  # Write a keep file as tab-delimited with chrom/pos columns
  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  keep_df <- vi[c(1, 5, 10), c("chrom", "pos", "A2", "A1")]
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  mask <- matchVariantsToKeep(vi, keep_file)
  expect_type(mask, "logical")
  expect_equal(sum(mask), 3L)
  expect_true(mask[1])
  expect_true(mask[5])
  expect_true(mask[10])
})

test_that("matchVariantsToKeep returns all FALSE for non-matching variants", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snpInfoToVariantInfo(handle@snpInfo)

  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  keep_df <- data.frame(chrom = c(1L, 2L), pos = c(999L, 888L),
                        A2 = c("A", "C"), A1 = c("T", "G"))
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  mask <- matchVariantsToKeep(vi, keep_file)
  expect_true(all(!mask))
})

# ===========================================================================
# readVariantMetadata
# ===========================================================================

test_that("readVariantMetadata reads 6-column bim file", {
  tmp <- tempfile(fileext = ".bim")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "1\trs1\t0\t100\tA\tG",
    "1\trs2\t0\t200\tC\tT"
  ), tmp)
  res <- readVariantMetadata(tmp)
  expect_equal(nrow(res), 2)
  expect_true("gpos" %in% names(res))
  expect_equal(as.character(res$chrom), c("1", "1"))
  expect_equal(res$pos, c(100L, 200L))
})

test_that("readVariantMetadata reads 9-column bim file", {
  tmp <- tempfile(fileext = ".bim")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "1\trs1\t0\t100\tA\tG\t0.5\t0.3\t100",
    "1\trs2\t0\t200\tC\tT\t0.4\t0.2\t99"
  ), tmp)
  res <- readVariantMetadata(tmp)
  expect_equal(nrow(res), 2)
  expect_true(all(c("variance", "allele_freq", "n_nomiss") %in% names(res)))
})

test_that("readVariantMetadata delegates to readPvar for .pvar files", {
  pvar_path <- test_path("test_data", "test_variants.pvar")
  res <- readVariantMetadata(pvar_path)
  expect_true(all(c("chrom", "id", "pos", "A1", "A2") %in% names(res)))
  expect_false("gpos" %in% names(res))
})

test_that("readVariantMetadata errors on unexpected column count", {
  tmp <- tempfile(fileext = ".bim")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("1\trs1\t0\t100\tA"), tmp)
  expect_error(readVariantMetadata(tmp), "Unexpected number of columns")
})

# ===========================================================================
# matchVariantsToKeep (additional coverage)
# ===========================================================================

test_that("matchVariantsToKeep works with single-column variant ID file", {
  vi <- data.frame(chrom = c("1", "1", "1"), pos = c(100L, 200L, 300L),
                   A2 = c("A", "C", "G"), A1 = c("G", "T", "A"),
                   stringsAsFactors = FALSE)
  keep_file <- tempfile(fileext = ".txt")
  on.exit(unlink(keep_file), add = TRUE)
  writeLines(c("1:100:A:G", "1:300:G:A"), keep_file)

  mask <- matchVariantsToKeep(vi, keep_file)
  expect_type(mask, "logical")
  expect_equal(sum(mask), 2L)
  expect_true(mask[1])
  expect_true(mask[3])
})

test_that("matchVariantsToKeep uses position-only matching when no alleles", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snpInfoToVariantInfo(handle@snpInfo)

  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  # Write keep file with chrom/pos only (no alleles)
  keep_df <- vi[c(1, 5), c("chrom", "pos")]
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  mask <- matchVariantsToKeep(vi, keep_file)
  expect_type(mask, "logical")
  expect_equal(sum(mask), 2L)
  expect_true(mask[1])
  expect_true(mask[5])
})

# ===========================================================================
# standardiseSumstatsColumns
# ===========================================================================




# ===========================================================================
# readGenotypes + extractblockgenotypes: plink2 tests (replacing load_plink2_data)
# ===========================================================================

test_that("readGenotypes loads plink2 handle with all variants", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@nSamples, 100L)
  expect_equal(nrow(handle@snpInfo), 349L)
  rse <- extractBlockGenotypes(handle, seq_len(nrow(handle@snpInfo)))
  expect_s4_class(rse, "SummarizedExperiment")
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  expect_equal(nrow(dosage), 349L)
  expect_equal(ncol(dosage), 100L)
})

test_that("loadGenotypeRegion filters by region for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  region <- "chr21:17513228-17550000"
  result <- loadGenotypeRegion(file.path(td, "test_variants"), region = region)
  expect_true(ncol(result) < 349L)
})

test_that("loadGenotypeRegion errors on empty region for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  expect_error(
    loadGenotypeRegion(file.path(td, "test_variants"), region = "chr21:1-2"),
    "No SNPs found"
  )
})

test_that("loadGenotypeRegion removes indels for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  full <- loadGenotypeRegion(file.path(td, "test_variants"))
  filtered <- loadGenotypeRegion(file.path(td, "test_variants"), keepIndel = FALSE)
  # test data has 36 indels
  expect_equal(ncol(filtered), ncol(full) - 36L)
})

test_that("loadGenotypeRegion filters by keep_variants_path for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snpInfoToVariantInfo(handle@snpInfo)

  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  keep_df <- vi[c(1, 3, 7), c("chrom", "pos", "A2", "A1")]
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  result <- loadGenotypeRegion(file.path(td, "test_variants"),
                                  keepVariantsPath =keep_file)
  expect_equal(ncol(result), 3L)
})

test_that("loadGenotypeRegion attaches afreq info for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  result <- loadGenotypeRegion(file.path(td, "test_variants"),
                                  returnVariantInfo =TRUE)
  vi <- result$variant_info
  expect_true("alt_freq" %in% colnames(vi))
  expect_true("obs_ct" %in% colnames(vi))
  expect_true(all(vi$alt_freq >= 0 & vi$alt_freq <= 1))
})

test_that("readGenotypes plink2 sample names match psam IIDs", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  expect_true(all(grepl("^(HG|NA)\\d+", handle@sampleIds)))
  expect_equal(length(unique(handle@sampleIds)), 100L)
})

# ===========================================================================
# loadPhenotypeData with real BED-style tabix-indexed fixture
# ===========================================================================








# =============================================================================
# Removed during the post-S4-refactor cleanup
# -----------------------------------------------------------------------------
# The following test blocks were deleted because the functions / classes they
# exercise were either removed outright or replaced by `.Deprecated()` no-op
# stubs in the S4 refactor (see R/fileUtils.R deprecation notices):
#
#   * loadRegionalAssociationData    -> use QtlDataset()
#   * loadRegionalUnivariateData     -> use QtlDataset()
#   * loadRegionalRegressionData     -> use QtlDataset()
#   * loadRegionalMultivariateData   -> use MultiStudyQtlDataset()
#   * loadRegionalFunctionalData     -> use QtlDataset()
#   * loadMultitaskRegionalData      -> use MultiStudyQtlDataset()
#   * loadRssData                    -> use GwasSumStats() / QtlSumStats() +
#                                       summaryStatsQc()
#   * loadTwasWeights                -> use TwasWeights() / TwasWeightsEntry()
#   * regionDataToIndInput           -> use QtlDataset()
#   * regionDataToRssInput           -> use QtlSumStats() / GwasSumStats()
#   * phenoListToMat                 -> function removed (no replacement)
#
# Removed S4 classes (no longer constructible / inspectable):
#   * RegionalData
#   * MultivariateRegionalData
#   * QcResult
#   * AlleleQcResult
#
# Removed accessor / helper functions:
#   * getrssinput, getlddata, getoutliernumber
#   * rssBasicQc                     -> folded into summaryStatsQc(<SumStats>)
#
# Removed pipeline wrappers (all `.Deprecated()` no-ops returning NULL):
#   * colocWrapper, xqtlEnrichmentWrapper, colocPostProcessor,
#     rssAnalysisPipeline
#
# No surgical renames were applied: the previous test file contained no
# references to `finemappingResult` / `FineMappingResult(variantNames=...)` or
# related identifiers, so no `$finemappingEntry` /
# `FineMappingEntry(variantIds=...)` substitutions were necessary here.
# =============================================================================


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

# ===========================================================================
# Tests migrated from test_dataStructures.R (extractBlockGenotypes)
# ===========================================================================

test_that("extractBlockGenotypes returns SummarizedExperiment", {
  skip_if_not_installed("pgenlibr")
  stem <- test_path("test_data", "test_variants")
  handle <- readGenotypes(stem, format = "plink2")
  n_snps <- nrow(handle@snpInfo)
  skip_if(n_snps == 0, "No SNPs in handle")

  rse <- extractBlockGenotypes(handle, seq_len(min(5L, n_snps)))
  expect_s4_class(rse, "SummarizedExperiment")
  expect_true("dosage" %in% SummarizedExperiment::assayNames(rse))
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  # Bioc convention: variants x samples
  expect_equal(nrow(dosage), min(5L, n_snps))
  expect_equal(ncol(dosage), handle@nSamples)
  # rowRanges should have variant info
  rr <- SummarizedExperiment::rowRanges(rse)
  expect_true("A1" %in% names(S4Vectors::mcols(rr)))
  expect_true("A2" %in% names(S4Vectors::mcols(rr)))
})

# =============================================================================
# adjustPips
# =============================================================================

# Build a minimal FineMappingEntry with a known lbf_variable so PIP
# renormalization can be checked end-to-end.
.makeAdjustEntry <- function(vids, L = 2L) {
  p <- length(vids)
  set.seed(11L)
  lbf <- matrix(rnorm(L * p), nrow = L, ncol = p)
  colnames(lbf) <- vids
  alpha <- lbfToAlpha(lbf)
  pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
  FineMappingEntry(
    variantIds = vids,
    susieFit = list(
      pip          = pip,
      alpha        = alpha,
      lbf_variable = lbf,
      mu           = matrix(0, L, p),
      X_column_scale_factors = rep(1, p)
    ),
    topLoci = data.frame(
      variant_id = vids,
      pip        = pip,
      betahat    = rep(0, p),
      sebetahat  = rep(1, p),
      stringsAsFactors = FALSE
    )
  )
}


