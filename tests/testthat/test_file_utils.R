context("file_utils")
library(tidyverse)

test_that("read_pvar reads real pvar file", {
    skip_if_not_installed("pgenlibr")
    pvar_path <- file.path(test_path("test_data"), "test_variants.pvar")
    res <- pecotmr:::read_pvar(pvar_path)
    expect_equal(colnames(res), c("chrom", "id", "pos", "A2", "A1"))
    expect_equal(nrow(res), 349L)
    expect_true(all(res$chrom == "21"))
})

test_that("read_bim dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- read_bim(example_path)
    expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
    expect_equal(nrow(res), 100)
})

test_that("read_fam dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- read_fam(example_path)
    expect_equal(nrow(res), 100)
})

test_that("open_bed dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- open_bed(example_path)
    expect_equal(res$class, "pgen")
})

test_that("find_valid_file_path works",{
    ref_path <- "test_data/protocol_example.genotype.bed"
    expect_error(
        find_valid_file_path(paste0(ref_path, "s"), "protocol_example.genotype.bamf"),
        "Both reference and target file paths do not work. Tried paths: 'test_data/protocol_example.genotype.beds' and 'test_data/protocol_example.genotype.bamf'")
    expect_equal(
        find_valid_file_path(ref_path, "abc"),
        ref_path)
    expect_equal(
        find_valid_file_path(ref_path, "protocol_example.genotype.bim"),
        "test_data/protocol_example.genotype.bim")
    expect_equal(
        find_valid_file_path(ref_path, "test_data/protocol_example.genotype.bim"),
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


test_that("Test load_genotype_region",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype")
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
})

test_that("Test load_genotype_region no indels",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype", keep_indel = F)
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

test_that("Test load_genotype_region with region",{
  res <- load_genotype_region(
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

test_that("Test load_genotype_region with region and no indels",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype",
    region = "chr22:20689453-20845958", keep_indel = F)
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

test_that("load_genotype_region errors on missing genotype files", {
  expect_error(
    load_genotype_region("/nonexistent/geno"),
    "Genotype files not found"
  )
})

# --- find_stochastic_meta tests ---

test_that("find_stochastic_meta finds generic sidecar from PLINK1 prefix", {
  td <- test_path("test_data")
  # test_harmonize_regions has .stochastic_meta.tsv alongside it
  result <- pecotmr:::find_stochastic_meta(file.path(td, "test_harmonize_regions"))
  expect_true(!is.null(result))
  expect_true(grepl("\\.(afreq|stochastic_meta\\.tsv)$", result))
})

test_that("find_stochastic_meta finds sidecar from VCF path", {
  td <- test_path("test_data")
  result <- pecotmr:::find_stochastic_meta(file.path(td, "test_harmonize_regions.vcf.gz"))
  expect_true(!is.null(result))
  expect_true(grepl("\\.(afreq|stochastic_meta\\.tsv)$", result))
})

test_that("find_stochastic_meta finds sidecar from GDS path", {
  td <- test_path("test_data")
  result <- pecotmr:::find_stochastic_meta(file.path(td, "test_harmonize_regions.gds"))
  expect_true(!is.null(result))
  expect_true(grepl("\\.(afreq|stochastic_meta\\.tsv)$", result))
})

test_that("find_stochastic_meta returns NULL when no sidecar exists", {
  td <- test_path("test_data")
  result <- pecotmr:::find_stochastic_meta(file.path(td, "protocol_example.genotype"))
  expect_null(result)
})

# --- read_stochastic_meta tests ---

test_that("read_stochastic_meta reads generic format", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
  result <- pecotmr:::read_stochastic_meta(path)
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
  expect_equal(nrow(result), 8L)
  expect_true(is.numeric(result$u_min))
  expect_true(is.numeric(result$u_max))
})

test_that("read_stochastic_meta reads afreq format", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.afreq")
  result <- pecotmr:::read_stochastic_meta(path)
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
  expect_equal(nrow(result), 8L)
  expect_true(all(grepl("^chr21_", result$id)))
})

test_that("read_stochastic_meta reads afreq.zst format", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.afreq.zst")
  result <- pecotmr:::read_stochastic_meta(path)
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
  expect_equal(nrow(result), 8L)
  # Should produce identical results to the plain afreq
  plain <- pecotmr:::read_stochastic_meta(file.path(td, "test_harmonize_regions.afreq"))
  expect_equal(result, plain)
})

test_that("find_stochastic_meta prefers afreq over afreq.zst", {
  td <- test_path("test_data")
  # Both .afreq and .afreq.zst exist; find_stochastic_meta should return .afreq first
  result <- pecotmr:::find_stochastic_meta(file.path(td, "test_harmonize_regions"))
  expect_true(grepl("\\.afreq$", result))
})

test_that("read_stochastic_meta auto-detects format from extension", {
  td <- test_path("test_data")
  # .afreq extension -> afreq parser
  afreq_result <- pecotmr:::read_stochastic_meta(file.path(td, "test_harmonize_regions.afreq"))
  # .tsv extension -> generic parser
  generic_result <- pecotmr:::read_stochastic_meta(
    file.path(td, "test_harmonize_regions.stochastic_meta.tsv"))
  # Both should return the same u_min/u_max values
  expect_equal(afreq_result$u_min, generic_result$u_min)
  expect_equal(afreq_result$u_max, generic_result$u_max)
  expect_equal(afreq_result$id, generic_result$id)
})

test_that("read_stochastic_meta respects format override", {
  td <- test_path("test_data")
  path <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
  # Explicit generic format should work
  result <- pecotmr:::read_stochastic_meta(path, format = "generic")
  expect_equal(nrow(result), 8L)
  expect_equal(colnames(result), c("id", "u_min", "u_max"))
})

test_that("read_stochastic_meta returns NULL for afreq without U_MIN/U_MAX", {
  td <- test_path("test_data")
  # test_variants.afreq has no U_MIN/U_MAX columns
  path <- file.path(td, "test_variants.afreq")
  result <- pecotmr:::read_stochastic_meta(path)
  expect_null(result)
})

test_that("read_stochastic_meta returns NULL for nonexistent file", {
  result <- pecotmr:::read_stochastic_meta("/nonexistent/file.tsv")
  expect_null(result)
})

# --- load_genotype_region stochastic inversion test ---

test_that("load_genotype_region applies stochastic inversion with explicit sidecar", {
  td <- test_path("test_data")
  meta_path <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
  smeta <- pecotmr:::read_stochastic_meta(meta_path)

  # Load with explicit sidecar - inversion transforms the integer dosages
  res <- load_genotype_region(
    file.path(td, "test_harmonize_regions"),
    return_variant_info = TRUE,
    stochastic_meta_path = meta_path
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
  raw <- load_genotype_region(
    file.path(td, "protocol_example.genotype"),
    region = "chr22:20689453-20845958"
  )
  # protocol_example has no sidecar, so raw values are unchanged (integer dosages)
  expect_true(all(raw == round(raw), na.rm = TRUE))

  # The inverted matrix should NOT be all integers (u_min != 0 or u_max != 2)
  expect_false(all(res$X == round(res$X), na.rm = TRUE))
})

test_that("Test load_covariate_data reads tab-delimited file", {
  # Create a temp covariate file: first column is sample ID, rest are numeric
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("SampleID\tPC1\tPC2", "S1\t0.1\t0.2", "S2\t0.3\t0.4"), tmp)
  result <- load_covariate_data(tmp)
  expect_type(result, "list")
  expect_length(result, 1)
  # Result should be transposed matrix (covariates x samples)
  expect_true(is.matrix(result[[1]]))
  file.remove(tmp)
})

test_that("load_covariate_data errors on non-numeric columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("SampleID\tPC1\tLabel", "S1\t0.1\tabc", "S2\t0.3\tdef"), tmp)
  expect_error(
    load_covariate_data(tmp),
    "Non-numeric columns found in covariate file.*Label.*must be numeric"
  )
  file.remove(tmp)
})

test_that("load_covariate_data errors on missing file", {
  expect_error(
    load_covariate_data("/nonexistent/covar.tsv"),
    "Covariate file.*not found"
  )
})

test_that("Test load_phenotype_data errors on invalid extract_region_name", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("ID\tgene1\tgene2", "S1\t1.0\t2.0"), tmp)
  expect_error(
    load_phenotype_data(tmp, region = NULL, extract_region_name = "not_a_list"),
    "must be NULL or a list"
  )
  file.remove(tmp)
})

test_that("load_phenotype_data errors when extract_region_name length mismatch", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("#chr\tstart\tend\tS1\tS2", "chr1\t100\t200\t1.0\t2.0"), tmp)
  expect_error(
    load_phenotype_data(c(tmp, tmp), region = NULL,
                        extract_region_name = list("gene1")),
    "same length as phenotype_path"
  )
  file.remove(tmp)
})

test_that("load_phenotype_data errors when all phenotype files are empty", {
  local_mocked_bindings(
    tabix_region = function(...) tibble::tibble()
  )
  expect_error(
    load_phenotype_data("fake.gz", region = "chr1:1-100"),
    class = "NoPhenotypeError"
  )
})

test_that("load_phenotype_data with region_name_col out of bounds errors", {
  mock_df <- data.frame(
    chr = "chr1", start = 100, end = 200,
    S1 = 1.0,
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    tabix_region = function(...) mock_df
  )
  expect_error(
    load_phenotype_data("fake.gz", region = "chr1:1-500",
                        extract_region_name = list("gene1"),
                        region_name_col = 99),
    "out of bounds"
  )
})

test_that("load_phenotype_data with extract_region_name and region_name_col filters properly", {
  mock_df <- data.frame(
    chr = c("chr1", "chr1"),
    gene = c("BRCA1", "TP53"),
    start = c(100, 200),
    end = c(150, 250),
    S1 = c(1.0, 2.0),
    S2 = c(3.0, 4.0),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    tabix_region = function(...) mock_df
  )
  result <- load_phenotype_data(
    "fake.gz", region = "chr1:1-500",
    extract_region_name = list("BRCA1"),
    region_name_col = 2
  )
  expect_true(length(result) >= 1)
})

test_that("load_phenotype_data stores kept_indices attribute", {
  mock_df1 <- data.frame(
    chr = "chr1", start = 100, end = 200, S1 = 1.0,
    stringsAsFactors = FALSE
  )
  call_count <- 0
  local_mocked_bindings(
    tabix_region = function(...) {
      call_count <<- call_count + 1
      if (call_count == 1) mock_df1 else tibble::tibble()
    }
  )
  result <- load_phenotype_data(c("f1.gz", "f2.gz"), region = "chr1:1-500")
  expect_true(!is.null(attr(result, "kept_indices")))
  expect_equal(attr(result, "kept_indices"), 1L)
})

test_that("load_phenotype_data assigns colnames from region_name_col without extract_region_name", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5",
               "ENSG002\t300\t400\t3.5\t4.5"), tmp)

  result <- load_phenotype_data(tmp, region = NULL, region_name_col = 1)
  expect_type(result, "list")
  expect_length(result, 1)
  expect_true(all(c("ENSG001", "ENSG002") %in% colnames(result[[1]])))
  file.remove(tmp)
})

test_that("load_phenotype_data errors on empty phenotype file", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2"), tmp)

  expect_error(
    load_phenotype_data(tmp, region = NULL),
    "empty"
  )
  file.remove(tmp)
})

test_that("load_phenotype_data kept_indices reflects filtering", {
  tmp1 <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5"), tmp1)

  tmp2 <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2"), tmp2)

  result <- tryCatch(
    load_phenotype_data(c(tmp1, tmp2), region = NULL),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    idx <- attr(result, "kept_indices")
    expect_true(1 %in% idx)
  }

  file.remove(tmp1, tmp2)
})

test_that("Test filter_by_common_samples",{
    common_samples <- c("Sample_1", "Sample_2", "Sample_3")
    dat <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8), nrow=4, ncol=2))
    rownames(dat) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dat) <- c("chr1:122:G:C", "chr1:123:G:C")
    expect_equal(nrow(filter_by_common_samples(dat, common_samples)), 3)
    expect_equal(rownames(filter_by_common_samples(dat, common_samples)), common_samples)
})

test_that("Test prepare_data_list multiple pheno",{
    # Create dummy data
    ## Prepare Genotype Data
    dummy_geno_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_14")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    ## Prepare Phenotype Data
    dummy_pheno_data_one <- matrix(c("chr1", "222", "223", "1","1","2",NA), nrow=7, ncol=1)
    rownames(dummy_pheno_data_one) <- c("#chr", "start", "end", "Sample_3", "Sample_1", "Sample_2", "Sample_10")
    dummy_pheno_data_two <- matrix(c("chr1", "222", "223", "2","1","2",NA), nrow=7, ncol=1)
    rownames(dummy_pheno_data_two) <- c("#chr", "start", "end", "Sample_3", "Sample_1", "Sample_2", "Sample_10")
    ## Prepare Covariate Data
    dummy_covar_data <- matrix(c(70,71,72,73, 28,30,15,20, 1,2,3,4), nrow=4, ncol=3)
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.025
    mac_cutoff <- 1.0
    xvar_cutoff <- 0.3
    keep_samples <- c("Sample_1", "Sample_2", "Sample_3")
    res <- prepare_data_list(
        dummy_geno_data, list(dummy_pheno_data_one, dummy_pheno_data_two), list(dummy_covar_data, dummy_covar_data),
        imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, phenotype_header = 3, keep_samples=keep_samples)
    # Check that Covar, X, and Y have the same number of rows
    expect_equal(nrow(res$covar[[1]]), 3)
    expect_equal(nrow(res$X[[1]]), 3)
    expect_equal(length(res$Y[[1]]), 3)
    # Check that filter_X occured
    # expect_equal(ncol(res$X[[1]]), 2)
    # Check that Covar, X, and Y have the same samples
    expect_equal(rownames(res$covar[[1]]), rownames(res$X[[1]]))
    expect_equal(rownames(res$covar[[1]]), rownames(res$Y[[1]]))
    expect_equal(rownames(res$X[[1]]), rownames(res$Y[[1]]))
})

test_that("Test prepare_data_list",{
    # Create dummy data
    ## Prepare Genotype Data
    dummy_geno_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_14")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    ## Prepare Phenotype Data
    dummy_pheno_data <- matrix(
        c(
            rep("chr1", 4),
            rep(10, 4),
            rep(11, 4),
            1, NA, NA, NA,
            1, 1, 2, NA,
            2, 1, 2, NA
        ), ncol = 6, nrow = 4
    )
    rownames(dummy_pheno_data) <- c("Pheno_1", "Pheno_2", "Pheno_3", "Pheno_4")
    colnames(dummy_pheno_data) <- c("chrom", "start", "end", "Sample_1", "Sample_2", "Sample_3")
    dummy_pheno_data <- t(dummy_pheno_data)
    ## Prepare Covariate Data
    dummy_covar_data <- matrix(c(70,71,72,73, 28,30,15,20, 1,2,3,4), nrow=4, ncol=3)
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.1
    mac_cutoff <- 1.8
    xvar_cutoff <- 0.3
    keep_samples <- c("Sample_1", "Sample_2", "Sample_3")
    res <- prepare_data_list(
        dummy_geno_data, list(dummy_pheno_data), list(dummy_covar_data), imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff,
        phenotype_header = 3, keep_samples=keep_samples)
    # Check that Covar, X, and Y have the same number of rows
    expect_equal(nrow(res$covar[[1]]), 3)
    expect_equal(nrow(res$X[[1]]), 3)
    expect_equal(nrow(res$Y[[1]]), 3)
    # Check that filter_X occured
    expect_equal(ncol(res$X[[1]]), 2)
    # Check that Covar, X, and Y have the same samples
    expect_equal(rownames(res$covar[[1]]), rownames(res$X[[1]]))
    expect_equal(rownames(res$covar[[1]]), rownames(res$Y[[1]]))
    expect_equal(rownames(res$X[[1]]), rownames(res$Y[[1]]))
})

test_that("Test prepare_X_matrix",{
    dummy_geno_data <- matrix(
        c(1,NA,NA,NA,2, 0,0,1,1,0, 2,2,2,2,2, 1,1,1,2,1, 2,2,0,1,2, 0,1,1,2,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=5, ncol=6)
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    dummy_covar_data <- matrix(
        c(70,71,72,73,74, 28,30,15,20,22, 1,2,3,4,5),
        nrow=5, ncol=3)
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    dummy_data_list <- tibble(
        covar = list(dummy_covar_data))
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.3
    mac_cutoff <- 1.8
    xvar_cutoff <- 0.3
    res <- prepare_X_matrix(dummy_geno_data, dummy_data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
    target <- matrix(c(2,2,0,1,2, 0,1,1,2,2), nrow=5, ncol=2)
    rownames(target) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(target) <- c("chr1:126:G:C", "chr1:127:G:C")
    expect_equal(res, target)
})

test_that("Test add_X_residuals",{
    dummy_geno_data <- matrix(
        c(2,2,0,1, 0,1,1,2),
        nrow=4, ncol=2)
    dummy_covar_data <- matrix(
        c(70,71,72,73, 28,30,15,20, 1,2,3,4),
        nrow=4, ncol=3)
    dummy_data_list <- tibble(
        X = list(dummy_geno_data),
        covar = list(dummy_covar_data))
    res <- add_X_residuals(dummy_data_list)
    res_X <- .lm.fit(x = cbind(1, dummy_covar_data), y = dummy_geno_data)$residuals %>% as.matrix()
    res_X_mean <- apply(res_X, 2, mean)
    res_X_sd <- apply(res_X, 2, sd)
    expect_equal(res$lm_res_X[[1]], res_X)
    expect_equal(res$X_resid_mean[[1]], res_X_mean)
    expect_equal(res$X_resid_sd[[1]], res_X_sd)
})

test_that("add_X_residuals with scale_residuals=TRUE scales output", {
  dummy_X <- matrix(c(2, 2, 0, 1, 0, 1, 1, 2), nrow = 4, ncol = 2)
  dummy_covar <- matrix(c(70, 71, 72, 73, 28, 30, 15, 20), nrow = 4, ncol = 2)
  data_list <- tibble::tibble(
    X = list(dummy_X),
    covar = list(dummy_covar)
  )
  result <- add_X_residuals(data_list, scale_residuals = TRUE)
  resid_mat <- result$X_resid[[1]]
  expect_true(is.matrix(resid_mat))
  col_means <- apply(resid_mat, 2, mean, na.rm = TRUE)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("Test add_Y_residuals",{
    dummy_pheno_data <- rnorm(4)
    dummy_covar_data <- matrix(
        c(70,71,72,73, 28,30,15,20, 1,2,3,4),
        nrow=4, ncol=3)
    dummy_data_list <- tibble(
        Y = list(dummy_pheno_data),
        covar = list(dummy_covar_data))
    conditions <- c("cond_1")
    res_Y <- .lm.fit(x = cbind(1, dummy_covar_data), y = dummy_pheno_data)$residuals %>% as.matrix()
    res_Y_mean <- apply(res_Y, 2, mean)
    res_Y_sd <- apply(res_Y, 2, sd)
    res <- add_Y_residuals(dummy_data_list, conditions)
    expect_equal(res$lm_res[[1]], res_Y)
    expect_equal(res$Y_resid_mean[[1]], res_Y_mean)
    expect_equal(res$Y_resid_sd[[1]], res_Y_sd)
})

test_that("add_Y_residuals with scale_residuals=TRUE scales output", {
  set.seed(42)
  dummy_Y <- rnorm(5)
  names(dummy_Y) <- paste0("S", 1:5)
  dummy_covar <- matrix(rnorm(15), nrow = 5, ncol = 3)
  rownames(dummy_covar) <- paste0("S", 1:5)
  data_list <- tibble::tibble(
    Y = list(dummy_Y),
    covar = list(dummy_covar)
  )
  result <- add_Y_residuals(data_list, conditions = "cond1", scale_residuals = TRUE)
  resid_mat <- result$Y_resid[[1]]
  expect_true(is.matrix(resid_mat))
})

# ===========================================================================
# load_regional_association_data tests
# ===========================================================================

test_that("Test load_regional_association_data complete overlap",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res@genotype_matrix), 10)
    expect_equal(ncol(res@genotype_matrix), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res@genotype_matrix[order(as.numeric(gsub("Sample_", "", rownames(res@genotype_matrix)))), , drop = FALSE], geno_data)
    expect_equal(length(res@phenotypes[[1]]), 10)
    expect_equal(
        as.vector(res@phenotypes[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))]),
        as.numeric(as.vector(asplit(pheno_data[[1]], 2)[[1]])[4:13]))
    expect_equal(nrow(res@covariates[[1]]), 10)
    expect_equal(ncol(res@covariates[[1]]), 5)
    expect_equal(res@covariates[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@covariates[[1]])))), , drop = FALSE], covar_data)
})

test_that("Test load_regional_association_data fewer covar samples",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 3)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res@genotype_matrix), 8)
    expect_equal(ncol(res@genotype_matrix), 9)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(
        res@genotype_matrix[order(as.numeric(gsub("Sample_", "", rownames(res@genotype_matrix)))), , drop = FALSE],
        geno_data[3:10,-6])
    expect_equal(length(res@phenotypes[[1]]), 8)
    expect_equal(
        setNames(res@phenotypes[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))],
rownames(res@phenotypes[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][6:13,]),
            names(pheno_data[[1]][6:13,])))
    expect_equal(nrow(res@covariates[[1]]), 8)
    expect_equal(ncol(res@covariates[[1]]), 5)
    expect_equal(
        res@covariates[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@covariates[[1]])))), , drop = FALSE],
        covar_data[1:8,])
})

test_that("Test load_regional_association_data slight overlap across geno, pheno, covar",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 3)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 7)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res@genotype_matrix), 4)
    expect_equal(ncol(res@genotype_matrix), 3)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(
        res@genotype_matrix[order(as.numeric(gsub("Sample_", "", rownames(res@genotype_matrix)))), , drop = FALSE],
        geno_data[7:10,c(2,4,7)])
    expect_equal(length(res@phenotypes[[1]]), 4)
    expect_equal(
        setNames(res@phenotypes[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))],
rownames(res@phenotypes[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][4:7,]),
            names(pheno_data[[1]][4:7,])))
    expect_equal(nrow(res@covariates[[1]]), 4)
    expect_equal(ncol(res@covariates[[1]]), 5)
    expect_equal(
        res@covariates[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@covariates[[1]])))), , drop = FALSE],
        covar_data[5:8,])
})

test_that("Test load_regional_association_data no overlap",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 11)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 21)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    expect_error(
        load_regional_association_data(
            "dummy_geno.bed.gz",
            "dummy_pheno.bed.gz",
            "dummy_covar.txt.gz",
            "chr1:1000-2000",
            "cond_1",
            imiss_cutoff = 0.70,
            maf_cutoff = 0.1,
            mac_cutoff = (0.1*10*2),
            xvar_cutoff = 0.2,
            phenotype_header = 3,
            keep_samples = NULL),
        "No common complete samples between genotype and phenotype/covariate data")
})

test_that("Test load_regional_association_data unordered samples",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = TRUE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = TRUE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        c("cond_1"),
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res@genotype_matrix), 10)
    expect_equal(ncol(res@genotype_matrix), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res@genotype_matrix[order(as.numeric(gsub("Sample_", "", rownames(res@genotype_matrix)))), , drop = FALSE], geno_data)
    expect_equal(length(res@phenotypes[[1]]), 10)
    expect_equal(
        setNames(res@phenotypes[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))],
rownames(res@phenotypes[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res@phenotypes[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][4:13,]),
            names(pheno_data[[1]][4:13,])))
    expect_equal(nrow(res@covariates[[1]]), 10)
    expect_equal(ncol(res@covariates[[1]]), 5)
    expect_equal(
        res@covariates[[1]][order(as.numeric(gsub("Sample_", "", rownames(res@covariates[[1]])))), , drop = FALSE],
        covar_data[order(as.numeric(gsub("Sample_", "", rownames(covar_data)))), , drop = FALSE])
})

test_that("load_regional_association_data aligns covariates when phenotypes are filtered", {
  geno_data <- dummy_geno_data(
    number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
    number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
  covar_data <- dummy_covar_data(
    number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
  pheno_data1 <- dummy_pheno_data(
    number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)

  local_mocked_bindings(
    load_genotype_region = function(...) geno_data,
    load_covariate_data = function(...) list(covar_data, covar_data),
    load_phenotype_data = function(...) {
      result <- pheno_data1
      attr(result, "kept_indices") <- 1L
      result
    }
  )
  result <- load_regional_association_data(
    "dummy_geno.bed.gz",
    c("dummy_pheno1.bed.gz", "dummy_pheno2.bed.gz"),
    c("dummy_covar1.txt.gz", "dummy_covar2.txt.gz"),
    "chr1:1000-2000",
    c("cond_1", "cond_2"),
    imiss_cutoff = 0.70,
    maf_cutoff = 0.1,
    mac_cutoff = (0.1 * 10 * 2),
    xvar_cutoff = 0.2,
    phenotype_header = 3,
    keep_samples = NULL
  )
  expect_true(!is.null(result@genotype_matrix))
})

test_that("load_regional_association_data returns scalar info when scale_residuals=TRUE", {
  geno_data <- dummy_geno_data(
    number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
    number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
  covar_data <- dummy_covar_data(
    number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
  pheno_data <- dummy_pheno_data(
    number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
  local_mocked_bindings(
    load_genotype_region = function(...) geno_data,
    load_covariate_data = function(...) list(covar_data),
    load_phenotype_data = function(...) pheno_data
  )
  result <- load_regional_association_data(
    "dummy_geno.bed.gz", "dummy_pheno.bed.gz", "dummy_covar.txt.gz",
    "chr1:1000-2000", "cond_1",
    imiss_cutoff = 0.70, maf_cutoff = 0.1, mac_cutoff = (0.1 * 10 * 2),
    xvar_cutoff = 0.2, phenotype_header = 3, keep_samples = NULL,
    scale_residuals = TRUE
  )
  ry_scalar <- getResidualYScalar(result, 1L)
  expect_true(is.numeric(ry_scalar))
  expect_true(all(ry_scalar > 0))
  rx_scalar <- getResidualXScalar(result, 1L)
  expect_true(is.numeric(rx_scalar))
})

test_that("Test load_regional_univariate_data",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_univariate_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_true(is(res, "RegionalData"))
    expect_true(is.matrix(getResidualX(res, 1L)))
    expect_true(is.matrix(getResidualY(res, 1L)))
})

test_that("Test load_regional_regression_data",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_regression_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_true(is(res, "RegionalData"))
    X_mat <- res@genotype_matrix
    expect_equal(nrow(X_mat), 10)
    expect_equal(ncol(X_mat), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(X_mat[order(as.numeric(gsub("Sample_", "", rownames(X_mat)))), , drop = FALSE], geno_data)
    Y_mat <- res@phenotypes[[1]]
    expect_equal(nrow(Y_mat), 10)
    expect_equal(
        setNames(as.numeric(Y_mat[order(as.numeric(gsub("Sample_", "", rownames(Y_mat)))), 1]),
                 rownames(Y_mat)[order(as.numeric(gsub("Sample_", "", rownames(Y_mat))))]),
        setNames(
            as.numeric(pheno_data[[1]][4:13,]),
            names(pheno_data[[1]][4:13,])))
    covar_mat <- res@covariates[[1]]
    expect_equal(nrow(covar_mat), 10)
    expect_equal(ncol(covar_mat), 5)
    expect_equal(covar_mat[order(as.numeric(gsub("Sample_", "", rownames(covar_mat)))), , drop = FALSE], covar_data)
})

test_that("load_regional_multivariate_data filters Y by min completeness", {
  geno_data <- dummy_geno_data(
    number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
    number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
  covar_data <- dummy_covar_data(
    number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
  pheno_data <- dummy_pheno_data(
    number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
  local_mocked_bindings(
    load_genotype_region = function(...) geno_data,
    load_covariate_data = function(...) list(covar_data),
    load_phenotype_data = function(...) pheno_data
  )
  result <- load_regional_multivariate_data(
    matrix_y_min_complete = 5,
    genotype = "dummy_geno.bed.gz",
    phenotype = "dummy_pheno.bed.gz",
    covariate = "dummy_covar.txt.gz",
    region = "chr1:1000-2000",
    conditions = "cond_1",
    imiss_cutoff = 0.70, maf_cutoff = 0.1, mac_cutoff = (0.1 * 10 * 2),
    xvar_cutoff = 0.2, phenotype_header = 3, keep_samples = NULL
  )
  expect_true(is(result, "MultivariateRegionalData"))
  expect_true(is.matrix(result@genotype_matrix))
  expect_true(is.numeric(getMaf(result)))
  expect_true(is.numeric(getXVariance(result)))
})

test_that("load_regional_functional_data returns full association data", {
  geno_data <- dummy_geno_data(
    number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
    number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
  covar_data <- dummy_covar_data(
    number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
  pheno_data <- dummy_pheno_data(
    number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
  local_mocked_bindings(
    load_genotype_region = function(...) geno_data,
    load_covariate_data = function(...) list(covar_data),
    load_phenotype_data = function(...) pheno_data
  )
  result <- load_regional_functional_data(
    genotype = "dummy.bed", phenotype = "dummy.bed.gz",
    covariate = "dummy.txt.gz", region = "chr1:1000-2000",
    conditions = "cond_1",
    imiss_cutoff = 0.70, maf_cutoff = 0.1, mac_cutoff = (0.1 * 10 * 2),
    xvar_cutoff = 0.2, phenotype_header = 3, keep_samples = NULL
  )
  expect_true(is(result, "RegionalData"))
  expect_true(is.matrix(getResidualY(result, 1L)))
  expect_true(is.matrix(result@genotype_matrix))
})

# ===========================================================================
# read_bim vroom-based tests
# ===========================================================================

test_that("read_bim returns correct columns and types", {
  bim_path <- tempfile(fileext = ".bim")
  cat("22\trs100\t0\t50000\tA\tG\n", file = bim_path)
  cat("22\trs200\t0\t60000\tT\tC\n", file = bim_path, append = TRUE)
  cat("22\trs300\t0\t70000\tC\tA\n", file = bim_path, append = TRUE)

  bed_path <- sub("\\.bim$", ".bed", bim_path)
  file.copy(bim_path, bim_path)
  res <- read_bim(bed_path)
  expect_equal(nrow(res), 3)
  expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
  expect_equal(res$id, c("rs100", "rs200", "rs300"))
  expect_equal(res$pos, c(50000, 60000, 70000))
  file.remove(bim_path)
})

# ===========================================================================
# tabix_region
# ===========================================================================

test_that("tabix_region stops when file does not exist", {
  expect_error(
    tabix_region("/nonexistent/path.tsv.gz", "chr1:1-100"),
    "Input file does not exist"
  )
})

test_that("tabix_region returns empty tibble on NULL cmd_output (error path)", {
  tmp <- tempfile()
  writeLines("dummy", tmp)
  local_mocked_bindings(
    read_tabix_region = function(...) stop("mock error")
  )
  result <- tabix_region(tmp, "chr1:1-100")
  expect_true(nrow(result) == 0)
  file.remove(tmp)
})

test_that("tabix_region filters with target and target_column_index", {
  mock_df <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    pos = c(100, 200, 300),
    gene = c("BRCA1", "TP53", "BRCA1"),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile()
  writeLines("dummy", tmp)
  local_mocked_bindings(
    read_tabix_region = function(...) mock_df
  )
  result <- tabix_region(tmp, "chr1:1-500", target = "BRCA1", target_column_index = 3)
  expect_equal(nrow(result), 2)
  file.remove(tmp)
})

test_that("tabix_region filters with target but no target_column_index (text path)", {
  mock_df <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(100, 200),
    name = c("ABC", "DEF"),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile()
  writeLines("dummy", tmp)
  local_mocked_bindings(
    read_tabix_region = function(...) mock_df
  )
  result <- tabix_region(tmp, "chr1:1-500", target = "ABC")
  expect_equal(nrow(result), 1)
  file.remove(tmp)
})

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

test_that("NoPhenotypeError creates proper error condition", {
  err <- NoPhenotypeError("no pheno")
  expect_true(inherits(err, "NoPhenotypeError"))
  expect_true(inherits(err, "error"))
  expect_equal(err$message, "no pheno")
})

# ===========================================================================
# extract_phenotype_coordinates
# ===========================================================================

test_that("extract_phenotype_coordinates returns correct structure", {
  pheno <- list(
    matrix(
      c("chr1", "100", "200", "1.0", "2.0"),
      nrow = 5, ncol = 1,
      dimnames = list(c("#chr", "start", "end", "S1", "S2"), NULL)
    )
  )
  result <- extract_phenotype_coordinates(pheno)
  expect_true(is.list(result))
  expect_true("start" %in% colnames(result[[1]]))
  expect_true("end" %in% colnames(result[[1]]))
  expect_true(is.numeric(result[[1]]$start))
})

# ===========================================================================
# clean_context_names
# ===========================================================================

test_that("clean_context_names removes gene suffix from context", {
  context <- c("tissue1_ENSG00001", "tissue2_ENSG00001", "tissue3_ENSG00002")
  gene <- c("ENSG00001", "ENSG00002")
  result <- clean_context_names(context, gene)
  expect_equal(result, c("tissue1", "tissue2", "tissue3"))
})

test_that("clean_context_names handles multiple gene IDs, longest match first", {
  context <- c("ctx_GENE_LONG", "ctx_GENE")
  gene <- c("GENE", "GENE_LONG")
  result <- clean_context_names(context, gene)
  expect_equal(result, c("ctx", "ctx"))
})

# ===========================================================================
# pheno_list_to_mat
# ===========================================================================

test_that("pheno_list_to_mat converts phenotype list to matrix", {
  data_list <- list(
    residual_Y = list(
      cond1 = matrix(c(1, 2), nrow = 2, ncol = 1, dimnames = list(c("S1", "S2"), "V1")),
      cond2 = matrix(c(5, 6), nrow = 2, ncol = 1, dimnames = list(c("S1", "S3"), "V2"))
    )
  )
  result <- pheno_list_to_mat(data_list)
  expect_true(is.matrix(result$residual_Y))
  expect_equal(sort(rownames(result$residual_Y)), c("S1", "S2", "S3"))
  expect_equal(ncol(result$residual_Y), 2)
})

test_that("pheno_list_to_mat fills NA for missing samples", {
  data_list <- list(
    residual_Y = list(
      cond1 = matrix(1:3, nrow = 3, dimnames = list(c("A", "B", "C"), "V1")),
      cond2 = matrix(4:5, nrow = 2, dimnames = list(c("A", "D"), "V2"))
    )
  )
  result <- pheno_list_to_mat(data_list)
  expect_true(is.na(result$residual_Y["D", 1]))
  expect_true(is.na(result$residual_Y["B", 2]))
})

# ===========================================================================
# load_tsv_region
# ===========================================================================

test_that("load_tsv_region reads plain tsv file", {
  tsv_path <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 300),
    value = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tsv_path)

  res <- suppressWarnings(load_tsv_region(tsv_path))
  expect_equal(nrow(res), 3)
  expect_equal(colnames(res), c("chrom", "pos", "value"))
  expect_equal(res$pos, c(100, 200, 300))
  file.remove(tsv_path)
})

test_that("load_tsv_region reads plain file with region_name filter", {
  tsv_path <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    gene = c("BRCA1", "TP53", "BRCA1"),
    value = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tsv_path)

  res <- suppressWarnings(load_tsv_region(tsv_path,
    extract_region_name = "BRCA1", region_name_col = 2))
  expect_equal(nrow(res), 2)
  expect_equal(res$gene, c("BRCA1", "BRCA1"))
  file.remove(tsv_path)
})

# ===========================================================================
# batch_load_twas_weights
# ===========================================================================

test_that("batch_load_twas_weights returns empty list for empty input", {
  result <- batch_load_twas_weights(list(), data.frame())
  expect_equal(result, list())
})

test_that("batch_load_twas_weights does not split when within memory limit", {
  mock_results <- list(
    gene1 = list(weights = matrix(1:10, nrow = 5)),
    gene2 = list(weights = matrix(1:10, nrow = 5))
  )
  meta_df <- data.frame(
    region_id = c("gene1", "gene2"),
    TSS = c(100, 200),
    stringsAsFactors = FALSE
  )
  result <- batch_load_twas_weights(mock_results, meta_df, max_memory_per_batch = 1000)
  expect_equal(names(result), "all_genes")
  expect_equal(names(result$all_genes), c("gene1", "gene2"))
})

test_that("batch_load_twas_weights splits when exceeding memory limit", {
  mock_results <- list(
    gene1 = list(weights = matrix(rnorm(10000), nrow = 100)),
    gene2 = list(weights = matrix(rnorm(10000), nrow = 100)),
    gene3 = list(weights = matrix(rnorm(10000), nrow = 100))
  )
  meta_df <- data.frame(
    region_id = c("gene1", "gene2", "gene3"),
    TSS = c(100, 200, 300),
    stringsAsFactors = FALSE
  )
  result <- batch_load_twas_weights(mock_results, meta_df, max_memory_per_batch = 0.0001)
  expect_true(length(result) >= 2)
})

# ===========================================================================
# load_rss_data
# ===========================================================================

test_that("load_rss_data errors on missing sumstat file", {
  expect_error(
    load_rss_data("/nonexistent/sumstat.tsv", "/nonexistent/col.txt"),
    "Summary statistics file not found"
  )
})

test_that("load_rss_data errors on missing column file", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  writeLines("dummy", tmp_sumstat)
  expect_error(
    load_rss_data(tmp_sumstat, "/nonexistent/col.txt"),
    "Column mapping file not found"
  )
  file.remove(tmp_sumstat)
})

test_that("load_rss_data computes z from beta and se when z is missing", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(100, 200),
    effect = c(0.5, -0.3),
    stderr = c(0.1, 0.15),
    n = c(1000, 1000),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)

  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:effect", "se:stderr", "n_sample:n"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col))
  expect_true("z" %in% colnames(result$sumstats))
  expect_equal(result$sumstats$z[1], 0.5 / 0.1, tolerance = 1e-10)
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data creates beta from z when beta is missing", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1"),
    pos = c(100),
    zscore = c(4.5),
    n = c(1000),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)

  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("z:zscore", "n_sample:n"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col))
  expect_equal(result$sumstats$beta[1], 4.5)
  expect_equal(result$sumstats$se[1], 1)
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data errors when both n_sample and n_case+n_control are provided", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = "chr1", pos = 100, b = 0.5, s = 0.1, stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s"), tmp_col)

  expect_error(
    suppressWarnings(load_rss_data(tmp_sumstat, tmp_col, n_sample = 100, n_case = 50, n_control = 50)),
    "not both"
  )
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data computes var_y from case/control counts", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = "chr1", pos = 100, b = 0.5, s = 0.1, stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col, n_case = 500, n_control = 500))
  expect_equal(result$n, 1000)
  phi <- 500 / 1000
  expect_equal(result$var_y, 1 / (phi * (1 - phi)))
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data returns NULL n when no sample size info available", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = "chr1", pos = 100, b = 0.5, s = 0.1, stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col))
  expect_null(result$n)
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data returns empty sumstats message for zero-row region", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  writeLines("chrom\tpos\tb\ts", tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col))
  expect_equal(nrow(result$sumstats), 0)
  expect_null(result$n)
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data handles tabix regions with no records", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  writeLines("chrom\tpos\tb\ts", tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s"), tmp_col)

  local_mocked_bindings(
    load_tsv_region = function(...) NULL
  )
  expect_message(
    result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col, region = "chr21:1-2")),
    "No variants in region chr21:1-2."
  )
  expect_true(is.data.frame(result$sumstats))
  expect_equal(nrow(result$sumstats), 0)
  expect_null(result$n)
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data extracts n from n_sample column in sumstats", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1"), pos = c(100, 200),
    b = c(0.5, 0.3), s = c(0.1, 0.1),
    ns = c(1000, 1200),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s", "n_sample:ns"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col))
  expect_equal(result$n, median(c(1000, 1200)))
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data treats NA sample-size arguments as unspecified", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1"), pos = c(100, 200),
    b = c(0.5, 0.3), s = c(0.1, 0.1),
    ns = c(1000, 1200),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s", "n_sample:ns"), tmp_col)

  result <- suppressWarnings(load_rss_data(
    tmp_sumstat, tmp_col, n_sample = NA_real_, n_case = NA_real_,
    n_control = NA_real_))
  expect_equal(result$n, median(c(1000, 1200)))
  file.remove(tmp_sumstat, tmp_col)
})

test_that("load_rss_data extracts n from n_case and n_control columns", {
  tmp_sumstat <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1"), pos = c(100),
    b = c(0.5), s = c(0.1),
    nc = c(500), nco = c(500),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tmp_sumstat)
  tmp_col <- tempfile(fileext = ".txt")
  writeLines(c("beta:b", "se:s", "n_case:nc", "n_control:nco"), tmp_col)

  result <- suppressWarnings(load_rss_data(tmp_sumstat, tmp_col))
  expect_equal(result$n, 1000)
  expect_true(!is.null(result$var_y))
  file.remove(tmp_sumstat, tmp_col)
})

# ===========================================================================
# get_filter_lbf_index
# ===========================================================================

test_that("get_filter_lbf_index returns numeric index vector", {
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

  result <- get_filter_lbf_index(mock_susie, coverage = 0.5, size_factor = 0.5)
  expect_true(is.numeric(result))
})

# ===========================================================================
# get_ref_variant_info
# ===========================================================================

test_that("get_ref_variant_info processes precomputed bim with 6 columns", {
  td <- test_path("test_data")
  meta_file <- file.path(td, "ld_meta_refinfo_6col_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- get_ref_variant_info(meta_file, "chr1:1000-1190")
  expect_true(is.data.frame(result))
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% colnames(result)))
  expect_equal(nrow(result), 5L)
})

test_that("get_ref_variant_info processes precomputed bim with 9 columns", {
  td <- test_path("test_data")
  meta_file <- file.path(td, "ld_meta_refinfo_9col_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("1", "1000", "1200",
            "LD_block_1.chr1_1000_1200.float16.txt.xz,LD_block_1.chr1_1000_1200.float16.9col.bim",
            sep = "\t"), "\n", file = meta_file, append = TRUE)
  result <- get_ref_variant_info(meta_file, "chr1:1000-1190")
  expect_true(all(c("chrom", "id", "pos", "A2", "A1", "variance", "allele_freq", "n_nomiss") %in% colnames(result)))
  expect_equal(nrow(result), 5L)
  expect_equal(result$allele_freq, c(0.3, 0.4, 0.2, 0.5, 0.15))
})

# ===========================================================================
# load_multitask_regional_data
# ===========================================================================

test_that("load_multitask_regional_data errors when no data sources provided", {
  expect_error(
    load_multitask_regional_data(region = "chr1:1-1000"),
    "Data load error"
  )
})

test_that("load_multitask_regional_data errors with multiple genotypes and no match_geno_pheno", {
  expect_error(
    load_multitask_regional_data(
      region = "chr1:1-1000",
      genotype_list = c("geno1.bed", "geno2.bed"),
      phenotype_list = c("pheno1.gz"),
      covariate_list = c("covar1.gz"),
      conditions_list_individual = "cond1"
    ),
    "match_geno_pheno"
  )
})

test_that("load_multitask_regional_data individual-level path returns expected structure", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_multitask_regional_data(
    region = "chr21:17513043-17593579",
    genotype_list = file.path(td, "test_variants"),
    phenotype_list = file.path(td, "test_phenotypes.tsv.gz"),
    covariate_list = file.path(td, "test_covariates.tsv"),
    conditions_list_individual = "cond1"
  )
  expect_true(is.list(result))
  expect_named(result, c("individual_data", "sumstat_data"))
  expect_false(is.null(result$individual_data))
  expect_true(is.null(result$sumstat_data))
  # Individual data should be a RegionalData
  expect_true(is(result$individual_data, "RegionalData"))
  expect_true(is.matrix(getResidualY(result$individual_data, 1L)))
  expect_true(is.matrix(result$individual_data@genotype_matrix))
  expect_true(!is.null(getChrom(result$individual_data)))
})

test_that("load_multitask_regional_data loads and merges multiple genotype groups", {
  calls <- list()
  make_mock_dat <- function(conditions) {
    x <- setNames(lapply(conditions, function(nm) {
      matrix(1, nrow = 2, ncol = 2,
             dimnames = list(c("s1", "s2"), paste0(nm, "_v", 1:2)))
    }), conditions)
    y <- setNames(lapply(conditions, function(nm) {
      matrix(1, nrow = 2, ncol = 1,
             dimnames = list(c("s1", "s2"), nm))
    }), conditions)
    covar_list <- setNames(lapply(conditions, function(nm) {
      matrix(numeric(0), nrow = 2, ncol = 0, dimnames = list(c("s1", "s2"), NULL))
    }), conditions)
    maf_list <- setNames(lapply(x, function(mat) rep(0.1, ncol(mat))), conditions)
    RegionalData(
      genotype_matrix = do.call(cbind, x),
      phenotypes = y,
      covariates = covar_list,
      scale_residuals = FALSE,
      maf = maf_list,
      region = NULL,
      dropped_samples = list(X = list(), Y = list(), covar = list()),
      Y_coordinates = NULL
    )
  }

  local_mocked_bindings(
    load_regional_univariate_data = function(genotype, phenotype, covariate,
                                             conditions, extract_region_name, ...) {
      calls[[length(calls) + 1]] <<- list(
        genotype = genotype,
        phenotype = phenotype,
        covariate = covariate,
        conditions = conditions,
        extract_region_name = extract_region_name
      )
      make_mock_dat(conditions)
    }
  )

  result <- expect_warning(
    load_multitask_regional_data(
      region = "chr1:1-100",
      genotype_list = c("geno1", "geno2"),
      phenotype_list = paste0("pheno", 1:4),
      covariate_list = paste0("covar", 1:4),
      conditions_list_individual = paste0("cond", 1:4),
      match_geno_pheno = c(1, 1, 2, 2),
      extract_region_name = as.list(paste0("gene", 1:4)),
      region_name_col = 4
    ),
    NA
  )

  expect_equal(length(calls), 2L)
  expect_equal(vapply(calls, `[[`, character(1), "genotype"), c("geno1", "geno2"))
  expect_equal(calls[[1]]$phenotype, paste0("pheno", 1:2))
  expect_equal(calls[[2]]$phenotype, paste0("pheno", 3:4))
  expect_equal(calls[[1]]$extract_region_name, as.list(paste0("gene", 1:2)))
  expect_equal(calls[[2]]$extract_region_name, as.list(paste0("gene", 3:4)))
  expect_true(is(result$individual_data, "RegionalData"))
  expect_equal(names(result$individual_data@phenotypes), paste0("cond", 1:4))
})

test_that("load_multitask_regional_data defaults missing individual condition names", {
  calls <- list()
  make_mock_dat <- function(conditions) {
    x <- setNames(lapply(conditions, function(nm) {
      matrix(1, nrow = 2, ncol = 2,
             dimnames = list(c("s1", "s2"), paste0(nm, "_v", 1:2)))
    }), conditions)
    y <- setNames(lapply(conditions, function(nm) {
      matrix(1, nrow = 2, ncol = 1,
             dimnames = list(c("s1", "s2"), nm))
    }), conditions)
    covar_list <- setNames(lapply(conditions, function(nm) {
      matrix(numeric(0), nrow = 2, ncol = 0, dimnames = list(c("s1", "s2"), NULL))
    }), conditions)
    maf_list <- setNames(lapply(x, function(mat) rep(0.1, ncol(mat))), conditions)
    RegionalData(
      genotype_matrix = do.call(cbind, x),
      phenotypes = y,
      covariates = covar_list,
      scale_residuals = FALSE,
      maf = maf_list,
      region = NULL,
      dropped_samples = list(X = list(), Y = list(), covar = list()),
      Y_coordinates = NULL
    )
  }

  local_mocked_bindings(
    load_regional_univariate_data = function(genotype, phenotype, covariate,
                                             conditions, ...) {
      calls[[length(calls) + 1]] <<- list(
        genotype = genotype,
        phenotype = phenotype,
        covariate = covariate,
        conditions = conditions
      )
      make_mock_dat(conditions)
    }
  )

  result <- expect_warning(
    load_multitask_regional_data(
      region = "chr1:1-100",
      genotype_list = c("geno1", "geno2"),
      phenotype_list = paste0("pheno", 1:4),
      covariate_list = paste0("covar", 1:4),
      match_geno_pheno = c(1, 1, 2, 2)
    ),
    "conditions_list_individual"
  )

  expect_equal(length(calls), 2L)
  expect_equal(calls[[1]]$conditions, paste0("condition", 1:2))
  expect_equal(calls[[2]]$conditions, paste0("condition", 3:4))
  expect_equal(names(result$individual_data@phenotypes), paste0("condition", 1:4))
})

test_that("load_multitask_regional_data validates individual input vector lengths", {
  expect_error(
    load_multitask_regional_data(
      region = "chr1:1-100",
      genotype_list = "geno",
      phenotype_list = paste0("pheno", 1:2),
      covariate_list = "covar",
      conditions_list_individual = paste0("cond", 1:2)
    ),
    "phenotype_list"
  )
  expect_error(
    load_multitask_regional_data(
      region = "chr1:1-100",
      genotype_list = "geno",
      phenotype_list = paste0("pheno", 1:2),
      covariate_list = paste0("covar", 1:2),
      conditions_list_individual = "cond"
    ),
    "conditions_list_individual"
  )
  expect_error(
    load_multitask_regional_data(
      region = "chr1:1-100",
      genotype_list = c("geno1", "geno2"),
      phenotype_list = paste0("pheno", 1:2),
      covariate_list = paste0("covar", 1:2),
      conditions_list_individual = paste0("cond", 1:2),
      match_geno_pheno = c(1, 3)
    ),
    "match_geno_pheno"
  )
})

test_that("load_multitask_regional_data sumstat path returns expected structure", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  # Create LD metadata file pointing to harmonize_regions genotype
  meta_file <- file.path(td, "ld_meta_harmonize_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("21", "0", "0", "test_harmonize_regions", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  sumstat_path <- file.path(td, "test_sumstats.tsv.gz")
  result <- suppressMessages(suppressWarnings(
    load_multitask_regional_data(
      region = "chr21:17014042-45433269",
      association_window = "chr21:17014042-45433269",
      sumstat_path_list = sumstat_path,
      LD_meta_file_path_list = meta_file,
      conditions_list_sumstat = "sumstat_cond1",
      n_samples = 1000,
      n_cases = 0,
      n_controls = 0
    )
  ))
  expect_true(is.list(result))
  expect_named(result, c("individual_data", "sumstat_data"))
  expect_true(is.null(result$individual_data))
  expect_false(is.null(result$sumstat_data))
  # Sumstat data should have sumstats and LD_info lists
  expect_true("sumstats" %in% names(result$sumstat_data))
  expect_true("LD_info" %in% names(result$sumstat_data))
  expect_equal(length(result$sumstat_data$sumstats), 1L)
  expect_equal(length(result$sumstat_data$LD_info), 1L)
  # The inner sumstats should be a named list with the condition
  inner <- result$sumstat_data$sumstats[[1]]
  expect_true("sumstat_cond1" %in% names(inner))
  ss <- inner[["sumstat_cond1"]]
  expect_true(is.data.frame(ss$sumstats))
  expect_true(nrow(ss$sumstats) > 0)
  expect_true("z" %in% names(ss$sumstats))
  expect_equal(ss$n, 1000)
})

test_that("load_multitask_regional_data both paths simultaneously", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  # Create LD metadata for sumstat path
  meta_file <- file.path(td, "ld_meta_both_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("21", "0", "0", "test_harmonize_regions", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  result <- suppressMessages(suppressWarnings(
    load_multitask_regional_data(
      region = "chr21:17513043-17593579",
      # Individual-level data
      genotype_list = file.path(td, "test_variants"),
      phenotype_list = file.path(td, "test_phenotypes.tsv.gz"),
      covariate_list = file.path(td, "test_covariates.tsv"),
      conditions_list_individual = "ind_cond1",
      # Summary statistics
      association_window = "chr21:17014042-45433269",
      sumstat_path_list = file.path(td, "test_sumstats.tsv.gz"),
      LD_meta_file_path_list = meta_file,
      conditions_list_sumstat = "ss_cond1",
      n_samples = 500,
      n_cases = 0,
      n_controls = 0
    )
  ))
  expect_false(is.null(result$individual_data))
  expect_false(is.null(result$sumstat_data))
  # Both paths should produce valid data
  expect_true(is(result$individual_data, "RegionalData"))
  expect_true(is.matrix(result$individual_data@genotype_matrix))
  expect_true(is.data.frame(result$sumstat_data$sumstats[[1]][["ss_cond1"]]$sumstats))
})

test_that("load_multitask_regional_data sumstat path errors on mismatched match_LD_sumstat", {
  td <- test_path("test_data")
  expect_error(
    load_multitask_regional_data(
      region = "chr21:17014042-45433269",
      sumstat_path_list = file.path(td, "test_sumstats.tsv.gz"),
      LD_meta_file_path_list = c("meta1.tsv", "meta2.tsv"),
      match_LD_sumstat = list("cond1"),
      conditions_list_sumstat = "cond1",
      n_samples = 100,
      n_cases = 0,
      n_controls = 0
    ),
    "match_LD_sumstat"
  )
})

# ---- invert_minmax_scaling ----

test_that("invert_minmax_scaling exactly recovers original U", {
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
  U_recovered <- invert_minmax_scaling(U_scaled, u_min, u_max)

  # Must be exactly the original (up to floating point)
  expect_equal(U_recovered, U_original, tolerance = 1e-12)
})

test_that("invert_minmax_scaling preserves correlation structure", {
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
  U_recovered <- invert_minmax_scaling(U_scaled, u_min, u_max)

  # Exact recovery
  expect_equal(U_recovered, U_original, tolerance = 1e-12)
})

test_that("invert_minmax_scaling handles monomorphic variant", {
  X <- matrix(c(1.0, 1.0, 1.0, 0.5, 1.0, 1.5), ncol = 2)
  u_min <- c(0.5, 0.0)
  u_max <- c(0.5, 1.0)  # first column is monomorphic
  result <- invert_minmax_scaling(X, u_min, u_max)
  expect_equal(ncol(result), 2)
})

test_that("invert_minmax_scaling errors on mismatched lengths", {
  X <- matrix(1:6, ncol = 2)
  expect_error(invert_minmax_scaling(X, c(0, 0, 0), c(1, 1, 1)),
               "Length of u_min")
})

# ===========================================================================
# batch_load_twas_weights
# ===========================================================================

test_that("batch_load_twas_weights returns empty list for empty input", {
  expect_message(
    result <- batch_load_twas_weights(list(), data.frame(region_id = character(), TSS = integer())),
    "No genes"
  )
  expect_equal(length(result), 0)
})

test_that("batch_load_twas_weights returns single batch when total memory fits", {
  twas <- list(
    gene1 = list(a = 1:10),
    gene2 = list(a = 1:10)
  )
  meta <- data.frame(region_id = c("gene1", "gene2"), TSS = c(100, 200))
  expect_message(
    result <- batch_load_twas_weights(twas, meta, max_memory_per_batch = 1000),
    "No need to split"
  )
  expect_equal(length(result), 1)
  expect_true(all(c("gene1", "gene2") %in% names(result[[1]])))
})

test_that("batch_load_twas_weights splits into multiple batches", {
  # Create data large enough to require splitting
  twas <- list(
    gene1 = rnorm(1e5),
    gene2 = rnorm(1e5),
    gene3 = rnorm(1e5)
  )
  # Each gene is ~0.76 MB
  gene_size_mb <- as.numeric(object.size(twas[[1]])) / (1024^2)
  # Set limit so at most 2 genes fit per batch
  max_mb <- gene_size_mb * 1.5
  meta <- data.frame(region_id = c("gene1", "gene2", "gene3"), TSS = c(100, 200, 300))
  result <- batch_load_twas_weights(twas, meta, max_memory_per_batch = max_mb)
  expect_true(length(result) >= 2)
  # All genes should be present across batches
  all_genes <- unlist(lapply(result, names))
  expect_true(all(c("gene1", "gene2", "gene3") %in% all_genes))
})

test_that("batch_load_twas_weights puts oversized gene in its own batch", {
  twas <- list(
    gene_small = list(a = 1:10),
    gene_big = rnorm(1e6)
  )
  big_size_mb <- as.numeric(object.size(twas$gene_big)) / (1024^2)
  small_size_mb <- as.numeric(object.size(twas$gene_small)) / (1024^2)
  # Set limit between small and big
  max_mb <- big_size_mb * 0.5
  meta <- data.frame(region_id = c("gene_small", "gene_big"), TSS = c(100, 200))
  result <- batch_load_twas_weights(twas, meta, max_memory_per_batch = max_mb)
  # Big gene should be in its own batch
  expect_true(length(result) >= 2)
})

# ===========================================================================
# load_covariate_data with real fixture
# ===========================================================================

test_that("load_covariate_data reads and transposes covariate file", {
  covar_path <- file.path(test_path("test_data"), "test_covariates.tsv")
  result <- pecotmr:::load_covariate_data(covar_path)
  expect_true(is.list(result))
  expect_equal(length(result), 1L)
  mat <- result[[1]]
  expect_true(is.matrix(mat))
  # Original: 8 rows (PCs) x 101 cols (variable + 100 samples)
  # After drop col 1 + transpose: 100 rows (samples) x 8 cols (PCs)
  expect_equal(nrow(mat), 100L)
  expect_equal(ncol(mat), 8L)
  expect_true(is.numeric(mat))
  expect_false(any(is.na(mat)))
})

test_that("load_covariate_data errors on missing file", {
  expect_error(
    pecotmr:::load_covariate_data("/nonexistent/covariate.tsv"),
    "not found"
  )
})

# ===========================================================================
# load_tsv_region with real tabix-indexed fixture
# ===========================================================================

test_that("load_tsv_region reads full gz file without region", {
  skip_if_not_installed("Rsamtools")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- load_tsv_region(sumstat_path)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 8L)
  expect_true("BETA" %in% names(result) || "beta" %in% names(result))
})

test_that("load_tsv_region queries region via tabix", {
  skip_if_not_installed("Rsamtools")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- load_tsv_region(sumstat_path, region = "chr21:17014042-45433269")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 8L)
})

test_that("load_tsv_region errors for non-overlapping region", {
  skip_if_not_installed("Rsamtools")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  expect_error(load_tsv_region(sumstat_path, region = "chr1:1-2"), "tabix-indexed")
})

test_that("load_tsv_region queries subregion correctly", {
  skip_if_not_installed("Rsamtools")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  # Only first 2 variants: pos 17014042 and 18759786
  result <- load_tsv_region(sumstat_path, region = "chr21:17014042-18759786")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2L)
})

# ===========================================================================
# load_rss_data with real tabix-indexed fixture
# ===========================================================================

test_that("load_rss_data reads summary statistics without region", {
  skip_if_not_installed("MungeSumstats")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- suppressMessages(load_rss_data(sumstat_path))
  expect_true(is.list(result))
  expect_true(is.data.frame(result$sumstats))
  expect_equal(nrow(result$sumstats), 8L)
  # Should have standardized column names including z
  expect_true("beta" %in% names(result$sumstats))
  expect_true("se" %in% names(result$sumstats))
  expect_true("z" %in% names(result$sumstats))
  # z should be beta / se
  expect_equal(result$sumstats$z, result$sumstats$beta / result$sumstats$se)
})

test_that("load_rss_data reads summary statistics with region", {
  skip_if_not_installed("MungeSumstats")
  skip_if_not_installed("Rsamtools")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- suppressMessages(
    load_rss_data(sumstat_path, region = "chr21:17014042-45433269")
  )
  expect_true(is.data.frame(result$sumstats))
  expect_equal(nrow(result$sumstats), 8L)
})

test_that("load_rss_data with n_sample returns sample size", {
  skip_if_not_installed("MungeSumstats")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- suppressMessages(load_rss_data(sumstat_path, n_sample = 500))
  expect_equal(result$n, 500)
  expect_null(result$var_y)
})

test_that("load_rss_data with n_case and n_control computes var_y", {
  skip_if_not_installed("MungeSumstats")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- suppressMessages(load_rss_data(sumstat_path, n_case = 200, n_control = 300))
  expect_equal(result$n, 500)
  expect_true(!is.null(result$var_y))
  # var_y = 1 / (phi * (1 - phi)) where phi = 200/500 = 0.4
  expect_equal(result$var_y, 1 / (0.4 * 0.6))
})

test_that("load_rss_data extracts n from sumstats N column", {
  skip_if_not_installed("MungeSumstats")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  # n_sample=0 means "get from file"
  result <- suppressMessages(load_rss_data(sumstat_path))
  # The fixture has N column with per-variant sample sizes; median should be used
  expect_true(!is.null(result$n))
  expect_true(result$n > 0)
})

test_that("load_rss_data errors on missing file", {
  expect_error(load_rss_data("/nonexistent/sumstats.tsv.gz"), "not found")
})

test_that("load_rss_data errors for non-overlapping region", {
  skip_if_not_installed("MungeSumstats")
  skip_if_not_installed("Rsamtools")
  sumstat_path <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  expect_error(
    suppressMessages(load_rss_data(sumstat_path, region = "chr1:1-2")),
    "tabix-indexed"
  )
})

# ===========================================================================
# get_ref_variant_info with PLINK2 fixture
# ===========================================================================

test_that("get_ref_variant_info returns variant info for PLINK2 source", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- get_ref_variant_info(meta_file, region = "chr21:17513228-17592874")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 349L)
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result)))
  # .afreq is present, so allele_freq should be populated
  expect_true("allele_freq" %in% names(result))
  expect_true(all(result$allele_freq > 0 & result$allele_freq < 1))
})

test_that("get_ref_variant_info filters by subregion", {
  skip_if_not_installed("pgenlibr")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_sub_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- get_ref_variant_info(meta_file, region = "chr21:17513228-17550000")
  expect_true(nrow(result) < 349L)
  expect_true(all(result$pos >= 17513228 & result$pos <= 17550000))
})

test_that("get_ref_variant_info returns variant info for VCF source", {
  skip_if_not_installed("VariantAnnotation")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_vcf_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- suppressWarnings(
    get_ref_variant_info(meta_file, region = "chr21:17513228-17592874")
  )
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 349L)
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result)))
  expect_true("allele_freq" %in% names(result))
  expect_true(all(result$allele_freq > 0 & result$allele_freq < 1))
})

test_that("get_ref_variant_info returns variant info for GDS source", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_gds_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.gds", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- get_ref_variant_info(meta_file, region = "chr21:17513228-17592874")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 349L)
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(result)))
  expect_true("allele_freq" %in% names(result))
})

test_that("get_ref_variant_info VCF filters by subregion", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("Rsamtools")
  meta_file <- file.path(test_path("test_data"), "ld_meta_refinfo_vcf_sub_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
  cat(paste("21", "0", "0", "test_variants.vcf.gz", sep = "\t"), "\n",
      file = meta_file, append = TRUE)
  result <- suppressWarnings(
    get_ref_variant_info(meta_file, region = "chr21:17513228-17550000")
  )
  expect_true(nrow(result) < 349L)
  expect_true(all(result$pos >= 17513228 & result$pos <= 17550000))
})

test_that("get_ref_variant_info returns consistent results across formats", {
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

  info_plink <- get_ref_variant_info(meta_plink, region = region)
  info_gds <- get_ref_variant_info(meta_gds, region = region)

  expect_equal(nrow(info_plink), nrow(info_gds))
  expect_equal(info_plink$pos, info_gds$pos)
})

# ===========================================================================
# read_afreq
# ===========================================================================

test_that("read_afreq returns correct structure from .afreq file", {
  td <- test_path("test_data")
  af <- read_afreq(file.path(td, "test_variants"))
  expect_true(is.data.frame(af))
  expect_equal(nrow(af), 349L)
  expect_true(all(c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct") %in% colnames(af)))
})

test_that("read_afreq returns correct types", {
  td <- test_path("test_data")
  af <- read_afreq(file.path(td, "test_variants"))
  expect_type(af$alt_freq, "double")
  expect_true(all(af$alt_freq >= 0 & af$alt_freq <= 1))
  expect_true(all(af$obs_ct > 0))
})

test_that("read_afreq returns NULL when no afreq file exists", {
  af <- read_afreq(file.path(tempdir(), "nonexistent_prefix"))
  expect_null(af)
})

test_that("read_afreq reads .afreq.zst file", {
  td <- test_path("test_data")
  # test_harmonize_regions has both .afreq and .afreq.zst; read_afreq prefers .zst
  af <- read_afreq(file.path(td, "test_harmonize_regions"))
  expect_true(is.data.frame(af))
  expect_true(all(c("id", "A2", "A1", "alt_freq", "obs_ct") %in% colnames(af)))
  # This afreq has U_MIN/U_MAX columns
  expect_true(all(c("u_min", "u_max") %in% colnames(af)))
  expect_equal(nrow(af), 8L)
})

test_that("read_afreq reads plain .afreq with U_MIN/U_MAX", {
  td <- test_path("test_data")
  # Temporarily hide the .zst so read_afreq falls through to plain .afreq
  zst_path <- file.path(td, "test_harmonize_regions.afreq.zst")
  tmp_path <- paste0(zst_path, ".bak")
  file.rename(zst_path, tmp_path)
  on.exit(file.rename(tmp_path, zst_path), add = TRUE)

  af <- read_afreq(file.path(td, "test_harmonize_regions"))
  expect_true(is.data.frame(af))
  expect_true(all(c("u_min", "u_max") %in% colnames(af)))
  expect_equal(nrow(af), 8L)
})

test_that("read_afreq IDs match pvar IDs", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  af <- read_afreq(file.path(td, "test_variants"))
  pvar <- read_pvar(file.path(td, "test_variants.pvar"))
  expect_equal(af$id, pvar$id)
})

# ===========================================================================
# match_variants_to_keep
# ===========================================================================

test_that("match_variants_to_keep filters to specified variants", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snp_info_to_variant_info(handle@snp_info)

  # Write a keep file as tab-delimited with chrom/pos columns
  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  keep_df <- vi[c(1, 5, 10), c("chrom", "pos", "A2", "A1")]
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  mask <- match_variants_to_keep(vi, keep_file)
  expect_type(mask, "logical")
  expect_equal(sum(mask), 3L)
  expect_true(mask[1])
  expect_true(mask[5])
  expect_true(mask[10])
})

test_that("match_variants_to_keep returns all FALSE for non-matching variants", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snp_info_to_variant_info(handle@snp_info)

  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  keep_df <- data.frame(chrom = c(1L, 2L), pos = c(999L, 888L),
                        A2 = c("A", "C"), A1 = c("T", "G"))
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  mask <- match_variants_to_keep(vi, keep_file)
  expect_true(all(!mask))
})

# ===========================================================================
# read_variant_metadata
# ===========================================================================

test_that("read_variant_metadata reads 6-column bim file", {
  tmp <- tempfile(fileext = ".bim")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "1\trs1\t0\t100\tA\tG",
    "1\trs2\t0\t200\tC\tT"
  ), tmp)
  res <- read_variant_metadata(tmp)
  expect_equal(nrow(res), 2)
  expect_true("gpos" %in% names(res))
  expect_equal(as.character(res$chrom), c("1", "1"))
  expect_equal(res$pos, c(100L, 200L))
})

test_that("read_variant_metadata reads 9-column bim file", {
  tmp <- tempfile(fileext = ".bim")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "1\trs1\t0\t100\tA\tG\t0.5\t0.3\t100",
    "1\trs2\t0\t200\tC\tT\t0.4\t0.2\t99"
  ), tmp)
  res <- read_variant_metadata(tmp)
  expect_equal(nrow(res), 2)
  expect_true(all(c("variance", "allele_freq", "n_nomiss") %in% names(res)))
})

test_that("read_variant_metadata delegates to read_pvar for .pvar files", {
  pvar_path <- test_path("test_data", "test_variants.pvar")
  res <- read_variant_metadata(pvar_path)
  expect_true(all(c("chrom", "id", "pos", "A1", "A2") %in% names(res)))
  expect_false("gpos" %in% names(res))
})

test_that("read_variant_metadata errors on unexpected column count", {
  tmp <- tempfile(fileext = ".bim")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("1\trs1\t0\t100\tA"), tmp)
  expect_error(read_variant_metadata(tmp), "Unexpected number of columns")
})

# ===========================================================================
# match_variants_to_keep (additional coverage)
# ===========================================================================

test_that("match_variants_to_keep works with single-column variant ID file", {
  vi <- data.frame(chrom = c("1", "1", "1"), pos = c(100L, 200L, 300L),
                   A2 = c("A", "C", "G"), A1 = c("G", "T", "A"),
                   stringsAsFactors = FALSE)
  keep_file <- tempfile(fileext = ".txt")
  on.exit(unlink(keep_file), add = TRUE)
  writeLines(c("1:100:A:G", "1:300:G:A"), keep_file)

  mask <- match_variants_to_keep(vi, keep_file)
  expect_type(mask, "logical")
  expect_equal(sum(mask), 2L)
  expect_true(mask[1])
  expect_true(mask[3])
})

test_that("match_variants_to_keep uses position-only matching when no alleles", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snp_info_to_variant_info(handle@snp_info)

  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  # Write keep file with chrom/pos only (no alleles)
  keep_df <- vi[c(1, 5), c("chrom", "pos")]
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  mask <- match_variants_to_keep(vi, keep_file)
  expect_type(mask, "logical")
  expect_equal(sum(mask), 2L)
  expect_true(mask[1])
  expect_true(mask[5])
})

# ===========================================================================
# standardise_sumstats_columns
# ===========================================================================

test_that("standardise_sumstats_columns renames standard headers", {
  skip_if_not_installed("MungeSumstats")
  df <- data.frame(
    SNPID = "rs1", CHR = 1, POS = 100,
    EFFECT_ALLELE = "A", OTHER_ALLELE = "G",
    BETA = 0.5, SE = 0.1, P = 0.01,
    stringsAsFactors = FALSE
  )
  result <- standardise_sumstats_columns(df)
  expect_true("chrom" %in% colnames(result))
  expect_true("pos" %in% colnames(result))
  expect_true("beta" %in% colnames(result))
  expect_true("se" %in% colnames(result))
  expect_true("p" %in% colnames(result))
})

test_that("standardise_sumstats_columns applies custom column mapping", {
  skip_if_not_installed("MungeSumstats")
  df <- data.frame(
    SNP = "rs1", CHR = 1, BP = 100,
    A1 = "A", A2 = "G",
    BETA = 0.5, SE = 0.1, P = 0.01,
    MY_FREQ = 0.3,
    stringsAsFactors = FALSE
  )
  col_file <- tempfile(fileext = ".txt")
  on.exit(unlink(col_file), add = TRUE)
  writeLines("maf:MY_FREQ", col_file)

  result <- standardise_sumstats_columns(df, column_file_path = col_file)
  expect_true("maf" %in% colnames(result))
})

test_that("standardise_sumstats_columns errors on missing column file", {
  skip_if_not_installed("MungeSumstats")
  df <- data.frame(SNP = "rs1", CHR = 1, BP = 100, A1 = "A", A2 = "G",
                   BETA = 0.5, SE = 0.1, P = 0.01, stringsAsFactors = FALSE)
  expect_error(
    standardise_sumstats_columns(df, column_file_path = "/no/such/file.txt"),
    "Column mapping file not found"
  )
})

# ===========================================================================
# readGenotypes + extractBlockGenotypes: plink2 tests (replacing load_plink2_data)
# ===========================================================================

test_that("readGenotypes loads plink2 handle with all variants", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  expect_s4_class(handle, "GenotypeHandle")
  expect_equal(handle@n_samples, 100L)
  expect_equal(nrow(handle@snp_info), 349L)
  rse <- extractBlockGenotypes(handle, seq_len(nrow(handle@snp_info)))
  expect_s4_class(rse, "SummarizedExperiment")
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  expect_equal(nrow(dosage), 349L)
  expect_equal(ncol(dosage), 100L)
})

test_that("load_genotype_region filters by region for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  region <- "chr21:17513228-17550000"
  result <- load_genotype_region(file.path(td, "test_variants"), region = region)
  expect_true(ncol(result) < 349L)
})

test_that("load_genotype_region errors on empty region for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  expect_error(
    load_genotype_region(file.path(td, "test_variants"), region = "chr21:1-2"),
    "No SNPs found"
  )
})

test_that("load_genotype_region removes indels for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  full <- load_genotype_region(file.path(td, "test_variants"))
  filtered <- load_genotype_region(file.path(td, "test_variants"), keep_indel = FALSE)
  # test data has 36 indels
  expect_equal(ncol(filtered), ncol(full) - 36L)
})

test_that("load_genotype_region filters by keep_variants_path for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  vi <- pecotmr:::.snp_info_to_variant_info(handle@snp_info)

  keep_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(keep_file), add = TRUE)
  keep_df <- vi[c(1, 3, 7), c("chrom", "pos", "A2", "A1")]
  vroom::vroom_write(keep_df, keep_file, delim = "\t")

  result <- load_genotype_region(file.path(td, "test_variants"),
                                  keep_variants_path = keep_file)
  expect_equal(ncol(result), 3L)
})

test_that("load_genotype_region attaches afreq info for plink2", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  result <- load_genotype_region(file.path(td, "test_variants"),
                                  return_variant_info = TRUE)
  vi <- result$variant_info
  expect_true("alt_freq" %in% colnames(vi))
  expect_true("obs_ct" %in% colnames(vi))
  expect_true(all(vi$alt_freq >= 0 & vi$alt_freq <= 1))
})

test_that("readGenotypes plink2 sample names match psam IIDs", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  handle <- readGenotypes(file.path(td, "test_variants"), format = "plink2")
  expect_true(all(grepl("^(HG|NA)\\d+", handle@sample_ids)))
  expect_equal(length(unique(handle@sample_ids)), 100L)
})

# ===========================================================================
# load_phenotype_data with real BED-style tabix-indexed fixture
# ===========================================================================

test_that("load_phenotype_data reads compressed file with tabix region", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- load_phenotype_data(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579"
  )
  expect_true(is.list(pheno))
  expect_equal(length(pheno), 1L)
  mat <- pheno[[1]]
  expect_true(is.matrix(mat))
  # 4 header rows (seqid, start, end, gene_id) + 100 samples = 104 rows, 1 gene column
  expect_equal(nrow(mat), 104L)
  expect_equal(ncol(mat), 1L)
  # Sample IDs start at row 5
  expect_equal(rownames(mat)[5], "HG02461")
})

test_that("load_phenotype_data filters by extract_region_name and region_name_col", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- load_phenotype_data(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579",
    extract_region_name = list(c("ENSG00000154639")),
    region_name_col = 4
  )
  expect_equal(length(pheno), 1L)
  expect_true("ENSG00000154639" %in% colnames(pheno[[1]]))
})

test_that("load_phenotype_data assigns gene names with region_name_col", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- load_phenotype_data(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579",
    region_name_col = 4
  )
  expect_equal(colnames(pheno[[1]]), "ENSG00000154639")
})

test_that("load_phenotype_data returns multiple genes for broad region", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- load_phenotype_data(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:14000000-18000000",
    region_name_col = 4
  )
  expect_equal(length(pheno), 1L)
  expect_true(ncol(pheno[[1]]) > 1)
  expect_true("ENSG00000154639" %in% colnames(pheno[[1]]))
})

test_that("load_phenotype_data errors on non-overlapping region", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  expect_error(
    load_phenotype_data(
      file.path(td, "test_phenotypes.tsv.gz"),
      region = "chr21:1-100"
    ),
    "empty"
  )
})

test_that("load_phenotype_data reads uncompressed file without region", {
  td <- test_path("test_data")
  pheno <- load_phenotype_data(
    file.path(td, "test_phenotypes.tsv"),
    region = NULL,
    region_name_col = 4
  )
  expect_equal(length(pheno), 1L)
  # All 93 genes in the uncompressed file
  expect_equal(ncol(pheno[[1]]), 93L)
})

test_that("load_phenotype_data stores kept_indices attribute", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- load_phenotype_data(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579"
  )
  expect_equal(attr(pheno, "kept_indices"), 1L)
})

# ===========================================================================
# load_regional_association_data with real fixtures (full pipeline)
# ===========================================================================

test_that("load_regional_association_data returns expected structure", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_association_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1"
  )
  expect_s4_class(result, "RegionalData")
  # 100 samples, 349 variants
  expect_equal(nrow(result@genotype_matrix), 100L)
  expect_equal(ncol(result@genotype_matrix), 349L)
  expect_equal(as.character(GenomicRanges::seqnames(result@region))[1], "chr21")
  expect_equal(length(result@phenotypes), 1L)
  # residual_Y should be a 100-sample x 1-gene matrix
  ry <- getResidualY(result, 1L)
  expect_equal(nrow(ry), 100L)
  expect_equal(ncol(ry), 1L)
  # Y_coordinates should have gene coordinates
  expect_true(is.data.frame(result@Y_coordinates[[1]]))
})

test_that("load_regional_association_data with scale_residuals returns scalars", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_association_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1",
    scale_residuals = TRUE
  )
  # With scale_residuals, residual_Y_scalar should be non-trivial
  ry_scalar <- getResidualYScalar(result, 1L)
  expect_true(is.numeric(ry_scalar))
  expect_true(all(ry_scalar > 0))
  rx_scalar <- getResidualXScalar(result, 1L)
  expect_true(is.numeric(rx_scalar))
})

test_that("load_regional_association_data with keep_indel=FALSE reduces variants", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_association_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1",
    keep_indel = FALSE
  )
  # 349 total - 36 indels = 313 SNPs
  expect_equal(ncol(result@genotype_matrix), 313L)
})

test_that("load_regional_association_data covariate residuals affect Y", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_association_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1"
  )
  # Raw Y and residual Y should differ (covariates regressed out)
  raw_y <- as.numeric(result@phenotypes[[1]])
  resid_y <- as.numeric(getResidualY(result, 1L))
  expect_false(isTRUE(all.equal(raw_y, resid_y)))
})

# ===========================================================================
# load_regional_univariate_data with real fixtures
# ===========================================================================

test_that("load_regional_univariate_data returns correct fields", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_univariate_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1"
  )
  expect_true(is(result, "RegionalData"))
  expect_equal(nrow(result@genotype_matrix), 100L)
  # Per-condition accessors should return valid data
  expect_true(is.matrix(getResidualX(result, 1L)))
  expect_true(is.matrix(getResidualY(result, 1L)))
  x_var <- getXVariance(result, 1L)
  expect_true(is.numeric(x_var))
  expect_equal(length(x_var), ncol(result@genotype_matrix))
  # Chrom and grange accessors should work
  expect_true(!is.null(getChrom(result)))
  expect_true(!is.null(getGrange(result)))
})

# ===========================================================================
# load_regional_regression_data with real fixtures
# ===========================================================================

test_that("load_regional_regression_data returns correct fields", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_regression_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1"
  )
  expect_true(is(result, "RegionalData"))
  expect_true(is.list(result@phenotypes))
  expect_true(is.matrix(result@phenotypes[[1]]))
  expect_true(is.matrix(result@genotype_matrix))
  expect_true(is.list(result@covariates))
  expect_true(is.matrix(result@covariates[[1]]))
})

# ===========================================================================
# load_regional_multivariate_data with real fixtures
# ===========================================================================

test_that("load_regional_multivariate_data returns correct fields", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  result <- load_regional_multivariate_data(
    genotype = file.path(td, "test_variants"),
    phenotype = file.path(td, "test_phenotypes.tsv.gz"),
    covariate = file.path(td, "test_covariates.tsv"),
    region = "chr21:17513043-17593579",
    conditions = "cond1"
  )
  expect_true(is(result, "MultivariateRegionalData"))
  expect_true(is.matrix(getYMatrix(result)))
  expect_equal(nrow(result@genotype_matrix), 100L)
  expect_true(is.numeric(getYScalar(result)))
  expect_true(is.numeric(getMaf(result)))
  expect_true(is.numeric(getXVariance(result)))
})

# =============================================================================
# load_twas_weights file pre-validation
# =============================================================================

test_that("load_twas_weights skips non-existent files with warning", {
  expect_warning(
    tryCatch(
      load_twas_weights(c("/nonexistent/path/fake.rds")),
      error = function(e) NULL
    ),
    "does not exist"
  )
})

test_that("load_twas_weights skips too-small files with warning", {
  tmp <- tempfile(fileext = ".rds")
  writeLines("x", tmp)  # tiny file, not valid RDS
  on.exit(unlink(tmp))
  expect_warning(
    tryCatch(
      load_twas_weights(tmp),
      error = function(e) NULL
    ),
    "too small"
  )
})

test_that("load_twas_weights skips corrupt RDS files with warning", {
  tmp <- tempfile(fileext = ".rds")
  writeBin(as.raw(rep(0L, 500)), tmp)  # 500 bytes of garbage
  on.exit(unlink(tmp))
  expect_warning(
    tryCatch(
      load_twas_weights(tmp),
      error = function(e) NULL
    ),
    "failed to read RDS"
  )
})

test_that("load_twas_weights skips non-list RDS with warning", {
  tmp <- tempfile(fileext = ".rds")
  saveRDS(paste0("x", seq_len(10000)), tmp)  # valid RDS but not a list; large enough to pass size check
  on.exit(unlink(tmp))
  expect_warning(
    tryCatch(
      load_twas_weights(tmp),
      error = function(e) NULL
    ),
    "unexpected structure"
  )
})

# =============================================================================
# load_rss_data sample size validation
# =============================================================================

test_that("load_rss_data rejects negative sample size", {
  skip_if_not_installed("MungeSumstats")
  sumstat_file <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  skip_if_not(file.exists(sumstat_file), "test sumstat file not found")
  expect_error(
    suppressMessages(load_rss_data(sumstat_file, n_sample = -100)),
    "Invalid sample size"
  )
})

test_that("load_rss_data rejects Inf sample size", {
  skip_if_not_installed("MungeSumstats")
  sumstat_file <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  skip_if_not(file.exists(sumstat_file), "test sumstat file not found")
  expect_error(
    suppressMessages(load_rss_data(sumstat_file, n_sample = Inf)),
    "Invalid sample size"
  )
})
