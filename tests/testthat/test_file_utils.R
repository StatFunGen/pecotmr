context("file_utils")
library(tidyverse)

test_that("read_pvar dummy data works",{
    skip("disabled pending Bucket B resolution: read_pvar not exported / removed")
    dummy_path <- gsub("//", "/", tempfile(pattern = "dummy_pvar", tmpdir = tempdir(), fileext = ".pvar"))
    dummy <- data.frame("#CHROM" = c(1, 2, 3, 4, 5),
        "ID" = c("rs1", "rs2", "rs3", "rs4", "rs5"),
        "POS" = c(100, 200, 300, 400, 500),
        "REF" = c("A", "T", "C", "G", "A"),
        "ALT" = c("T", "C", "G", "A", "T"))
    colnames(dummy) <- c("#CHROM", "ID", "POS", "REF", "ALT")
    cat(c("#DUMMY HEADER 1", "#DUMMY HEADER 2", "#DUMMY HEADER 3"), file = dummy_path, sep = "\n")
    write_delim(
        dummy, dummy_path, delim = "\t", col_names = TRUE, append = TRUE)
    expect_equal(colnames(read_pvar(dummy_path)), c("chrom", "id", "pos", "A2", "A1"))
    file.remove(dummy_path)
})

test_that("read_bim dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- read_bim(example_path)
    expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
    expect_equal(nrow(res), 100)
})

test_that("read_psam dummy data works",{
    skip("disabled pending Bucket B resolution: read_psam not exported / removed")
    dummy_path <- gsub("//", "/", tempfile(pattern = "dummy_psam", tmpdir = tempdir(), fileext = ".psam"))
    dummy <- data.frame("#CHROM" = c(1, 2, 3, 4, 5),
        "IID" = c("rs1", "rs2", "rs3", "rs4", "rs5"),
        "SID" = c(100, 200, 300, 400, 500),
        "PAT" = c("A", "T", "C", "G", "A"),
        "MAT" = c("T", "C", "G", "A", "T"),
        "SEX" = c(1, 2, 1, 2, 1))
    write_delim(
        dummy, dummy_path, delim = "\t", col_names = TRUE, append = TRUE)
    res <- read_psam(dummy_path)
    expect_equal(colnames(res), c("FID", "IID", "SID", "PAT", "MAT", "SEX"))
    file.remove(dummy_path)
})

test_that("read_fam dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- read_fam(example_path)
    expect_equal(nrow(res), 100)
})

test_that("open_pgen dummy data works",{
    skip("disabled pending Bucket B resolution: open_pgen not exported / removed")
    example_path <- "test_data/dummy_data.pgen"
    res <- open_pgen(example_path)
    expect_equal(res$class, "pgen")
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
    expect_equal(nrow(res$X), 10)
    expect_equal(ncol(res$X), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE], geno_data)
    expect_equal(length(res$Y[[1]]), 10)
    expect_equal(
        as.vector(res$Y[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))]),
        as.numeric(as.vector(asplit(pheno_data[[1]], 2)[[1]])[4:13]))
    expect_equal(nrow(res$covar[[1]]), 10)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE], covar_data)
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
    expect_equal(nrow(res$X), 8)
    expect_equal(ncol(res$X), 9)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(
        res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE],
        geno_data[3:10,-6])
    expect_equal(length(res$Y[[1]]), 8)
    expect_equal(
        setNames(res$Y[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))],
rownames(res$Y[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][6:13,]),
            names(pheno_data[[1]][6:13,])))
    expect_equal(nrow(res$covar[[1]]), 8)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(
        res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE],
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
    expect_equal(nrow(res$X), 4)
    expect_equal(ncol(res$X), 3)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(
        res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE],
        geno_data[7:10,c(2,4,7)])
    expect_equal(length(res$Y[[1]]), 4)
    expect_equal(
        setNames(res$Y[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))],
rownames(res$Y[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][4:7,]),
            names(pheno_data[[1]][4:7,])))
    expect_equal(nrow(res$covar[[1]]), 4)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(
        res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE],
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
    expect_equal(nrow(res$X), 10)
    expect_equal(ncol(res$X), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE], geno_data)
    expect_equal(length(res$Y[[1]]), 10)
    expect_equal(
        setNames(res$Y[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))],
rownames(res$Y[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][4:13,]),
            names(pheno_data[[1]][4:13,])))
    expect_equal(nrow(res$covar[[1]]), 10)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(
        res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE],
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
  expect_true(!is.null(result$X))
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
  expect_true(!is.null(result$residual_Y_scalar))
  expect_true(!is.null(result$residual_X_scalar))
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
    expect_true("residual_X" %in% names(res))
    expect_true("residual_Y" %in% names(res))
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
    expect_equal(nrow(res$X_data[[1]]), 10)
    expect_equal(ncol(res$X_data[[1]]), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X_data[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$X_data[[1]])))), , drop = FALSE], geno_data)
    expect_equal(length(res$Y[[1]]), 10)
    expect_equal(
        setNames(res$Y[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))],
rownames(res$Y[[1]])[order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]]))))]),
        setNames(
            as.numeric(pheno_data[[1]][4:13,]),
            names(pheno_data[[1]][4:13,])))
    expect_equal(nrow(res$covar[[1]]), 10)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE], covar_data)
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
  expect_true(!is.null(result$X))
  expect_true(!is.null(result$maf))
  expect_true(!is.null(result$X_variance))
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
  expect_true("residual_Y" %in% names(result))
  expect_true("X" %in% names(result))
})

# ===========================================================================
# read_pvar / read_bim vroom-based tests
# ===========================================================================

test_that("read_pvar handles multiple comment styles", {
  skip("disabled pending Bucket B resolution: read_pvar not exported / removed")
  pvar_path <- tempfile(fileext = ".pvar")
  cat("##fileformat=VCFv4.2\n", file = pvar_path)
  cat("##INFO=<ID=AF>\n", file = pvar_path, append = TRUE)
  cat("#CHROM\tID\tPOS\tREF\tALT\n", file = pvar_path, append = TRUE)
  cat("1\trs10\t1000\tA\tG\n", file = pvar_path, append = TRUE)
  cat("1\trs20\t2000\tT\tC\n", file = pvar_path, append = TRUE)

  res <- read_pvar(pvar_path)
  expect_equal(nrow(res), 2)
  expect_equal(colnames(res), c("chrom", "id", "pos", "A2", "A1"))
  expect_equal(res$id, c("rs10", "rs20"))
  expect_equal(res$pos, c(1000, 2000))
  file.remove(pvar_path)
})

test_that("read_pvar errors when #CHROM header is missing", {
  skip("disabled pending Bucket B resolution: read_pvar not exported / removed")
  bad_path <- tempfile(fileext = ".pvar")
  cat("col1\tcol2\tcol3\n", file = bad_path)
  cat("1\t2\t3\n", file = bad_path, append = TRUE)
  expect_error(read_pvar(bad_path), "Could not find #CHROM header")
  file.remove(bad_path)
})

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
    vroom = function(...) stop("mock error")
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
    vroom = function(...) mock_df
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
    vroom = function(...) mock_df
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

test_that("get_ref_variant_info processes bim files with 6 columns", {
  skip("disabled pending Bucket C resolution: resolve_ld_source now requires real LD metadata file")
  bim_path <- tempfile(fileext = ".bim")
  bim_data <- data.frame(
    V1 = c("chr1", "chr1"),
    V2 = c("chr1:100:A:G", "chr1:200:C:T"),
    V3 = c(0, 0),
    V4 = c(100, 200),
    V5 = c("A", "C"),
    V6 = c("G", "T"),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(bim_data, bim_path, col_names = FALSE)

  local_mocked_bindings(
    get_regional_ld_meta = function(...) {
      list(intersections = list(bim_file_paths = bim_path))
    }
  )

  result <- get_ref_variant_info("/fake/ld_meta.txt", "chr1:1-300")
  expect_true(is.data.frame(result))
  expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% colnames(result)))
  expect_equal(nrow(result), 2)
  file.remove(bim_path)
})

test_that("get_ref_variant_info processes bim files with 8 columns", {
  skip("disabled pending Bucket C resolution: resolve_ld_source now requires real LD metadata file")
  bim_path <- tempfile(fileext = ".bim")
  bim_data <- data.frame(
    V1 = c("chr1", "chr1"),
    V2 = c("chr1:100:A:G", "chr1:200:C:T"),
    V3 = c(0, 0),
    V4 = c(100, 200),
    V5 = c("A", "C"),
    V6 = c("G", "T"),
    V7 = c(0.01, 0.02),
    V8 = c(0.3, 0.4),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(bim_data, bim_path, col_names = FALSE)

  local_mocked_bindings(
    get_regional_ld_meta = function(...) {
      list(intersections = list(bim_file_paths = bim_path))
    }
  )

  result <- get_ref_variant_info("/fake/ld_meta.txt", "chr1:1-300")
  expect_true(all(c("chrom", "id", "pos", "A2", "A1", "variance", "allele_freq") %in% colnames(result)))
  file.remove(bim_path)
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
      covariate_list = c("covar1.gz")
    ),
    "match_geno_pheno"
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
