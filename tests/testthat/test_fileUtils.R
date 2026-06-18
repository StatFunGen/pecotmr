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
        findValidFilePath(paste0(ref_path, "s"), "protocol_example.genotype.bamf"),
        "Both reference and target file paths do not work. Tried paths: 'test_data/protocol_example.genotype.beds' and 'test_data/protocol_example.genotype.bamf'")
    expect_equal(
        findValidFilePath(ref_path, "abc"),
        ref_path)
    expect_equal(
        findValidFilePath(ref_path, "protocol_example.genotype.bim"),
        "test_data/protocol_example.genotype.bim")
    expect_equal(
        findValidFilePath(ref_path, "test_data/protocol_example.genotype.bim"),
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

test_that("Test loadCovariateData reads tab-delimited file", {
  # Create a temp covariate file: first column is sample ID, rest are numeric
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("SampleID\tPC1\tPC2", "S1\t0.1\t0.2", "S2\t0.3\t0.4"), tmp)
  result <- loadCovariateData(tmp)
  expect_type(result, "list")
  expect_length(result, 1)
  # Result should be transposed matrix (covariates x samples)
  expect_true(is.matrix(result[[1]]))
  file.remove(tmp)
})

test_that("loadCovariateData errors on non-numeric columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("SampleID\tPC1\tLabel", "S1\t0.1\tabc", "S2\t0.3\tdef"), tmp)
  expect_error(
    loadCovariateData(tmp),
    "Non-numeric columns found in covariate file.*Label.*must be numeric"
  )
  file.remove(tmp)
})

test_that("loadCovariateData errors on missing file", {
  expect_error(
    loadCovariateData("/nonexistent/covar.tsv"),
    "Covariate file.*not found"
  )
})

test_that("Test loadPhenotypeData errors on invalid extract_region_name", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("ID\tgene1\tgene2", "S1\t1.0\t2.0"), tmp)
  expect_error(
    loadPhenotypeData(tmp, region = NULL, extractRegionName ="not_a_list"),
    "must be NULL or a list"
  )
  file.remove(tmp)
})

test_that("loadPhenotypeData errors when extract_region_name length mismatch", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("#chr\tstart\tend\tS1\tS2", "chr1\t100\t200\t1.0\t2.0"), tmp)
  expect_error(
    loadPhenotypeData(c(tmp, tmp), region = NULL,
                        extractRegionName =list("gene1")),
    "same length as phenotype_path"
  )
  file.remove(tmp)
})

test_that("loadPhenotypeData errors when all phenotype files are empty", {
  local_mocked_bindings(
    tabixRegion = function(...) tibble::tibble()
  )
  expect_error(
    loadPhenotypeData("fake.gz", region = "chr1:1-100"),
    class = "NoPhenotypeError"
  )
})

test_that("loadPhenotypeData with region_name_col out of bounds errors", {
  mock_df <- data.frame(
    chr = "chr1", start = 100, end = 200,
    S1 = 1.0,
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    tabixRegion = function(...) mock_df
  )
  expect_error(
    loadPhenotypeData("fake.gz", region = "chr1:1-500",
                        extractRegionName =list("gene1"),
                        regionNameCol =99),
    "out of bounds"
  )
})

test_that("loadPhenotypeData with extract_region_name and region_name_col filters properly", {
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
    tabixRegion = function(...) mock_df
  )
  result <- loadPhenotypeData(
    "fake.gz", region = "chr1:1-500",
    extractRegionName =list("BRCA1"),
    regionNameCol =2
  )
  expect_true(length(result) >= 1)
})

test_that("loadPhenotypeData stores kept_indices attribute", {
  mock_df1 <- data.frame(
    chr = "chr1", start = 100, end = 200, S1 = 1.0,
    stringsAsFactors = FALSE
  )
  call_count <- 0
  local_mocked_bindings(
    tabixRegion = function(...) {
      call_count <<- call_count + 1
      if (call_count == 1) mock_df1 else tibble::tibble()
    }
  )
  result <- loadPhenotypeData(c("f1.gz", "f2.gz"), region = "chr1:1-500")
  expect_true(!is.null(attr(result, "kept_indices")))
  expect_equal(attr(result, "kept_indices"), 1L)
})

test_that("loadPhenotypeData assigns colnames from region_name_col without extract_region_name", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5",
               "ENSG002\t300\t400\t3.5\t4.5"), tmp)

  result <- loadPhenotypeData(tmp, region = NULL, regionNameCol =1)
  expect_type(result, "list")
  expect_length(result, 1)
  expect_true(all(c("ENSG001", "ENSG002") %in% colnames(result[[1]])))
  file.remove(tmp)
})

test_that("loadPhenotypeData errors on empty phenotype file", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2"), tmp)

  expect_error(
    loadPhenotypeData(tmp, region = NULL),
    "empty"
  )
  file.remove(tmp)
})

test_that("loadPhenotypeData kept_indices reflects filtering", {
  tmp1 <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5"), tmp1)

  tmp2 <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2"), tmp2)

  result <- tryCatch(
    loadPhenotypeData(c(tmp1, tmp2), region = NULL),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    idx <- attr(result, "kept_indices")
    expect_true(1 %in% idx)
  }

  file.remove(tmp1, tmp2)
})

test_that("Test filterByCommonSamples",{
    common_samples <- c("Sample_1", "Sample_2", "Sample_3")
    dat <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8), nrow=4, ncol=2))
    rownames(dat) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dat) <- c("chr1:122:G:C", "chr1:123:G:C")
    expect_equal(nrow(filterByCommonSamples(dat, common_samples)), 3)
    expect_equal(rownames(filterByCommonSamples(dat, common_samples)), common_samples)
})

test_that("Test prepareDataList multiple pheno",{
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
    res <- prepareDataList(
        dummy_geno_data, list(dummy_pheno_data_one, dummy_pheno_data_two), list(dummy_covar_data, dummy_covar_data),
        imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, phenotypeHeader =3, keepSamples=keep_samples)
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

test_that("Test prepareDataList",{
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
    res <- prepareDataList(
        dummy_geno_data, list(dummy_pheno_data), list(dummy_covar_data), imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff,
        phenotypeHeader =3, keepSamples=keep_samples)
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

test_that("Test prepareXMatrix",{
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
    res <- prepareXMatrix(dummy_geno_data, dummy_data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
    target <- matrix(c(2,2,0,1,2, 0,1,1,2,2), nrow=5, ncol=2)
    rownames(target) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(target) <- c("chr1:126:G:C", "chr1:127:G:C")
    expect_equal(res, target)
})

test_that("Test addXResiduals",{
    dummy_geno_data <- matrix(
        c(2,2,0,1, 0,1,1,2),
        nrow=4, ncol=2)
    dummy_covar_data <- matrix(
        c(70,71,72,73, 28,30,15,20, 1,2,3,4),
        nrow=4, ncol=3)
    dummy_data_list <- tibble(
        X = list(dummy_geno_data),
        covar = list(dummy_covar_data))
    res <- addXResiduals(dummy_data_list)
    res_X <- .lm.fit(x = cbind(1, dummy_covar_data), y = dummy_geno_data)$residuals %>% as.matrix()
    res_X_mean <- apply(res_X, 2, mean)
    res_X_sd <- apply(res_X, 2, sd)
    expect_equal(res$lm_res_X[[1]], res_X)
    expect_equal(res$X_resid_mean[[1]], res_X_mean)
    expect_equal(res$X_resid_sd[[1]], res_X_sd)
})

test_that("addXResiduals with scale_residuals=TRUE scales output", {
  dummy_X <- matrix(c(2, 2, 0, 1, 0, 1, 1, 2), nrow = 4, ncol = 2)
  dummy_covar <- matrix(c(70, 71, 72, 73, 28, 30, 15, 20), nrow = 4, ncol = 2)
  data_list <- tibble::tibble(
    X = list(dummy_X),
    covar = list(dummy_covar)
  )
  result <- addXResiduals(data_list, scaleResiduals = TRUE)
  resid_mat <- result$X_resid[[1]]
  expect_true(is.matrix(resid_mat))
  col_means <- apply(resid_mat, 2, mean, na.rm = TRUE)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("Test addYResiduals",{
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
    res <- addYResiduals(dummy_data_list, conditions)
    expect_equal(res$lm_res[[1]], res_Y)
    expect_equal(res$Y_resid_mean[[1]], res_Y_mean)
    expect_equal(res$Y_resid_sd[[1]], res_Y_sd)
})

test_that("addYResiduals with scale_residuals=TRUE scales output", {
  set.seed(42)
  dummy_Y <- rnorm(5)
  names(dummy_Y) <- paste0("S", 1:5)
  dummy_covar <- matrix(rnorm(15), nrow = 5, ncol = 3)
  rownames(dummy_covar) <- paste0("S", 1:5)
  data_list <- tibble::tibble(
    Y = list(dummy_Y),
    covar = list(dummy_covar)
  )
  result <- addYResiduals(data_list, conditions = "cond1", scaleResiduals = TRUE)
  resid_mat <- result$Y_resid[[1]]
  expect_true(is.matrix(resid_mat))
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

test_that("tabixRegion stops when file does not exist", {
  expect_error(
    tabixRegion("/nonexistent/path.tsv.gz", "chr1:1-100"),
    "Input file does not exist"
  )
})

test_that("tabixRegion returns empty tibble on NULL cmd_output (error path)", {
  tmp <- tempfile()
  writeLines("dummy", tmp)
  local_mocked_bindings(
    readTabixRegion = function(...) stop("mock error")
  )
  result <- tabixRegion(tmp, "chr1:1-100")
  expect_true(nrow(result) == 0)
  file.remove(tmp)
})

test_that("tabixRegion filters with target and target_column_index", {
  mock_df <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    pos = c(100, 200, 300),
    gene = c("BRCA1", "TP53", "BRCA1"),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile()
  writeLines("dummy", tmp)
  local_mocked_bindings(
    readTabixRegion = function(...) mock_df
  )
  result <- tabixRegion(tmp, "chr1:1-500", target = "BRCA1", targetColumnIndex =3)
  expect_equal(nrow(result), 2)
  file.remove(tmp)
})

test_that("tabixRegion filters with target but no target_column_index (text path)", {
  mock_df <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(100, 200),
    name = c("ABC", "DEF"),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile()
  writeLines("dummy", tmp)
  local_mocked_bindings(
    readTabixRegion = function(...) mock_df
  )
  result <- tabixRegion(tmp, "chr1:1-500", target = "ABC")
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
# extractPhenotypeCoordinates
# ===========================================================================

test_that("extractPhenotypeCoordinates returns correct structure", {
  pheno <- list(
    matrix(
      c("chr1", "100", "200", "1.0", "2.0"),
      nrow = 5, ncol = 1,
      dimnames = list(c("#chr", "start", "end", "S1", "S2"), NULL)
    )
  )
  result <- extractPhenotypeCoordinates(pheno)
  expect_true(is.list(result))
  expect_true("start" %in% colnames(result[[1]]))
  expect_true("end" %in% colnames(result[[1]]))
  expect_true(is.numeric(result[[1]]$start))
})

# ===========================================================================
# cleanContextNames
# ===========================================================================

test_that("cleanContextNames removes gene suffix from context", {
  context <- c("tissue1_ENSG00001", "tissue2_ENSG00001", "tissue3_ENSG00002")
  gene <- c("ENSG00001", "ENSG00002")
  result <- cleanContextNames(context, gene)
  expect_equal(result, c("tissue1", "tissue2", "tissue3"))
})

test_that("cleanContextNames handles multiple gene IDs, longest match first", {
  context <- c("ctx_GENE_LONG", "ctx_GENE")
  gene <- c("GENE", "GENE_LONG")
  result <- cleanContextNames(context, gene)
  expect_equal(result, c("ctx", "ctx"))
})

# ===========================================================================
# loadTsvRegion
# ===========================================================================

test_that("loadTsvRegion reads plain tsv file", {
  tsv_path <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 300),
    value = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tsv_path)

  res <- suppressWarnings(loadTsvRegion(tsv_path))
  expect_equal(nrow(res), 3)
  expect_equal(colnames(res), c("chrom", "pos", "value"))
  expect_equal(res$pos, c(100, 200, 300))
  file.remove(tsv_path)
})

test_that("loadTsvRegion reads plain file with region_name filter", {
  tsv_path <- tempfile(fileext = ".tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    gene = c("BRCA1", "TP53", "BRCA1"),
    value = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(df, tsv_path)

  res <- suppressWarnings(loadTsvRegion(tsv_path,
    extractRegionName ="BRCA1", regionNameCol =2))
  expect_equal(nrow(res), 2)
  expect_equal(res$gene, c("BRCA1", "BRCA1"))
  file.remove(tsv_path)
})

# ===========================================================================
# batchLoadTwasWeights
# ===========================================================================

test_that("batchLoadTwasWeights returns empty list for empty input", {
  result <- batchLoadTwasWeights(list(), data.frame())
  expect_equal(result, list())
})

test_that("batchLoadTwasWeights does not split when within memory limit", {
  mock_results <- list(
    gene1 = list(weights = matrix(1:10, nrow = 5)),
    gene2 = list(weights = matrix(1:10, nrow = 5))
  )
  meta_df <- data.frame(
    region_id = c("gene1", "gene2"),
    TSS = c(100, 200),
    stringsAsFactors = FALSE
  )
  result <- batchLoadTwasWeights(mock_results, meta_df, maxMemoryPerBatch =1000)
  expect_equal(names(result), "allGenes")
  expect_equal(names(result$allGenes), c("gene1", "gene2"))
})

test_that("batchLoadTwasWeights splits when exceeding memory limit", {
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
  result <- batchLoadTwasWeights(mock_results, meta_df, maxMemoryPerBatch =0.0001)
  expect_true(length(result) >= 2)
})

# ===========================================================================
# getFilterLbfIndex
# ===========================================================================

test_that("getFilterLbfIndex returns numeric index vector", {
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

  result <- getFilterLbfIndex(mock_susie, coverage = 0.5, sizeFactor =0.5)
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

test_that("batchLoadTwasWeights returns empty list for empty input (with message)", {
  expect_message(
    result <- batchLoadTwasWeights(list(), data.frame(region_id = character(), TSS = integer())),
    "No genes"
  )
  expect_equal(length(result), 0)
})

test_that("batchLoadTwasWeights returns single batch when total memory fits", {
  twas <- list(
    gene1 = list(a = 1:10),
    gene2 = list(a = 1:10)
  )
  meta <- data.frame(region_id = c("gene1", "gene2"), TSS = c(100, 200))
  expect_message(
    result <- batchLoadTwasWeights(twas, meta, maxMemoryPerBatch =1000),
    "No need to split"
  )
  expect_equal(length(result), 1)
  expect_true(all(c("gene1", "gene2") %in% names(result[[1]])))
})

test_that("batchLoadTwasWeights splits into multiple batches", {
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
  result <- batchLoadTwasWeights(twas, meta, maxMemoryPerBatch =max_mb)
  expect_true(length(result) >= 2)
  # All genes should be present across batches
  allGenes <- unlist(lapply(result, names))
  expect_true(all(c("gene1", "gene2", "gene3") %in% allGenes))
})

test_that("batchLoadTwasWeights puts oversized gene in its own batch", {
  twas <- list(
    gene_small = list(a = 1:10),
    gene_big = rnorm(1e6)
  )
  big_size_mb <- as.numeric(object.size(twas$gene_big)) / (1024^2)
  small_size_mb <- as.numeric(object.size(twas$gene_small)) / (1024^2)
  # Set limit between small and big
  max_mb <- big_size_mb * 0.5
  meta <- data.frame(region_id = c("gene_small", "gene_big"), TSS = c(100, 200))
  result <- batchLoadTwasWeights(twas, meta, maxMemoryPerBatch =max_mb)
  # Big gene should be in its own batch
  expect_true(length(result) >= 2)
})

# ===========================================================================
# loadCovariateData with real fixture
# ===========================================================================

test_that("loadCovariateData reads and transposes covariate file", {
  covar_path <- file.path(test_path("test_data"), "test_covariates.tsv")
  result <- pecotmr:::loadCovariateData(covar_path)
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

test_that("loadCovariateData errors on missing file", {
  expect_error(
    pecotmr:::loadCovariateData("/nonexistent/covariate.tsv"),
    "not found"
  )
})

# ===========================================================================
# loadTsvRegion with real tabix-indexed fixture
# ===========================================================================

test_that("loadTsvRegion reads full gz file without region", {
  skip_if_not_installed("Rsamtools")
  sumstatPath <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- loadTsvRegion(sumstatPath)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 8L)
  expect_true("BETA" %in% names(result) || "beta" %in% names(result))
})

test_that("loadTsvRegion queries region via tabix", {
  skip_if_not_installed("Rsamtools")
  sumstatPath <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  result <- loadTsvRegion(sumstatPath, region = "chr21:17014042-45433269")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 8L)
})

test_that("loadTsvRegion errors for non-overlapping region", {
  skip_if_not_installed("Rsamtools")
  sumstatPath <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  expect_error(loadTsvRegion(sumstatPath, region = "chr1:1-2"), "tabix-indexed")
})

test_that("loadTsvRegion queries subregion correctly", {
  skip_if_not_installed("Rsamtools")
  sumstatPath <- file.path(test_path("test_data"), "test_sumstats.tsv.gz")
  # Only first 2 variants: pos 17014042 and 18759786
  result <- loadTsvRegion(sumstatPath, region = "chr21:17014042-18759786")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2L)
})

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

test_that("standardiseSumstatsColumns renames standard headers", {
  skip_if_not_installed("MungeSumstats")
  df <- data.frame(
    SNPID = "rs1", CHR = 1, POS = 100,
    EFFECT_ALLELE = "A", OTHER_ALLELE = "G",
    BETA = 0.5, SE = 0.1, P = 0.01,
    stringsAsFactors = FALSE
  )
  result <- standardiseSumstatsColumns(df)
  expect_true("chrom" %in% colnames(result))
  expect_true("pos" %in% colnames(result))
  expect_true("beta" %in% colnames(result))
  expect_true("se" %in% colnames(result))
  expect_true("p" %in% colnames(result))
})

test_that("standardiseSumstatsColumns applies custom column mapping", {
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

  result <- standardiseSumstatsColumns(df, columnFilePath =col_file)
  expect_true("maf" %in% colnames(result))
})

test_that("standardiseSumstatsColumns errors on missing column file", {
  skip_if_not_installed("MungeSumstats")
  df <- data.frame(SNP = "rs1", CHR = 1, BP = 100, A1 = "A", A2 = "G",
                   BETA = 0.5, SE = 0.1, P = 0.01, stringsAsFactors = FALSE)
  expect_error(
    standardiseSumstatsColumns(df, columnFilePath ="/no/such/file.txt"),
    "Column mapping file not found"
  )
})

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

test_that("loadPhenotypeData reads compressed file with tabix region", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- loadPhenotypeData(
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

test_that("loadPhenotypeData filters by extract_region_name and region_name_col", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- loadPhenotypeData(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579",
    extractRegionName =list(c("ENSG00000154639")),
    regionNameCol =4
  )
  expect_equal(length(pheno), 1L)
  expect_true("ENSG00000154639" %in% colnames(pheno[[1]]))
})

test_that("loadPhenotypeData assigns gene names with region_name_col", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- loadPhenotypeData(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579",
    regionNameCol =4
  )
  expect_equal(colnames(pheno[[1]]), "ENSG00000154639")
})

test_that("loadPhenotypeData returns multiple genes for broad region", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- loadPhenotypeData(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:14000000-18000000",
    regionNameCol =4
  )
  expect_equal(length(pheno), 1L)
  expect_true(ncol(pheno[[1]]) > 1)
  expect_true("ENSG00000154639" %in% colnames(pheno[[1]]))
})

test_that("loadPhenotypeData errors on non-overlapping region", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  expect_error(
    loadPhenotypeData(
      file.path(td, "test_phenotypes.tsv.gz"),
      region = "chr21:1-100"
    ),
    "empty"
  )
})

test_that("loadPhenotypeData reads uncompressed file without region", {
  td <- test_path("test_data")
  pheno <- loadPhenotypeData(
    file.path(td, "test_phenotypes.tsv"),
    region = NULL,
    regionNameCol =4
  )
  expect_equal(length(pheno), 1L)
  # All 93 genes in the uncompressed file
  expect_equal(ncol(pheno[[1]]), 93L)
})

test_that("loadPhenotypeData stores kept_indices attribute", {
  skip_if_not_installed("Rsamtools")
  td <- test_path("test_data")
  pheno <- loadPhenotypeData(
    file.path(td, "test_phenotypes.tsv.gz"),
    region = "chr21:17513043-17593579"
  )
  expect_equal(attr(pheno, "kept_indices"), 1L)
})

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
#   * loadRegionalMultivariateData   -> use MultiTaskQtlDataset()
#   * loadRegionalFunctionalData     -> use QtlDataset()
#   * loadMultitaskRegionalData      -> use MultiTaskQtlDataset()
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
