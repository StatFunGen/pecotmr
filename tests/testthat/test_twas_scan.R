context("twas_scan")

generate_mock_data <- function(seed=1, num_snps=100, empty_sets = F, gwas_mismatch = F, LD_mismatch = F) {
  set.seed(seed)

  random_weights <- function(num_true_effects = 10, num_snps = num_snps, seed = 1) {
    # Generate true effect sizes for a subset of variants
    true_effects = rnorm(num_true_effects, mean = 0, sd = 1)  # assuming a normal distribution

    # Set zero effects for the remaining variants
    zero_effects = rep(0, num_snps - num_true_effects)

    # Combine and shuffle the effect sizes
    set.seed(seed)
    effect_sizes = sample(c(true_effects, zero_effects))
    noise = rnorm(num_snps, mean = 0, sd = 0.1)
    effect_sizes = effect_sizes + noise
    return(effect_sizes)
  }

  weights_all_matrix <- matrix(
    c(random_weights(num_true_effects = 10, num_snps = num_snps, seed = 1),
    random_weights(num_true_effects = 10, num_snps = num_snps, seed = 2),
    random_weights(num_true_effects = 10, num_snps = num_snps, seed = 3),
    random_weights(num_true_effects = 10, num_snps = num_snps, seed = 4)),
    ncol = 4, nrow=100)
  
  sumstat_A1 <- sample(c("A", "T", "G", "C"), num_snps, replace = TRUE)
  sumstat_A2 <- unlist(lapply(sumstat_A1, function(x) {
    if (x == "A") {
      return(sample(c("G", "C"), 1))
    } else if (x == "T") {
      return(sample(c("G", "C"), 1))
    } else if (x == "G") {
      return(sample(c("A", "T"), 1))
    } else if (x == "C") {
      return(sample(c("A", "T"), 1))
    }
  }))
  gwas_sumstats_db <- data.frame(
    chr = rep(22, num_snps),
    pos = sample(1e6:5e6, num_snps, replace = TRUE),
    A1.sumstats = sumstat_A1,
    A2.sumstats = sumstat_A2,
    beta = rnorm(num_snps),
    se = runif(num_snps, 0.01, 0.1),
    z = rnorm(num_snps)
  ) %>% arrange(chr, pos)
  
  LD_matrix <- matrix(runif(num_snps^2), nrow = num_snps, ncol = num_snps)
  colnames(LD_matrix) <- rownames(LD_matrix) <- variants_id_all <- paste(
    gwas_sumstats_db$chr, gwas_sumstats_db$pos, gwas_sumstats_db$A1.sumstats, gwas_sumstats_db$A2.sumstats, sep = ":")
  
  # Ensure row names of weights_all_matrix match variant IDs
  rownames(weights_all_matrix) <- variants_id_all
   
  extract_variants_objs <- variants_id_all

  gwas_sumstats_db$variant_id <- if (!gwas_mismatch) variants_id_all else paste(
    12, gwas_sumstats_db$pos, gwas_sumstats_db$A1.sumstats, gwas_sumstats_db$A2.sumstats, sep = ":")

  colnames(LD_matrix) <- rownames(LD_matrix) <- if (!LD_mismatch) colnames(LD_matrix) else paste(
    3, gwas_sumstats_db$pos, gwas_sumstats_db$A1.sumstats, gwas_sumstats_db$A2.sumstats, sep = ":")

  return(list(
    weights_all_matrix = weights_all_matrix,
    gwas_sumstats_db = gwas_sumstats_db,
    LD_matrix = LD_matrix,
    extract_variants_objs = extract_variants_objs))
}

mock_weights_db <- function(seed = 1, region = "chr1:100-200", condition = "monocytes", has_variable_names = TRUE, var_row_lengths = FALSE, same_gene = TRUE) {
  if (same_gene) set.seed(1) else set.seed(seed)
  gene <- paste0("ENSG000000000", sample(seq(100,999), 1))
  if (var_row_lengths) {
    if (same_gene) set.seed(seed) else set.seed(1)
    n_variants <- sample(100:150, 1)
    r_variants <- sort(sample(0:200, n_variants))
  } else {
    r_variants <- sort(sample(0:200, 100))
    if (same_gene) set.seed(seed) else set.seed(1)
  }

  # Define mock data
  #weights <- matrix(runif(100), ncol = 10) # Adjust size as needed
  variant_names <- if (has_variable_names) paste0("variant_", r_variants) else NULL
  #colnames(weights) <- rownames(weights) <- variant_names
  susie_result_trimmed <- runif(10) # Example data
  top_loci <- sample(r_variants, 5) # Example data
  region_info <- list(region = region, condition = condition)
  
  # Combine data into a list
  weights_db_data <- list(
    region_info = region_info,
    preset_variants_result = list(
      variant_names = variant_names,
      susie_result_trimmed = susie_result_trimmed,
      top_loci = top_loci
    ),
    twas_weights = list(
      model_one_weights = runif(length(variant_names)),
      model_two_weights = runif(length(variant_names)),
      model_three_weights = runif(length(variant_names)),
      variant_names = variant_names
    ),
    twas_cv_result = list(
      performance=data.frame(corr=0.5, rsq=0.2, adj_rsq=0.2, pval=8e-20, RMSE=0.4, MAE=0.3)
    )
  )
  weights_db <- list(random_gene = list(condition = weights_db_data))
  names(weights_db$random_gene) <- condition
  names(weights_db) <- gene
  
  # Save the list as an RDS file
  return(weights_db)
}

test_that("Confirm twas_scan works with simulated data",{
  data <- generate_mock_data()
  res <- twas_analysis(data$weights_all_matrix, data$gwas_sumstats_db, data$LD_matrix, data$extract_variants_objs)
  expect_equal(length(res), 4)
  expect_true(all(unlist(lapply(res, function(x) {"z" %in% names(x)}))))
  expect_true(all(unlist(lapply(res, function(x) {"pval" %in% names(x)}))))
})

# test_that("twas_analysis raises error if empty gwas",{
#   data <- generate_mock_data(gwas_mismatch = T, LD_mismatch = F)
#   expect_error(
#     twas_analysis(data$weights_all_matrix, data$gwas_sumstats_db, data$LD_matrix, data$extract_variants_objs),
#     "No GWAS summary statistics found for the specified variants.")
# })
#
# test_that("twas_analysis raises error if specified variants are not in LD_matrix", {
#   data <- generate_mock_data(gwas_mismatch = FALSE, LD_mismatch = TRUE)
#   expect_error(
#     twas_analysis(data$weights_matrix, data$gwas_sumstats_db, data$LD_matrix, data$extract_variants_objs),
#     "None of the specified variants are present in the LD matrix."
#   )
# })

setup_weight_db_vector <- function(seed = 1, n_rds = 2, n_cond = 4, condition = NA, same_condition = FALSE, same_gene = TRUE, var_row_lengths = FALSE) {
  set.seed(seed)
  weight_db_vector <- lapply(1:n_rds, function(i) {
    cond <- if (same_condition) condition else gsub(", ", "", toString(sample(LETTERS, 3)))
    mock_weights_db(seed = i, condition = cond, same_gene = same_gene, var_row_lengths = var_row_lengths)
  })
  
  weight_db_paths <- lapply(1:n_rds, function(i) {
    weight_db_path <- gsub("//", "/", tempfile(pattern = paste0("weights_db_", i), tmpdir = tempdir(), fileext = ".RDS"))
    saveRDS(weight_db_vector[[i]], weight_db_path)
    return(weight_db_path)
  })

  return(list(weight_vec = weight_db_vector, weight_paths = weight_db_paths))
}

cleanup_weight_db_vector <- function(weight_db_paths) {
  lapply(weight_db_paths, file.remove)
}

# Test unique regions
test_that("load_twas_weights raises error if different regions specified", {
  weight_db <- setup_weight_db_vector(n_rds = 2, n_cond = 4, same_gene = FALSE)
  expect_true(
    inherits(
      load_twas_weights(weight_db$weight_paths, conditions = NULL, variable_name_obj = c("preset_variants_result", "variant_names"),
                       susie_obj = c("preset_variants_result", "susie_result_trimmed"),twas_weights_table = "twas_weights"),
      "try-error"))
  cleanup_weight_db_vector(weight_db$weight_paths)
})

# Test unique regions
test_that("load_twas_weights raises error if different number of conditions per rds file", {
  weight_db <- setup_weight_db_vector(n_rds = 2, n_cond = 4)
  expect_true(
    inherits(
      load_twas_weights(weight_db$weight_paths, conditions = "not_found", variable_name_obj = c("preset_variants_result", "variant_names"),
                       susie_obj = c("preset_variants_result", "susie_result_trimmed"),twas_weights_table = "twas_weights"),
      "try-error"))
  cleanup_weight_db_vector(weight_db$weight_paths)
})

# Test null conditions
test_that("load_twas_weights works with null condition", {
  weight_db <- setup_weight_db_vector(n_rds = 2, n_cond = 4)
  res <- load_twas_weights(weight_db$weight_paths, conditions = NULL, variable_name_obj = c("preset_variants_result", "variant_names"),
                          susie_obj =c("preset_variants_result", "susie_result_trimmed"),twas_weights_table = "twas_weights")
  expect_true(all(c("susie_results", "weights") %in% names(res)))
  cleanup_weight_db_vector(weight_db$weight_paths)
})

# Specify conditions
test_that("load_twas_weights works with specified condition", {
  weight_db <- setup_weight_db_vector(n_rds = 2, n_cond = 4, same_condition = TRUE, condition = "cond_1_joe_eQTL")
  res <- load_twas_weights(weight_db$weight_paths, conditions = "cond_1_joe_eQTL", variable_name_obj = c("preset_variants_result", "variant_names"),
                          susie_obj =c("preset_variants_result", "susie_result_trimmed"),twas_weights_table = "twas_weights")
  expect_true(all(c("susie_results", "weights") %in% names(res)))
  cleanup_weight_db_vector(weight_db$weight_paths)
})

# # Test null variable_name_obj
# variable_name_obj cannot be NULL
# test_that("load_twas_weights works with null variable_name_obj", {
#   weight_db <- setup_weight_db_vector(n_rds = 2, n_cond = 4)
#   res <- load_twas_weights(weight_db$weight_paths, conditions = NULL, variable_name_obj = NULL)
#   expect_true(all(c("susie_results", "weights") %in% names(res)))
#   cleanup_weight_db_vector(weight_db$weight_paths)
# })

# Test different number of rows per condition
# different number of variants from different context will not affect the loading
test_that("load_twas_weights raises error with null variable_name_obj and variable row lengths", {
  weight_db <- setup_weight_db_vector(n_rds = 2, n_cond = 4, var_row_lengths = TRUE)
  expect_false( # it should not return error 
    inherits(
      load_twas_weights(weight_db$weight_paths, conditions = NULL, variable_name_obj = c("preset_variants_result", "variant_names"),
                       susie_obj = c("preset_variants_result", "susie_result_trimmed"),twas_weights_table = "twas_weights"),
      "try-error"))
  cleanup_weight_db_vector(weight_db$weight_paths)
}) 
