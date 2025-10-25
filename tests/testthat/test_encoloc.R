context("encoloc")
library(tidyverse)
library(coloc)
library(data.table)

generate_mock_ld_files <- function(seed = 1, num_blocks = 5) {
    data(coloc_test_data)
    attach(coloc_test_data)
    set.seed(seed)

    # Generate mock LD files
    blocks <- seq(100, num_blocks*100, 100)
    ld_blocks <- lapply(
        blocks,
        function(i) {
            variants <- paste0("s", (i-100+1):i)
            ld_block <- as.data.frame(D1$LD[variants, variants])
            return(ld_block)
    })

    bim_files <- lapply(
        blocks,
        function(i) {
            bim_df <- data.frame(
                chrom = "chr1",
                id = paste0("chr1:", (i-100+1):i, "_A_G"),
                rand = 0,
                pos = (i-100+1):i,
                ref = "A",
                alt = "G"
            )
            return(bim_df)
        }
    )

    bim_paths <- lapply(blocks, function(i) {
        gsub("//", "/", tempfile(pattern = paste0("LD_block_", i/100, ".chr1_", i-100+1, "_", i, ".float16"), tmpdir = tempdir(), fileext = ".bim"))
    })

    ld_paths <- lapply(blocks, function(i) {
        gsub("//", "/", tempfile(pattern = paste0("LD_block_", i/100, ".chr1_", i-100+1, "_", i, ".float16"), tmpdir = tempdir(), fileext = ".txt.xz"))
    })

    lapply(
        1:num_blocks,
        function(i) {
            xzfile <- xzfile(ld_paths[[i]], "wb")
            write_delim(ld_blocks[[i]], xzfile, delim = "\t", col_names = FALSE)
            close(xzfile)
        })

    lapply(
        1:num_blocks,
        function(i) {
            write_delim(bim_files[[i]], bim_paths[[i]], delim = "\t", col_names = FALSE)
        })
    
    meta_df <- data.frame(
        chrom = "chr1",
        start = seq(1, 401, 100),
        end = seq(100, 500, 100),
        path = unlist(lapply(1:num_blocks, function(i) paste0(ld_paths[[i]], ",", bim_paths[[i]]))))
    
    meta_path <- gsub("//", "/", tempfile(pattern = paste0("ld_meta_file_path"), tmpdir = tempdir(), fileext = ".txt.gz"))
    write_delim(meta_df, meta_path, delim = "\t")

    return(list(ld_paths = ld_paths, bim_paths = bim_paths, meta_path = meta_path))
}

generate_mock_data_for_enrichment <- function(seed=1, num_files=2) {
    gwas_finemapped_data <- as.vector(paste0("gwas_file", 1:num_files))
    xqtl_finemapped_data <- "xqtl_file.rds"
    return(list(gwas_finemapped_data = gwas_finemapped_data,
                xqtl_finemapped_data = xqtl_finemapped_data))
}

generate_mock_susie_fit <- function(seed=1, num_samples = 10, num_features=10, unique_names=T) {
    set.seed(seed)
    start_index <- sample(1:100000, 1)
    alpha_raw <- matrix(runif(num_samples * num_features), nrow = num_samples)
    alpha_normalized <- t(apply(alpha_raw, 1, function(x) x / sum(x)))
    if (unique_names) {
        variant_names <- paste0("chr22:", start_index:(num_features + start_index -1), ":A:C")
    } else {
        variant_names <- rep("chr22:1:A:C", num_features + start_index - 1)
    }
    susie_fit <- list(
        pip = setNames(runif(num_features), paste0("rs", start_index:(num_features+start_index-1))),
        variant_names = variant_names,
        lbf_variable = matrix(runif(num_features * 10), nrow = 10, ncol=num_features),
        alpha = alpha_normalized,
        V = runif(10), 
        prior_variance = runif(num_features))
    return(susie_fit)
}

test_that("xqtl_enrichment_wrapper works with dummy input single threaded",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = NULL, gwas_varname_obj = c("variant_names"),
        xqtl_finemapping_obj = "susie_fit", xqtl_varname_obj = c("susie_fit", "variant_names"),
        num_gwas = 5000, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 1)
    expect_length(res,n = 2)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("xqtl_enrichment_wrapper fails non-unique variant names",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i, unique_names=F)), input_data$gwas_finemapped_data[i])
    }
    expect_error(xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = NULL, gwas_varname_obj = c("variant_names"),
        xqtl_finemapping_obj = "susie_fit", xqtl_varname_obj = c("susie_fit", "variant_names"),
        num_gwas = 5000, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 1))
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})


test_that("xqtl_enrichment_wrapper works with dummy input multi threaded",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = NULL, gwas_varname_obj = c("variant_names"),
        xqtl_finemapping_obj = "susie_fit", xqtl_varname_obj = c("susie_fit", "variant_names"),
        num_gwas = 5000, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 2)
    expect_length(res,n = 2)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("xqtl_enrichment_wrapper works with dummy input single and multi threaded",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res_single <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = NULL, gwas_varname_obj = c("variant_names"),
        xqtl_finemapping_obj = "susie_fit", xqtl_varname_obj = c("susie_fit", "variant_names"),
        num_gwas = 5000, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 1)
    res_multi <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = NULL, gwas_varname_obj = c("variant_names"),
        xqtl_finemapping_obj = "susie_fit", xqtl_varname_obj = c("susie_fit", "variant_names"),
        num_gwas = 5000, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 2)
    expect_equal(res_single, res_multi)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("coloc_wrapper works with dummy input",{
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- coloc_wrapper(input_data$xqtl_finemapped_data, input_data$gwas_finemapped_data, 
                     xqtl_finemapping_obj = "susie_fit", gwas_finemapping_obj=  NULL, 
                     xqtl_varname_obj= c("susie_fit", "variant_names"), gwas_varname_obj = c("variant_names"))
    expect_true(all(names(res) %in% c("summary","results","priors","analysis_region")))
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})


test_that("filter_and_order_coloc_results raises error with insufficient columns",{
    expect_error(filter_and_order_coloc_results(data.frame()))
})

test_that("filter_and_order_coloc_results works with dummy data",{
    data(coloc_test_data)
    attach(coloc_test_data)
    data <- generate_mock_ld_files()
    region <- "chr1:1-500"
    B1 <- D1
    B2 <- D2
    B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
    mock_coloc_results <- coloc.signals(B1, B2, p12 = 1e-5)
    # Mimic the path to using filter_and_order_coloc_results
    coloc_summary <- as.data.frame(mock_coloc_results$summary)
    coloc_pip <- coloc_summary[, grepl("PP", colnames(coloc_summary))]
    PPH4_thres <- 0.8
    coloc_index <- "PP.H4.abf"
    coloc_results_df <- as.data.frame(mock_coloc_results$results)
    coloc_filter <- apply(coloc_pip, 1, function(row) {
        max_index <- which.max(row)
        max_value <- row[max_index]
        return(max_value > PPH4_thres && colnames(coloc_pip)[max_index] == coloc_index)
    })
    coloc_results_fil <- coloc_results_df[, c(1, which(coloc_filter) + 1), drop = FALSE]
    coloc_summary_fil <- coloc_summary[which(coloc_filter),, drop = FALSE]
    ordered_results <- filter_and_order_coloc_results(coloc_results_fil)
    expect_equal(length(ordered_results), 1)
    lapply(unlist(data), function(x) {
        file.remove(x)
    })
})

test_that("calculate_cumsum works with dummy data", {
    data(coloc_test_data)
    attach(coloc_test_data)
    data <- generate_mock_ld_files()
    region <- "chr1:1-500"
    B1 <- D1
    B2 <- D2
    B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
    mock_coloc_results <- coloc.signals(B1, B2, p12 = 1e-5)
    # Mimic the path to using filter_and_order_coloc_results
    coloc_summary <- as.data.frame(mock_coloc_results$summary)
    coloc_pip <- coloc_summary[, grepl("PP", colnames(coloc_summary))]
    PPH4_thres <- 0.8
    coloc_index <- "PP.H4.abf"
    coloc_results_df <- as.data.frame(mock_coloc_results$results)
    coloc_filter <- apply(coloc_pip, 1, function(row) {
        max_index <- which.max(row)
        max_value <- row[max_index]
        return(max_value > PPH4_thres && colnames(coloc_pip)[max_index] == coloc_index)
    })
    coloc_results_fil <- coloc_results_df[, c(1, which(coloc_filter) + 1), drop = FALSE]
    coloc_summary_fil <- coloc_summary[which(coloc_filter),, drop = FALSE]
    ordered_results <- filter_and_order_coloc_results(coloc_results_fil)
    cs <- list()
    
    for (n in 1:length(ordered_results)) {
      tmp_coloc_results_fil <- ordered_results[[n]]
      tmp_coloc_results_fil_csm <- calculate_cumsum(tmp_coloc_results_fil)
      expect_equal(tmp_coloc_results_fil_csm, cumsum(tmp_coloc_results_fil[,2]))  
    }
    lapply(unlist(data), function(x) {
        file.remove(x)
    })
})

# test_that("load_and_extract_ld_matrix works with dummy data", {
#     data(coloc_test_data)
#     attach(coloc_test_data)
#     data <- generate_mock_ld_files()
#     region <- "chr1:1-5"
#     B1 <- D1
#     B2 <- D2
#     B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
#     variants <- paste0("1:", 1:5, ":A:G")
#     res <- load_and_extract_ld_matrix(data$meta_path, region, variants)
#     expect_equal(nrow(res), 5)
#     expect_equal(ncol(res), 5)
#     lapply(unlist(data), function(x) {
#         file.remove(x)
#     })
# })

# test_that("calculate_purity works with dummy data", {
#     data(coloc_test_data)
#     attach(coloc_test_data)
#     data <- generate_mock_ld_files()
#     region <- "chr1:1-5"
#     B1 <- D1
#     B2 <- D2
#     B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
#     variants <- paste0("1:", 1:5, ":A:G")
#     ext_ld <- load_and_extract_ld_matrix(data$meta_path, region, variants)
#     res <- calculate_purity(variants, ext_ld, squared = TRUE)
#     expect_equal(ncol(res), 3)
#     lapply(unlist(data), function(x) {
#         file.remove(x)
#     })
# })

# test_that("process_coloc_results works with dummy data", {
#     data(coloc_test_data)
#     attach(coloc_test_data)
#     data <- generate_mock_ld_files()
#     region <- "chr1:1-500"
#     B1 <- D1
#     B2 <- D2
#     B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
#     mock_coloc_results <- coloc.signals(B1, B2, p12 = 1e-5)
#     res <- process_coloc_results(mock_coloc_results, data$meta_path, region)
#     expect_equal(length(res$sets$cs), 1)
#     lapply(unlist(data), function(x) {
#         file.remove(x)
#     })
# })
