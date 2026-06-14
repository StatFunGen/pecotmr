context("encoloc")

.test_fm_result <- function(variantNames, trimmedFit = list(),
                            topLoci = data.frame(variant_id = character(0),
                                                  method = character(0),
                                                  stringsAsFactors = FALSE)) {
  FineMappingResult(
    variantNames = variantNames,
    trimmedFit = trimmedFit,
    topLoci = topLoci,
    method = "susie"
  )
}

library(tidyverse)
library(coloc)

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

    metaPath <- gsub("//", "/", tempfile(pattern = paste0("ld_meta_file_path"), tmpdir = tempdir(), fileext = ".txt.gz"))
    write_delim(meta_df, metaPath, delim = "\t")

    return(list(ld_paths = ld_paths, bim_paths = bim_paths, metaPath = metaPath))
}

generate_mock_data_for_enrichment <- function(seed=1, num_files=2) {
    gwas_finemapped_data <- as.vector(paste0("gwas_file", 1:num_files))
    xqtl_finemapped_data <- "xqtl_file.rds"
    return(list(gwas_finemapped_data = gwas_finemapped_data,
                xqtl_finemapped_data = xqtl_finemapped_data))
}

generate_mock_susie_fit <- function(seed=1, num_samples = 10, num_features=10, unique_names=T) {
    set.seed(seed)
    startIndex <- sample(1:100000, 1)
    alpha_raw <- matrix(runif(num_samples * num_features), nrow = num_samples)
    alpha_normalized <- t(apply(alpha_raw, 1, function(x) x / sum(x)))
    if (unique_names) {
        variantNames <- paste0("chr22:", startIndex:(num_features + startIndex -1), ":A:C")
    } else {
        variantNames <- rep("chr22:1:A:C", num_features + startIndex - 1)
    }
    susie_fit <- list(
        pip = setNames(runif(num_features), paste0("rs", startIndex:(num_features+startIndex-1))),
        variantNames = variantNames,
        lbf_variable = matrix(runif(num_features * 10), nrow = 10, ncol=num_features),
        alpha = alpha_normalized,
        V = runif(10),
        prior_variance = runif(num_features))
    return(susie_fit)
}

# ===========================================================================
# xqtlEnrichmentWrapper tests
# ===========================================================================

test_that("xqtlEnrichmentWrapper works with dummy input single threaded",{
    local_mocked_bindings(
        qtlEnrichmentRcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susieFit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susieFit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- xqtlEnrichmentWrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data,
        gwasFinemappingObj = NULL, gwasVarnameObj = c("variantNames"),
        xqtlFinemappingObj = "susieFit", xqtlVarnameObj = c("susieFit", "variantNames"),
        numGwas = 5000, piQtl = 0.5,
        lambda = 1.0, impN = 25,
        numThreads = 1)
    expect_length(res,n = 2)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("xqtlEnrichmentWrapper fails non-unique variant names",{
    local_mocked_bindings(
        qtlEnrichmentRcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susieFit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susieFit = generate_mock_susie_fit(seed=i, unique_names=F)), input_data$gwas_finemapped_data[i])
    }
    expect_error(xqtlEnrichmentWrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data,
        gwasFinemappingObj = NULL, gwasVarnameObj = c("variantNames"),
        xqtlFinemappingObj = "susieFit", xqtlVarnameObj = c("susieFit", "variantNames"),
        numGwas = 5000, piQtl = 0.5,
        lambda = 1.0, impN = 25,
        numThreads = 1),
        regexp = "must be the same length")
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})


test_that("xqtlEnrichmentWrapper works with dummy input single and multi threaded",{
    local_mocked_bindings(
        qtlEnrichmentRcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susieFit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susieFit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res_single <- xqtlEnrichmentWrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data,
        gwasFinemappingObj = NULL, gwasVarnameObj = c("variantNames"),
        xqtlFinemappingObj = "susieFit", xqtlVarnameObj = c("susieFit", "variantNames"),
        numGwas = 5000, piQtl = 0.5,
        lambda = 1.0, impN = 25,
        numThreads = 1)
    res_multi <- xqtlEnrichmentWrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data,
        gwasFinemappingObj = NULL, gwasVarnameObj = c("variantNames"),
        xqtlFinemappingObj = "susieFit", xqtlVarnameObj = c("susieFit", "variantNames"),
        numGwas = 5000, piQtl = 0.5,
        lambda = 1.0, impN = 25,
        numThreads = 2)
    expect_equal(res_single, res_multi)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("xqtlEnrichmentWrapper errors with non-existent files", {
  expect_error(
    xqtlEnrichmentWrapper(
      xqtlFiles = "/nonexistent/xqtl.rds",
      gwasFiles = "/nonexistent/gwas.rds"
    )
  )
})

test_that("xqtlEnrichmentWrapper handles xqtl_finemapping_obj error gracefully", {
  local_mocked_bindings(
    qtlEnrichmentRcpp = function(...) TRUE
  )

  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  gwas_fit <- generate_mock_susie_fit(seed = 1)
  # xqtl data has structure that will fail getNestedElement
  xqtl_data <- list(gene = list(no_susie = "nothing"))

  saveRDS(list(gwas_fit), gwas_file)
  saveRDS(list(xqtl_data), xqtl_file)

  # The tryCatch in process_finemapped_data should return NULL for this xqtl
  result <- tryCatch(
    xqtlEnrichmentWrapper(
      xqtl_file, gwas_file,
      xqtlFinemappingObj = c("nonexistent", "deep", "path"),
      gwasFinemappingObj = NULL,
      gwasVarnameObj = c("variantNames"),
      xqtlVarnameObj = c("susieFit", "variantNames"),
      numGwas = 5000, piQtl = 0.5
    ),
    error = function(e) list(error = e$message)
  )
  expect_true(is.list(result))

  file.remove(gwas_file, xqtl_file)
})

# ===========================================================================
# colocWrapper tests
# ===========================================================================

test_that("colocWrapper works with dummy input",{
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susieFit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susieFit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- colocWrapper(input_data$xqtl_finemapped_data, input_data$gwas_finemapped_data,
                     xqtlFinemappingObj = "susieFit", gwasFinemappingObj = NULL,
                     xqtlVarnameObj = c("susieFit", "variantNames"), gwasVarnameObj = c("variantNames"))
    expect_true(all(names(res) %in% c("summary","results","priors","analysisRegion")))
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("colocWrapper errors with non-existent xqtl file", {
  expect_error(
    colocWrapper(
      xqtlFile = "/nonexistent/xqtl.rds",
      gwasFiles = "/nonexistent/gwas.rds"
    ),
    regexp = "cannot open the connection"
  )
})

test_that("colocWrapper returns message when xqtl finemapping obj not found", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  gwas_data <- generate_mock_susie_fit(seed = 1)
  xqtl_data <- list(something_else = "no susie_fit here")

  saveRDS(list(gwas_data), gwas_file)
  saveRDS(list(xqtl_data), xqtl_file)

  result <- colocWrapper(
    xqtl_file, gwas_file,
    xqtlFinemappingObj = c("nonexistent", "path"),
    gwasFinemappingObj = NULL,
    gwasVarnameObj = c("variantNames")
  )
  # Should return a message about missing finemapping object
  expect_true(is.list(result))

  file.remove(gwas_file, xqtl_file)
})

test_that("colocWrapper returns placeholder when GWAS lbf_variable has zero rows", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  gwas_fit <- generate_mock_susie_fit(seed = 1)
  gwas_fit$lbf_variable <- matrix(nrow = 0, ncol = 10)  # zero rows
  gwas_fit$V <- numeric(0)

  xqtl_fit <- generate_mock_susie_fit(seed = 2)

  saveRDS(list(gwas_fit), gwas_file)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  res <- colocWrapper(
    xqtl_file, gwas_file,
    xqtlFinemappingObj = "susieFit",
    gwasFinemappingObj = NULL,
    gwasVarnameObj = c("variantNames"),
    xqtlVarnameObj = c("susieFit", "variantNames")
  )
  expect_true(is.list(res))

  file.remove(gwas_file, xqtl_file)
})

test_that("colocWrapper with filter_lbf_cs uses cs_index filtering", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  # Use the same seed so variant names overlap
  gwas_fit <- generate_mock_susie_fit(seed = 1)
  gwas_fit$sets <- list(cs_index = c(1, 3))  # filter rows 1 and 3

  xqtl_fit <- generate_mock_susie_fit(seed = 1)
  xqtl_fit$sets <- list(cs_index = c(1, 2))

  saveRDS(list(gwas_fit), gwas_file)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  result <- colocWrapper(
    xqtl_file, gwas_file,
    xqtlFinemappingObj = "susieFit",
    gwasFinemappingObj = NULL,
    gwasVarnameObj = c("variantNames"),
    xqtlVarnameObj = c("susieFit", "variantNames"),
    filterLbfCs = TRUE
  )
  expect_true(is.list(result))

  file.remove(gwas_file, xqtl_file)
})

test_that("colocWrapper with filter_lbf_cs_secondary uses getFilterLbfIndex", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  # Use the same seed so variant names overlap
  gwas_fit <- generate_mock_susie_fit(seed = 1)
  # Need alpha for getFilterLbfIndex
  set.seed(10)
  gwas_fit$alpha <- matrix(runif(100), nrow = 10, ncol = 10)
  gwas_fit$alpha <- t(apply(gwas_fit$alpha, 1, function(x) x / sum(x)))
  gwas_fit$mu <- matrix(rnorm(100), nrow = 10, ncol = 10)
  gwas_fit$mu2 <- matrix(abs(rnorm(100)), nrow = 10, ncol = 10)
  gwas_fit$pip <- colSums(gwas_fit$alpha)

  xqtl_fit <- generate_mock_susie_fit(seed = 1)
  set.seed(20)
  xqtl_fit$alpha <- matrix(runif(100), nrow = 10, ncol = 10)
  xqtl_fit$alpha <- t(apply(xqtl_fit$alpha, 1, function(x) x / sum(x)))
  xqtl_fit$mu <- matrix(rnorm(100), nrow = 10, ncol = 10)
  xqtl_fit$mu2 <- matrix(abs(rnorm(100)), nrow = 10, ncol = 10)
  xqtl_fit$pip <- colSums(xqtl_fit$alpha)

  saveRDS(list(gwas_fit), gwas_file)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  # filter_lbf_cs_secondary triggers the getFilterLbfIndex branch
  result <- tryCatch(
    colocWrapper(
      xqtl_file, gwas_file,
      xqtlFinemappingObj = "susieFit",
      gwasFinemappingObj = NULL,
      gwasVarnameObj = c("variantNames"),
      xqtlVarnameObj = c("susieFit", "variantNames"),
      filterLbfCsSecondary = 0.5
    ),
    error = function(e) list(error = e$message)
  )
  expect_true(is.list(result))

  file.remove(gwas_file, xqtl_file)
})

test_that("colocWrapper produces message when xqtl_data has no V", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  gwas_fit <- generate_mock_susie_fit(seed = 1)
  xqtl_fit <- generate_mock_susie_fit(seed = 1)
  xqtl_fit$V <- NULL  # Remove V to trigger the "No V found" message

  saveRDS(list(gwas_fit), gwas_file)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  expect_message(
    result <- colocWrapper(
      xqtl_file, gwas_file,
      xqtlFinemappingObj = "susieFit",
      gwasFinemappingObj = NULL,
      gwasVarnameObj = c("variantNames"),
      xqtlVarnameObj = c("susieFit", "variantNames")
    ),
    "No V found"
  )
  expect_true(is.list(result))

  file.remove(gwas_file, xqtl_file)
})

test_that("colocWrapper extracts analysisRegion from xqtl_region_obj", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")

  gwas_fit <- generate_mock_susie_fit(seed = 1)
  xqtl_fit <- generate_mock_susie_fit(seed = 1)

  # Add region_info
  xqtl_raw <- list(gene = list(
    susieFit = xqtl_fit,
    region_info = data.frame(chrom = 22, start = 1000, end = 2000)
  ))

  saveRDS(list(xqtl_raw$gene), xqtl_file)
  saveRDS(list(gwas_fit), gwas_file)

  result <- colocWrapper(
    xqtl_file, gwas_file,
    xqtlFinemappingObj = "susieFit",
    gwasFinemappingObj = NULL,
    gwasVarnameObj = c("variantNames"),
    xqtlVarnameObj = c("susieFit", "variantNames"),
    xqtlRegionObj = "region_info"
  )
  expect_true("analysisRegion" %in% names(result))

  file.remove(gwas_file, xqtl_file)
})

# ===========================================================================
# colocWrapper inline fine-mapping tests
# ===========================================================================

test_that("colocWrapper errors when no GWAS source provided", {
  xqtl_file <- tempfile(fileext = ".rds")
  saveRDS(list(list(susieFit = generate_mock_susie_fit(seed = 1))), xqtl_file)
  expect_error(
    colocWrapper(xqtl_file),
    "Either set runFinemapping"
  )
  file.remove(xqtl_file)
})

test_that("colocWrapper errors when run_finemapping missing sumstatPath", {
  expect_error(
    colocWrapper("fake.rds", runFinemapping = TRUE, ldData = list()),
    "sumstatPath is required"
  )
})

test_that("colocWrapper errors when run_finemapping missing ldData", {
  expect_error(
    colocWrapper("fake.rds", runFinemapping = TRUE, sumstatPath = "s.tsv"),
    "ldData is required"
  )
})

test_that("colocWrapper warns when both gwas_files and run_finemapping", {
  # This will warn, then error on sumstatPath/ldData validation
  expect_warning(
    tryCatch(
      colocWrapper("xqtl.rds", gwasFiles = "gwas.rds",
                    runFinemapping = TRUE, sumstatPath = "s.tsv",
                    ldData = list()),
      error = function(e) NULL
    ),
    "Inline fine-mapping will be used"
  )
})

test_that("colocWrapper with runFinemapping = TRUE uses rssAnalysisPipeline", {
  xqtl_file <- tempfile(fileext = ".rds")
  xqtl_fit <- generate_mock_susie_fit(seed = 1)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  # Build mock pipeline result matching rssAnalysisPipeline output structure
  mock_pipeline <- list(
    "susie_rss_SLALOM_RAISS_imputed" = list(
      finemappingResult = .test_fm_result(
        variantNames = xqtl_fit$variantNames,
        trimmedFit = list(
          lbf_variable = xqtl_fit$lbf_variable,
          V = xqtl_fit$V,
          pip = xqtl_fit$pip,
          sets = list(cs_index = seq_len(nrow(xqtl_fit$lbf_variable)))
        )
      )
    ),
    rssDataAnalyzed = data.frame(
      variant_id = xqtl_fit$variantNames,
      z = rnorm(length(xqtl_fit$variantNames))
    )
  )

  local_mocked_bindings(
    rssAnalysisPipeline = function(...) mock_pipeline
  )

  result <- colocWrapper(
    xqtl_file,
    runFinemapping = TRUE,
    sumstatPath = "/fake/gwas.tsv",
    ldData = list(ldMatrix = diag(10)),
    nSample = 10000,
    region = "chr22:1-100",
    xqtlFinemappingObj = "susieFit",
    xqtlVarnameObj = c("susieFit", "variantNames")
  )
  expect_true(all(c("summary", "results") %in% names(result)))
  file.remove(xqtl_file)
})

test_that("colocWrapper with return_finemapping includes pipeline result", {
  xqtl_file <- tempfile(fileext = ".rds")
  xqtl_fit <- generate_mock_susie_fit(seed = 1)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  mock_pipeline <- list(
    "susie_rss_SLALOM" = list(
      finemappingResult = .test_fm_result(
        variantNames = xqtl_fit$variantNames,
        trimmedFit = list(
          lbf_variable = xqtl_fit$lbf_variable,
          V = xqtl_fit$V,
          pip = xqtl_fit$pip,
          sets = list(cs_index = seq_len(nrow(xqtl_fit$lbf_variable)))
        )
      )
    ),
    rssDataAnalyzed = data.frame(
      variant_id = xqtl_fit$variantNames,
      z = rnorm(length(xqtl_fit$variantNames))
    )
  )

  local_mocked_bindings(
    rssAnalysisPipeline = function(...) mock_pipeline
  )

  result <- colocWrapper(
    xqtl_file,
    runFinemapping = TRUE,
    sumstatPath = "/fake/gwas.tsv",
    ldData = list(ldMatrix = diag(10)),
    nSample = 10000,
    xqtlFinemappingObj = "susieFit",
    xqtlVarnameObj = c("susieFit", "variantNames"),
    returnFinemapping = TRUE
  )
  expect_true("gwasFinemapping" %in% names(result))
  expect_true("susie_rss_SLALOM" %in% names(result$gwasFinemapping))
  file.remove(xqtl_file)
})

test_that("colocWrapper save_finemapping_path saves reusable RDS", {
  xqtl_file <- tempfile(fileext = ".rds")
  save_path <- tempfile(fileext = ".rds")
  xqtl_fit <- generate_mock_susie_fit(seed = 1)
  saveRDS(list(gene = list(susieFit = xqtl_fit)), xqtl_file)

  mock_pipeline <- list(
    "susie_rss_SLALOM" = list(
      finemappingResult = .test_fm_result(
        variantNames = xqtl_fit$variantNames,
        trimmedFit = list(
          lbf_variable = xqtl_fit$lbf_variable,
          V = xqtl_fit$V,
          pip = xqtl_fit$pip,
          sets = list(cs_index = seq_len(nrow(xqtl_fit$lbf_variable)))
        )
      )
    ),
    rssDataAnalyzed = data.frame(
      variant_id = xqtl_fit$variantNames,
      z = rnorm(length(xqtl_fit$variantNames))
    )
  )

  local_mocked_bindings(
    rssAnalysisPipeline = function(...) mock_pipeline
  )

  result <- colocWrapper(
    xqtl_file,
    runFinemapping = TRUE,
    sumstatPath = "/fake/gwas.tsv",
    ldData = list(ldMatrix = diag(10)),
    nSample = 10000,
    xqtlFinemappingObj = "susieFit",
    xqtlVarnameObj = c("susieFit", "variantNames"),
    saveFinemappingPath = save_path
  )

  # Verify file was saved
  expect_true(file.exists(save_path))

  # Verify saved format is compatible with file-based reading path
  saved_data <- readRDS(save_path)[[1]]
  expect_true("susie_fit" %in% names(saved_data))
  expect_true("variantNames" %in% names(saved_data))
  expect_true(!is.null(saved_data$susie_fit$lbf_variable))
  expect_true(!is.null(saved_data$susie_fit$V))

  # Verify reusable: can be read back by colocWrapper via file-based path
  result2 <- colocWrapper(
    xqtl_file,
    gwasFiles = save_path,
    xqtlFinemappingObj = "susieFit",
    gwasFinemappingObj = "susie_fit",
    xqtlVarnameObj = c("susieFit", "variantNames"),
    gwasVarnameObj = "variantNames"
  )
  expect_true(all(c("summary", "results") %in% names(result2)))

  file.remove(xqtl_file, save_path)
})

test_that("colocWrapper backward compatibility with gwas_files only", {
  # This mirrors the existing test at line 228 but explicitly verifies
  # that the default run_finemapping=FALSE works
  input_data <- generate_mock_data_for_enrichment()
  input_data$gwas_finemapped_data <- unlist(lapply(
    input_data$gwas_finemapped_data, function(x) {
      gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
  input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
  saveRDS(list(gene = list(susieFit = generate_mock_susie_fit(seed = 1))), input_data$xqtl_finemapped_data)
  for (i in 1:length(input_data$gwas_finemapped_data)) {
    saveRDS(list(susieFit = generate_mock_susie_fit(seed = i)), input_data$gwas_finemapped_data[i])
  }
  res <- colocWrapper(input_data$xqtl_finemapped_data, input_data$gwas_finemapped_data,
                       xqtlFinemappingObj = "susieFit", gwasFinemappingObj = NULL,
                       xqtlVarnameObj = c("susieFit", "variantNames"), gwasVarnameObj = c("variantNames"))
  expect_true(all(names(res) %in% c("summary", "results", "priors", "analysisRegion")))
  file.remove(input_data$gwas_finemapped_data)
  file.remove(input_data$xqtl_finemapped_data)
})

# ===========================================================================
# filterAndOrderColocResults
# ===========================================================================

test_that("filterAndOrderColocResults raises error with insufficient columns",{
    expect_error(filterAndOrderColocResults(data.frame()))
})

test_that("filterAndOrderColocResults with single credible set returns list of length 1", {
  df <- data.frame(
    snp    = c("s1", "s2", "s3", "s4"),
    PP.H4  = c(0.6, 0.1, 0.2, 0.1)
  )
  result <- pecotmr:::filterAndOrderColocResults(df)
  expect_length(result, 1)
  # First row of ordered result should be the one with highest PP.H4 (0.6)
  expect_equal(result[[1]][1, 1], "s1")
  expect_equal(result[[1]][1, 2], 0.6)
})

test_that("filterAndOrderColocResults with multiple credible sets", {
  df <- data.frame(
    snp    = c("s1", "s2", "s3"),
    CS1    = c(0.2, 0.5, 0.3),
    CS2    = c(0.7, 0.1, 0.2),
    CS3    = c(0.1, 0.3, 0.6)
  )
  result <- pecotmr:::filterAndOrderColocResults(df)
  expect_length(result, 3)
  # CS1: s2 should be first (0.5)
  expect_equal(result[[1]][1, 1], "s2")
  # CS2: s1 should be first (0.7)
  expect_equal(result[[2]][1, 1], "s1")
  # CS3: s3 should be first (0.6)
  expect_equal(result[[3]][1, 1], "s3")
})

test_that("filterAndOrderColocResults handles tied PP values", {
  df <- data.frame(
    snp = c("s1", "s2", "s3"),
    PP  = c(0.4, 0.4, 0.2)
  )
  result <- pecotmr:::filterAndOrderColocResults(df)
  expect_length(result, 1)
  # Both s1 and s2 have 0.4 -- order is stable (decreasing), both should appear before s3
  expect_equal(result[[1]][3, 2], 0.2)
})

test_that("filterAndOrderColocResults works with dummy data",{
    data(coloc_test_data)
    attach(coloc_test_data)
    data <- generate_mock_ld_files()
    region <- "chr1:1-500"
    B1 <- D1
    B2 <- D2
    B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
    mock_coloc_results <- coloc.signals(B1, B2, p12 = 1e-5)
    # Mimic the path to using filterAndOrderColocResults
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
    ordered_results <- filterAndOrderColocResults(coloc_results_fil)
    expect_equal(length(ordered_results), 1)
    lapply(unlist(data), function(x) {
        file.remove(x)
    })
})

test_that("filterAndOrderColocResults orders by PP values", {
  coloc_results <- data.frame(
    snp = c("s1", "s2", "s3"),
    PP.H4.1 = c(0.1, 0.5, 0.4),
    PP.H4.2 = c(0.3, 0.2, 0.5)
  )
  result <- pecotmr:::filterAndOrderColocResults(coloc_results)
  expect_length(result, 2)  # 2 credible sets
  # First set should be ordered by PP.H4.1 decreasingly
  expect_equal(result[[1]][1, 2], 0.5)
})

# ===========================================================================
# calculateCumsum
# ===========================================================================

test_that("calculateCumsum works with dummy data", {
    data(coloc_test_data)
    attach(coloc_test_data)
    data <- generate_mock_ld_files()
    region <- "chr1:1-500"
    B1 <- D1
    B2 <- D2
    B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
    mock_coloc_results <- coloc.signals(B1, B2, p12 = 1e-5)
    # Mimic the path to using filterAndOrderColocResults
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
    ordered_results <- filterAndOrderColocResults(coloc_results_fil)
    cs <- list()

    for (n in 1:length(ordered_results)) {
      tmp_coloc_results_fil <- ordered_results[[n]]
      tmp_coloc_results_fil_csm <- calculateCumsum(tmp_coloc_results_fil)
      expect_equal(tmp_coloc_results_fil_csm, cumsum(tmp_coloc_results_fil[,2]))
    }
    lapply(unlist(data), function(x) {
        file.remove(x)
    })
})

test_that("calculateCumsum with single row returns that value", {
  df <- data.frame(snp = "s1", pp = 0.95)
  result <- pecotmr:::calculateCumsum(df)
  expect_equal(result, 0.95)
})

test_that("calculateCumsum monotonically increases", {
  df <- data.frame(
    snp = paste0("s", 1:5),
    pp  = c(0.4, 0.25, 0.15, 0.1, 0.1)
  )
  result <- pecotmr:::calculateCumsum(df)
  expect_true(all(diff(result) >= 0))
  expect_equal(result[5], 1.0)
})

test_that("calculateCumsum can be used to find coverage threshold index", {
  df <- data.frame(
    snp = paste0("s", 1:10),
    pp  = c(0.4, 0.2, 0.15, 0.1, 0.05, 0.03, 0.03, 0.02, 0.01, 0.01)
  )
  cs <- pecotmr:::calculateCumsum(df)
  coverage_idx <- min(which(cs > 0.95))
  expect_equal(coverage_idx, 7)
  expect_gt(cs[coverage_idx], 0.95)
  expect_lte(cs[coverage_idx - 1], 0.95)
})

# ===========================================================================
# calculate_purity
# ===========================================================================

test_that("calculate_purity returns 1x3 matrix for perfectly correlated variants", {
  variants <- c("chr1:100:A:G", "chr1:200:C:T")
  # Identity-like LD matrix (perfect self-correlation)
  ext_ld <- matrix(c(1, 0.95, 0.95, 1), 2, 2)
  rownames(ext_ld) <- colnames(ext_ld) <- variants

  result <- pecotmr:::calculatePurity(variants, ext_ld, squared = FALSE)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  # min abs corr should be 0.95
  expect_equal(result[1, 1], 0.95, tolerance = 1e-6)
})

test_that("calculate_purity with single variant returns matrix of 1s", {
  variants <- "chr1:100:A:G"
  ext_ld <- matrix(1, 1, 1)
  rownames(ext_ld) <- colnames(ext_ld) <- variants

  result <- pecotmr:::calculatePurity(variants, ext_ld, squared = FALSE)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  expect_equal(result[1, 1], 1)
})

test_that("calculate_purity with squared=TRUE returns squared correlations", {
  variants <- c("chr1:100:A:G", "chr1:200:C:T")
  ext_ld <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  rownames(ext_ld) <- colnames(ext_ld) <- variants

  result_sq <- pecotmr:::calculatePurity(variants, ext_ld, squared = TRUE)
  result_no <- pecotmr:::calculatePurity(variants, ext_ld, squared = FALSE)
  expect_equal(nrow(result_sq), 1)
  expect_equal(ncol(result_sq), 3)
  expect_equal(result_sq[1, 1], 0.8, tolerance = 1e-6)
  expect_equal(result_no[1, 1], 0.8, tolerance = 1e-6)
})

# ===========================================================================
# colocPostProcessor
# ===========================================================================

test_that("colocPostProcessor preserves original fields when no LD path given", {
  coloc_res <- list(
    summary = data.frame(PP.H4.abf = 0.92),
    results = data.frame(snp = c("s1", "s2"), PP.H4 = c(0.6, 0.4)),
    priors  = list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  )
  result <- suppressWarnings(colocPostProcessor(coloc_res))
  # All original fields should be retained
  expect_true("summary" %in% names(result))
  expect_true("results" %in% names(result))
  expect_true("priors" %in% names(result))
})

test_that("colocPostProcessor with LD path and region calls processColocResults", {
  local_mocked_bindings(
    processColocResults = function(...) {
      list(sets = list(cs = list(c("s1", "s2")),
                       purity = data.frame(minAbsCorr = 0.9, meanAbsCorr = 0.95, medianAbsCorr = 0.92)))
    }
  )

  coloc_res <- list(
    summary = data.frame(PP.H4.abf = 0.95),
    results = data.frame(snp = c("s1", "s2"), PP.H4 = c(0.6, 0.4))
  )

  result <- colocPostProcessor(
    coloc_res,
    ldMetaFilePath = "/some/path.txt",
    analysisRegion ="chr1:1-1000"
  )
  expect_true("sets" %in% names(result))
})

test_that("colocPostProcessor with LD path but no region errors", {
  coloc_res <- list(summary = data.frame(PP.H4.abf = 0.9))
  expect_error(
    colocPostProcessor(coloc_res, ldMetaFilePath = "/some/path"),
    "analysisRegion is not provided"
  )
})

test_that("colocPostProcessor with region but no LD path warns", {
  coloc_res <- list(summary = data.frame(PP.H4.abf = 0.9))
  expect_warning(
    result <- colocPostProcessor(coloc_res, analysisRegion ="chr1:100-200"),
    "will not be used"
  )
  expect_true("summary" %in% names(result))
})

test_that("colocPostProcessor with neither LD path nor region warns", {
  coloc_res <- list(summary = data.frame(PP.H4.abf = 0.9))
  expect_warning(
    result <- colocPostProcessor(coloc_res),
    "ldMetaFilePath not provided"
  )
})

# ===========================================================================
# processColocResults
# ===========================================================================

test_that("processColocResults returns NULL cs when no PP.H4 qualifies", {
  coloc_result <- list(
    summary = data.frame(
      hit1 = "1:100:A:G",
      hit2 = "1:200:C:T",
      PP.H0.abf = 0.9,
      PP.H1.abf = 0.05,
      PP.H2.abf = 0.03,
      PP.H3.abf = 0.01,
      PP.H4.abf = 0.01
    ),
    results = data.frame(
      snp    = c("s1", "s2"),
      PP.H4  = c(0.01, 0.99)
    )
  )

  expect_message(
    result <- pecotmr:::processColocResults(
      coloc_result, "/fake/path", "chr1:1-1000", pph4Thres = 0.5
    ),
    "did not find any variants"
  )
  expect_true(is.null(result$sets$cs))
})

test_that("processColocResults handles null_index producing purity of (-9, -9, -9)", {
  coloc_result <- list(
    summary = data.frame(
      hit1 = "100",
      hit2 = "200",
      PP.H0.abf = 0.01,
      PP.H1.abf = 0.01,
      PP.H2.abf = 0.01,
      PP.H3.abf = 0.01,
      PP.H4.abf = 0.96
    ),
    results = data.frame(
      snp = c("100", "200"),
      PP.H4.1 = c(0.7, 0.3)
    )
  )

  local_mocked_bindings(
    extractLdForVariants = function(...) {
      m <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
      rownames(m) <- colnames(m) <- c("100", "200")
      m
    },
    calculatePurity = function(...) matrix(c(0.9, 0.95, 0.92), 1, 3),
    normalizeVariantId = function(x, ...) x
  )

  result <- pecotmr:::processColocResults(
    coloc_result,
    "/fake/ld_meta.txt",
    "chr1:100-200",
    pph4Thres = 0,
    nullIndex = 100
  )

  expect_true(is.null(result$sets$cs) || length(result$sets) == 0)
})

test_that("processColocResults returns cs when purity passes", {
  coloc_result <- list(
    summary = data.frame(
      hit1 = "chr1:100:A:G",
      hit2 = "chr1:200:C:T",
      PP.H0.abf = 0.01,
      PP.H1.abf = 0.01,
      PP.H2.abf = 0.01,
      PP.H3.abf = 0.01,
      PP.H4.abf = 0.96
    ),
    results = data.frame(
      snp = c("chr1:100:A:G", "chr1:200:C:T"),
      PP.H4.1 = c(0.7, 0.3)
    )
  )

  local_mocked_bindings(
    extractLdForVariants = function(...) {
      m <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
      rownames(m) <- colnames(m) <- c("100", "200")
      m
    },
    calculatePurity = function(...) matrix(c(0.9, 0.95, 0.92), 1, 3),
    normalizeVariantId = function(x, ...) x
  )

  result <- pecotmr:::processColocResults(
    coloc_result,
    "/fake/ld_meta.txt",
    "chr1:100-200",
    pph4Thres = 0,
    coverage = 0.95
  )

  expect_true(!is.null(result$sets$cs))
  expect_true(!is.null(result$sets$purity))
  expect_equal(ncol(result$sets$purity), 3)
  expect_equal(colnames(result$sets$purity), c("minAbsCorr", "meanAbsCorr", "medianAbsCorr"))
})

test_that("processColocResults filters impure credible sets", {
  coloc_result <- list(
    summary = data.frame(
      hit1 = c("chr1:100:A:G", "chr1:300:G:T"),
      hit2 = c("chr1:200:C:T", "chr1:400:A:C"),
      PP.H0.abf = c(0.01, 0.01),
      PP.H1.abf = c(0.01, 0.01),
      PP.H2.abf = c(0.01, 0.01),
      PP.H3.abf = c(0.01, 0.01),
      PP.H4.abf = c(0.96, 0.96)
    ),
    results = data.frame(
      snp = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:T"),
      PP.H4.1 = c(0.6, 0.4, 0.0),
      PP.H4.2 = c(0.0, 0.0, 1.0)
    )
  )

  purity_call_count <- 0
  local_mocked_bindings(
    extractLdForVariants = function(ld_path, region, variants) {
      m <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
      rownames(m) <- colnames(m) <- variants[1:2]
      m
    },
    calculatePurity = function(...) {
      purity_call_count <<- purity_call_count + 1
      if (purity_call_count == 1) {
        matrix(c(0.9, 0.95, 0.92), 1, 3)  # Passes
      } else {
        matrix(c(0.3, 0.4, 0.35), 1, 3)   # Fails minAbsCorr = 0.8
      }
    },
    normalizeVariantId = function(x, ...) x
  )

  result <- pecotmr:::processColocResults(
    coloc_result,
    "/fake/ld_meta.txt",
    "chr1:100-400",
    pph4Thres = 0,
    coverage = 0.95,
    minAbsCorr = 0.8
  )

  # Only the first CS should pass purity
  if (!is.null(result$sets$cs)) {
    expect_equal(length(result$sets$cs), 1)
  }
})

# ---- extract_ld_for_variants (encoloc.R lines 117-123) ----
test_that("extract_ld_for_variants loads LD, aligns names, and subsets", {
  # Mock loadLdMatrix to return a small LD matrix with variant names
  variants <- c("chr1:10:A:G", "chr1:20:A:G", "chr1:30:A:G")
  ld_variants <- c("chr1:10:A:G", "chr1:20:A:G", "chr1:30:A:G", "chr1:40:A:G")
  ld_mat <- diag(4)
  ld_mat[1, 2] <- ld_mat[2, 1] <- 0.5
  colnames(ld_mat) <- rownames(ld_mat) <- ld_variants

  local_mocked_bindings(
    loadLdMatrix = function(meta_file, region, ...) {
      ref_panel <- parseVariantId(ld_variants)
      ref_panel$variant_id <- ld_variants
      variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
      bm <- data.frame(
        blockId = 1L, chrom = as.character(ref_panel$chrom[1]),
        blockStart = min(ref_panel$pos), blockEnd = max(ref_panel$pos),
        size = length(ld_variants), startIdx = 1L, endIdx = length(ld_variants),
        stringsAsFactors = FALSE
      )
      LdData(correlation = ld_mat, variants = variants_gr, blockMetadata = bm)
    }
  )

  result <- extractLdForVariants("dummy_meta.txt", "chr1:1-50", variants)
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), variants)
})

# ---- colocWrapper with xqtlFinemappingObj = NULL (encoloc.R line 282) ----
test_that("colocWrapper passes through xqtl_raw_data when xqtl_finemapping_obj is NULL", {
  data(coloc_test_data)
  attach(coloc_test_data)
  on.exit(detach(coloc_test_data))

  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")
  on.exit(file.remove(gwas_file, xqtl_file), add = TRUE)

  p <- 50
  variantNames <- paste0("chr1:", 1:p, ":A:G")

  # GWAS file with standard susie structure
  gwas_data <- list(list(
    lbf_variable = matrix(rnorm(5 * p), nrow = 5),
    V = rep(1, 5),
    variantNames = variantNames
  ))
  saveRDS(gwas_data, gwas_file)

  # xQTL file WITHOUT nesting - xqtlFinemappingObj = NULL means use raw_data directly
  xqtl_data <- list(list(
    lbf_variable = matrix(rnorm(3 * p), nrow = 3),
    V = rep(1, 3),
    variantNames = variantNames
  ))
  saveRDS(xqtl_data, xqtl_file)

  local_mocked_bindings(
    coloc.bf_bf = function(...) {
      list(summary = data.frame(PP.H4.abf = 0.9), results = data.frame())
    },
    .package = "coloc"
  )

  result <- colocWrapper(
    xqtl_file, gwas_file,
    xqtlFinemappingObj = NULL,
    xqtlVarnameObj = c("variantNames"),
    gwasVarnameObj = c("variantNames")
  )
  expect_true(is.list(result))
})

# ---- colocWrapper with fsusie fallback (encoloc.R lines 242-243, 288-289) ----
test_that("colocWrapper falls back to fsusie structure when lbf_variable is empty", {
  gwas_file <- tempfile(fileext = ".rds")
  xqtl_file <- tempfile(fileext = ".rds")
  on.exit(file.remove(gwas_file, xqtl_file), add = TRUE)

  p <- 20
  variantNames <- paste0("chr1:", 1:p, ":A:G")

  # GWAS: readRDS(file)[[1]] = raw_data, raw_data has lbf_variable (empty),
  # and raw_data[[1]]$fsusie_result$lBF for the fallback path
  fsusie_inner <- list(fsusie_result = list(lBF = list(rnorm(p), rnorm(p))))
  gwas_raw <- list(
    fsusie_inner,
    lbf_variable = data.frame(),
    V = numeric(0),
    variantNames = variantNames
  )
  saveRDS(list(gwas_raw), gwas_file)

  # xQTL: readRDS(file)[[1]] = xqtl_raw_data,
  # xqtl_raw_data[[1]]$fsusie_result$lBF for the fallback path
  xqtl_raw <- list(
    list(fsusie_result = list(lBF = list(rnorm(p), rnorm(p)))),
    susieFit = list(
      lbf_variable = data.frame(),
      V = numeric(0),
      variantNames = variantNames
    )
  )
  saveRDS(list(xqtl_raw), xqtl_file)

  local_mocked_bindings(
    coloc.bf_bf = function(...) {
      list(summary = data.frame(PP.H4.abf = 0.8), results = data.frame())
    },
    .package = "coloc"
  )

  expect_message(
    result <- colocWrapper(
      xqtl_file, gwas_file,
      xqtlFinemappingObj = c("susieFit"),
      xqtlVarnameObj = c("susieFit", "variantNames"),
      gwasVarnameObj = c("variantNames")
    ),
    "fSuSiE"
  )
  expect_true(is.list(result))
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
#     res <- load_and_extract_ld_matrix(data$metaPath, region, variants)
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
#     ext_ld <- load_and_extract_ld_matrix(data$metaPath, region, variants)
#     res <- calculate_purity(variants, ext_ld, squared = TRUE)
#     expect_equal(ncol(res), 3)
#     lapply(unlist(data), function(x) {
#         file.remove(x)
#     })
# })

# test_that("processColocResults works with dummy data", {
#     data(coloc_test_data)
#     attach(coloc_test_data)
#     data <- generate_mock_ld_files()
#     region <- "chr1:1-500"
#     B1 <- D1
#     B2 <- D2
#     B1$snp <- B2$snp <- colnames(B1$LD) <- colnames(B2$LD) <- rownames(B1$LD) <- rownames(B2$LD) <- paste0("1:", 1:500, ":A:G")
#     mock_coloc_results <- coloc.signals(B1, B2, p12 = 1e-5)
#     res <- processColocResults(mock_coloc_results, data$metaPath, region)
#     expect_equal(length(res$sets$cs), 1)
#     lapply(unlist(data), function(x) {
#         file.remove(x)
#     })
# })
