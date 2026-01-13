# read PLINK files
#' @importFrom dplyr rename
#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_pvar <- function(pgen, num_threads = 1) {
  pvarf <- paste0(file_path_sans_ext(pgen), ".pvar")
  pvardt <- fread(pvarf, skip = "#CHROM", nThread = num_threads)
  pvardt <- rename(pvardt,
    "chrom" = "#CHROM", "pos" = "POS",
    "alt" = "ALT", "ref" = "REF", "id" = "ID"
  )
  pvardt <- select(pvardt, chrom, id, pos, alt, ref)
  return(pvardt)
}

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
read_bim <- function(bed, num_threads = 1) {
  bimf <- paste0(file_path_sans_ext(bed), ".bim")
  bim <- vroom(bimf, col_names = FALSE, num_threads = num_threads)
  colnames(bim) <- c("chrom", "id", "gpos", "pos", "a1", "a0")
  return(bim)
}

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
read_psam <- function(pgen, num_threads = 1) {
  psamf <- paste0(file_path_sans_ext(pgen), ".psam")
  psam <- vroom(psamf, num_threads = num_threads)
  colnames(psam)[1:2] <- c("FID", "IID")
  return(psam)
}

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
read_fam <- function(bed, num_threads = 1) {
  famf <- paste0(file_path_sans_ext(bed), ".fam")
  return(vroom(famf, col_names = FALSE, num_threads = num_threads))
}

# open pgen/pvar PLINK 2 data format
open_pgen <- function(pgenf) {
  # Make sure pgenlibr is installed
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("To use this function, please install pgenlibr: https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  return(pgenlibr::NewPgen(pgenf))
}

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
open_bed <- function(bed, num_threads = 1) {
  # Make sure pgenlibr is installed
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("To use this function, please install pgenlibr: https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  raw_s_ct <- nrow(read_fam(bed, num_threads = num_threads))
  return(pgenlibr::NewPgen(bed, raw_sample_ct = raw_s_ct))
}

read_pgen <- function(pgen, variantidx = NULL, meanimpute = F) {
  # Make sure pgenlibr is installed
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("To use this function, please install pgenlibr: https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  if (is.null(variantidx)) {
    variantidx <- 1:pgenlibr::GetVariantCt(pgen)
  }

  pgenlibr::ReadList(pgen,
    variant_subset = variantidx,
    meanimpute = meanimpute
  )
}

#' @importFrom data.table fread
#' @importFrom dplyr as_tibble mutate filter
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect

tabix_region <- function(file, region, tabix_header = "auto", target = "", target_column_index = "", num_threads = 1) {
  cmd_output <- tryCatch(
    {
      fread(cmd = paste0("tabix -h ", file, " ", region), sep = "auto", header = tabix_header, nThread = num_threads)
    },
    error = function(e) NULL
  )

  if (!is.null(cmd_output) && target != "" && target_column_index != "") {
    cmd_output <- cmd_output %>%
      filter(str_detect(.[[target_column_index]], target))
  } else if (!is.null(cmd_output) && target != "") {
    cmd_output <- cmd_output %>%
      mutate(text = apply(., 1, function(row) paste(row, collapse = "_"))) %>%
      filter(str_detect(text, target)) %>%
      select(-text)
  }

  if (is.null(cmd_output) || nrow(cmd_output) == 0) {
    return(tibble())
  }

  cmd_output %>%
    as_tibble() %>%
    mutate(
      !!names(.)[1] := as.character(.[[1]]),
      !!names(.)[2] := as.numeric(.[[2]])
    )
}


NoSNPsError <- function(message) {
  structure(list(message = message), class = c("NoSNPsError", "error", "condition"))
}

#' Load genotype data for a specific region using vroom for efficiency
#'
#' By default, plink usage dosage of the *major* allele, since "effect allele" A1 is
#' usually the minor allele and the code "1" refers to the "other allele" A2,
#' so that "11" is A2/A2 or major/major. We always use effect allele dosage, to
#' be more consistent with the minor allele based convention ie, plink --recodeA which used minor allele
#' dosage by default.
#'
#' @param genotype Path to the genotype data file (without extension).
#' @param region The target region in the format "chr:start-end".
#' @param keep_indel Whether to keep indel SNPs.
#' @return A vector of SNP IDs in the specified region.
#'
#' @importFrom vroom vroom
#' @importFrom magrittr %>%
#' @export
load_genotype_region <- function(genotype, region = NULL, keep_indel = TRUE, keep_variants_path = NULL, num_threads = 1) {
  # Make sure snpStats is installed
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("To use this function, please install snpStats: https://bioconductor.org/packages/release/bioc/html/snpStats.html")
  }
  if (!is.null(region)) {
    # Get SNP IDs from bim file
    parsed_region <- parse_region(region)
    chrom <- parsed_region$chrom
    start <- parsed_region$start
    end <- parsed_region$end
    # 6 columns for bim file
    col_types <- list(col_character(), col_character(), col_guess(), col_integer(), col_guess(), col_guess())
    # Read a few lines of the bim file to check for 'chr' prefix
    bim_sample <- vroom(paste0(genotype, ".bim"), n_max = 5, col_names = FALSE, col_types = col_types, num_threads = num_threads)
    chr_prefix_present <- any(grepl("^chr", bim_sample$X1))
    # Read the bim file and remove 'chr' prefix if present
    bim_data <- vroom(paste0(genotype, ".bim"), col_names = FALSE, col_types = col_types, num_threads = num_threads)
    if (chr_prefix_present) {
      bim_data$X1 <- gsub("^chr", "", bim_data$X1)
    }
    snp_ids <- filter(bim_data, X1 == chrom & start <= X4 & X4 <= end) %>% pull(X2)
    if (length(snp_ids) == 0) {
      stop(NoSNPsError(paste("No SNPs found in the specified region", region)))
    }
  } else {
    snp_ids <- NULL
  }
  # Read genotype data using snpStats read.plink
  geno <- snpStats::read.plink(genotype, select.snps = snp_ids)

  # Remove indels if specified
  if (!keep_indel) {
    is_indel <- with(geno$map, grepl("[^ATCG]", allele.1) | grepl("[^ATCG]", allele.2) | nchar(allele.1) > 1 | nchar(allele.2) > 1)
    geno_bed <- geno$genotypes[, !is_indel]
    geno_map <- geno$map[!is_indel, ]
  } else {
    geno_bed <- geno$genotypes
    geno_map <- geno$map
  }
  if (!is.null(keep_variants_path)) {
    keep_variants <- vroom(keep_variants_path, num_threads = num_threads)
    if (!("chrom" %in% names(keep_variants)) | !("pos" %in% names(keep_variants))) {
      keep_variants <- do.call(rbind, lapply(strsplit(format_variant_id(keep_variants[[1]]), ":", fixed = TRUE), function(x) {
        data.frame(
          chrom = x[1],
          pos = as.integer(x[2]),
          ref = x[3],
          alt = x[4]
        )
      }))
    }
    if (any(grepl("^chr", keep_variants$chrom))) {
      keep_variants <- keep_variants %>% mutate(chrom = gsub("^chr", "", chrom))
    }
    keep_variants_index <- paste0(geno_map$chromosome, geno_map$position, sep = ":") %in% paste0(keep_variants$chrom, keep_variants$pos, sep = ":")
    geno_bed <- geno_bed[, keep_variants_index]
  } else {
    geno_bed <- geno_bed
  }
  return(2 - as(geno_bed, "numeric"))
}

#' @importFrom purrr map
#' @importFrom readr read_delim cols
#' @importFrom dplyr select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
load_covariate_data <- function(covariate_path, num_threads = 1) {
  return(map(covariate_path, ~ read_delim(.x, "\t", col_types = cols(), num_threads = num_threads) %>%
    select(-1) %>%
    mutate(across(everything(), as.numeric)) %>%
    t()))
}

NoPhenotypeError <- function(message) {
  structure(list(message = message), class = c("NoPhenotypeError", "error", "condition"))
}

#' @importFrom purrr map2
#' @importFrom readr read_delim cols
#' @importFrom dplyr filter select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
load_phenotype_data <- function(phenotype_path, region, extract_region_name = NULL, region_name_col = NULL, tabix_header = TRUE, num_threads = 1) {
  if (is.null(extract_region_name)) {
    extract_region_name <- rep(list(NULL), length(phenotype_path))
  } else if (is.list(extract_region_name) && length(extract_region_name) != length(phenotype_path)) {
    stop("extract_region_name must be NULL or a list with the same length as phenotype_path.")
  } else if (!is.null(extract_region_name) && !is.list(extract_region_name)) {
    stop("extract_region_name must be NULL or a list.")
  }

  # Use `map2` to iterate over `phenotype_path` and `extract_region_name` simultaneously
  phenotype_data <- map2(phenotype_path, extract_region_name, ~ {
    tabix_data <- if (!is.null(region)) tabix_region(.x, region, tabix_header = tabix_header) else read_delim(.x, "\t", col_types = cols(), num_threads = num_threads)
    if (nrow(tabix_data) == 0) {
      message(paste("Phenotype file ", .x, " is empty for the specified region", if (!is.null(region)) "" else region))
      return(NULL)
    }
    # Transpose phenotype data: rows become columns (samples), columns become rows (genes + metadata)
    tabix_data_transposed <- if (!is.null(.y) && is.vector(.y) && !is.null(region_name_col) && (region_name_col %% 1 == 0)) {
      # Filter by extract_region_name if provided
      if (region_name_col <= ncol(tabix_data)) {
        region_col_name <- colnames(tabix_data)[region_name_col]
        tabix_data %>%
          filter(.data[[region_col_name]] %in% .y) %>%
          t()
      } else {
        stop("region_name_col is out of bounds for the number of columns in tabix_data.")
      }
    } else {
      # No filtering, just transpose
      tabix_data %>% t()
    }
    
    # Set colnames from region_name_col row (contains gene IDs) if region_name_col is provided
    # After transposing, the gene ID column becomes a row at position region_name_col
    if (!is.null(region_name_col) && (region_name_col %% 1 == 0) && region_name_col <= nrow(tabix_data_transposed)) {
      colnames(tabix_data_transposed) <- tabix_data_transposed[region_name_col, ]
    }
    
    return(tabix_data_transposed)
  })

  # Check if all phenotype files are empty (all elements NULL)
  if (all(vapply(phenotype_data, is.null, logical(1)))) {
    stop(NoPhenotypeError(paste("All phenotype files are empty for the specified region", if (!is.null(region)) "" else region)))
  }
  return(phenotype_data)
}

#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @noRd
extract_phenotype_coordinates <- function(phenotype_list) {
  return(map(phenotype_list, ~ t(.x[1:3, ]) %>%
    as_tibble() %>%
    mutate(start = as.numeric(start), end = as.numeric(end))))
}

#' @importFrom magrittr %>%
#' @noRd
filter_by_common_samples <- function(dat, common_samples) {
  dat[common_samples, , drop = FALSE] %>% .[order(rownames(.)), , drop = FALSE]
}

#' @importFrom tibble tibble
#' @importFrom dplyr mutate select
#' @importFrom purrr map map2
#' @importFrom magrittr %>%
#' @noRd
prepare_data_list <- function(geno_bed, phenotype, covariate, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, phenotype_header = 4, keep_samples = NULL) {
  data_list <- tibble(
    covar = covariate,
    Y = lapply(phenotype, function(x) apply(x[-c(1:phenotype_header), , drop = F], c(1, 2), as.numeric))
  ) %>%
    mutate(
      # Determine common complete samples across Y, covar, and geno_bed, considering missing values
      common_complete_samples = map2(covar, Y, ~ {
        covar_non_na <- rownames(.x)[!apply(.x, 1, function(row) all(is.na(row)))]
        y_non_na <- rownames(.y)[!apply(.y, 1, function(row) all(is.na(row)))]
        if (length(intersect(intersect(covar_non_na, y_non_na), rownames(geno_bed))) == 0) {
          stop("No common complete samples between genotype and phenotype/covariate data")
        }
        intersect(intersect(covar_non_na, y_non_na), rownames(geno_bed))
      }),
      # Further intersect with keep_samples if provided
      common_complete_samples = if (!is.null(keep_samples) && length(keep_samples) > 0) {
        map(common_complete_samples, ~ intersect(.x, keep_samples))
      } else {
        common_complete_samples
      },
      # Determine dropped samples before filtering
      dropped_samples_covar = map2(covar, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_Y = map2(Y, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_X = map(common_complete_samples, ~ setdiff(rownames(geno_bed), .x)),
      # Filter data based on common complete samples
      Y = map2(Y, common_complete_samples, ~ filter_by_common_samples(.x, .y)),
      covar = map2(covar, common_complete_samples, ~ filter_by_common_samples(.x, .y)),
      # Apply filter_X on the geno_bed data filtered by common complete samples and then format column names
      X = map(common_complete_samples, ~ {
        filtered_geno_bed <- filter_by_common_samples(geno_bed, .x)
        mac_val <- if (nrow(filtered_geno_bed) == 0) 0 else (mac_cutoff / (2 * nrow(filtered_geno_bed)))
        maf_val <- max(maf_cutoff, mac_val)
        filtered_data <- filter_X(filtered_geno_bed, imiss_cutoff, maf_val, var_thresh = xvar_cutoff)
        colnames(filtered_data) <- format_variant_id(colnames(filtered_data)) # Format column names right after filtering
        colnames(filtered_data) <- normalize_variant_id(colnames(filtered_data)) # Normalize to chr{chrom}:{pos}:{A2}:{A1} format
        filtered_data
      })
    ) %>%
    select(covar, Y, X, dropped_samples_Y, dropped_samples_X, dropped_samples_covar)
  return(data_list)
}

#' @importFrom purrr map
#' @importFrom dplyr intersect
#' @importFrom magrittr %>%
#' @noRd
prepare_X_matrix <- function(geno_bed, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff) {
  # Calculate the union of all samples from data_list: any of X, covar and Y would do
  all_samples_union <- map(data_list$covar, ~ rownames(.x)) %>%
    unlist() %>%
    unique()
  # Find the intersection of these samples with the samples in geno_bed
  common_samples <- intersect(all_samples_union, rownames(geno_bed))
  # Filter geno_bed using common_samples
  X_filtered <- filter_by_common_samples(geno_bed, common_samples)
  # Calculate MAF cutoff considering the number of common samples
  maf_val <- max(maf_cutoff, mac_cutoff / (2 * length(common_samples)))
  # Apply further filtering on X
  X_filtered <- filter_X(X_filtered, imiss_cutoff, maf_val, xvar_cutoff)
  colnames(X_filtered) <- format_variant_id(colnames(X_filtered))
  # Normalize variant IDs to chr{chrom}:{pos}:{A2}:{A1} format
  colnames(X_filtered) <- normalize_variant_id(colnames(X_filtered))

  # To keep a log message
  variants <- as.data.frame(do.call(rbind, lapply(format_variant_id(colnames(X_filtered)), function(x) strsplit(x, ":")[[1]][1:2])), stringsAsFactors = FALSE)
  message(paste0("Dimension of input genotype data is ", nrow(X_filtered), " rows and ", ncol(X_filtered), " columns for genomic region of ", variants[1, 1], ":", min(as.integer(variants[, 2])), "-", max(as.integer(variants[, 2]))))
  return(X_filtered)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
add_X_residuals <- function(data_list, scale_residuals = FALSE) {
  # Compute residuals for X and add them to data_list
  n_conditions <- length(data_list$X)
  lm_res_X_list <- list()
  for (idx in seq_along(data_list$X)) {
    lm_res_X_list[[idx]] <- .lm.fit(x = cbind(1, data_list$covar[[idx]]), y = data_list$X[[idx]])$residuals %>% as.matrix()
  }
  data_list <- data_list %>%
    mutate(
      lm_res_X = lm_res_X_list,
      X_resid_mean = map(lm_res_X, ~ apply(.x, 2, mean)),
      X_resid_sd = map(lm_res_X, ~ apply(.x, 2, sd)),
      X_resid = map(lm_res_X, ~ {
        if (scale_residuals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  return(data_list)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
add_Y_residuals <- function(data_list, conditions, scale_residuals = FALSE) {
  # Compute residuals, their mean, and standard deviation, and add them to data_list
  # Preserve colnames from original Y (gene IDs) through the residual computation
  n_conditions <- length(data_list$Y)
  lm_res_list <- list()
  for (idx in seq_along(data_list$Y)) {
    y <- data_list$Y[[idx]]
    res <- .lm.fit(x = cbind(1, data_list$covar[[idx]]), y = y)$residuals %>% as.matrix()
    # Preserve colnames from original Y (gene IDs)
    if (!is.null(colnames(y))) {
      colnames(res) <- colnames(y)
    }
    lm_res_list[[idx]] <- res
  }
  data_list <- data_list %>%
    mutate(
      lm_res = lm_res_list,
      Y_resid_mean = map(lm_res, ~ apply(.x, 2, mean)),
      Y_resid_sd = map(lm_res, ~ apply(.x, 2, sd)),
      Y_resid = map2(lm_res, Y, function(lm, y_orig) {
        if (scale_residuals) {
          res <- scale(lm)
        } else {
          res <- lm
        }
        # Preserve colnames from original Y (gene IDs)
        if (!is.null(colnames(y_orig))) {
          colnames(res) <- colnames(y_orig)
        }
        return(res)
      })
    )

  names(data_list$Y_resid) <- conditions

  return(data_list)
}

#' Load regional association data
#'
#' This function loads genotype, phenotype, and covariate data for a specific region and performs data preprocessing.
#'
#' @param genotype PLINK bed file containing genotype data.
#' @param phenotype A vector of phenotype file names.
#' @param covariate A vector of covariate file names corresponding to the phenotype file vector.
#' @param region A string of chr:start-end for the phenotype region.
#' @param conditions A vector of strings representing different conditions or groups.
#' @param maf_cutoff Minimum minor allele frequency (MAF) cutoff. Default is 0.
#' @param mac_cutoff Minimum minor allele count (MAC) cutoff. Default is 0.
#' @param xvar_cutoff Minimum variance cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
#' @param region_name_col Column name containing the region name. Default is NULL.
#' @param keep_indel Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.
#' @param keep_samples A vector of sample names to keep. Default is NULL.
#' @param phenotype_header Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).
#' @param scale_residuals Logical indicating whether to scale residuals. Default is FALSE.
#' @param tabix_header Logical indicating whether the tabix file has a header. Default is TRUE.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
#'   \item residual_X: A list of residualized genotype matrices for each condition.
#'   \item residual_Y_scalar: Scaling factor for residualized phenotype values.
#'   \item residual_X_scalar: Scaling factor for residualized genotype values.
#'   \item dropped_sample: A list of dropped samples for X, Y, and covariates.
#'   \item covar: Covariate data.
#'   \item Y: Original phenotype data.
#'   \item X_data: Original genotype data.
#'   \item X: Filtered genotype matrix.
#'   \item maf: Minor allele frequency (MAF) for each variant.
#'   \item chrom: Chromosome of the region.
#'   \item grange: Genomic range of the region (start and end positions).
#'   \item Y_coordinates: Phenotype coordinates if a region is specified.
#' }
#'
#' @export
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end for phenotype region
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           xvar_cutoff = 0,
                                           imiss_cutoff = 0,
                                           association_window = NULL,
                                           extract_region_name = NULL,
                                           region_name_col = NULL,
                                           keep_indel = TRUE,
                                           keep_samples = NULL,
                                           keep_variants = NULL,
                                           phenotype_header = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                           scale_residuals = FALSE,
                                           tabix_header = TRUE, 
                                           num_threads = 1) {
  ## Load genotype
  geno <- load_genotype_region(genotype, association_window, keep_indel, keep_variants_path = keep_variants, num_threads = num_threads)
  ## Load phenotype and covariates and perform some pre-processing
  covar <- load_covariate_data(covariate, num_threads = num_threads)
  pheno <- load_phenotype_data(phenotype, region, extract_region_name = extract_region_name, region_name_col = region_name_col, tabix_header = tabix_header, num_threads = num_threads)
  # Keep covariates / conditions aligned with successfully loaded phenotypes.
  # load_phenotype_data() returns NULL for phenotype files that are empty in
  # the requested region; we must drop the corresponding covariate / condition
  # entries before constructing data_list.
  keep_idx <- !vapply(pheno, is.null, logical(1))
  if (!any(keep_idx)) {
    stop(NoPhenotypeError(paste("All phenotype files are empty for the specified region", if (!is.null(region)) "" else region)))
  }
  if (any(!keep_idx)) {
    dropped_conditions <- conditions[!keep_idx]
    message(paste(
      "Dropping phenotypes with no data in region", region, ":",
      paste(dropped_conditions, collapse = ", ")
    ))
  }
  pheno <- pheno[keep_idx]
  covar <- covar[keep_idx]
  conditions <- conditions[keep_idx]
  ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
  data_list <- prepare_data_list(geno, pheno, covar, imiss_cutoff,
    maf_cutoff, mac_cutoff, xvar_cutoff,
    phenotype_header = phenotype_header, keep_samples = keep_samples
  )
  maf_list <- setNames(lapply(data_list$X, function(x) apply(x, 2, compute_maf)), colnames(data_list$X))
  
  ## Get residue Y for each of condition and its mean and sd
  data_list <- add_Y_residuals(data_list, conditions, scale_residuals)
  
  ## Get residue X for each of condition and its mean and sd
  data_list <- add_X_residuals(data_list, scale_residuals)
  
  # Get X matrix for union of samples
  X <- prepare_X_matrix(geno, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
  region <- if (!is.null(region)) unlist(strsplit(region, ":", fixed = TRUE))
  ## residual_Y: a list of y either vector or matrix (CpG for example), and they need to match with residual_X in terms of which samples are missing.
  ## residual_X: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
  ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y matrix form (not list); but the matrices inside residual_X would be subsets of sample name of residual_Y matrix form (not list).
  return(list(
    residual_Y = data_list$Y_resid,
    residual_X = data_list$X_resid,
    residual_Y_scalar = if (scale_residuals) data_list$Y_resid_sd else rep(1, length(data_list$Y_resid)),
    residual_X_scalar = if (scale_residuals) data_list$X_resid_sd else rep(1, length(data_list$X_resid)),
    dropped_sample = list(X = data_list$dropped_samples_X, Y = data_list$dropped_samples_Y, covar = data_list$dropped_samples_covar),
    covar = data_list$covar,
    Y = data_list$Y,
    X_data = data_list$X,
    X = X,
    maf = maf_list,
    chrom = region[1],
    grange = if (!is.null(region)) unlist(strsplit(region[2], "-", fixed = TRUE)) else NULL,
    Y_coordinates = if (!is.null(region)) extract_phenotype_coordinates(pheno) else NULL
  ))
}

#' Load Regional Univariate Association Data
#'
#' This function loads regional association data for univariate analysis.
#' It includes residual matrices, original genotype data, and additional metadata.
#'
#' @importFrom matrixStats colVars
#' @return A list
#' @export
load_regional_univariate_data <- function(..., num_threads = 1) {
  dat <- load_regional_association_data(..., num_threads = num_threads)
  return(list(
    residual_Y = dat$residual_Y,
    residual_X = dat$residual_X,
    residual_Y_scalar = dat$residual_Y_scalar,
    residual_X_scalar = dat$residual_X_scalar,
    dropped_sample = dat$dropped_sample,
    maf = dat$maf,
    X = dat$X, # X unadjusted by covariate
    chrom = dat$chrom,
    grange = dat$grange,
    X_variance = lapply(dat$residual_X, function(x) colVars(x))
  ))
}

#' Load Regional Data for Regression Modeling
#'
#' This function loads regional association data formatted for regression modeling.
#' It includes phenotype, genotype, and covariate matrices along with metadata.
#'
#' @return A list
#' @export
load_regional_regression_data <- function(..., num_threads = 1) {
  dat <- load_regional_association_data(..., num_threads = num_threads)
  return(list(
    Y = dat$Y,
    X_data = dat$X_data,
    covar = dat$covar,
    dropped_sample = dat$dropped_sample,
    maf = dat$maf,
    chrom = dat$chrom,
    grange = dat$grange
  ))
}

# return matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column.
#' @noRd
pheno_list_to_mat <- function(data_list) {
  all_row_names <- unique(unlist(lapply(data_list$residual_Y, rownames)))
  # Step 2: Align matrices and fill with NA where necessary
  aligned_mats <- lapply(data_list$residual_Y, function(mat) {
    ### change the ncol of each matrix
    expanded_mat <- matrix(NA, nrow = length(all_row_names), ncol = ncol(mat), dimnames = list(all_row_names, colnames(mat)))
    common_rows <- intersect(rownames(mat), all_row_names)
    expanded_mat[common_rows, ] <- mat[common_rows, ]
    return(expanded_mat)
  })
  Y_resid_matrix <- do.call(cbind, aligned_mats)
  if (!is.null(names(data_list$residual_Y))) {
    colnames(Y_resid_matrix) <- names(data_list$residual_Y)
  }
  data_list$residual_Y <- Y_resid_matrix
  return(data_list)
}

#' Load and Preprocess Regional Multivariate Data
#'
#' This function loads regional association data and processes it into a multivariate format.
#' It optionally filters out samples based on missingness thresholds in the response matrix.
#'
#' @importFrom matrixStats colVars
#' @return A list
#' @export
load_regional_multivariate_data <- function(matrix_y_min_complete = NULL, # when Y is saved as matrix, remove those with non-missing counts less than this cutoff
                                            ..., 
                                            num_threads = 1) {
  dat <- pheno_list_to_mat(load_regional_association_data(..., num_threads = num_threads))
  if (!is.null(matrix_y_min_complete)) {
    Y <- filter_Y(dat$residual_Y, matrix_y_min_complete)
    if (length(Y$rm_rows) > 0) {
      X <- dat$X[-Y$rm_rows, ]
      Y_scalar <- dat$residual_Y_scalar[-Y$rm_rows]
      dropped_sample <- rownames(dat$residual_Y)[Y$rm_rows]
    } else {
      X <- dat$X
      Y_scalar <- dat$residual_Y_scalar
      dropped_sample <- dat$dropped_sample
    }
  } else {
    Y <- dat$residual_Y
    X <- dat$X
    Y_scalar <- dat$residual_Y_scalar
    dropped_sample <- dat$dropped_sample
  }
  return(list(
    residual_Y = Y,
    residual_Y_scalar = Y_scalar,
    dropped_sample = dropped_sample,
    X = X,
    maf = apply(X, 2, compute_maf),
    chrom = dat$chrom,
    grange = dat$grange,
    X_variance = colVars(X)
  ))
}

#' Load Regional Functional Association Data
#'
#' This function loads precomputed regional functional association data.
#'
#' @return A list
#' @export
load_regional_functional_data <- function(..., num_threads = 1) {
  dat <- load_regional_association_data(..., num_threads = num_threads)
  return(dat)
}



# Function to remove gene name at the end of context name
#' @export
clean_context_names <- function(context, gene) {
  # Remove gene name if it matches the last part of the context
  gene <- gene[order(-nchar(unique(gene)))]
  for (gene_id in gene) {
    context <- gsub(paste0("_", gene_id), "", context)
  }
  return(context)
}

#' Load, Validate, and Consolidate TWAS Weights from Multiple RDS Files
#'
#' This function loads TWAS weight data from multiple RDS files, checks for the presence
#' of specified region and condition. If variable_name_obj is provided, it aligns and
#' consolidates weight matrices based on the object's variant names, filling missing data
#' with zeros. If variable_name_obj is NULL, it checks that all files have the same row
#' numbers for the condition and consolidates weights accordingly.
#'
#' @param weight_db_file weight_db_files Vector of file paths for RDS files containing TWAS weights..
#' Each element organized as region/condition/weights
#' @param condition The specific condition to be checked and consolidated across all files.
#' @param variable_name_obj The name of the variable/object to fetch from each file, if not NULL.
#' @return A consolidated list of weights for the specified condition and a list of SuSiE results.
#' @examples
#' # Example usage (replace with actual file paths, condition, region, and variable_name_obj):
#' weight_db_files <- c("path/to/file1.rds", "path/to/file2.rds")
#' condition <- "example_condition"
#' region <- "example_region"
#' variable_name_obj <- "example_variable" # or NULL for standard processing
#' consolidated_weights <- load_twas_weights(weight_db_files, condition, region, variable_name_obj)
#' print(consolidated_weights)
#' @export
load_twas_weights <- function(weight_db_files, conditions = NULL,
                              variable_name_obj = c("preset_variants_result", "variant_names"),
                              susie_obj = c("preset_variants_result", "susie_result_trimmed"),
                              twas_weights_table = "twas_weights") {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, conditions, variable_name_obj) {
    all_data <- do.call(c, lapply(unname(weight_db_files), function(rds_file) {
      db <- readRDS(rds_file)
      gene <- names(db)
      if (any(unique(names(find_data(db, c(3, "twas_weights")))) %in% c("mrmash_weights", "mvsusie_weights"))) {
        names(db[[1]]) <- clean_context_names(names(db[[1]]), gene = gene)
        db <- list(mnm_rs = db[[1]])
      } else {
        # Check if region from all RDS files are the same
        if (length(gene) != 1) {
          stop("More than one region provided in the RDS file. ")
        } else {
          names(db[[gene]]) <- clean_context_names(names(db[[gene]]), gene = gene)
        }
      }
      return(db)
    }))

    # Combine the lists with the same region name
    gene <- unique(names(all_data)[!names(all_data) %in% "mnm_rs"])
    if (length(gene) > 1) stop("More than one region of twas weights data provided. ")
    combined_all_data <- lapply(split(all_data, names(all_data)), function(lst) {
      if (length(lst) > 1) {
        lst <- do.call(c, unname(lst))
      }
      if (isTRUE(names(lst) == "mnm_rs")) lst <- lst[[1]]
      if (gene %in% names(lst)) lst <- do.call(c, lapply(unname(lst), function(x) x))
      return(lst)
    })

    # merge univariate and multivariate results for same gene-context pair
    if ("mnm_rs" %in% names(combined_all_data)) {
      # gene <- names(combined_all_data)[!names(combined_all_data) %in% "mnm_rs"]
      overl_contexts <- names(combined_all_data[["mnm_rs"]])[names(combined_all_data[["mnm_rs"]]) %in% names(combined_all_data[[gene]])]
      multi_variants <- unique(find_data(combined_all_data$mnm_rs, c(2, variable_name_obj)))
      for (context in overl_contexts) {
        uni_variants <- get_nested_element(combined_all_data[[gene]][[context]], variable_name_obj)
        multi_weights <- setNames(rep(0, length(uni_variants)), uni_variants)
        multi_weights <- lapply(combined_all_data[["mnm_rs"]][[context]]$twas_weights, function(weight_list) {
          aligned_weights <- setNames(rep(0, length(uni_variants)), uni_variants)
          method_weight_variants <- names(unlist(weight_list))
          overlap_variants <- method_weight_variants[method_weight_variants %in% multi_variants[multi_variants %in% uni_variants]] # overlapping variants from method, multivariate, univariate
          aligned_weights[overlap_variants] <- unlist(weight_list)[overlap_variants]
          aligned_weights <- as.matrix(aligned_weights)
        })
        combined_all_data[[gene]][[context]]$twas_weights <- c(combined_all_data[[gene]][[context]]$twas_weights, multi_weights)
        combined_all_data[[gene]][[context]]$twas_cv_result$performance <- c(
          combined_all_data[[gene]][[context]]$twas_cv_result$performance,
          combined_all_data[["mnm_rs"]][[context]]$twas_cv_result$performance
        )
      }
      combined_all_data[["mnm_rs"]] <- NULL
    }
    if (gene %in% names(combined_all_data)) combined_all_data <- do.call(c, unname(combined_all_data))
    if (gene %in% names(combined_all_data)) combined_all_data <- combined_all_data[[1]]

    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data)
    }
    ## Check if the specified condition and variable_name_obj are available in all files
    if (!all(conditions %in% names(combined_all_data))) {
      stop("The specified condition is not available in all RDS files.")
    }
    return(combined_all_data)
  }

  # Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    if (length(weights_list) != length(variable_objs)) {
      stop("The length of the weights_list and variable_objs must be the same.")
    }
    # Loop through each weight matrix and assign variant names as rownames
    for (i in seq_along(weights_list)) {
      # Ensure dimensions match
      if (nrow(weights_list[[i]]) != length(variable_objs[[i]])) {
        stop(paste("Number of rows in weights_list[[", i, "]] does not match the length of variable_objs[[", i, "]]", sep = ""))
      }
      # Apply variant names to the row names of the weight matrix
      rownames(weights_list[[i]]) <- variable_objs[[i]]
    }
    return(weights_list)
  }

  # Internal function to consolidate weights for given condition
  consolidate_weights_list <- function(combined_all_data, conditions, variable_name_obj, twas_weights_table) {
    combined_weights_by_condition <- lapply(conditions, function(condition) {
      temp_list <- get_nested_element(combined_all_data, c(condition, twas_weights_table))
      sapply(temp_list, cbind)
    })
    names(combined_weights_by_condition) <- conditions
    if (is.null(variable_name_obj)) {
      # Standard processing: Check for identical row numbers and consolidate
      row_numbers <- sapply(combined_weights_by_condition, function(data) nrow(data))
      if (length(unique(row_numbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- combined_weights_by_condition
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variable_objs <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, variable_name_obj))
      })
      weights <- align_and_merge(combined_weights_by_condition, variable_objs)
    }
    names(weights) <- conditions
    return(weights)
  }

  ## Load, validate, and consolidate data
  try(
    {
      combined_all_data <- load_and_validate_data(weight_db_files, conditions, variable_name_obj)
      if (is.null(combined_all_data)) {
        return(NULL)
      }
      # update condition in case of merging rds files
      conditions <- names(combined_all_data)
      weights <- consolidate_weights_list(combined_all_data, conditions, variable_name_obj, twas_weights_table)
      combined_susie_result <- lapply(combined_all_data, function(context) get_nested_element(context, susie_obj))
      performance_tables <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, "twas_cv_result", "performance"))
      })
      names(performance_tables) <- conditions
      return(list(susie_results = combined_susie_result, weights = weights, twas_cv_performance = performance_tables))
    },
    silent = TRUE
  )
}

#' Load summary statistic data
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. The mapping is performed through the specified column file.
#' Additionally, it extracts sample size, case number, control number, and variance of Y.
#' Missing values in n_sample, n_case, and n_control are backfilled with median values.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column file for mapping.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#' @param region The region where tabix use to subset the input dataset.
#' @param extract_region_name User-specified gene/phenotype name used to further subset the phenotype data.
#' @param region_name_col Filter this specific column for the extract_region_name.
#' @param comment_string Comment sign in the column_mapping file, default is #
#' @return A list of rss_input, including the column-name-formatted summary statistics,
#' sample size (n), and var_y.
#'
#' @importFrom dplyr mutate group_by summarise
#' @importFrom magrittr %>%
#' @export
load_rss_data <- function(sumstat_path, column_file_path, n_sample = 0, n_case = 0, n_control = 0, region = NULL,
                          extract_region_name = NULL, region_name_col = NULL, comment_string = "#", num_threads = 1) {
  message("Loading RSS data from ", sumstat_path, " with n_sample = ", n_sample, " n_case = ", n_case, " n_control = ", n_control)
  # Read and preprocess column mapping
  if (is.null(comment_string)) {
    column_data <- read.table(column_file_path,
      header = FALSE, sep = ":",
      comment.char = "", # This tells R not to treat any character as comment
      stringsAsFactors = FALSE
    ) %>%
      rename(standard = V1, original = V2)
  } else {
    column_data <- read.table(column_file_path,
      header = FALSE, sep = ":",
      comment.char = comment_string,
      stringsAsFactors = FALSE
    ) %>%
      rename(standard = V1, original = V2)
  }
  # Initialize sumstats variable
  sumstats <- NULL
  var_y <- NULL
  sumstats <- load_tsv_region(file_path = sumstat_path, region = region, extract_region_name = extract_region_name, region_name_col = region_name_col, num_threads = num_threads)
  
  # To keep a log message
  n_variants <- nrow(sumstats)
  if (n_variants == 0){
      message(paste0("No variants in region ", region, "."))
      return(list(sumstats = sumstats, n = NULL, var_y = NULL))
  } else {
      message(paste0("Region ", region, " include ", n_variants, " in input sumstats."))
  }
  
  # Standardize column names based on mapping
  for (name in colnames(sumstats)) {
    if (name %in% column_data$original) {
      index <- which(column_data$original == name)
      colnames(sumstats)[colnames(sumstats) == name] <- column_data$standard[index]
    }
  }
  if (!"z" %in% colnames(sumstats) && all(c("beta", "se") %in%
    colnames(sumstats))) {
    sumstats$z <- sumstats$beta / sumstats$se
  }
  if (!"beta" %in% colnames(sumstats) && "z" %in% colnames(sumstats)) {
    sumstats$beta <- sumstats$z
    sumstats$se <- 1
  }
  for (col in c("n_sample", "n_case", "n_control")) {
    if (col %in% colnames(sumstats)) {
      sumstats[[col]][is.na(sumstats[[col]])] <- median(sumstats[[col]],
        na.rm = TRUE
      )
    }
  }
  if (n_sample != 0 && (n_case + n_control) != 0) {
    stop("Please provide sample size, or case number with control number, but not both")
  } else if (n_sample != 0) {
    n <- n_sample
  } else if ((n_case + n_control) != 0) {
    n <- n_case + n_control
    phi <- n_case / n
    var_y <- 1 / (phi * (1 - phi))
  } else {
    if ("n_sample" %in% colnames(sumstats)) {
      n <- median(sumstats$n_sample)
    } else if (all(c("n_case", "n_control") %in% colnames(sumstats))) {
      n <- median(sumstats$n_case + sumstats$n_control)
      phi <- median(sumstats$n_case / n)
      var_y <- 1 / (phi * (1 - phi))
    } else {
      warning("Sample size and variance of Y could not be determined from the summary statistics.")
      n <- NULL
    }
  }
  return(list(sumstats = sumstats, n = n, var_y = var_y))
}

#' Load and filter tabular data with optional region subsetting
#'
#' This function loads summary statistics data from tabular files (TSV, TXT).
#' For compressed (.gz) and tabix-indexed files, it can subset data by genomic region.
#' Additionally, it can filter results by a specified target value in a designated column.
#'
#' @param file_path Path to the summary statistics file.
#' @param region Genomic region for subsetting tabix-indexed files. Format: chr:start-end (e.g., "9:10000-50000").
#' @param extract_region_name Value to filter for in the specified filter column.
#' @param region_name_col Index of the column to apply the extract_region_name against.
#'
#' @return A dataframe containing the filtered summary statistics.
#'
#' @importFrom data.table fread
#' @importFrom vroom vroom
#' @export
load_tsv_region <- function(file_path, region = NULL, extract_region_name = NULL, region_name_col = NULL, num_threads = 1) {
  sumstats <- NULL
  cmd <- NULL

  if (!is.null(region)) {
    if (grepl("^chr", region)) {
      region <- sub("^chr", "", region)
    }
  }

  if (grepl("\\.gz$", file_path)) {
    if (is.null(sumstats) || nrow(sumstats) == 0) {
      # Determine the appropriate command based on provided parameters
      if (!is.null(extract_region_name) && !is.null(region) && !is.null(region_name_col)) {
        # Both region and filter specified
        cmd <- paste0(
          "zcat ", file_path, " | head -1 && tabix ", file_path, " ", region,
          " | awk '$", region_name_col, " ~ /", extract_region_name, "/'"
        )
      } else if (!is.null(extract_region_name) && is.null(region) && !is.null(region_name_col)) {
        # Only filter specified, no region
        cmd <- paste0("zcat ", file_path, " | awk '$", region_name_col, " ~ /", extract_region_name, "/'")
      } else if (!is.null(region) && (is.null(region_name_col) || is.null(extract_region_name))) {
        # Only region specified, no filter
        cmd <- paste0("zcat ", file_path, " | head -1 && tabix ", file_path, " ", region)
      } else {
        # Neither region nor filter specified
        cmd <- file_path
      }

      sumstats <- tryCatch(
        {
          fread(cmd)
        },
        error = function(e) {
          stop("Data read error. Please make sure this gz file is tabix-indexed and the specified filter column exists.")
        }
      )
    }
  } else {
    warning("Not a tabix-indexed gz file, loading the entire dataset.")
    sumstats <- vroom(file_path, num_threads = num_threads)
    # Apply filter if specified
    if (!is.null(extract_region_name) && !is.null(region_name_col)) {
      keep_index <- which(str_detect(sumstats[[region_name_col]], extract_region_name))
      sumstats <- sumstats[keep_index, ]
    }
  }

  return(sumstats)
}

#' Split loaded twas_weights_results into batches based on maximum memory usage
#'
#' @param twas_weights_results List of loaded gene data by load_twas_weights()
#' @param meta_data_df Dataframe containing gene metadata with region_id and TSS columns
#' @param max_memory_per_batch Maximum memory per batch in MB (default: 750)
#' @return List of batches, where each batch contains a subset of twas_weights_results
#' @export
batch_load_twas_weights <- function(twas_weights_results, meta_data_df, max_memory_per_batch = 750) {
  gene_names <- names(twas_weights_results)
  if (length(gene_names) == 0) {
    message("No genes in twas_weights_results.")
    return(list())
  }

  gene_memory_df <- data.frame(
    gene_name = gene_names, memory_mb = sapply(gene_names, function(gene) {
      as.numeric(object.size(twas_weights_results[[gene]])) / (1024^2) # Get object size in bytes and convert to MB
    })
  )

  # Merge with meta_data_df to get TSS information
  gene_memory_df <- merge(gene_memory_df, meta_data_df[, c("region_id", "TSS")],
    by.x = "gene_name",
    by.y = "region_id", all.x = TRUE
  )
  gene_memory_df <- gene_memory_df[order(gene_memory_df$TSS), ]

  # Check if we need to split into batches
  total_memory_mb <- sum(gene_memory_df$memory_mb)
  message("Total memory usage: ", round(total_memory_mb, 2), " MB")
  if (total_memory_mb <= max_memory_per_batch) {
    message("All genes fit within the memory limit. No need to split into batches.")
    return(list(all_genes = twas_weights_results))
  }

  # Create batches by adding genes until we reach the memory limit
  batches <- list()
  current_batch_genes <- character(0)
  current_batch_memory <- 0
  batch_index <- 1

  for (i in 1:nrow(gene_memory_df)) {
    gene <- gene_memory_df$gene_name[i]
    gene_memory <- gene_memory_df$memory_mb[i]
    # If a single gene exceeds the memory limit, include it in its own batch
    if (gene_memory > max_memory_per_batch) {
      batches[[paste0("batch_", batch_index)]] <- twas_weights_results[gene]
      batch_index <- batch_index + 1
      next
    }
    # If adding this gene would exceed the memory limit, start a new batch
    if (current_batch_memory + gene_memory > max_memory_per_batch && length(current_batch_genes) > 0) {
      batches[[paste0("batch_", batch_index)]] <- twas_weights_results[current_batch_genes]
      current_batch_genes <- character(0)
      current_batch_memory <- 0
      batch_index <- batch_index + 1
    }
    current_batch_genes <- c(current_batch_genes, gene)
    current_batch_memory <- current_batch_memory + gene_memory
  }
  # Add the last batch if not empty
  if (length(current_batch_genes) > 0) {
    batches[[batch_index]] <- twas_weights_results[current_batch_genes]
  }
  message("Split into ", length(batches), " batches")
  names(batches) <- NULL
  return(batches)
}

#' Compute genotype correlation 
#' @export
get_cormat <- function(X, intercepte = TRUE) {
  X <- t(X)
  # Center each variable
  if (intercepte) {
    X <- X - rowMeans(X)
  }
  # Standardize each variable
  X <- X / sqrt(rowSums(X^2))
  # Calculate correlations
  cr <- tcrossprod(X)
  return(cr)
}

# Function to filter a single credible set based on coverage and purity
#' @importFrom susieR susie_get_cs
#' @importFrom purrr map_lgl
#' @export      
get_filter_lbf_index <- function(susie_obj, coverage = 0.5, size_factor = 0.5) {
  susie_obj$V <- NULL  # ensure no filtering by estimated prior

  # Get CS list with coverage
  cs_list <- susie_get_cs(susie_obj, coverage = coverage, dedup = FALSE)
  
  # Total number of variants
  total_variants <- ncol(susie_obj$alpha)
  
  # Maximum allowed CS size to be considered 'concentrated'
  max_size <- total_variants * coverage * size_factor
  
  # Identify which CSs are 'concentrated enough'
  keep_idx <- map_lgl(cs_list$cs, ~ length(.x) < max_size)
  
  # Extract the CS indices that pass the filter
  cs_index <- which(keep_idx) %>% names %>% gsub("L","", .) %>% as.numeric

  # Return filtered lbf_variable rows (one per CS)
  return(cs_index)
}

#' Function to load LD reference data variants
#' @export
#' @noRd
load_ld_snp_info <- function(ld_meta_file_path, region_of_interest, num_threads = 1) {
  bim_file_path <- get_regional_ld_meta(ld_meta_file_path, region_of_interest)$intersections$bim_file_paths
  bim_data <- lapply(bim_file_path, function(bim_file) as.data.frame(vroom(bim_file, col_names = FALSE, num_threads = num_threads)))
  snp_info <- setNames(lapply(bim_data, function(info_table) {
    # for TWAS and MR, the variance and allele_freq are not necessary
    if (ncol(info_table) >= 8) {
      info_table <- info_table[, c(1, 2, 4:8)]
      colnames(info_table) <- c("chrom", "id", "pos", "alt", "ref", "variance", "allele_freq")
    } else if (ncol(info_table) == 6) {
      info_table <- info_table[, c(1, 2, 4:6)]
      colnames(info_table) <- c("chrom", "id", "pos", "alt", "ref")
    } else {
      warning("Unexpected number of columns; skipping this element.")
      return(NULL)
    }
    info_table$id <- gsub("chr", "", gsub("_", ":", info_table$id))
    return(info_table)
  }), sapply(names(bim_data), function(x) gsub("chr", "", paste(strsplit(basename(x), "[_:/.]")[[1]][1:3], collapse = "_"))))
  return(snp_info)
}