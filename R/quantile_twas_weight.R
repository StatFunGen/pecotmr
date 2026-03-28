#' @title Quantile TWAS Weight Calculation and QTL Analysis
#'
#' @description
#' This file contains functions for performing Quantile Transcriptome-Wide
#' Association Studies (TWAS) weight calculations and Quantile QTL analysis.
#' It provides tools for screening quantile regression results, performing
#' LD clumping and pruning, and calculating TWAS weights.
#'
#' @details
#' The main function in this file is `quantile_twas_weight_pipeline`, which
#' orchestrates the entire analysis process. Other functions are helper
#' functions used within the main pipeline.
#'

#' Qrank Score Test Screen
#' @param X Matrix of predictors
#' @param Y Matrix or vector of response variables
#' @param Z Matrix of covariates (optional)
#' @param tau.list Vector of quantiles to be analyzed
#' @param screen_threshold Significance threshold for adjusted p-values
#' @param screen_method Method for p-value adjustment ('fdr' or 'qvalue')
#' @param top_count Number of top SNPs to select
#' @param top_percent Percentage of top SNPs to select
#' @return A list containing various results from the QR screen
#' @importFrom tidyr separate
#' @importFrom dplyr %>% mutate select
#' @export
qr_screen <- function(
    X, Y, Z = NULL, tau.list = seq(0.05, 0.95, by = 0.05),
    screen_threshold = 0.05, screen_method = "qvalue", top_count = 10, top_percent = 15) {
  # Make sure quantreg is installed
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("To use this function, please install quantreg: https://cran.r-project.org/web/packages/quantreg/index.html")
  }
  p <- ncol(X)
  pvec <- rep(NA, p)
  ltau <- length(tau.list)
  quantile.pvalue <- matrix(NA, nrow = p, ncol = ltau, dimnames = list(colnames(X), paste("p_qr", tau.list, sep = "_")))
  quantile.zscore <- matrix(NA, nrow = p, ncol = ltau, dimnames = list(colnames(X), paste("zscore_qr", tau.list, sep = "_")))
  y <- as.matrix(Y)

  if (!is.null(Z)) {
    zz <- cbind(rep(1, nrow(y)), Z)
  } else {
    zz <- matrix(1, nrow = nrow(y), ncol = 1)
  }

  ranks_list <- lapply(tau.list, function(tau) {
    suppressWarnings(quantreg::rq.fit.br(zz, y, tau = tau)$dual - (1 - tau))
  })

  for (ip in 1:p) {
    x <- as.matrix(X[, ip])
    VN <- matrix(0, nrow = ltau, ncol = ltau)
    for (i in 1:ltau) {
      for (j in 1:ltau) {
        VN[i, j] <- min(tau.list[i], tau.list[j]) - tau.list[i] * tau.list[j]
      }
    }

    if (!is.null(Z)) {
      # xstar = lm(x ~ zz - 1)$residual
      xstar <- .lm.fit(zz, x)$residual
    } else {
      xstar <- x
    }

    SN <- NULL
    for (itau in 1:ltau) {
      Sn <- as.matrix(t(xstar) %*% ranks_list[[itau]])
      SN <- c(SN, Sn)
    }
    VN2 <- matrix(outer(VN, t(xstar) %*% xstar, "*"), nrow = ltau)
    z_score <- SN / sqrt(diag(VN2))
    pvalue1 <- pchisq(SN^2 / diag(VN2), 1, lower.tail = FALSE)
    names(pvalue1) <- tau.list
    quantile.pvalue[ip, ] <- pvalue1
    quantile.zscore[ip, ] <- z_score
    e <- solve(chol(VN2))
    SN2 <- t(e) %*% SN
    pvalue <- pchisq(sum(SN2^2), ltau, lower.tail = FALSE)
    pvec[ip] <- pvalue
  }

  pvec <- apply(quantile.pvalue, 1, pval_cauchy)

  if (screen_method == "fdr") {
    adjusted_pvalues <- p.adjust(pvec)
    method_col_name <- "fdr_p_qr"
    method_quantile_names <- paste0("fdr_p_qr_", tau.list)
    quantile_adjusted_pvalues <- apply(quantile.pvalue, 2, p.adjust)
  } else if (screen_method == "qvalue") {
    adjusted_pvalues <- compute_qvalues(pvec)
    method_col_name <- "qvalue_qr"
    method_quantile_names <- paste0("qvalue_qr_", tau.list)
    quantile_adjusted_pvalues <- apply(quantile.pvalue, 2, compute_qvalues)
  } else {
    stop("Invalid screen_method. Choose 'fdr' or 'qvalue'.")
  }

  # Ensure quantile_adjusted_pvalues is always a matrix (handles single-variant case)
  if (!is.matrix(quantile_adjusted_pvalues)) {
    quantile_adjusted_pvalues <- matrix(quantile_adjusted_pvalues, nrow = 1,
                                        dimnames = list(colnames(X), names(quantile_adjusted_pvalues)))
  }

  sig_SNP_threshold <- which(adjusted_pvalues < screen_threshold)
  sig_SNP_top_count <- order(adjusted_pvalues)[1:top_count]
  sig_SNP_top_percent <- order(adjusted_pvalues)[1:max(1, round(length(adjusted_pvalues) * top_percent / 100))]

  sig.SNPs_names <- colnames(X)[sig_SNP_threshold]
  sig.SNPs_names_top_count <- colnames(X)[sig_SNP_top_count]
  sig.SNPs_names_top_percent <- colnames(X)[sig_SNP_top_percent]
  phenotype_id <- colnames(y)[1]

  df_result <- data.frame(
    phenotype_id = phenotype_id,
    variant_id = colnames(X),
    p_qr = pvec
  )

  # Add quantile-specific p-values
  for (tau in tau.list) {
    df_result[[paste0("p_qr_", tau)]] <- quantile.pvalue[, paste0("p_qr_", tau)]
  }

  # Add overall q-value
  df_result[[method_col_name]] <- adjusted_pvalues

  # Add quantile-specific q-values
  for (tau in tau.list) {
    df_result[[paste0(method_col_name, "_", tau)]] <- quantile_adjusted_pvalues[, paste0("p_qr_", tau)]
  }

  # Add quantile-specific z-scores
  for (tau in tau.list) {
    df_result[[paste0("zscore_qr_", tau)]] <- quantile.zscore[, paste0("zscore_qr_", tau)]
  }

  # Split variant_id and reorder columns
  parsed <- parse_variant_id(df_result$variant_id)
  df_result <- df_result %>%
    mutate(chr = parsed$chrom, pos = parsed$pos, A2 = parsed$A2, A1 = parsed$A1)

  # Define the column order
  col_order <- c("chr", "pos", "A2", "A1", "phenotype_id", "variant_id", "p_qr")
  col_order <- c(col_order, paste0("p_qr_", tau.list))
  col_order <- c(col_order, method_col_name)
  col_order <- c(col_order, paste0(method_col_name, "_", tau.list))
  col_order <- c(col_order, paste0("zscore_qr_", tau.list))

  # Reorder the columns
  df_result <- df_result %>% select(all_of(col_order))

  return(list(
    df_result = df_result,
    tau_list = tau.list,
    quantile_pvalue = quantile.pvalue,
    quantile_zscore = quantile.zscore,
    pvec = pvec,
    adjusted_pvalues = adjusted_pvalues,
    sig_SNP_threshold = sig_SNP_threshold,
    sig.SNPs_names = sig.SNPs_names,
    sig_SNP_top_count = sig_SNP_top_count,
    sig_SNP_top_percent = sig_SNP_top_percent
  ))
}

#' Perform Clumping and Pruning
#' @param X Matrix of genotypes
#' @param qr_results Results from QR_screen
#' @param maf_list List of minor allele frequencies (optional)
#' @param ld_clump_r2 R-squared threshold for initial LD clumping based on pvalue
#' @param final_clump_r2 R-squared threshold for final LD clumping based on MAF
#' @return A list containing final SNPs and clumped SNPs
#' @importFrom bigstatsr FBM.code256
#' @importFrom bigsnpr snp_clumping
#' @export
multicontext_ld_clumping <- function(X, qr_results, maf_list = NULL, ld_clump_r2 = 0.2, final_clump_r2 = 0.8) {
  # Extract significant SNP names
  sig_SNPs_names <- qr_results$sig.SNPs_names

  # If no significant SNPs, return empty result
  if (length(sig_SNPs_names) == 0) {
    return(list(final_SNPs = NULL, clumped_SNPs = NULL, message = "No significant SNPs found"))
  }

  # Check if X only contains one SNP (one column)
  if (ncol(X) == 1) {
    print("Only one SNP in X. Skipping LD clumping.")
    final_SNPs <- colnames(X)
    return(list(final_SNPs = final_SNPs, clumped_SNPs = final_SNPs, message = "Only one SNP, no LD clumping performed"))
  }

  # Extract genotype matrix for significant SNPs
  G_all <- X # Only retain columns for significant SNPs

  # Convert the genotype matrix into FBM format
  code_vec <- c(0, 1, 2, rep(NA, 256 - 3))
  G_all <- FBM.code256(
    nrow = nrow(G_all),
    ncol = ncol(G_all),
    init = G_all,
    code = code_vec
  )

  # Parse SNP names to extract chromosome and position information
  parsed_snp_info <- do.call(rbind, strsplit(sig_SNPs_names, ":"))
  chr <- as.numeric(strip_chr_prefix(parsed_snp_info[, 1]))
  pos <- as.numeric(parsed_snp_info[, 2]) # Extract position

  # Step 1: Perform LD clumping for each tau based on p-values
  clumped_snp_list <- list()

  for (itau in seq_along(qr_results$tau_list)) {
    tau_name <- paste0("p_qr_", qr_results$tau_list[itau])

    # Extract p-values for the given quantile
    p_values_quantile <- qr_results$quantile_pvalue[, tau_name][qr_results$sig_SNP_threshold]
    log_p_values <- -log10(p_values_quantile) # Calculate log10 p-values

    # Perform LD clumping using p-values as S
    ind_clumped <- snp_clumping(
      G = G_all,
      infos.chr = chr,
      infos.pos = pos,
      S = log_p_values, # Use log10 p-values for each quantile
      thr.r2 = ld_clump_r2,
      size = 100 / ld_clump_r2
    ) # Window size of 500kb

    # Store clumping results for each quantile
    clumped_snp_list[[tau_name]] <- ind_clumped
  }

  # Step 2: Take the union of clumping results across all quantiles
  clumped_snp_union <- unique(unlist(clumped_snp_list)) # This is the SNP index, not the name
  clumped_SNPs_name <- sig_SNPs_names[clumped_snp_union]
  print(paste("Number of SNPs after union of clumping:", length(clumped_snp_union)))

  if (length(clumped_snp_union) == 1) {
    message("Only one SNP found in the union. Skipping LD pruning and returning the single SNP directly.")
    final_SNPs <- sig_SNPs_names[clumped_snp_union]
    return(list(final_SNPs = final_SNPs, clumped_SNPs = clumped_SNPs_name))
  }


  # Step 3: Sort results from union
  sorted_indices <- order(chr[clumped_snp_union], pos[clumped_snp_union])
  chr_sorted <- chr[clumped_snp_union][sorted_indices]
  pos_sorted <- pos[clumped_snp_union][sorted_indices]

  # Step 4: Initialize maf_values to NULL, and update only if maf_list is provided
  maf_values <- NULL
  if (!is.null(maf_list)) {
    maf_values <- maf_list[qr_results$sig_SNP_threshold][clumped_snp_union][sorted_indices]
  }
  G_union <- X[, clumped_snp_union, drop = FALSE][, sorted_indices]
  if (!inherits(G_union, "FBM")) {
    G_union <- FBM.code256(
      nrow = nrow(G_union),
      ncol = ncol(G_union),
      init = G_union,
      code = code_vec
    )
  }
  # Step 5: Perform final clumping using maf_values (if available), otherwise proceed with NULL S
  final_clumped <- snp_clumping(
    G = G_union,
    infos.chr = chr_sorted,
    infos.pos = pos_sorted,
    S = maf_values, # Use MAF values if provided, otherwise NULL
    thr.r2 = final_clump_r2,
    size = 100 / final_clump_r2
  ) # Final clumping

  # Restrict the genotype matrix to only the final clumped SNPs
  G_final_clumped <- G_union[, final_clumped, drop = FALSE] # Limit G to the final clumped SNPs
  # Get the final SNP names
  final_SNPs <- sig_SNPs_names[clumped_snp_union][sorted_indices][final_clumped]
  print(paste("Number of final SNPs after MAF-based clumping (if applied):", length(final_SNPs)))

  return(list(final_SNPs = final_SNPs, clumped_SNPs = clumped_SNPs_name))
}

#' Perform Quantile Regression Analysis to get beta
#' @param X Matrix of predictors
#' @param Y Matrix or vector of response variables
#' @param Z Matrix of covariates (optional)
#' @param tau_values Vector of quantiles to be analyzed
#' @return A data frame with QR coefficients for each quantile
#' @importFrom tidyr pivot_wider separate
#' @importFrom dplyr %>% mutate select
#' @export
perform_qr_analysis <- function(X, Y, Z = NULL, tau_values = seq(0.05, 0.95, by = 0.05)) {
  # Make sure quantreg is installed
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("To use this function, please install quantreg: https://cran.r-project.org/web/packages/quantreg/index.html")
  }
  # Convert Y and X to matrices if they aren't already
  pheno.mat <- as.matrix(Y)
  geno.mat <- as.matrix(X)

  # Initialize an empty result table to store results
  result_table <- data.frame(
    phenotype_id = character(),
    variant_id = character(),
    tau = numeric(),
    predictor_coef = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over each tau value to perform quantile regression
  for (tau in tau_values) {
    # Loop over each SNP/variant in geno.mat (X)
    for (n in 1:ncol(geno.mat)) {
      response <- pheno.mat # Y
      predictor <- geno.mat[, n] # X
      phenotype_id <- colnames(pheno.mat)
      variant_id <- colnames(geno.mat)[n]

      # Construct the design matrix based on whether Z is provided
      if (is.null(Z)) {
        # If no covariates, include intercept and predictor
        X_design <- cbind(1, predictor)
      } else {
        # If covariates are provided, include them in the design matrix
        X_design <- cbind(1, predictor, as.matrix(Z))
      }

      # Fit the quantile regression model using rq.fit.br
      mod <- suppressWarnings(quantreg::rq.fit.br(X_design, response, tau = tau))

      # Extract the coefficient for the predictor (second coefficient)
      predictor_coef <- mod$coefficients[2] # Coefficient for predictor

      # Create a row with the results and append to the result table
      row <- data.frame(
        phenotype_id = phenotype_id,
        variant_id = variant_id,
        tau = tau,
        predictor_coef = predictor_coef,
        stringsAsFactors = FALSE
      )
      result_table <- rbind(result_table, row)
    }
  }

  # Reshape result_table to a wide format, so each tau's results are in separate columns
  result_table_wide <- result_table %>%
    pivot_wider(
      id_cols = c(phenotype_id, variant_id),
      names_from = tau,
      values_from = predictor_coef,
      names_prefix = "coef_qr_"
    )
  parsed_ids <- parse_variant_id(result_table_wide$variant_id)
  result_table_wide <- result_table_wide %>%
    mutate(chr = parsed_ids$chrom, pos = parsed_ids$pos, A2 = parsed_ids$A2, A1 = parsed_ids$A1) %>%
    select(chr, pos, A2, A1, everything())

  # Return the wide format result table
  return(result_table_wide)
}

#' Filter Highly Correlated SNPs
#' @param X Matrix of genotypes
#' @param cor_thres Correlation threshold for filtering
#' @return A list containing filtered X matrix and filter IDs
corr_filter <- function(X, cor_thres = 0.8) {
  p <- ncol(X)

  # Use Rfast for faster correlation calculation if available
  if (requireNamespace("Rfast", quietly = TRUE)) {
    cor.X <- Rfast::cora(X, large = TRUE)
  } else {
    cor.X <- cor(X)
  }
  Sigma.distance <- as.dist(1 - abs(cor.X))
  fit <- hclust(Sigma.distance, method = "single")
  clusters <- cutree(fit, h = 1 - cor_thres)
  groups <- unique(clusters)
  ind.delete <- NULL
  X.new <- X
  filter.id <- c(1:p)
  for (ig in seq_along(groups)) {
    temp.group <- which(clusters == groups[ig])
    if (length(temp.group) > 1) {
      ind.delete <- c(ind.delete, temp.group[-1])
    }
  }
  ind.delete <- unique(ind.delete)
  if ((length(ind.delete)) > 0) {
    X.new <- as.matrix(X[, -ind.delete])
    filter.id <- filter.id[-ind.delete]
  }

  # Check if X.new has only one column and ensure column names are preserved
  if (ncol(X.new) == 1) {
    colnames(X.new) <- colnames(X)[-ind.delete]
  }

  return(list(X.new = X.new, filter.id = filter.id))
}



#' Check and Remove Problematic Columns to Ensure Full Rank
#'
#' This function checks for problematic columns in the design matrix that cause it to be not full rank,
#' and iteratively removes them based on the chosen strategy until the matrix is full rank.
#'
#' @param X Matrix of SNPs
#' @param C Matrix of covariates (unnamed)
#' @param strategy The strategy for removing problematic columns ("variance", "correlation", or "response_correlation")
#' @param response Optional response vector for "response_correlation" strategy
#' @param max_iterations Maximum number of iterations to attempt removing problematic columns
#' @return Cleaned matrix X with problematic columns removed
#' @noRd
check_remove_highcorr_snp <- function(X = X, C = C, strategy = c("correlation", "variance", "response_correlation"), response = NULL, max_iterations = 300, corr_thresholds = seq(0.75, 0.5, by = -0.05)) {
  strategy <- match.arg(strategy)
  original_colnames <- colnames(X)
  initial_ncol <- ncol(X) # Store the initial number of columns in X
  iteration <- 0
  # Combine the design matrix with X (SNPs) and C (covariates), keeping C without column names
  X_design <- cbind(1, X, C) # Add an intercept column (1)
  colnames_X_design <- c("Intercept", colnames(X)) # Assign column names only to X (SNPs) part

  # Assign column names only to the X part, leaving C without names
  colnames(X_design)[1:(length(colnames_X_design))] <- colnames_X_design

  # Check the initial rank of the design matrix
  matrix_rank <- qr(X_design)$rank
  message("Initial rank of the design matrix: ", matrix_rank, " / ", ncol(X_design), " columns.")

  # Skip remove_highcorr_snp if removing all problematic columns doesn't achieve full rank
  skip_remove_highcorr <- FALSE

  # First check: Try removing all problematic columns at once
  if (matrix_rank < ncol(X_design)) {
    message("Design matrix is not full rank, identifying all problematic columns...")

    # QR decomposition to identify linearly dependent columns
    qr_decomp <- qr(X_design)
    R <- qr_decomp$rank
    Q <- qr_decomp$pivot

    # Get all problematic columns
    problematic_cols <- Q[(R + 1):ncol(X_design)]
    problematic_colnames <- colnames(X_design)[problematic_cols]
    problematic_colnames <- problematic_colnames[problematic_colnames %in% colnames(X)]

    if (length(problematic_colnames) > 0) {
      message("Attempting to remove all problematic columns at once: ", paste(problematic_colnames, collapse = ", "))

      # Remove all problematic columns at once
      X_temp <- X[, !(colnames(X) %in% problematic_colnames), drop = FALSE]
      X_design_temp <- cbind(1, X_temp, C)
      colnames_X_design_temp <- c("Intercept", colnames(X_temp))
      colnames(X_design_temp)[1:length(colnames_X_design_temp)] <- colnames_X_design_temp

      # Check if removing all problematic columns achieves full rank
      matrix_rank_temp <- qr(X_design_temp)$rank

      if (matrix_rank_temp == ncol(X_design_temp)) {
        message("Achieved full rank by removing all problematic columns at once. Proceeding with original logic...")
      } else {
        message("Removing all problematic columns did not achieve full rank. Skipping to corr_filter...")
        skip_remove_highcorr <- TRUE
      }
    }
  }

  # Only proceed with remove_highcorr_snp if not skipping
  if (!skip_remove_highcorr) {
    while (matrix_rank < ncol(X_design) && iteration < max_iterations) {
      message("Design matrix is not full rank, identifying problematic columns...")

      qr_decomp <- qr(X_design)
      R <- qr_decomp$rank
      Q <- qr_decomp$pivot
      problematic_cols <- Q[(R + 1):ncol(X_design)]
      problematic_colnames <- colnames(X_design)[problematic_cols]
      problematic_colnames <- problematic_colnames[problematic_colnames %in% colnames(X)]

      if (length(problematic_colnames) == 0) {
        message("No more problematic SNP columns found in X. Breaking the loop.")
        break
      }

      message("Problematic SNP columns identified: ", paste(problematic_colnames, collapse = ", "))
      X <- remove_highcorr_snp(X, problematic_colnames, strategy = strategy, response = response)

      X_design <- cbind(1, X, C)
      colnames_X_design <- c("Intercept", colnames(X))
      colnames(X_design)[1:length(colnames_X_design)] <- colnames_X_design

      matrix_rank <- qr(X_design)$rank
      message("New rank of the design matrix: ", matrix_rank, " / ", ncol(X_design), " columns.")
      iteration <- iteration + 1
    }

    if (iteration == max_iterations) {
      warning("Maximum iterations reached. The design matrix may still not be full rank.")
    }
  }

  # Final check and corr_filter if needed
  matrix_rank <- qr(cbind(1, X, C))$rank
  if (matrix_rank < ncol(cbind(1, X, C))) {
    message("Applying corr_filter to ensure the design matrix is full rank...")
    for (threshold in corr_thresholds) {
      filter_result <- corr_filter(X, cor_thres = threshold)
      X <- filter_result$X.new
      X_design <- cbind(1, X, C)
      colnames_X_design <- c("Intercept", colnames(X))
      colnames(X_design)[1:length(colnames_X_design)] <- colnames_X_design

      matrix_rank <- qr(X_design)$rank
      message("Rank after corr_filter with threshold ", threshold, ": ", matrix_rank, " / ", ncol(X_design), " columns.")
      if (matrix_rank == ncol(X_design)) {
        break
      }
    }
  }


  if (iteration == max_iterations) {
    warning("Maximum iterations reached. The design matrix may still not be full rank.")
  }

  if (ncol(X) == 1 && initial_ncol == 1) {
    colnames(X) <- original_colnames
  }
  return(X)
}

#' Remove Problematic Columns Based on a Given Strategy
#'
#' This function removes problematic columns from a matrix based on different strategies, such as smallest variance,
#' highest correlation, or lowest correlation with the response variable.
#'
#' @param X Matrix of SNPs
#' @param problematic_cols A vector of problematic columns to be removed
#' @param strategy The strategy for removing problematic columns ("variance", "correlation", or "response_correlation")
#' @param response Optional response vector for "response_correlation" strategy
#' @return Cleaned matrix X with the selected column removed
#' @importFrom stats var cor
#' @noRd
remove_highcorr_snp <- function(X, problematic_cols, strategy = c("correlation", "variance", "response_correlation"), response = NULL) {
  # Set default strategy
  strategy <- match.arg(strategy)

  message("Identified problematic columns: ", paste(problematic_cols, collapse = ", "))

  if (length(problematic_cols) == 0) {
    return(X) # If there are no problematic columns, return as is
  }

  if (length(problematic_cols) == 1) {
    message("Only one problematic column: ", problematic_cols)
    col_to_remove <- problematic_cols[1]
    message("Removing column: ", col_to_remove)
    X <- X[, !(colnames(X) %in% col_to_remove), drop = FALSE]
    # If X only has one column left after removal, ensure its column name is preserved
    if (ncol(X) == 1) {
      colnames(X) <- colnames(X)[colnames(X) != col_to_remove] # Preserve remaining SNP name
    }
    return(X)
  }

  # Choose columns to remove based on the strategy
  if (strategy == "variance") {
    # Strategy 1: Remove the column with the smallest variance
    variances <- apply(X[, problematic_cols, drop = FALSE], 2, var)
    col_to_remove <- problematic_cols[which.min(variances)]
    message("Removing column with the smallest variance: ", col_to_remove)
  } else if (strategy == "correlation") {
    # Strategy 2: Remove the column with the highest sum of absolute correlations
    cor_matrix <- abs(cor(X[, problematic_cols, drop = FALSE])) # Calculate absolute correlation matrix
    diag(cor_matrix) <- 0 # Ignore the diagonal (self-correlation)

    if (length(problematic_cols) == 2) {
      # If there are only two problematic columns, randomly remove one
      col_to_remove <- sample(problematic_cols, 1)
      message("Only two problematic columns, randomly removing: ", col_to_remove)
    } else {
      # Calculate sum of absolute correlations for each column
      cor_sums <- colSums(cor_matrix)
      col_to_remove <- problematic_cols[which.max(cor_sums)] # Remove the column with the largest sum of correlations
      message("Removing column with highest sum of absolute correlations: ", col_to_remove)
    }
  } else if (strategy == "response_correlation" && !is.null(response)) {
    # Strategy 3: Remove the column with the lowest correlation with the response variable
    # FIXME: This strategy is potentially biased based on corr of response and variants
    cor_with_response <- apply(X[, problematic_cols, drop = FALSE], 2, function(col) cor(col, response))
    col_to_remove <- problematic_cols[which.min(abs(cor_with_response))]
    message("Removing column with lowest correlation with the response: ", col_to_remove)
  } else {
    stop("Invalid strategy or missing response variable for 'response_correlation' strategy.")
  }

  # Remove the selected column from X
  X <- X[, !(colnames(X) %in% col_to_remove), drop = FALSE]
  if (ncol(X) == 1) {
    colnames(X) <- colnames(X)[colnames(X) != col_to_remove] # Preserve remaining SNP name
  }
  return(X)
}

#' Calculate QR Coefficients and Pseudo R-squared Across Multiple Quantiles
#'
#' This function calculates quantile regression coefficients and pseudo R-squared values across multiple quantiles,
#' while handling problematic columns that might affect the rank of the design matrix.
#'
#' @param AssocData List containing X, Y, C, and X.filter
#' @param tau.list Vector of quantiles to be analyzed
#' @param strategy The strategy for removing problematic columns ("variance", "correlation", or "response_correlation")
#' @return A list containing the cleaned X matrix, beta matrix as twas weight, and pseudo R-squared values
#' @noRd
calculate_qr_and_pseudo_R2 <- function(AssocData, tau.list, strategy = c("correlation", "variance", "response_correlation"),
                                       corr_thresholds = seq(0.75, 0.5, by = -0.05)) {
  # Make sure quantreg is installed
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("To use this function, please install quantreg: https://cran.r-project.org/web/packages/quantreg/index.html")
  }
  strategy <- match.arg(strategy)
  # Check and handle problematic columns affecting the full rank of the design matrix
  AssocData$X.filter <- check_remove_highcorr_snp(X = AssocData$X.filter, C = AssocData$C, strategy = strategy, response = AssocData$Y, corr_thresholds = corr_thresholds)
  snp_names <- colnames(AssocData$X.filter)
  # Build the cleaned design matrix using the filtered X and unnamed C

  # Fit the models for all tau values
  message("Start fitting full model for all taus...")
  fit_full <- suppressWarnings(quantreg::rq(Y ~ X.filter + C, tau = tau.list, data = AssocData))
  message("Finished fitting full model. Start fitting intercept-only model for all taus...")
  fit_intercept <- suppressWarnings(quantreg::rq(AssocData$Y ~ 1, tau = tau.list, data = AssocData))
  message("Finished fitting intercept-only model.")
  # Define the rho function for pseudo R2 calculation
  rho <- function(u, tau) {
    u * (tau - (u < 0))
  }

  # Prepare to store the pseudo R2 results
  pseudo_R2 <- numeric(length(tau.list))
  names(pseudo_R2) <- tau.list

  # Calculate pseudo R2 for each tau
  for (i in seq_along(tau.list)) {
    tau <- tau.list[i]

    # Get residuals for the intercept-only and full models
    residuals0 <- residuals(fit_intercept, subset = i)
    residuals1 <- residuals(fit_full, subset = i)

    # Calculate and store pseudo R2 for each tau
    rho0 <- sum(rho(residuals0, tau))
    rho1 <- sum(rho(residuals1, tau))
    pseudo_R2[i] <- 1 - rho1 / rho0
  }

  # Extract the coefficients for the SNPs
  num_filter_vars <- ncol(AssocData$X.filter)
  beta_mat <- coef(fit_full)[2:(1 + num_filter_vars), , drop = FALSE]
  rownames_beta <- rownames(beta_mat)
  if (ncol(AssocData$X.filter) == 1) {
    rownames(beta_mat) <- snp_names
  } else {
    rownames_beta <- rownames(beta_mat)
    rownames(beta_mat) <- gsub("^X.filter", "", rownames_beta)
  }
  return(list(X.filter = AssocData$X.filter, beta_mat = beta_mat, pseudo_R2 = pseudo_R2))
}

#' Calculate Heterogeneity of Beta Coefficients Across Quantiles
#'
#' This function calculates the heterogeneity of beta coefficients across multiple quantiles for each variant_id.
#' Heterogeneity is computed as log(sd(beta) / abs(mean(beta))).
#'
#' @param rq_coef_result Data frame containing variant_id and QR coefficient columns
#' @return A data frame with variant_id and heterogeneity values
#' @noRd
calculate_coef_heterogeneity <- function(rq_coef_result) {
  # Identify all the columns starting with "coef_qr_" (quantile regression coefficient columns)
  coef_cols <- grep("^coef_qr_", colnames(rq_coef_result), value = TRUE)

  # Create a new data frame with variant_id and heterogeneity
  heterogeneity_result <- data.frame(
    variant_id = rq_coef_result$variant_id,
    coef_heter = apply(rq_coef_result[, coef_cols], 1, function(beta) {
      # Compute the mean and standard deviation, ignoring NAs
      beta_mean <- mean(beta, na.rm = TRUE)
      beta_sd <- sd(beta, na.rm = TRUE)

      # Handle the case where mean(beta) is 0 to avoid division by zero
      if (abs(beta_mean) == 0) {
        return(NA) # Return NA if mean is zero
      }

      # Compute the heterogeneity: log(sd(beta) / abs(mean(beta)))
      heterogeneity <- log(beta_sd / abs(beta_mean))
      return(heterogeneity)
    }),
    stringsAsFactors = FALSE
  )

  # Return only variant_id and heterogeneity
  return(heterogeneity_result)
}

#' Calculate Xi Correlation for QR Coefficients
#'
#' This function calculates the xi correlation coefficient and p-value for each variant,
#' measuring the functional dependence between tau values and QR coefficients.
#' Uses coef_qr_0.1 to coef_qr_0.9 (17 values, excluding 0.05 and 0.95).
#'
#' @param rq_coef_result Data frame containing variant_id and coef_qr_* columns
#' @param tau_range Numeric vector of tau values to use (default: seq(0.1, 0.9, by = 0.05))
#' @param min_valid Minimum number of valid (non-NA) coefficients required (default: 10)
#' @return A data frame with variant_id, xi, and xi_pval columns
#' @export
calculate_xi_correlation <- function(rq_coef_result, tau_range = seq(0.1, 0.9, by = 0.05), min_valid = 10) {
  if (!requireNamespace("XICOR", quietly = TRUE)) {
    stop("Package 'XICOR' is required for xi correlation calculation. Please install it.")
  }
  # Build column names for the specified tau range
  coef_col_names <- paste0("coef_qr_", tau_range)

  # Check which columns exist
  existing_cols <- coef_col_names[coef_col_names %in% colnames(rq_coef_result)]

  if (length(existing_cols) == 0) {
    warning("No coef_qr columns found in the specified tau range")
    return(data.frame(
      variant_id = rq_coef_result$variant_id,
      xi = NA_real_,
      xi_pval = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  # Extract tau values from existing columns
  existing_tau <- as.numeric(gsub("coef_qr_", "", existing_cols))

  # Calculate xi for each variant
  xi_results <- apply(rq_coef_result[, existing_cols, drop = FALSE], 1, function(coef_values) {
    # Get valid (non-NA) coefficients
    valid_indices <- !is.na(coef_values)
    valid_coefs <- as.numeric(coef_values[valid_indices])
    valid_tau <- existing_tau[valid_indices]

    # Check if enough valid values
    if (length(valid_coefs) < min_valid) {
      return(c(xi = NA_real_, xi_pval = NA_real_))
    }

    # Calculate xi correlation
    tryCatch({
      xicor_result <- XICOR::xicor(valid_tau, y = valid_coefs, pvalue = TRUE, method = "asymptotic")
      return(c(xi = xicor_result$xi, xi_pval = xicor_result$pval))
    }, error = function(e) {
      return(c(xi = NA_real_, xi_pval = NA_real_))
    })
  })

  # Convert to data frame
  xi_df <- data.frame(
    variant_id = rq_coef_result$variant_id,
    xi = xi_results["xi", ],
    xi_pval = xi_results["xi_pval", ],
    stringsAsFactors = FALSE
  )

  return(xi_df)
}

#' Quantile TWAS Weight Pipeline
#'
#' @param X Matrix of genotypes
#' @param Y Matrix or vector of phenotypes
#' @param Z Matrix of covariates (optional)
#' @param maf Vector of minor allele frequencies (optional)
#' @param region_id Name of the region being analyzed
#' @param quantile_qtl_tau_list Vector of quantiles for QTL analysis
#' @param quantile_twas_tau_list Vector of quantiles for TWAS analysis
#'
#'
#' @return A list containing various results from the TWAS weight pipeline:
#' \itemize{
#'   \item qr_screen_pvalue_df: Data frame with QR screening results: pavlue, qvalue and zscore.
#'   \item message: Any informational or warning messages.
#'   \item twas_variant_names: Names of variants used in TWAS weight calculation.
#'   \item rq_coef_df: Data frame with quantile regression coefficients.
#'   \item twas_weight: Matrix of TWAS weights.
#'   \item pseudo_R2: Vector of pseudo R-squared values.
#'   \item quantile_twas_prediction: Matrix of TWAS predictions.
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. QR screening to identify significant SNPs.
#' 2. Filtering of highly correlated SNPs.
#' 3. LD clumping and pruning(use filtered SNPs from step 1).
#' 4. Calculation of QR coefficients for selected SNPs(use filtered SNPs from step 3).
#' 5. Calculation of TWAS weights and pseudo R-squared values(use filtered SNPs from step 2).
#'
#' @examples
#' # Example usage:
#' # X <- matrix of genotypes
#' # Y <- vector of phenotypes
#' # Z <- matrix of covariates
#' # results <- quantile_twas_weight_pipeline(X, Y, Z, region_id = "GeneA")
#'
#' @export
quantile_twas_weight_pipeline <- function(X, Y, Z = NULL, maf = NULL, region_id = "",
                                          ld_reference_meta_file = NULL, twas_maf_cutoff = 0.01,
                                          ld_clumping = FALSE, ld_pruning = FALSE,
                                          screen_significant = FALSE,
                                          quantile_qtl_tau_list = seq(0.05, 0.95, by = 0.05),
                                          quantile_twas_tau_list = seq(0.01, 0.99, by = 0.01),
                                          screen_threshold = 0.05,
                                          xi_tau_range = seq(0.1, 0.9, by = 0.05),
                                          keep_variants = NULL,
                                          marginal_beta_calculate = TRUE,
                                          twas_weight_calculate = TRUE,
                                          qrank_screen_calculate = TRUE) {
  # Step 1: vQTL
  # Step 1-1: Calculate vQTL rank scores
  message("Step 0: Calculating vQTL rank scores for region ", region_id)
  num_tau_levels <- length(quantile_qtl_tau_list) # Convert tau.list to numeric count
  rank_score <- QUAIL_rank_score_pipeline(
    phenotype = Y,
    covariates = Z,
    num_tau_levels = num_tau_levels,
    method = "equal",
    num_cores = 1
  )
  message("vQTL rank scores calculated.")

  # Step 1-2: Run vQTL pipeline
  message("Step 0.5: Running vQTL analysis for rank scores in region ", region_id)
  vqtl_results <- QUAIL_pipeline(
    genotype = X,
    rank_score = rank_score,
    covariates = Z,
    phenotype_id = colnames(Y)[1]
  )
  message("vQTL analysis completed. Proceeding to QR screen.")

  # Initialize results list
  results <- list(
    vqtl_results = vqtl_results
  )

  if (qrank_screen_calculate) {
    # Step 2: QR screen
    message("Starting QR screen for region ", region_id)
    p.screen <- qr_screen(X = X, Y = Y, Z = Z, tau.list = quantile_qtl_tau_list, screen_threshold = screen_threshold, screen_method = "qvalue", top_count = 10, top_percent = 15)
    message(paste0("Number of SNPs after QR screening: ", length(p.screen$sig_SNP_threshold)))
    message("QR screen completed. Screening significant SNPs")
    results$qr_screen_pvalue_df <- p.screen$df_result

    if (screen_significant && length(p.screen$sig_SNP_threshold) == 0) {
      results$message <- paste0("No significant SNPs detected in region ", region_id)
      return(results)
    }

    if (screen_significant) {
      X_filtered <- X[, p.screen$sig_SNP_threshold, drop = FALSE]
    } else {
      X_filtered <- X
    }

    # # Step 3: Optional LD clumping and pruning from results of QR_screen (using original QR screen results)
    if (ld_clumping) {
      message("Performing LD clumping and pruning from QR screen results...")
      LD_SNPs <- multicontext_ld_clumping(X = X[, p.screen$sig_SNP_threshold, drop = FALSE], qr_results = p.screen, maf_list = NULL)
      selected_snps <- if (ld_pruning) LD_SNPs$final_SNPs else LD_SNPs$clumped_SNPs
      x_clumped <- X[, p.screen$sig_SNP_threshold, drop = FALSE][, selected_snps, drop = FALSE]
    } else {
      message("Skipping LD clumping.")
    }

  } else {
    message("Skipping QR screen.")
  }

  # Determine whether to skip marginal beta calculation:
  # - skip if marginal_beta_calculate = FALSE
  # - skip if keep_variants is provided but empty (length 0)
  # - skip if qrank_screen_calculate = FALSE and keep_variants is NULL (no variants to select)
  skip_marginal_beta <- !marginal_beta_calculate ||
    (!is.null(keep_variants) && length(keep_variants) == 0) ||
    (!qrank_screen_calculate && is.null(keep_variants))

  if (!skip_marginal_beta) {
    # Step 4: Fit marginal QR to get beta with SNPs for quantile_qtl_tau_list values
    message("Fitting marginal QR for selected SNPs...")
    if (qrank_screen_calculate) {
      X_for_qr <- if (ld_clumping) x_clumped else X_filtered
      if (!is.null(keep_variants)) {
        variants_to_keep <- intersect(keep_variants, colnames(X_for_qr))
        if (length(variants_to_keep) > 0) {
          X_for_qr <- X_for_qr[, variants_to_keep, drop = FALSE]
          message("Filtered to ", ncol(X_for_qr), " variants from keep_variants list for QR analysis")
        } else {
          message("Warning: No variants from keep_variants found in selected SNPs, using all selected SNPs")
        }
      }
    } else {
      # qrank_screen_calculate = FALSE but keep_variants provided
      variants_to_keep <- intersect(keep_variants, colnames(X))
      if (length(variants_to_keep) > 0) {
        X_for_qr <- X[, variants_to_keep, drop = FALSE]
        message("Using ", ncol(X_for_qr), " variants from keep_variants list for QR analysis (QR screen skipped)")
      } else {
        message("Warning: No variants from keep_variants found in X, skipping marginal beta calculation")
        skip_marginal_beta <- TRUE
      }
    }
  }

  if (!skip_marginal_beta) {
    rq_coef_result <- perform_qr_analysis(X = X_for_qr, Y = Y, Z = Z, tau_values = quantile_qtl_tau_list)

    # Step 5: Heterogeneity calculation
    # Step 5-1: beta_heterogeneity index in marginal model
    message("Marginal QR for selected SNPs completed. Calculating beta heterogeneity...")
    beta_heterogeneity <- calculate_coef_heterogeneity(rq_coef_result)
    message("Beta heterogeneity calculation completed.")

    # Step 5-2: Calculate xi correlation (Chatterjee correlation test)
    message("Calculating xi correlation for QR coefficients...")
    xi_correlation <- calculate_xi_correlation(rq_coef_result, tau_range = xi_tau_range, min_valid = 10)
    message("Xi correlation calculation completed.")

    # Merge xi and xi_pval into rq_coef_result (using left_join to preserve row order)
    rq_coef_result <- rq_coef_result %>%
      dplyr::left_join(xi_correlation, by = "variant_id")

    results$rq_coef_df <- rq_coef_result
    results$beta_heterogeneity <- beta_heterogeneity
    results$xi_correlation <- xi_correlation
  } else {
    message("Skipping marginal beta calculation and heterogeneity analysis.")
  }

  if (twas_weight_calculate && qrank_screen_calculate) {
  # Step 6: Optional LD panel filtering and MAF filtering from results of QR_screen
  if (!is.null(ld_reference_meta_file)) {
    message("Starting LD panel filtering...")
    ld_result <- tryCatch(
      {
        variants_kept <- filter_variants_by_ld_reference(colnames(X_filtered), ld_reference_meta_file)
        if (length(variants_kept$data) == 0) NULL else variants_kept
      },
      error = function(e) {
        message("Error in LD filtering for region ", region_id, ": ", e$message)
        NULL
      }
    )

    if (is.null(ld_result)) {
      results$message <- paste0("No SNPs left or error in LD filtering in region ", region_id)
      return(results)
    }

    X_filtered <- X_filtered[, ld_result$data, drop = FALSE]
    message(paste0("Number of SNPs after LD filtering: ", ncol(X_filtered)))

    # MAF filtering
    if (!is.null(maf)) {
      maf_filtered <- maf[colnames(X_filtered)] > twas_maf_cutoff
      X_filtered <- X_filtered[, maf_filtered, drop = FALSE]

      # Check if any SNPs are left after MAF filtering
      if (ncol(X_filtered) == 0) {
        results$message <- paste0("No SNPs left after MAF filtering in region ", region_id)
        return(results)
      }

      message(paste0("Number of SNPs after MAF filtering: ", ncol(X_filtered)))
    }
  }

    # Step 7-9: Pre-filter by p_qr, corr_filter, then fit QR model
    # First attempt: p_qr < 0.05, corr_filter down to 0.5
    # If still singular: retry with p_qr < 0.01, corr_filter down to 0.2
    pqr_cutoffs <- list(
      list(pval = 0.05, corr_thresholds = seq(0.75, 0.5, by = -0.05)),
      list(pval = 0.01, corr_thresholds = seq(0.75, 0.2, by = -0.05))
    )

    qr_beta_R2_results <- NULL
    for (attempt in seq_along(pqr_cutoffs)) {
      pqr_cutoff <- pqr_cutoffs[[attempt]]$pval
      corr_thresholds <- pqr_cutoffs[[attempt]]$corr_thresholds

      # Step 7: Pre-filter variants by raw p_qr
  raw_pvals <- p.screen$pvec[colnames(X_filtered)]
      sig_raw <- which(raw_pvals < pqr_cutoff)
      if (length(sig_raw) == 0) {
        message("No variants with raw p_qr < ", pqr_cutoff, " in region ", region_id)
        next
      }
      if (length(sig_raw) < ncol(X_filtered)) {
        message("Pre-filtering variants by raw p_qr < ", pqr_cutoff, ": keeping ", length(sig_raw), " out of ", ncol(X_filtered), " variants")
        X_pqr_filtered <- X_filtered[, sig_raw, drop = FALSE]
      } else {
        X_pqr_filtered <- X_filtered
  }

  # Step 8: Filter highly correlated SNPs
  message("Filtering highly correlated SNPs...")
      if (ncol(X_pqr_filtered) > 1) {
        filtered <- corr_filter(X_pqr_filtered, 0.8)
    X.filter <- filtered$X.new
  } else {
        X.filter <- X_pqr_filtered
  }

  # Step 9: Fit QR and get twas weight and R2 for all taus
      message("Fitting full QR to calculate TWAS weights and pseudo R-squared values...")
  AssocData <- list(X = X, Y = Y, C = Z, X.filter = X.filter)
      qr_beta_R2_results <- tryCatch(
        {
          res <- calculate_qr_and_pseudo_R2(AssocData, quantile_twas_tau_list, corr_thresholds = corr_thresholds)
          res
        },
        error = function(e) {
          message("Attempt ", attempt, " failed (p_qr < ", pqr_cutoff, "): ", e$message)
          NULL
        }
      )

      if (!is.null(qr_beta_R2_results)) break
    }

    if (is.null(qr_beta_R2_results)) {
      results$message <- paste0("Failed to fit QR model after all attempts in region ", region_id)
      return(results)
    }

  X.filter <- qr_beta_R2_results$X.filter
  message("TWAS weights and pseudo R-squared calculations completed.")

  # Add additional results
  results$twas_variant_names <- colnames(X.filter)
  results$twas_weight <- qr_beta_R2_results$beta_mat
  results$pseudo_R2 <- qr_beta_R2_results$pseudo_R2
  results$quantile_twas_prediction <- X.filter %*% results$twas_weight
  } else {
    message("Skipping TWAS weight calculation.")
  }

  return(results)
}
