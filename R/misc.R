compute_qvalues <- function(pvalues) {
  # Make sure qvalue is installed
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop("To use this function, please install qvalue: https://www.bioconductor.org/packages/release/bioc/html/qvalue.html")
  }
  if (all(is.na(pvalues))) {
    message("All p-values are NA. Returning NA vector.")
    return(rep(NA_real_, length(pvalues)))
  }      
  tryCatch(
    {
      if (length(pvalues) < 2) {
        return(pvalues)
      } else {
        return(qvalue::qvalue(pvalues)$qvalues)
      }
    },
    error = function(e) {
      message("Too few p-values to calculate qvalue, fall back to BH")
      qvalue::qvalue(pvalues, pi0 = 1)$qvalues
    }
  )
}

matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

compute_maf <- function(geno) {
  f <- mean(geno, na.rm = TRUE) / 2
  return(min(f, 1 - f))
}

compute_missing <- function(geno) {
  miss <- sum(is.na(geno)) / length(geno)
  return(miss)
}

compute_non_missing_y <- function(y) {
  nonmiss <- sum(!is.na(y))
  return(nonmiss)
}

compute_all_missing_y <- function(y) {
  allmiss <- all(is.na(y))
  return(allmiss)
}

mean_impute <- function(geno) {
  f <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (i in seq_along(f)) geno[, i][which(is.na(geno[, i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) length(unique(x)) == 1

#' Safe truncated SVD with numerical stability
#'
#' Computes a thin SVD and optionally truncates small singular values.
#' Useful for avoiding numerical issues when working with rank-deficient
#' or near-singular matrices.
#'
#' @param mat Input matrix (n x p).
#' @param tol Relative tolerance for filtering singular values.
#'   Singular values smaller than \code{tol * max(d)} are discarded.
#'   Set to 0 to keep all singular values.
#' @param max_rank Optional maximum number of singular values to retain.
#'   If NULL, all singular values passing the tolerance filter are kept.
#' @return A list with components:
#'   \describe{
#'     \item{u}{Left singular vectors (n x r matrix).}
#'     \item{d}{Singular values (length-r numeric vector).}
#'     \item{v}{Right singular vectors (p x r matrix).}
#'   }
#'   where r is the number of retained singular values.
#' @noRd
safe_svd <- function(mat, tol = 1e-8, max_rank = NULL) {
  if (max(abs(mat)) == 0) {
    stop("Cannot compute SVD of an all-zero matrix.")
  }
  # Compute thin SVD
  s <- svd(mat)
  d <- s$d
  # Filter by relative tolerance
  if (tol > 0 && length(d) > 0) {
    keep <- d / d[1] > tol
    if (!any(keep)) {
      stop("All singular values are below the tolerance threshold.")
    }
  } else {
    keep <- rep(TRUE, length(d))
  }
  # Apply max_rank cap
  if (!is.null(max_rank) && max_rank > 0) {
    n_keep <- min(sum(keep), max_rank)
    keep_idx <- which(keep)
    if (length(keep_idx) > n_keep) {
      keep[keep_idx[(n_keep + 1):length(keep_idx)]] <- FALSE
    }
  }
  r <- sum(keep)
  list(
    u = s$u[, keep, drop = FALSE],
    d = d[keep],
    v = s$v[, keep, drop = FALSE]
  )
}

#' Compute LD (Linkage Disequilibrium) Correlation Matrix from Genotypes
#'
#' Computes a pairwise Pearson correlation matrix from a genotype matrix.
#' Supports three variance conventions:
#' \describe{
#'   \item{\code{"sample"}}{Standard sample variance with N-1 denominator (default).
#'     Uses mean imputation for missing genotypes, then \code{Rfast::cora} (if available)
#'     or base \code{cor()}.}
#'   \item{\code{"population"}}{Population variance with N denominator, matching
#'     GCTA-style tools (e.g. DENTIST, GCTA --make-grm). Per-SNP means are computed
#'     from non-missing values; missing entries are set to zero after centering so they
#'     do not contribute to cross-products. Cross-products are normalized by the total
#'     sample count N, not by pairwise non-missing counts.}
#'   \item{\code{"gcta"}}{GCTA per-pair missing data correction. Like \code{"population"}
#'     but applies a correction term for each SNP pair based on the number of jointly
#'     non-missing samples. Matches the exact formula from the DENTIST C++ binary's
#'     \code{calcLDFromBfile_gcta}. Use this when missingness varies substantially
#'     across SNPs and accuracy of individual LD entries matters.}
#' }
#'
#' @param X Numeric genotype matrix (samples x SNPs). May contain \code{NA}
#'   for missing genotypes.
#' @param method Character, one of \code{"sample"} (default, N-1 denominator),
#'   \code{"population"} (N denominator, GCTA-style), or \code{"gcta"} (per-pair
#'   missing data correction). Partial matching is supported.
#' @param backend Character, one of \code{"internal"} (default), \code{"snprelate"},
#'   or \code{"snpstats"}. Controls which library computes the correlation matrix
#'   when \code{method = "sample"}:
#'   \describe{
#'     \item{\code{"internal"}}{Uses \code{Rfast::cora} if available, otherwise
#'       base \code{cor()}.}
#'     \item{\code{"snprelate"}}{Requires a temporary GDS file; uses
#'       \code{SNPRelate::snpgdsLDMat(method = "corr")}.}
#'     \item{\code{"snpstats"}}{Converts to \code{SnpMatrix}; uses
#'       \code{snpStats::ld(, stat = "R")}.}
#'   }
#'   The \code{"snprelate"} and \code{"snpstats"} backends are only supported
#'   with \code{method = "sample"}; combining them with other methods will
#'   raise an error.
#' @param trim_samples Logical. If \code{TRUE} and \code{method} is
#'   \code{"population"} or \code{"gcta"}, drops trailing samples so that
#'   \code{nrow(X)} is a multiple of 4, matching PLINK .bed file chunk processing.
#'   Ignored when \code{method = "sample"}. Default is \code{FALSE}.
#' @param shrinkage Numeric in (0, 1]. Shrink the LD matrix toward the identity:
#'   \code{R_s = (1 - shrinkage) * R + shrinkage * I}. Useful for regularizing
#'   LD for summary-statistics-based methods such as lassosum (Mak et al 2017).
#'   Default is 0 (no shrinkage).
#'
#' @return A symmetric correlation matrix with row and column names taken from
#'   \code{colnames(X)}.
#'
#' @details
#' \strong{Missing data handling.}
#' With \code{method = "sample"}, missing values are mean-imputed per SNP
#' before computing the full Pearson correlation matrix.
#' With \code{method = "population"}, per-SNP means are computed from
#' non-missing values, the matrix is centered, then \code{NA}s are set to 0
#' so that missing pairs contribute nothing to the cross-product.
#' The denominator is always the total sample count \code{N}
#' (after optional trimming), matching the original GCTA formula:
#' \deqn{\text{Var}(X_i) = E[X_i^2] - E[X_i]^2}
#' \deqn{\text{Cor}(X_i, X_j) = \frac{\text{Cov}(X_i, X_j)}{\sqrt{\text{Var}(X_i)\,\text{Var}(X_j)}}}
#'
#' \strong{Zero-variance SNPs.}
#' Any monomorphic SNP will have zero variance, producing \code{NaN}
#' correlations. These are set to 0 in the returned matrix; the diagonal
#' is forced to 1.
#'
#' @examples
#' \dontrun{
#' X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 50)
#' colnames(X) <- paste0("rs", 1:10)
#'
#' # Standard sample correlation (default)
#' R1 <- compute_LD(X)
#'
#' # GCTA-style population variance
#' R2 <- compute_LD(X, method = "population")
#'
#' # GCTA-style with per-pair missing data correction
#' R3 <- compute_LD(X, method = "gcta")
#' }
#'
#' @export
compute_LD <- function(X, method = c("sample", "population", "gcta"),
                       backend = c("internal", "snprelate", "snpstats"),
                       trim_samples = FALSE, shrinkage = 0) {
  if (is.null(X)) {
    stop("X must be provided.")
  }
  method <- match.arg(method)
  backend <- match.arg(backend)
  nms <- colnames(X)

  if (method == "sample") {
    # ---- Standard sample correlation (N-1 denominator) ----
    if (backend == "snprelate") {
      R <- .compute_ld_snprelate(X)
    } else if (backend == "snpstats") {
      R <- .compute_ld_snpstats(X)
    } else {
      # internal backend: Rfast::cora if available, else base cor()
      # Mean impute only if NAs exist (PLINK2 data typically has none)
      X_imp <- X
      if (anyNA(X_imp)) {
        col_means <- colMeans(X_imp, na.rm = TRUE)
        na_pos <- which(is.na(X_imp), arr.ind = TRUE)
        X_imp[na_pos] <- col_means[na_pos[, 2]]
      }
      if (requireNamespace("Rfast", quietly = TRUE)) {
        # large=FALSE uses tcrossprod internally, ~40x faster than large=TRUE
        R <- Rfast::cora(X_imp, large = FALSE)
      } else {
        R <- cor(X_imp)
      }
    }
  } else if (method == "population") {
    if (backend != "internal") {
      stop("backend '", backend, "' is only supported with method='sample'.")
    }
    # ---- Population variance (N denominator, GCTA-style) ----
    # Optionally trim trailing samples to a multiple of 4 (matches .bed processing)
    if (trim_samples) {
      N_kept <- (nrow(X) %/% 4L) * 4L
      if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
    }
    N <- nrow(X)
    # Per-SNP means from non-missing values
    col_means <- colMeans(X, na.rm = TRUE)
    # Population variance: E[X^2] - E[X]^2
    col_vars <- colMeans(X^2, na.rm = TRUE) - col_means^2
    # Center; set NA -> 0 so missing pairs don't contribute to cross-products.
    # NOTE: the covariance divides by total N (not pairwise non-missing counts),
    # which is an approximation that assumes uniform missingness across SNPs.
    # With heterogeneous missingness, correlations between high-missing and
    # low-missing columns will be slightly deflated. This matches the GCTA
    # convention and is standard for PLINK-style LD computation.
    if (anyNA(X)) {
      na_rates <- colMeans(is.na(X))
      if (max(na_rates) - min(na_rates) > 0.1) {
        warning("Population LD method with heterogeneous missingness ",
                "(max NA rate ", round(max(na_rates), 3),
                ", min ", round(min(na_rates), 3),
                "): correlations may be biased. Consider using method='sample' ",
                "which handles missingness via mean imputation.")
      }
    }
    X_c <- sweep(X, 2, col_means)
    X_c[is.na(X_c)] <- 0
    # Covariance with N denominator
    cov_mat <- crossprod(X_c) / N
    # Correlation
    sd_vec <- sqrt(col_vars)
    R <- cov_mat / outer(sd_vec, sd_vec)
  } else {
    if (backend != "internal") {
      stop("backend '", backend, "' is only supported with method='sample'.")
    }
    # ---- GCTA per-pair missing data correction ----
    # Matches the DENTIST binary's calcLDFromBfile_gcta formula exactly.
    # Unlike "population" which divides by total N, this method tracks
    # per-pair missing counts and applies a correction term.
    if (trim_samples) {
      N_kept <- (nrow(X) %/% 4L) * 4L
      if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
    }
    N <- nrow(X)
    p <- ncol(X)

    # Marginal statistics from non-missing values
    col_means <- colMeans(X, na.rm = TRUE)
    col_mean_sq <- colMeans(X^2, na.rm = TRUE)
    col_vars <- col_mean_sq - col_means^2

    # Build indicator matrix for non-missing values
    not_na <- !is.na(X)
    # Replace NA with 0 for cross-product computation
    X_zero <- X
    X_zero[is.na(X_zero)] <- 0

    # Per-pair non-missing counts: not_na'not_na gives count of jointly observed
    pair_counts <- crossprod(not_na * 1.0)
    n_missing <- N - pair_counts

    # Per-pair sums: sum of X_i over samples where both i and j are observed
    # For the correction term we need E_i2 = sum_i_pair / N (pair-specific mean)
    # X_zero' %*% not_na gives, for each (i,j), sum of X_i where j is not missing
    pair_sums <- crossprod(X_zero, not_na * 1.0)

    # Cross-product sum: sum(X_i * X_j) over jointly non-missing samples
    sum_XY <- crossprod(X_zero)

    # GCTA correction formula:
    # E_i2[i,j] = pair_sums[i,j] / N  (mean of SNP i restricted to non-missing-j samples, divided by N)
    # cov = sum_XY/N + E[i]*E[j]*(N-m)/N - E[i]*E_j2 - E_i2*E[j]
    E_i2 <- pair_sums / N  # p x p: row i, col j = sum of X_i where j non-missing, / N
    E_j2 <- t(E_i2)        # transposed version

    cov_mat <- sum_XY / N +
      outer(col_means, col_means) * (pair_counts / N) -
      col_means * E_j2 -
      E_i2 * rep(col_means, each = p)

    # Correlation
    sd_vec <- sqrt(col_vars)
    sd_outer <- outer(sd_vec, sd_vec)
    R <- matrix(0.001, p, p)
    valid <- sd_outer > 0
    R[valid] <- cov_mat[valid] / sd_outer[valid]
  }

  # Ensure clean output
  diag(R) <- 1.0
  R[is.na(R) | is.nan(R)] <- 0

  # Optional shrinkage toward identity: R_s = (1 - shrinkage) * R + shrinkage * I
  # Used e.g. by lassosum (Mak et al 2017) to regularize LD for RSS methods.
  if (shrinkage > 0 && shrinkage <= 1) {
    R <- (1 - shrinkage) * R + shrinkage * diag(nrow(R))
  }

  colnames(R) <- rownames(R) <- nms
  R
}

#' Compute LD via SNPRelate (creates a temporary GDS file from the dosage matrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix.
#' @noRd
.compute_ld_snprelate <- function(X) {
  if (!requireNamespace("SNPRelate", quietly = TRUE))
    stop("Package 'SNPRelate' is required for backend='snprelate'")
  if (!requireNamespace("gdsfmt", quietly = TRUE))
    stop("Package 'gdsfmt' is required for backend='snprelate'")

  tmp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(tmp_gds), add = TRUE)

  # Round to integer dosage for GDS (0/1/2)
  X_int <- round(X)
  storage.mode(X_int) <- "integer"
  X_int[is.na(X_int)] <- 3L  # GDS missing code

  snp_ids <- colnames(X) %||% seq_len(ncol(X))
  sample_ids <- rownames(X) %||% seq_len(nrow(X))

  SNPRelate::snpgdsCreateGeno(tmp_gds,
    genmat = X_int,
    sample.id = sample_ids,
    snp.id = snp_ids,
    snp.chromosome = rep(1L, ncol(X)),
    snp.position = seq_len(ncol(X)),
    snpfirstdim = FALSE
  )

  gds <- SNPRelate::snpgdsOpen(tmp_gds, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  ld_obj <- SNPRelate::snpgdsLDMat(gds, method = "corr",
                                    slide = -1, verbose = FALSE)
  ld_obj$LD
}

#' Compute LD via snpStats (converts dosage matrix to SnpMatrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix (r, not r²).
#' @noRd
.compute_ld_snpstats <- function(X) {
  if (!requireNamespace("snpStats", quietly = TRUE))
    stop("Package 'snpStats' is required for backend='snpstats'")

  # snpStats expects counts of the B allele as raw codes: 1=AA, 2=AB, 3=BB, 0=NA
  # pecotmr dosage is ALT count (0/1/2), so map: 0->1, 1->2, 2->3, NA->0
  X_raw <- round(X) + 1L
  X_raw[is.na(X) | X_raw < 1L] <- 0L
  X_raw[X_raw > 3L] <- 3L
  storage.mode(X_raw) <- "raw"
  sm <- new("SnpMatrix", X_raw)

  R <- as.matrix(snpStats::ld(sm, stats = "R", depth = ncol(X) - 1L))
  # snpStats::ld returns a sparse-like matrix; ensure full dense
  R[is.na(R)] <- 0
  diag(R) <- 1
  R
}

#' @importFrom matrixStats colVars
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh = 0, maf = NULL, X_variance = NULL) {
  tol_variants <- ncol(X)
  if (!is.null(missing_rate_thresh) && missing_rate_thresh < 1.0) {
    rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
    if (length(rm_col)) X <- X[, -rm_col, drop = FALSE]
  }

  # Check if non-NA values are valid genotypes before MAF filtering
  if (!is.null(maf_thresh) && maf_thresh > 0.0) {
    valid_genotypes <- all(sapply(1:ncol(X), function(i) {
      x <- X[!is.na(X[, i]), i]
      all(x %in% c(0, 1, 2))
    }))

    if (valid_genotypes || !is.null(maf)) {
      rm_col <- if (!is.null(maf)) which(maf <= maf_thresh) else which(apply(X, 2, compute_maf) <= maf_thresh)
      if (length(rm_col)) X <- X[, -rm_col, drop = FALSE]
    } else {
      message("Skipping MAF filtering as X does not appear to be 0/1/2 matrix, and no external MAF information is provided")
    }
  }

  rm_col <- which(apply(X, 2, is_zero_variance))
  if (length(rm_col)) X <- X[, -rm_col, drop = FALSE]
  X <- mean_impute(X)
  if (var_thresh > 0) {
    rm_col <- if (!is.null(X_variance)) which(X_variance < var_thresh) else which(colVars(X) < var_thresh)
    if (length(rm_col)) X <- X[, -rm_col, drop = FALSE]
  }
  n_dropped <- tol_variants - ncol(X)
  if (n_dropped > 0) {
    message(paste0(n_dropped, " out of ", tol_variants, " total variants dropped due to quality control on X matrix."))
  }
  return(X)
}

#' This function performing filters on X variants based on Y subjects for TWAS analysis. This function checks
#' whether the absence (NA) of certain subjects would lead to monomorphic in some variants in X after removing
#' of these subjects data from X.
#' @param missing_rate_thresh Maximum individual missingness cutoff.
#' @param maf_thresh Minimum minor allele frequency (MAF) cutoff.
#' @param var_thresh Minimum variance cutoff for a variant. Default is 0.
#' @param X_variance A vector of variance for X variants.
filter_X_with_Y <- function(X, Y, missing_rate_thresh, maf_thresh, var_thresh = 0, maf = NULL, X_variance = NULL) {
  tol_variants <- ncol(X)
  X <- filter_X(X, missing_rate_thresh, maf_thresh, var_thresh = var_thresh, maf = maf, X_variance = X_variance)
  drop_idx <- do.call(c, lapply(colnames(Y), function(context) {
    subjects_with_na_Y <- rownames(Y)[is.na(Y[, context])]
    X_temp <- X
    X_temp[subjects_with_na_Y, ] <- NA
    rm_col <- which(apply(X_temp, 2, function(x) is_zero_variance(na.omit(x))))
    return(unique(rm_col))
  }))
  drop_idx <- unique(sort(drop_idx))
  if (length(drop_idx)) X <- X[, -drop_idx, drop = FALSE]
  if (length(drop_idx) > 0) {
    message(paste0("Additional ", length(drop_idx), " variants dropped after considering missing data in Y matrix, with ", ncol(X), " variants left."))
  }
  return(X)
}

filter_Y <- function(Y, n_nonmiss) {
  rm_col <- which(apply(Y, 2, compute_non_missing_y) < n_nonmiss)
  if (length(rm_col)) Y <- Y[, -rm_col]
  rm_rows <- NULL
  if (is.matrix(Y)) {
    rm_rows <- which(apply(Y, 1, compute_all_missing_y))
    if (length(rm_rows)) Y <- Y[-rm_rows, ]
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(list(Y = Y, rm_rows = rm_rows))
}

# Retrieve a nested element from a list structure
#' @export
get_nested_element <- function(nested_list, name_vector) {
  if (is.null(name_vector)) {
    return(NULL)
  }
  name_vector <- name_vector[name_vector!='']
  current_element <- nested_list
  for (name in name_vector) {
    if (is.null(current_element[[name]])) {
      stop("Element not found in the list")
    }
    current_element <- current_element[[name]]
  }
  return(current_element)
}



#' Utility function to specify the path to access the target list item in a nested list, especially when some list layers
#' in between are dynamic or uncertain.
#' @export
find_data <- function(x, depth_obj, show_path = FALSE, rm_null = TRUE, rm_dup = FALSE, docall = c, last_obj = NULL) {
  depth <- as.integer(depth_obj[1])
  list_name <- if (length(depth_obj) > 1) depth_obj[2:length(depth_obj)] else NULL
  if (depth == 1 || depth == 0) {
    if (!is.null(list_name)) {
      if (list_name[1] %in% names(x)) {
        if (any(grepl("^[0-9]+$", list_name))) { # list names, indx name, list names
          second_depth <- which(grepl("^[0-9]+$", list_name))[1]
          data <- get_nested_element(x, list_name[1:second_depth[1] - 1])
          remaining_path <- list_name[second_depth:length(list_name)]
          return(find_data(data, remaining_path,
            show_path = show_path,
            rm_null = rm_null, rm_dup = rm_dup, last_obj = names(data)
          ))
        }
        return(get_nested_element(x, list_name))
      }
    } else {
      return(x)
    }
  } else if (is.list(x)) {
    result <- lapply(x, find_data,
      depth_obj = c(depth - 1, list_name), show_path = show_path,
      rm_null = rm_null, rm_dup = rm_dup, last_obj = names(x)
    )
    shared_list_names <- list()
    if (isTRUE(rm_null)) {
      result <- result[!sapply(result, is.null)]
      result <- result[!sapply(result, function(x) length(x) == 0)]
    }
    if (isTRUE(rm_dup)) {
      unique_result <- list()
      unique_counter <- 1
      for (i in seq_along(result)) {
        duplicate_found <- FALSE
        for (j in seq_along(unique_result)) {
          if (identical(result[[i]], unique_result[[j]])) {
            duplicate_found <- TRUE
            shared_list_names[[paste0("unique_list_", j)]] <- c(shared_list_names[[paste0("unique_list_", j)]], names(result)[i])
            break
          }
        }
        if (!duplicate_found) {
          unique_name <- paste0("unique_list_", unique_counter)
          unique_result[[names(result)[i]]] <- result[[i]]
          shared_list_names[[unique_name]] <- names(result)[i]
          unique_counter <- unique_counter + 1
        }
      }
      result <- unique_result
    }

    if (isTRUE(show_path)) {
      if (length(shared_list_names) > 0 & depth == 2) result$shared_list_names <- shared_list_names
      return(result) # Carry original list structure
    } else {
      flat_result <- do.call(docall, unname(result))
      if (length(shared_list_names) > 0 & depth == 2) {
        names(result) <- paste0("unique_list_", seq_along(result))
        result$shared_list_names <- shared_list_names
        return(result)
      } else {
        return(flat_result) # Only return values
      }
    }
  } else {
    message(paste0("list ", depth_obj[length(depth_obj)], " is not found in ", last_obj, ".  \n"))
  }
}


thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    ## Rscript
    path <- cmdArgs[match]
    path <- gsub("\\~\\+\\~", " ", path)
    return(normalizePath(sub(needle, "", path)))
  } else {
    ## 'source'd via R console
    return(sys.frames()[[1]]$ofile)
  }
}

load_script <- function() {
  fileName <- thisFile()
  return(ifelse(!is.null(fileName) && file.exists(fileName),
    readChar(fileName, file.info(fileName)$size), ""
  ))
}

#' Find Valid File Path
find_valid_file_path <- function(reference_file_path, target_file_path) {
  # Check if the reference file path exits
  try_reference <- function() {
    if (file.exists(reference_file_path)) {
      return(reference_file_path)
    } else {
      return(NULL)
    }
  }
  # Check if the target file path exists
  try_target <- function() {
    if (file.exists(target_file_path)) {
      return(target_file_path)
    } else {
      # If not, construct a new target path by combining the directory of the reference file path with the target file path
      target_full_path <- file.path(dirname(reference_file_path), target_file_path)
      if (file.exists(target_full_path)) {
        return(target_full_path)
      } else {
        return(NULL)
      }
    }
  }

  target_result <- try_target()
  if (!is.null(target_result)) {
    return(target_result)
  }

  reference_result <- try_reference()
  if (!is.null(reference_result)) {
    return(reference_result)
  }

  stop(sprintf(
    "Both reference and target file paths do not work. Tried paths: '%s' and '%s'",
    reference_file_path, file.path(dirname(reference_file_path), target_file_path)
  ))
}

find_valid_file_paths <- function(reference_file_path, target_file_paths) sapply(target_file_paths, function(x) find_valid_file_path(reference_file_path, x))

#' Filter a vector based on a correlation matrix
#'
#' This function filters a vector `z` based on a correlation matrix `LD` and a correlation threshold `rThreshold`.
#' It keeps only one element among those having an absolute correlation value greater than the threshold.
#'
#' @param z A numeric vector to be filtered.
#' @param LD A square correlation matrix with dimensions equal to the length of `z`.
#' @param rThreshold The correlation threshold for filtering.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{filteredZ}{The filtered vector `z` based on the correlation threshold.}
#'     \item{filteredLD}{The filtered matrix `LD` based on the correlation threshold.}
#'     \item{dupBearer}{A vector indicating the duplicate status of each element in `z`.}
#'     \item{corABS}{A vector storing the absolute correlation values of duplicates.}
#'     \item{sign}{A vector storing the sign of the correlation values (-1 for negative, 1 for positive).}
#'     \item{minValue}{The minimum absolute correlation value encountered.}
#'   }
#'
#' @examples
#' z <- c(1, 2, 3, 4, 5)
#' LD <- matrix(c(
#'   1.0, 0.8, 0.2, 0.1, 0.3,
#'   0.8, 1.0, 0.4, 0.2, 0.5,
#'   0.2, 0.4, 1.0, 0.6, 0.1,
#'   0.1, 0.2, 0.6, 1.0, 0.3,
#'   0.3, 0.5, 0.1, 0.3, 1.0
#' ), nrow = 5, ncol = 5)
#' rThreshold <- 0.5
#'
#' result <- find_duplicate_variants(z, LD, rThreshold)
#' print(result)
#'
#' @export
find_duplicate_variants <- function(z, LD, rThreshold) {
  p <- length(z)
  dupBearer <- rep(-1, p)
  corABS <- rep(0, p)
  sign <- rep(1, p)
  count <- 1
  minValue <- 1

  for (i in 1:(p - 1)) {
    if (dupBearer[i] != -1) next

    idx <- (i + 1):p
    corVec <- abs(LD[i, idx])
    dupIdx <- which(dupBearer[idx] == -1 & corVec > rThreshold)

    if (length(dupIdx) > 0) {
      j <- idx[dupIdx]
      sign[j] <- ifelse(LD[i, j] < 0, -1, sign[j])
      corABS[j] <- corVec[dupIdx]
      dupBearer[j] <- count
    }

    minValue <- min(minValue, min(corVec))
    count <- count + 1
  }

  # Filter z based on dupBearer
  filteredZ <- z[dupBearer == -1]
  filteredLD <- LD[dupBearer == -1, dupBearer == -1, drop = FALSE]

  return(list(filteredZ = filteredZ, filteredLD = filteredLD, dupBearer = dupBearer, corABS = corABS, sign = sign, minValue = minValue))
}

#' Convert Z-scores to Beta and Standard Error
#'
#' This function estimates the effect sizes (beta) and standard errors (SE) from
#' given z-scores, minor allele frequencies (MAF), and a sample size (n) in genetic studies.
#' It supports vector inputs for z-scores and MAFs to process multiple variants simultaneously.
#'
#' @param z Numeric vector. The z-scores of the genetic variants.
#' @param maf Numeric vector. The minor allele frequencies of the genetic variants (0 < maf <= 0.5).
#' @param n Integer. The sample size of the study (assumed to be the same for all variants).
#'
#' @return A data frame containing three columns:
#' \describe{
#'   \item{beta}{The estimated effect sizes.}
#'   \item{se}{The estimated standard errors.}
#'   \item{maf}{The input minor allele frequencies (possibly adjusted if > 0.5).}
#' }
#'
#' @details
#' The function uses the following formulas to estimate beta and SE:
#' Beta = z / sqrt(2p(1-p)(n + z^2))
#' SE = 1 / sqrt(2p(1-p)(n + z^2))
#' Where p is the minor allele frequency.
#'
#' @examples
#' z <- c(2.5, -1.8, 3.2, 0.7)
#' maf <- c(0.3, 0.1, 0.4, 0.05)
#' n <- 10000
#' result <- z_to_beta_se(z, maf, n)
#' print(result)
#' test_data_with_results <- cbind(test_data, results)
#' print(test_data_with_results)
#'
#' @note
#' This function assumes that the input z-scores are normally distributed and
#' that the genetic model is additive. It may not be accurate for rare variants
#' or in cases of imperfect imputation. The function automatically adjusts MAF > 0.5
#' to ensure it's always working with the minor allele.
#' @noRd
z_to_beta_se <- function(z, maf, n) {
  if (length(z) != length(maf)) {
    stop("z and maf must be vectors of the same length")
  }
  # Ensure MAF is the minor allele frequency
  p <- pmin(maf, 1 - maf)
  denominator <- sqrt(2 * p * (1 - p) * (n + z^2))
  beta <- z / denominator
  se <- 1 / denominator
  return(data.frame(beta = beta, se = se, maf = p))
}

#' Convert Z-scores to P-values
#'
#' This function calculates p-values from given z-scores using a two-tailed normal distribution.
#' It supports vector input to process multiple z-scores simultaneously.
#'
#' @param z Numeric vector. The z-scores to be converted to p-values.
#'
#' @return A numeric vector of p-values corresponding to the input z-scores.
#'
#' @details
#' The function uses the following formula to calculate p-values:
#' p-value = 2 * Phi(-|z|)
#' Where Phi is the cumulative distribution function of the standard normal distribution.
#'
#' @examples
#' z <- c(2.5, -1.8, 3.2, 0.7)
#' pvalues <- z_to_pvalue(z)
#' print(pvalues)
#'
#' @note
#' This function assumes that the input z-scores are from a two-tailed test and
#' are normally distributed. It calculates two-sided p-values.
#' For extremely large absolute z-scores, the resulting p-values may be computed as zero
#' due to floating-point limitations in R. This occurs when the absolute z-score > 37.
#'
#' @export
z_to_pvalue <- function(z) {
  2 * pnorm(-abs(z))
}
                                                                                 
#' Filter events based on provided context name pattern
#'       
#' @param events A character vector of event names 
#' @param filters A data frame with character column of type_pattern, valid_pattern, and exclude_pattern. 
#' @param condition Optional label context name 
#' @param remove_all_group Logical if \code{TRUE}, removes all events from the same group and character-defined context.
filter_molecular_events <- function(events, filters, condition = NULL, remove_all_group = FALSE) {
  # filters is a list of filter specifications
  # Each filter spec must have:
  #   type_pattern: pattern to identify event type
  #   And at least ONE of:
  #   valid_pattern: pattern that must exist in group
  #   exclude_pattern: pattern to exclude

  filtered_events <- events
  for (filter in filters) {
    if (is.null(filter$type_pattern) ||
      (is.null(filter$valid_pattern) && is.null(filter$exclude_pattern))) {
      stop("Each filter must specify type_pattern and at least one of valid_pattern or exclude_pattern")
    }
    # Get events of this type
    type_events <- filtered_events[grepl(filter$type_pattern, filtered_events)]
    type_events_all <- type_events
    if (length(type_events) == 0) next
    # Apply valid pattern if specified
    if (!is.null(filter$valid_pattern)) {
      filter$valid_pattern <- strsplit(filter$valid_pattern, ",")[[1]]
      valid_groups <- unique(gsub(
        filter$type_pattern, "\\1",
        type_events[grepl(paste(filter$valid_pattern, collapse = "|"), type_events)]
      ))
      if (length(valid_groups) > 0) {
        type_events <- type_events[grepl(paste(filter$valid_pattern, collapse = "|"), type_events)] # filter for valid pattern in type events
      } else {
        type_events <- character(0)
      }
    }
    # Apply exclusions if specified
    if (!is.null(filter$exclude_pattern)) {
      filter$exclude_pattern <- strsplit(filter$exclude_pattern, ",")[[1]]
      type_events <- type_events[!grepl(paste(filter$exclude_pattern, collapse = "|"), type_events)]
    }
    if (is.null(condition)) condition <- events
    if (length(type_events) == length(events)) {
      message(paste("All events matching", filter$type_pattern, "in", condition, "included in following analysis."))
    } else if (length(type_events) == 0) {
      message(paste("No events matching", filter$type_pattern, "in", condition, "pass the filtering."))
      return(NULL)
    } else {
      exclude_events <- paste0(setdiff(type_events_all, type_events), collapse = ";")
      message(paste("Some events,", exclude_events, "in", condition, "are removed. \n"))
      if (remove_all_group) {
        exclude_events <- setdiff(type_events_all, type_events)
        exclude_groups <- gsub(filter$type_pattern, "\\1", 
                               exclude_events[grepl(paste(filter$exclude_pattern, collapse = "|"), exclude_events)]
        )
        for (i in seq_along(exclude_events)) {
            #if (!any(grepl(exclude_groups[i], type_events))) next  # skip the event if the corresponding group is all removed
            for (x in filter$exclude_pattern) exclude_events[i] <- gsub(x, ".*", exclude_events[i]) # remove exclude pattern from the context
            context_key <- gsub("\\b\\d+\\b", "", exclude_events[i]) # remove stand alone numbers (strings such as "lf2" or "chr8" will be kept)
            # General pattern to match all events of same group ID and similar character structure
            pattern_to_remove <- paste0(".*", exclude_groups[i], ".*")
            # Identify all events that match both the context structure and group ID
            same_group_events <- type_events[grepl(pattern_to_remove, type_events) & grepl(gsub("\\d+", "", context_key), gsub("\\d+", "", type_events))]
            type_events <- setdiff(type_events, same_group_events)
        }
      }
    }
    # Update events list
    filtered_events <- unique(c(
      filtered_events[!grepl(filter$type_pattern, filtered_events)],
      type_events
    ))
  }

  return(filtered_events)
}


#' Robust Mahalanobis Distance
#'
#' Drop-in replacement for \code{\link[stats]{mahalanobis}} that handles
#' singular (rank-deficient) covariance matrices by falling back to the
#' Moore–Penrose pseudoinverse via \code{MASS::ginv}.
#'
#' @param x Numeric matrix (samples x features) or vector.
#' @param center Numeric vector of column means (length = number of features).
#'   If \code{NULL}, computed from \code{x}.
#' @param cov Covariance matrix. If \code{NULL}, computed from \code{x}.
#' @param inverted Logical; if \code{TRUE}, \code{cov} is already inverted.
#' @return Named numeric vector of Mahalanobis distances.
#' @importFrom MASS ginv
#' @importFrom stats cov quantile
#' @export
robust_mahalanobis <- function(x, center = NULL, cov = NULL,
                               inverted = FALSE) {
  x <- if (is.vector(x)) matrix(x, ncol = length(x)) else as.matrix(x)
  if (is.null(center)) center <- colMeans(x)
  if (is.null(cov)) cov <- cov(x)
  x <- sweep(x, 2L, center)
  if (!inverted) {
    cov <- tryCatch(solve(cov), error = function(cond) {
      ginv(cov)
    })
  }
  setNames(rowSums(x %*% cov * x), rownames(x))
}

#' Detect Outliers via Mahalanobis Distance
#'
#' Identifies outlier samples in a numeric matrix (e.g., PCA scores) using
#' Mahalanobis distance with chi-squared-based p-values. Useful for QC
#' in genotype PCA or expression PCA workflows.
#'
#' @param x Numeric matrix (samples x features). Rownames are used as
#'   sample IDs in the output.
#' @param prob Numeric in (0, 1); quantile threshold for the Mahalanobis
#'   distance cutoff (default 0.99).
#' @param pval_threshold P-value threshold for outlier classification
#'   (default 0.05). A sample is flagged only if its distance exceeds
#'   the quantile cutoff \emph{and} its p-value is below this threshold.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{sample_id}{Row names from \code{x}, or row indices if unnamed.}
#'     \item{mahal}{Mahalanobis distance.}
#'     \item{pvalue}{Chi-squared p-value (df = number of features).}
#'     \item{is_outlier}{Logical; TRUE if distance > quantile cutoff and
#'       p-value < \code{pval_threshold}.}
#'   }
#' @export
detect_outliers_mahalanobis <- function(x, prob = 0.99,
                                        pval_threshold = 0.05) {
  x <- as.matrix(x)
  sample_ids <- rownames(x) %||% as.character(seq_len(nrow(x)))
  center <- colMeans(x)
  cov_mat <- cov(x)
  d <- robust_mahalanobis(x, center, cov_mat)
  p <- ncol(x)
  pvals <- pchisq(d, df = p, lower.tail = FALSE)
  cutoff <- quantile(d, probs = prob)
  data.frame(
    sample_id = sample_ids,
    mahal = as.numeric(d),
    pvalue = pvals,
    is_outlier = (d > cutoff) & (pvals < pval_threshold),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
