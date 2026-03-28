#' Perform Univariate Regression for Genotype Data
#'
#' @param X Numeric matrix of genotypes (n x p), where n is the number of samples and p is the number of SNPs.
#' @param y Numeric vector of phenotypes (length n).
#' @param Z Optional numeric matrix of covariates (n x k), where k is the number of covariates.
#' @param center Logical, whether to center the data (default: TRUE).
#' @param scale Logical, whether to scale the data (default: FALSE).
#' @param return_residuals Logical, whether to return residuals (default: FALSE).
#' @return A list containing regression results: \code{betahat}, \code{sebetahat}, \code{z_scores}, \code{p_values}, and \code{q_values}.
#' @examples
#' X <- matrix(rnorm(1000), 100, 10)
#' y <- rnorm(100)
#' results <- univariate_regression(X, y)
#' @noRd
univariate_regression <- function(X, y, Z = NULL, center = TRUE,
                                  scale = FALSE, return_residuals = FALSE) {
  y_na <- which(is.na(y))
  if (length(y_na)) {
    X <- X[-y_na, ]
    y <- y[-y_na]
  }

  if (center) {
    y <- y - mean(y)
    X <- scale(X, center = TRUE, scale = scale)
  } else {
    X <- scale(X, center = FALSE, scale = scale)
  }

  X[is.nan(X)] <- 0

  if (!is.null(Z)) {
    if (center) {
      Z <- scale(Z, center = TRUE, scale = scale)
    }
    residual_x <- X # Store X before modifying y
    y <- .lm.fit(Z, y)$residuals
  } else {
    residual_x <- X
  }

  output <- try(do.call(
    rbind,
    lapply(1:ncol(X), function(i) {
      g <- .lm.fit(cbind(1, X[, i]), y)
      return(c(coef(g)[2], calc_stderr(cbind(1, X[, i]), g$residuals)[2]))
    })
  ), silent = TRUE)

  # Exception handling
  if (inherits(output, "try-error")) {
    output <- matrix(0, ncol(X), 2)
    for (i in 1:ncol(X)) {
      fit <- summary(lm(y ~ X[, i]))$coef
      if (nrow(fit) == 2) {
        output[i, ] <- as.vector(summary(lm(y ~ X[, i]))$coef[2, 1:2])
      } else {
        output[i, ] <- c(0, 0)
      }
    }
  }

  # Calculate z-scores (t-statistics)
  z_scores <- output[, 1] / output[, 2]

  # Calculate p-values from z-scores
  p_values <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)

  # Calculate q-values using the provided function
  q_values <- compute_qvalues(p_values)

  # Use residual_x's column names as the column names for the results
  rownames(output) <- colnames(residual_x)

  result_list <- list(
    betahat = setNames(output[, 1], rownames(output)),
    sebetahat = setNames(output[, 2], rownames(output)),
    z_scores = setNames(z_scores, rownames(output)),
    p_values = setNames(p_values, rownames(output)),
    q_values = setNames(q_values, rownames(output))
  )

  if (return_residuals && !is.null(Z)) {
    result_list$residuals <- y
  }

  return(result_list)
}


#' Perform Linear Regression for GWAS
#'
#' @param genotype Numeric matrix of genotypes (n x p), where n is the number of samples and p is the number of SNPs.
#' @param phenotype Numeric vector of phenotypes (length n).
#' @param covariates Optional numeric matrix of covariates (n x k), where k is the number of covariates.
#' @return A data frame containing regression results for each SNP, including \code{BETA}, \code{SE}, \code{Z}, \code{P}, and \code{Q}.
#' @examples
#' genotype <- matrix(rbinom(1000 * 10, 2, 0.3), 1000, 10)
#' phenotype <- rnorm(1000)
#' results <- run_linear_regression1(genotype, phenotype)
#' @noRd
run_linear_regression <- function(genotype, phenotype, covariates = NULL, phenotype_id = NULL) {
  if (!is.null(covariates)) {
    covariates <- as.data.frame(lapply(covariates, as.numeric))
  }

  reg_results <- univariate_regression(
    X = genotype,
    y = phenotype,
    Z = covariates,
    center = TRUE,
    scale = FALSE
  )

  snp_info <- parse_variant_id(colnames(genotype))

  data.frame(
    phenotype_id = if (!is.null(phenotype_id)) phenotype_id else NA,
    chr = snp_info$chrom,
    pos = snp_info$pos,
    A1 = snp_info$A1,
    A2 = snp_info$A2,
    variant_id = colnames(genotype),
    beta = reg_results$betahat,
    se = reg_results$sebetahat,
    z = reg_results$z_scores,
    p = reg_results$p_values,
    q = reg_results$q_values,
    N = nrow(genotype)
  )
}

#' Main QUAIL pipeline
#' QUAIL vQTL Analysis Pipeline
#'
#' @param genotype numeric matrix (n x p) of genotypes.
#' @param rank_score numeric vector (n x 1) of rank scores from Step 1.
#' @param phenotype optional numeric vector (n x 1) of original phenotype values.
#' @param covariates optional numeric matrix (n x k) of covariates.
#' @return A data frame containing vQTL results.
#' @export
#' @examples
#' \dontrun{
#' results <- QUAIL_pipeline(genotype, rank_score, covariates = covariates)
#' }
QUAIL_pipeline <- function(genotype, rank_score, phenotype = NULL,
                           covariates = NULL, phenotype_id = NULL) {
  start_time <- Sys.time()

  # Validate rank_score
  if (!is.numeric(rank_score)) {
    stop("rank_score must be a numeric vector.")
  }

  # Validate covariates
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) {
      covariates <- as.data.frame(covariates)
    }
    covariates <- as.data.frame(lapply(covariates, as.numeric))
  }

  # Perform vQTL analysis
  vqtl_results <- run_linear_regression(genotype, rank_score, covariates, phenotype_id = phenotype_id)

  end_time <- Sys.time()
  cat("\nTotal vQTL runtime:", round(difftime(end_time, start_time, units = "secs"), 2), " seconds\n")

  return(vqtl_results)
}
