#' Get Hierarchical Clusters with Modularity Optimization
#'
#' This function performs hierarchical clustering on a correlation matrix and identifies the optimal number of clusters using modularity optimization. It ensures continuous quantiles for each cluster and computes a modularity score.
#' FIXME: # The initial modularity (Q_modularity_initial) is calculated based on the clustering before we enforce continuous ranges in the clusters. After adjusting the cluster
#' assignments with `ensure_continuous_clusters`, we are not recalculating modularity for the new cluster structure. If required, you should compute the modularity
#' again based on the new cluster structure after ensuring continuous ranges.
#'
#' @param cormat A P x P correlation matrix (numeric matrix) where P is the number of items to be clustered.
#' @param between_cluster A numeric value specifying the correlation threshold for clustering. If all correlations are greater than this threshold, the function will form a single cluster. Default is 0.8.
#'
#' @return A list with the following elements:
#' \item{cluster}{A P x K binary matrix where K is the number of clusters. Each row represents an item, and a value of 1 in column k means that the item is assigned to cluster k.}
#' \item{Q_modularity_initial}{The modularity score based on the initial clustering, before enforcing continuous clusters.}
#' \item{cluster_ranges}{A list of ranges (quantiles) for each cluster.}
#'
#' @examples
#' # Example usage
#' cormat <- matrix(runif(100), nrow = 10)
#' cormat <- cormat %*% t(cormat) # Make it symmetric
#' res <- get_hierarchical_clusters(cormat, between_cluster = 0.8)
#' B <- res$cluster # The binary cluster matrix
#' Q <- res$Q_modularity_initial # The initial modularity score
#' ranges <- res$cluster_ranges # The quantile ranges for each cluster
#'
#' @importFrom stats hclust as.dist cutree
#'
#' @noRd
get_hierarchical_clusters <- function(cormat, between_cluster = 0.8) {
  # Perform hierarchical clustering on the correlation matrix
  hc <- hclust(as.dist(1 - cormat))

  # Get the optimal number of clusters
  opt_cluster <- get_n_cluster(hc, cormat, between_cluster = between_cluster)
  n_cluster <- opt_cluster$n_cluster
  Q_modularity_initial <- opt_cluster$Qmodularity # Save the initial modularity score

  # Obtain the initial clusters
  index <- cutree(hc, n_cluster)

  # Adjust the clusters to ensure each cluster has continuous quantiles
  index <- ensure_continuous_clusters(index)

  # Create the adjacency matrix B based on the adjusted continuous clusters
  B <- sapply(1:max(index), function(t) as.numeric(index == t))
  B <- as.matrix(B)

  # Get the range of quantiles for each cluster
  cluster_ranges <- get_cluster_ranges(index)

  # Return the cluster structure, initial modularity, and quantile ranges
  return(list(
    "cluster" = B,
    "Q_modularity_initial" = Q_modularity_initial,
    "cluster_ranges" = cluster_ranges
  ))
}

#' Ensure Clusters Have Continuous Quantiles
#'
#' This helper function ensures that clusters produced by hierarchical clustering are continuous in terms of quantiles.
#'
#' @param index A numeric vector representing the cluster assignments.
#' @return A numeric vector with adjusted cluster assignments ensuring continuity.
#' @examples
#' index <- c(1, 1, 2, 2, 3)
#' new_index <- ensure_continuous_clusters(index)
#' @noRd
ensure_continuous_clusters <- function(index) {
  n <- length(index)
  if (n <= 1) return(index)
  new_index <- integer(n)
  cluster_counter <- 1
  new_index[1] <- cluster_counter

  # Traverse the index to ensure continuous clusters
  for (i in 2:n) {
    # If the current quantile does not follow the previous one, start a new cluster
    if (index[i] != index[i - 1]) {
      cluster_counter <- cluster_counter + 1
    }
    new_index[i] <- cluster_counter
  }

  return(new_index)
}

#' Get Range of Quantiles for Each Cluster
#'
#' This function returns the range of quantiles for each identified cluster.
#'
#' @param index A numeric vector representing the cluster assignments.
#' @return A list of quantile ranges for each cluster.
#' @examples
#' index <- c(1, 1, 2, 2, 3)
#' cluster_ranges <- get_cluster_ranges(index)
#' @noRd
get_cluster_ranges <- function(index) {
  ranges <- list()
  cluster_numbers <- unique(index)

  # For each cluster, find the range of quantiles
  for (cluster in cluster_numbers) {
    members <- which(index == cluster)
    ranges[[cluster]] <- paste(min(members), "-", max(members))
  }

  return(ranges)
}

#' Get Optimal Number of Clusters Using Modularity
#'
#' This function computes the optimal number of clusters using modularity optimization. It iterates over possible numbers of clusters and evaluates each using a modularity score.
#'
#' @param hc An object of class \code{hclust}, resulting from hierarchical clustering.
#' @param Sigma A P x P correlation matrix.
#' @param m The maximum number of clusters to evaluate. Default is the number of columns in \code{Sigma}.
#' @param between_cluster A numeric value specifying the correlation threshold. Default is 0.8.
#' @return A list containing the number of clusters and the modularity score for each.
#' @examples
#' cormat <- matrix(runif(100), nrow = 10)
#' cormat <- cormat %*% t(cormat) # Make it symmetric
#' hc <- hclust(as.dist(1 - cormat))
#' result <- get_n_cluster(hc, cormat, between_cluster = 0.8)
#' n_cluster <- result$n_cluster
#' @noRd
get_n_cluster <- function(hc, Sigma, m = ncol(Sigma), between_cluster = 0.8) {
  # If the minimum correlation is greater than the threshold, use one cluster
  if (min(Sigma) > between_cluster) {
    IND <- 1
    Q <- 1
  } else {
    Q <- c()
    if (ncol(Sigma) < 10) {
      m <- ncol(Sigma)
    }
    # Iterate over possible numbers of clusters and compute modularity
    for (i in 1:m) {
      index <- cutree(hc, i)
      B <- sapply(1:i, function(t) as.numeric(index == t))
      Q[i] <- get_modularity(Sigma, B)
    }
    # Find the number of clusters with the maximum modularity
    IND <- which(Q == max(Q))
    L <- length(IND)
    if (L > 1) IND <- IND[1] # If multiple solutions, choose the first one
  }
  return(list("n_cluster" = IND, "Qmodularity" = Q))
}

#' Calculate Modularity
#'
#' This function computes the modularity score for a given cluster structure and weight matrix.
#'
#' @param Weight A P x P weight matrix.
#' @param B A binary matrix representing the cluster structure.
#' @return The modularity score (numeric).
#' @examples
#' W <- matrix(runif(100), nrow = 10)
#' B <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' Q <- get_modularity(W, B)
#' @noRd
get_modularity <- function(Weight, B) {
  if (dim(Weight)[1] == 1) {
    Q <- 0
  } else {
    W_pos <- Weight * (Weight > 0)
    W_neg <- Weight * (Weight < 0)
    N <- dim(Weight)[1]
    K_pos <- colSums(W_pos)
    K_neg <- colSums(W_neg)
    m_pos <- sum(K_pos)
    m_neg <- sum(K_neg)
    m <- m_pos + m_neg
    cate <- B %*% t(B)
    if (m_pos == 0 & m_neg == 0) {
      Q <- 0
    } else {
      if (m_pos == 0) {
        Q_positive <- 0
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      } else if (m_neg == 0) {
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- 0
      } else {
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      }
    }
    Q <- m_pos / m * Q_positive - m_neg / m * Q_negative
  }
  return(Q)
}


#' Perform Grouped Integration for TWAS Weights
#'
#' This function integrates TWAS weights using two methods: correlation-based clustering and fixed tau intervals. The weights are clustered and averaged within each cluster, and the corresponding performance metrics (e.g., pseudo R2) are calculated.
#'
#' @param twas_weight A matrix where rows represent variants and columns represent tau values. Each entry represents the TWAS weight for a specific variant and tau.
#' @param tau_values A numeric vector representing the tau values (quantiles) corresponding to the columns of \code{twas_weight}. Default is a sequence from 0.01 to 0.99.
#' @param pseudo_R2 A numeric vector of pseudo R2 values corresponding to the tau values.
#' @param between_cluster A numeric value specifying the correlation threshold for clustering. If all correlations are greater than this threshold, the function will form a single cluster. Default is 0.8.
#' @param num_intervals The number of fixed non-overlapping intervals to divide the tau values. Default is 3.
#'
#' @return A list containing:
#' \item{weights}{A matrix of integrated TWAS weights for each cluster or interval.}
#' \item{twas_weight_performance}{A list of performance metrics, including quantile ranges and average pseudo R2 values for each cluster/interval.}
#'
#' @examples
#' # Example usage
#' twas_weight <- matrix(runif(100), nrow = 10)
#' tau_values <- seq(0.01, 0.99, length.out = 10)
#' pseudo_R2 <- runif(10)
#' result <- perform_grouped_integration(twas_weight, tau_values, pseudo_R2)
#' weights <- result$weights
#' performance <- result$twas_weight_performance
#'
#' @noRd
perform_grouped_integration <- function(twas_weight, tau_values, pseudo_R2, between_cluster = 0.8, num_intervals = 3) {
  variant_id <- rownames(twas_weight)

  if (length(variant_id) > 1) {
    # Method 1: Correlation-based clustering (C1, C2, ..., Cn)
    cormat <- cor(twas_weight)
    res <- get_hierarchical_clusters(cormat, between_cluster = between_cluster)
    cluster_C <- res$cluster
    n_clusters_C <- ncol(cluster_C)

    integrated_twas_weight_C <- list()
    tau_ranges_C <- list()
    pseudo_R2_avg_C <- list()

    for (k in 1:n_clusters_C) {
      tau_indices_C <- which(cluster_C[, k] == 1)
      beta_values_C <- matrix(twas_weight[, tau_indices_C], nrow = nrow(twas_weight))
      rownames(beta_values_C) <- variant_id
      tau_subset_C <- tau_values[tau_indices_C]

      integrated_twas_weight_C[[k]] <- integrate_tau(tau_subset_C, beta_values_C)
      tau_ranges_C[[k]] <- c(min(tau_subset_C), max(tau_subset_C))
      pseudo_R2_avg_C[[k]] <- mean(pseudo_R2[tau_indices_C])
    }
  } else {
    # If only one variant_id, skip Method 1
    integrated_twas_weight_C <- list()
    tau_ranges_C <- list()
    pseudo_R2_avg_C <- list()
    n_clusters_C <- 0 # No clusters from Method 1
  }

  # Method 2: Fixed tau intervals (A1, A2, A3) with non-overlapping divisions
  tau_cuts <- cut(tau_values, breaks = num_intervals, labels = FALSE, include.lowest = TRUE)
  fixed_intervals <- tapply(tau_values, tau_cuts, range)

  integrated_twas_weight_A <- list()
  tau_ranges_A <- list()
  pseudo_R2_avg_A <- list()

  for (k in 1:num_intervals) {
    tau_indices_A <- which(tau_cuts == k)
    beta_values_A <- matrix(twas_weight[, tau_indices_A], nrow = nrow(twas_weight))
    rownames(beta_values_A) <- variant_id
    tau_subset_A <- tau_values[tau_indices_A]

    integrated_twas_weight_A[[k]] <- integrate_tau(tau_subset_A, beta_values_A)
    tau_ranges_A[[k]] <- c(fixed_intervals[[k]][1], fixed_intervals[[k]][2])
    pseudo_R2_avg_A[[k]] <- mean(pseudo_R2[tau_indices_A])
  }

  # Combine the results of Method 2 and Method 1 if Method 1 was applied
  combined_weights <- do.call(cbind, c(integrated_twas_weight_A, integrated_twas_weight_C))

  # Ensure the correct number of column names is provided
  total_columns <- ncol(combined_weights)
  expected_columns <- num_intervals + n_clusters_C

  if (total_columns == expected_columns) {
    colnames(combined_weights) <- c(
      paste0("A", 1:num_intervals, "_weights"),
      if (n_clusters_C > 0) paste0("C", 1:n_clusters_C, "_weights") else NULL
    )
  } else {
    stop("Mismatch between the number of columns and the number of intervals/clusters")
  }

  rownames(combined_weights) <- variant_id

  twas_weight_performance <- list()

  # Store results for Method 2 (A1, A2, A3) as data frames
  for (k in 1:num_intervals) {
    twas_weight_performance[[paste0("A", k, "_performance")]] <- data.frame(
      "quantile_start" = round(fixed_intervals[[k]][1], 2),
      "quantile_end" = round(fixed_intervals[[k]][2], 2),
      "pseudo_R2_avg" = pseudo_R2_avg_A[[k]]
    )
  }

  # Store results for Method 1 (C1, C2, ..., Cn) as data frames if it was applied
  if (n_clusters_C > 0) {
    for (k in 1:n_clusters_C) {
      twas_weight_performance[[paste0("C", k, "_performance")]] <- data.frame(
        "quantile_start" = round(tau_ranges_C[[k]][1], 2),
        "quantile_end" = round(tau_ranges_C[[k]][2], 2),
        "pseudo_R2_avg" = pseudo_R2_avg_C[[k]]
      )
    }
  }

  return(list("weights" = combined_weights, "twas_weight_performance" = twas_weight_performance))
}


#' Integrate Tau Across Quantiles
#'
#' This function integrates TWAS weights across a set of tau values (quantiles) using trapezoidal integration.
#'
#' @param tau.temp A numeric vector of tau values.
#' @param a_tau A numeric vector or matrix of TWAS weights corresponding to the tau values.
#'
#' @return A numeric vector of integrated TWAS weights.
#'
#' @examples
#' tau_values <- seq(0.01, 0.99, length.out = 10)
#' a_tau <- matrix(runif(100), nrow = 10)
#' integrated_weights <- integrate_tau(tau_values, a_tau)
#'
#' @noRd
integrate_tau <- function(tau.temp, a_tau) {
  tau.temp1 <- c(min(tau.temp) - 0.01, tau.temp, max(tau.temp) +
    0.01)
  if (is.matrix(a_tau)) {
    a_tau <- cbind(rep(0, nrow(a_tau)), a_tau, rep(0, nrow(a_tau)))
    out <- (0.5 * (a_tau[, -ncol(a_tau)] + a_tau[, -1]) %*%
      diff(tau.temp1))
  } else {
    a_tau <- c(0, a_tau, 0)
    out <- (0.5 * (a_tau[-length(a_tau)] + a_tau[-1]) %*%
      diff(tau.temp1))
  }
  return(out)
}

#' Load Quantile TWAS Weights
#'
#' This function loads TWAS weights from RDS files and performs grouped integration using both correlation-based clustering and fixed tau intervals.
#'
#' @param weight_db_files A character vector of file paths to the RDS files containing TWAS weights.
#' @param tau_values A numeric vector representing the tau values (quantiles) to use. Default is \code{seq(0.01, 0.99, 0.01)}.
#' @param between_cluster A numeric value specifying the correlation threshold for clustering. Default is 0.8.
#' @param num_intervals The number of fixed non-overlapping intervals to divide the tau values. Default is 3.
#'
#' @return A list containing:
#' \item{weights}{A list of integrated TWAS weights for each context.}
#' \item{twas_cv_performance}{A list of TWAS cross-validation performance metrics for each context.}
#'
#' @examples
#' weight_db_files <- c("file1.rds", "file2.rds")
#' tau_values <- seq(0.01, 0.99, by = 0.01)
#' result <- load_quantile_twas_weights(weight_db_files, tau_values)
#' weights <- result$weights
#' performance <- result$twas_cv_performance
#'
#' @export
load_quantile_twas_weights <- function(weight_db_files, tau_values = seq(0.01, 0.99, 0.01),
                                       between_cluster = 0.8, num_intervals = 3) {
  # Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files) {
    all_data <- lapply(weight_db_files, function(rds_file) {
      db <- readRDS(rds_file)

      # Extract gene data and check for valid entries
      lapply(names(db), function(gene) {
        gene_data <- db[[gene]]
        lapply(names(gene_data), function(context) {
          context_data <- gene_data[[context]]
          if (!is.null(context_data$twas_weight) && !is.null(context_data$pseudo_R2)) {
            list(
              gene = gene,
              context = clean_context_names(context, gene),
              data = context_data
            )
          } else {
            message(paste("Warning: Missing twas_weight or pseudo_R2 for", gene, "-", context))
            NULL
          }
        })
      })
    })

    result <- do.call(c, unlist(all_data, recursive = FALSE))
    result <- Filter(Negate(is.null), result) # Remove NULL elements
    return(result)
  }

  # Load and validate the data
  combined_data <- load_and_validate_data(weight_db_files)
  if (length(combined_data) == 0) {
    stop("No valid data found in the provided files.")
  }

  # Merge contexts for the same gene from different RDS files
  merged_data <- list()
  for (item in combined_data) {
    gene <- item$gene
    context <- item$context
    if (is.null(merged_data[[gene]])) {
      merged_data[[gene]] <- list()
    }
    # Ensure we are not duplicating gene entries unnecessarily
    if (!is.null(item$data)) {
      merged_data[[gene]][[context]] <- item$data
    }
  }

  # Perform grouped integration for each context and restructure results
  results <- list(weights = list(), twas_cv_performance = list())

  for (gene in names(merged_data)) {
    message(paste("Processing gene:", gene))

    for (context in names(merged_data[[gene]])) {
      context_data <- merged_data[[gene]][[context]]

      if (!is.null(context_data$twas_weight) && !is.null(context_data$pseudo_R2)) {
        tryCatch(
          {
            result <- perform_grouped_integration(
              context_data$twas_weight, tau_values,
              context_data$pseudo_R2, between_cluster, num_intervals
            )
            results$weights[[context]] <- result$weights
            results$twas_cv_performance[[context]] <- result$twas_weight_performance
          },
          error = function(e) {
            message(paste("Error in perform_grouped_integration for", gene, "-", context, ":", e$message))
          }
        )
      }
    }
  }

  # Output the final structure
  if (length(results$weights) == 0) {
    message("Warning: No valid results were generated.")
  } else {
    message(paste("Number of contexts with valid results:", length(results$weights)))
  }

  return(results)
}
