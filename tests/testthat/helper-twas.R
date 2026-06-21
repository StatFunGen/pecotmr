# Shared test helper: simulate a small (X, Y) pair for twas weights /
# causal-inference tests. Migrated from the old test_twas.R preamble.
generate_X_Y <- function(seed=1, num_samples=10, num_features=10, X_rownames=TRUE, y_rownames=TRUE) {
  set.seed(seed)
  X <- scale(
    matrix(rnorm(num_samples * num_features), nrow = num_samples),
    center = TRUE, scale = TRUE)

  if (X_rownames) {
    rownames(X) <- paste0("sample", 1:num_samples)
  } else {
    rownames(X) <- NULL
  }

  beta = rep(0, num_features)
  beta[1:4] = 1
  y <- X %*% beta + rnorm(num_samples)
  y <- matrix(y, nrow = num_samples, ncol = 1)
  if (y_rownames) {
    rownames(y) <- paste0("sample", 1:num_samples)
  } else {
    rownames(y) <- NULL
  }
  colnames(y) <- c("Outcome")

  return(list(X=X, Y=y))
}
