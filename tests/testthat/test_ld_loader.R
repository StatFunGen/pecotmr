context("ld_loader")

# ===========================================================================
# ld_loader: input validation
# ===========================================================================

test_that("ld_loader errors when no source is provided", {
  expect_error(ld_loader(), "Provide exactly one")
})

test_that("ld_loader errors when multiple sources are provided", {
  R <- list(matrix(1, 2, 2))
  X <- list(matrix(1, 3, 2))
  expect_error(ld_loader(R_list = R, X_list = X), "Provide exactly one")
})

# ===========================================================================
# ld_loader: R_list branch
# ===========================================================================

test_that("ld_loader with R_list returns a function", {
  R <- list(matrix(c(1, 0.5, 0.5, 1), 2, 2))
  loader <- ld_loader(R_list = R)
  expect_type(loader, "closure")
})

test_that("ld_loader R_list returns correct matrix", {
  R1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  R2 <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  loader <- ld_loader(R_list = list(R1, R2))
  expect_equal(loader(1), R1)
  expect_equal(loader(2), R2)
})

test_that("ld_loader R_list with max_variants downsamples", {
  set.seed(42)
  R <- matrix(0.1, 10, 10)
  diag(R) <- 1
  loader <- ld_loader(R_list = list(R), max_variants = 5)
  result <- loader(1)
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
})

test_that("ld_loader R_list without max_variants returns full matrix", {
  R <- matrix(0.1, 10, 10)
  diag(R) <- 1
  loader <- ld_loader(R_list = list(R))
  result <- loader(1)
  expect_equal(nrow(result), 10)
})

test_that("ld_loader R_list max_variants larger than matrix returns full matrix", {
  R <- matrix(0.1, 3, 3)
  diag(R) <- 1
  loader <- ld_loader(R_list = list(R), max_variants = 100)
  result <- loader(1)
  expect_equal(nrow(result), 3)
})

# ===========================================================================
# ld_loader: X_list branch
# ===========================================================================

test_that("ld_loader with X_list returns a function", {
  X <- list(matrix(rnorm(30), 10, 3))
  loader <- ld_loader(X_list = X)
  expect_type(loader, "closure")
})

test_that("ld_loader X_list returns correct matrix", {
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(1:8, 4, 2)
  loader <- ld_loader(X_list = list(X1, X2))
  expect_equal(loader(1), X1)
  expect_equal(loader(2), X2)
})

test_that("ld_loader X_list with max_variants downsamples columns", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  loader <- ld_loader(X_list = list(X), max_variants = 3)
  result <- loader(1)
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
})

test_that("ld_loader X_list max_variants larger than ncol returns full matrix", {
  X <- matrix(rnorm(12), 4, 3)
  loader <- ld_loader(X_list = list(X), max_variants = 100)
  result <- loader(1)
  expect_equal(ncol(result), 3)
})

# ===========================================================================
# ld_loader: ld_meta_path branch validation
# ===========================================================================

test_that("ld_loader with ld_meta_path but no regions errors", {
  expect_error(
    ld_loader(ld_meta_path = "/some/path"),
    "regions.*required"
  )
})

# ===========================================================================
# ld_loader: LD_info branch validation
# ===========================================================================

test_that("ld_loader with LD_info errors when not a data.frame", {
  expect_error(
    ld_loader(LD_info = "not_a_df"),
    "LD_info must be a data.frame"
  )
})

test_that("ld_loader with LD_info errors when missing LD_file column", {
  expect_error(
    ld_loader(LD_info = data.frame(col1 = "a")),
    "LD_info must be a data.frame with column 'LD_file'"
  )
})
