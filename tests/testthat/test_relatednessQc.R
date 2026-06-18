context("filterRelatedness")

# Helper: check no remaining pairs in the kept set exceed threshold
no_related_pairs_remain <- function(relatedness, excluded, threshold,
                                    iid1 = "IID1", iid2 = "IID2",
                                    value = "PI_HAT") {
  kept <- relatedness[
    !(relatedness[[iid1]] %in% excluded) &
    !(relatedness[[iid2]] %in% excluded), ]
  all(kept[[value]] < threshold)
}

test_that("maximize_unrelated removes related individuals and leaves clean set", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("plinkQC")

  # 10 individuals, several pairs above threshold 0.125
  rel <- data.frame(
    IID1 = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "A"),
    IID2 = c("B", "C", "D", "E", "F", "G", "H", "I", "J", "J"),
    PI_HAT = c(0.25, 0.15, 0.30, 0.08, 0.20, 0.05, 0.18, 0.03, 0.22, 0.14),
    stringsAsFactors = FALSE
  )
  threshold <- 0.125

  result <- filterRelatedness(
    relatedness = rel,
    relatednessThreshold = threshold,
    analysisType = "maximize_unrelated"
  )

  expect_type(result, "character")
  expect_true(length(result) > 0)
  expect_true(no_related_pairs_remain(rel, result, threshold))
})

test_that("no related pairs returns empty exclusion vector", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("plinkQC")

  rel <- data.frame(
    IID1 = c("A", "B", "C", "D"),
    IID2 = c("B", "C", "D", "E"),
    PI_HAT = c(0.01, 0.02, 0.03, 0.04),
    stringsAsFactors = FALSE
  )
  threshold <- 0.125

  result <- filterRelatedness(
    relatedness = rel,
    relatednessThreshold = threshold,
    analysisType = "maximize_unrelated"
  )

  expect_type(result, "character")
  expect_equal(length(result), 0)
})

test_that("large component pre-pruning removes individuals", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("plinkQC")

  # Build a chain of 30 individuals: 1-2, 2-3, ..., 29-30, all above threshold
  n <- 30
  ids <- paste0("IND", seq_len(n))
  rel <- data.frame(
    IID1 = ids[1:(n - 1)],
    IID2 = ids[2:n],
    PI_HAT = rep(0.20, n - 1),
    stringsAsFactors = FALSE
  )
  threshold <- 0.125

  result <- filterRelatedness(
    relatedness = rel,
    relatednessThreshold = threshold,
    analysisType = "maximize_unrelated",
    maxComponentSize = 10
  )

  expect_type(result, "character")
  expect_true(length(result) > 0)
  # The remaining kept individuals should have no related pairs
  expect_true(no_related_pairs_remain(rel, result, threshold))
})

test_that("maximize_cases preferentially retains cases", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("plinkQC")

  # Cases: C1, C2, C3; Controls: X1, X2, X3
  # Related pairs above threshold:
  #   C1-X1 (case-control), C2-X2 (case-control), C1-C2 (case-case),
  #   X1-X3 (control-control)
  rel <- data.frame(
    IID1 = c("C1", "C2", "C1", "X1", "C3", "X2"),
    IID2 = c("X1", "X2", "C2", "X3", "X3", "X3"),
    PI_HAT = c(0.25, 0.20, 0.15, 0.18, 0.05, 0.04),
    stringsAsFactors = FALSE
  )

  pheno <- data.frame(
    IID = c("C1", "C2", "C3", "X1", "X2", "X3"),
    pheno = c(1, 1, 1, 0, 0, 0),
    stringsAsFactors = FALSE
  )

  threshold <- 0.125
  result <- filterRelatedness(
    relatedness = rel,
    relatednessThreshold = threshold,
    analysisType = "maximize_cases",
    phenoData = pheno,
    phenoCol = "pheno"
  )

  expect_type(result, "character")
  # Controls related to kept cases should be excluded preferentially
  retained <- setdiff(pheno$IID, result)
  retained_cases <- intersect(retained, pheno$IID[pheno$pheno == 1])
  retained_controls <- intersect(retained, pheno$IID[pheno$pheno == 0])
  # We expect at least 2 of 3 cases kept, and controls sacrificed
  expect_true(length(retained_cases) >= 2)
  expect_true(no_related_pairs_remain(rel, result, threshold))
})

test_that("maximize_cases errors without pheno_data", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("plinkQC")

  rel <- data.frame(
    IID1 = c("A", "B"),
    IID2 = c("B", "C"),
    PI_HAT = c(0.25, 0.20),
    stringsAsFactors = FALSE
  )

  expect_error(
    filterRelatedness(
      relatedness = rel,
      relatednessThreshold = 0.125,
      analysisType = "maximize_cases"
    ),
    "Must provide phenoData"
  )
})
