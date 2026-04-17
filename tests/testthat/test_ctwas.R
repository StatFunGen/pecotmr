context("ctwas")

# ===========================================================================
#  ctwas wrapper tests
# ===========================================================================


# ---------- trim_ctwas_variants --------------------------------------------

# Helper: build a minimal region_data structure that trim_ctwas_variants expects
make_mock_region_data <- function() {
  # Variant IDs in canonical format (chr:pos:A2:A1)
  variant_ids <- c("chr1:1000:A:G", "chr1:2000:C:T", "chr1:3000:G:A", "chr1:4000:T:C")

  # Weight matrix (4 variants x 1 weight column)
  wgt <- matrix(c(0.5, 0.0001, 0.3, -0.2), nrow = 4, ncol = 1)
  rownames(wgt) <- variant_ids

  gene_id <- "GENE1|ctx1"
  context <- "ctx1"
  study <- "study1"

  weights <- list()
  weights[[gene_id]] <- list()
  weights[[gene_id]][[study]] <- list(
    wgt = wgt,
    context = context,
    p0 = 1000,
    p1 = 4000
  )

  # SuSiE intermediate info
  pip_vals <- c(0.8, 0.05, 0.6, 0.02)
  names(pip_vals) <- variant_ids

  susie_weights_intermediate <- list()
  susie_weights_intermediate[["GENE1"]] <- list()
  susie_weights_intermediate[["GENE1"]][[context]] <- list(
    pip = pip_vals,
    cs_variants = list(variant_ids[c(1, 3)]),
    cs_purity = list(min.abs.corr = 0.9)
  )

  list(
    weights = weights,
    susie_weights_intermediate = susie_weights_intermediate
  )
}

test_that("trim_ctwas_variants removes variants below weight cutoff", {
  rd <- make_mock_region_data()
  # Default cutoff 1e-5, variant 2 has weight 0.0001 (above), so all 4 should pass default cutoff
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5)
  expect_true(is.list(result))
  # With a higher cutoff, the near-zero variant should be removed
  result_strict <- trim_ctwas_variants(rd, twas_weight_cutoff = 0.001)
  # study1 should exist in result
  expect_true("study1" %in% names(result_strict))
  # Get the gene-level result
  gene_weights <- result_strict[["study1"]][["GENE1|ctx1"]]
  # Variant 2 has abs(weight) = 0.0001 < 0.001, so should be removed
  expect_false("chr1:2000:C:T" %in% rownames(gene_weights$wgt))
})

test_that("trim_ctwas_variants removes gene when all weights below cutoff", {
  rd <- make_mock_region_data()
  # Set cutoff so high that all variants are dropped
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 10)
  # Gene should be removed entirely since no weights pass the cutoff
  # Result should be an empty list
  expect_equal(length(result), 0)
})

test_that("trim_ctwas_variants returns result keyed by study", {
  rd <- make_mock_region_data()
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5)
  # merge_by_study reorganizes: weights[[study]][[group]]
  expect_true("study1" %in% names(result))
  expect_true("GENE1|ctx1" %in% names(result[["study1"]]))
})

test_that("trim_ctwas_variants updates p0 and p1 positions", {
  rd <- make_mock_region_data()
  # Use a weight cutoff that removes the variant at position 2000
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 0.001)
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  # p0 and p1 should reflect the range of remaining variant positions
  remaining_positions <- as.integer(sapply(
    rownames(gene_weights$wgt),
    function(v) strsplit(v, ":")[[1]][2]
  ))
  expect_equal(gene_weights$p0, min(remaining_positions))
  expect_equal(gene_weights$p1, max(remaining_positions))
})

test_that("trim_ctwas_variants respects max_num_variants", {
  rd <- make_mock_region_data()
  # Request max 2 variants; since nrow(wgt) == 4 >= max_num_variants == 2,
  # it triggers select_variants which picks by PIP priority
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5, max_num_variants = 2)
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  expect_true(nrow(gene_weights$wgt) <= 2)
})

test_that("trim_ctwas_variants handles NA weights by removing group", {
  rd <- make_mock_region_data()
  # Replace all weights with NA
  rd$weights[["GENE1|ctx1"]][["study1"]]$wgt[, 1] <- NA
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 0)
  # The group should be removed because all weights are NA
  expect_equal(length(result), 0)
})

test_that("trim_ctwas_variants handles multiple genes", {
  rd <- make_mock_region_data()

  # Add a second gene
  variant_ids2 <- c("chr1:5000:A:G", "chr1:6000:C:T")
  wgt2 <- matrix(c(0.4, -0.3), nrow = 2, ncol = 1)
  rownames(wgt2) <- variant_ids2

  rd$weights[["GENE2|ctx1"]] <- list()
  rd$weights[["GENE2|ctx1"]][["study1"]] <- list(
    wgt = wgt2,
    context = "ctx1",
    p0 = 5000,
    p1 = 6000
  )

  pip_vals2 <- c(0.7, 0.4)
  names(pip_vals2) <- variant_ids2
  rd$susie_weights_intermediate[["GENE2"]] <- list()
  rd$susie_weights_intermediate[["GENE2"]][["ctx1"]] <- list(
    pip = pip_vals2,
    cs_variants = list(variant_ids2[1]),
    cs_purity = list(min.abs.corr = 0.95)
  )

  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5)
  expect_true("GENE1|ctx1" %in% names(result[["study1"]]))
  expect_true("GENE2|ctx1" %in% names(result[["study1"]]))
})

test_that("trim_ctwas_variants select_variants uses cs_min_cor to include CS variants", {
  rd <- make_mock_region_data()
  # cs_purity min.abs.corr = 0.9, so with cs_min_cor = 0.8 the CS variants
  # (variant 1 and 3) should be included. Max 2 variants.
  result <- trim_ctwas_variants(rd,
    twas_weight_cutoff = 1e-5,
    cs_min_cor = 0.8,
    min_pip_cutoff = 0.0,
    max_num_variants = 2
  )
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  included <- rownames(gene_weights$wgt)
  # CS variants chr1:1000:A:G and chr1:3000:G:A have highest PIPs (0.8 and 0.6)
  # and are in the CS, so they should be prioritized
  expect_true("chr1:1000:A:G" %in% included)
  expect_true("chr1:3000:G:A" %in% included)
})
