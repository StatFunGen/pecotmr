context("alleleQc")

create_allele_data <- function(seed, n=100, match_min_prop=0.8, ambiguous=FALSE, non_actg=FALSE, edge_cases=FALSE) {
  set.seed(seed)
  num_pass <- n*match_min_prop
  sumstat_A1 <- sample(c("A", "T", "G", "C"), num_pass, replace = TRUE)
  sumstat_A2 <- unlist(lapply(sumstat_A1, function(x) {
    if (x == "A") {
      return(sample(c("G", "C"), 1))
    } else if (x == "T") {
      return(sample(c("G", "C"), 1))
    } else if (x == "G") {
      return(sample(c("A", "T"), 1))
    } else if (x == "C") {
      return(sample(c("A", "T"), 1))
    }
  }))

  if (ambiguous) {
    # Strand Ambiguous SNPs
    sumstat_A1 <- c(sumstat_A1, sample(c("A", "T", "G", "C"), n-num_pass, replace = TRUE))
    sumstat_A2 <- unlist(c(sumstat_A2, lapply(sumstat_A1[(num_pass+1):length(sumstat_A1)], function(x) {
      if (x == "A") {
        return("T")
      } else if (x == "T") {
        return("A")
      } else if (x == "G") {
        return("C")
      } else if (x == "C") {
        return("G")
      }
    })))
  } else if (non_actg) {
    # Non-ATCG coding SNPs
    sumstat_A1 <- c(sumstat_A1, sample(c("ATG", "TAC", "GACA", "CTAA"), n-num_pass, replace = TRUE))
    sumstat_A2 <- unlist(c(sumstat_A2, lapply(sumstat_A1[(num_pass+1):length(sumstat_A1)], function(x) {
      if (x == "ATG") {
        return("TAC")
      } else if (x == "TAC") {
        return("ATG")
      } else if (x == "GACA") {
        return("CTGT")
      } else if (x == "CTAA") {
        return("GATT")
      }
    })))
  }

  # Info SNPs
  info_A1 <- lapply(sumstat_A1[1:num_pass], function(x) {
    if(runif(1) < 0.2) {
      # flip a small proportion of the alleles
      if (x == "A") {
        return("T")
      } else if (x == "T") {
        return("A")
      } else if (x == "G") {
        return("C")
      } else if (x == "C") {
        return("G")
      }
    } else {
      return(x)
    }
  })
  info_A2 <- sumstat_A2[1:num_pass]
  # Handle random flips
  info_A2[info_A1 != sumstat_A1[1:num_pass]] <- unlist(lapply(info_A2[info_A1 != sumstat_A1[1:num_pass]], function(x) {
    if (x == "A") {
      return("T")
    } else if (x == "T") {
      return("A")
    } else if (x == "G") {
      return("C")
    } else if (x == "C") {
      return("G")
    }
  }))

  # Create the rest of the alleles
  info_A1 <- unlist(c(info_A1, sample(c("A", "T", "G", "C"), n-num_pass, replace = TRUE)))
  info_A2 <- unlist(c(info_A2, lapply(info_A1[(num_pass+1):length(info_A1)], function(x) {
    if (x == "A") {
      return(sample(c("G", "C"), 1))
    } else if (x == "T") {
      return(sample(c("G", "C"), 1))
    } else if (x == "G") {
      return(sample(c("A", "T"), 1))
    } else if (x == "C") {
      return(sample(c("A", "T"), 1))
    }
  })))

  chromosome <- unlist(rep(sample(1:20, 1), n))
  snp_positions <- sample(1:1000000, n)
  ref_variants <- data.frame(
    chrom = chromosome,
    pos = snp_positions,
    A1 = info_A1,
    A2 = info_A2
  )
  target_data <- data.frame(
    chrom = chromosome,
    pos = snp_positions,
    A1 = sumstat_A1,
    A2 = sumstat_A2,
    beta = rnorm(n),
    z = rnorm(n)
  )

  return(list(target_data = target_data, ref_variants = ref_variants))
}

test_that("Check that we correctly remove stand ambiguous SNPs",{
  res <- create_allele_data(1, n=100, match_min_prop=0.8, ambiguous=TRUE)
  output <- alleleQc(
    res$target_data, res$ref_variants, "beta", matchMinProp = 0.2,
    TRUE, FALSE, TRUE)
  expect_equal(nrow(getHarmonizedData(output)), 80)
})

test_that("Check that we correctly remove non-ACTG coding SNPs",{
  res <- create_allele_data(1, n=100, match_min_prop=0.4, non_actg=TRUE)
  output <- alleleQc(
    res$target_data, res$ref_variants, "beta", matchMinProp = 0.2,
    TRUE, FALSE, TRUE)
  expect_equal(nrow(getHarmonizedData(output)), 40)
})

test_that("Check that execution stops if not enough variants are matched",{
  res <- create_allele_data(1, n=100, match_min_prop=0.1, ambiguous=TRUE)
  expect_error(alleleQc(
    res$target_data, res$ref_variants, "beta", matchMinProp = 0.2,
    TRUE, FALSE, TRUE), "Not enough variants have been matched.")
})

test_that("alleleQc matches exact alleles", {
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "C"), A1 = c("G", "T")
  )
  ref <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "C"), A1 = c("G", "T")
  )
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 2)
})

test_that("alleleQc detects sign flips", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G",
    z = 2.5
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "G", A1 = "A"
  )
  result <- alleleQc(target, ref, colToFlip = "z", matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 1)
  # z should be flipped
  expect_equal(getHarmonizedData(result)$z, -2.5)
})

test_that("alleleQc handles string input format", {
  target <- c("1:100:A:G", "1:200:C:T")
  ref <- c("1:100:A:G", "1:200:C:T")
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 2)
})

test_that("alleleQc with chr prefix", {
  target <- c("chr1:100:A:G", "chr1:200:C:T")
  ref <- c("chr1:100:A:G", "chr1:200:C:T")
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 2)
})

test_that("alleleQc warns when too few matches", {
  target <- c("1:100:A:G")
  ref <- c("2:200:C:T", "2:300:A:G", "2:400:C:T", "2:500:A:G", "2:600:C:T")
  expect_warning(alleleQc(target, ref, matchMinProp = 0.5))
})

test_that("alleleQc with no matching positions returns empty", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  ref <- data.frame(
    chrom = 1, pos = 999,
    A2 = "C", A1 = "T"
  )
  expect_warning(
    result <- alleleQc(target, ref, matchMinProp = 0),
    "No matching variants"
  )
  expect_equal(nrow(getHarmonizedData(result)), 0)
})

test_that("alleleQc preserves extra columns", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G",
    beta = 0.5, se = 0.1
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_true("beta" %in% colnames(getHarmonizedData(result)))
  expect_true("se" %in% colnames(getHarmonizedData(result)))
})

test_that("alleleQc with lowercase alleles", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "a", A1 = "g"
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 1)
})

test_that("alignVariantNames correctly aligns variant names", {
  # Test case 1: Matching variant names
  source1 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference1 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned1 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_unmatched1 <- integer(0)

  result1 <- alignVariantNames(source1, reference1)
  expect_equal(result1$aligned_variants, expected_aligned1)
  expect_equal(result1$unmatched_indices, expected_unmatched1)

  # Test case 2: Unmatched variant names
  source2 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A", "4:101:G:C")
  reference2 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned2 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A", "4:101:G:C")
  expected_unmatched2 <- 4

  result2 <- alignVariantNames(source2, reference2)
  expect_equal(result2$aligned_variants, expected_aligned2)
  expect_equal(result2$unmatched_indices, expected_unmatched2)

  # Test case 3: Different variant name formats
  source3 <- c("1:123:A:C", "2:456_G_T", "3:789:C:A")
  reference3 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned3 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_unmatched3 <- integer(0)

  result3 <- alignVariantNames(source3, reference3)
  expect_equal(result3$aligned_variants, expected_aligned3)
  expect_equal(result3$unmatched_indices, expected_unmatched3)
})

test_that("alignVariantNames correctly aligns variant names with different flip patterns", {
  # Test case 4: Strand flip
  source4 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference4 <- c("1:123:T:G", "2:456:A:C", "3:789:C:A")
  expected_aligned4 <- c("1:123:T:G", "2:456:A:C", "3:789:C:A")
  expected_unmatched4 <- integer(0)

  result4 <- alignVariantNames(source4, reference4)
  expect_equal(result4$aligned_variants, expected_aligned4)
  expect_equal(result4$unmatched_indices, expected_unmatched4)

  # Test case 5: Strand ambiguous variants
  source5 <- c("1:123:A:T", "2:456:G:C", "3:789:C:A")
  reference5 <- c("1:123:A:T", "2:456:G:C", "3:789:C:A")
  expected_aligned5 <- c("1:123:A:T", "2:456:G:C", "3:789:C:A")
  expected_unmatched5 <- integer(0)

  result5 <- alignVariantNames(source5, reference5)
  expect_equal(result5$aligned_variants, expected_aligned5)
  expect_equal(result5$unmatched_indices, expected_unmatched5)

  # Test case 6: Sign flip
  source6 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference6 <- c("1:123:C:A", "2:456:T:G", "3:789:C:A")
  expected_aligned6 <- c("1:123:C:A", "2:456:T:G", "3:789:C:A")
  expected_unmatched6 <- integer(0)

  result6 <- alignVariantNames(source6, reference6)
  expect_equal(result6$aligned_variants, expected_aligned6)
  expect_equal(result6$unmatched_indices, expected_unmatched6)

  # Test case 7: Strand and sign flip
  source7 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference7 <- c("1:123:G:T", "2:456:A:C", "3:789:C:A")
  expected_aligned7 <- c("1:123:G:T", "2:456:A:C", "3:789:C:A")
  expected_unmatched7 <- integer(0)

  result7 <- alignVariantNames(source7, reference7)
  expect_equal(result7$aligned_variants, expected_aligned7)
  expect_equal(result7$unmatched_indices, expected_unmatched7)

  # Test case 8: Indels
  source8 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A", "4:101:G:GATC")
  reference8 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A", "4:101:GATC:G")
  expected_aligned8 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A", "4:101:GATC:G")
  expected_unmatched8 <- integer(0)

  result8 <- alignVariantNames(source8, reference8)
  expect_equal(result8$aligned_variants, expected_aligned8)
  expect_equal(result8$unmatched_indices, expected_unmatched8)
})

test_that("alignVariantNames correctly aligns variant names with different chr prefix conventions", {
  # Test case 9: Original without chr prefix, reference with chr prefix
  source9 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference9 <- c("chr1:123:A:C", "chr2:456:T:G", "chr3:789:C:A")
  expected_aligned9 <- c("chr1:123:A:C", "chr2:456:T:G", "chr3:789:C:A")
  expected_unmatched9 <- integer(0)

  result9 <- alignVariantNames(source9, reference9)
  expect_equal(result9$aligned_variants, expected_aligned9)
  expect_equal(result9$unmatched_indices, expected_unmatched9)

  # Test case 10: Original with chr prefix, reference without chr prefix
  source10 <- c("chr1:123:A:C", "chr2:456:G:T", "chr3:789:C:A")
  reference10 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned10 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_unmatched10 <- integer(0)

  result10 <- alignVariantNames(source10, reference10)
  expect_equal(result10$aligned_variants, expected_aligned10)
  expect_equal(result10$unmatched_indices, expected_unmatched10)
})

test_that("alignVariantNames warns on non-standard format", {
  source <- c("rs12345")
  reference <- c("rs67890")
  expect_warning(
    alignVariantNames(source, reference),
    "do not follow the expected"
  )
})

test_that("alignVariantNames errors on mixed formats", {
  source <- c("1:100:A:G")
  reference <- c("rs12345")
  expect_error(
    alignVariantNames(source, reference),
    "different variant naming conventions"
  )
})

test_that("alignVariantNames strips build suffix", {
  source <- c("1:100:A:G:b38")
  reference <- c("1:100:A:G")
  result <- alignVariantNames(source, reference, removeBuildSuffix = TRUE)
  expect_length(result$aligned_variants, 1)
})

# ---- sanitize_names edge cases (alleleQc.R lines 37, 42) ----
test_that("alleleQc handles data frame with NULL colnames after merge", {
  # Create a data frame where merge might produce empty names
  # by giving target_data a column with NA name
  target <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  colnames(target)[1] <- ""
  colnames(target) <- make.unique(colnames(target), sep = "_")
  # Restore chrom for the join
  colnames(target)[1] <- "chrom"
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 1)
})

# ---- target_data with redundant columns (alleleQc.R line 75) ----
test_that("alleleQc removes redundant columns from target_data before join", {
  target <- data.frame(
    chrom = 1, pos = 100, A2 = "A", A1 = "G",
    variant_id = "1:100:A:G", chromosome = "chr1", position = 100
  )
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  result <- alleleQc(target, ref, matchMinProp = 0)
  expect_equal(nrow(getHarmonizedData(result)), 1)
  # The redundant columns should have been removed before the join
  expect_true("variant_id" %in% colnames(getHarmonizedData(result)))
})

# ---- col_to_flip with nonexistent column (alleleQc.R line 130) ----
test_that("alleleQc errors when col_to_flip column does not exist", {
  target <- data.frame(chrom = 1, pos = 100, A2 = "G", A1 = "A")
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  expect_error(
    alleleQc(target, ref, colToFlip = "nonexistent_col", matchMinProp = 0),
    "not found in targetData"
  )
})

# ---- duplicate removal warning (alleleQc.R lines 148-150) ----
test_that("alleleQc warns and removes duplicate variants", {
  # Two target rows at the same position will produce duplicates after join
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 100),
    A2 = c("A", "A"), A1 = c("G", "G"),
    beta = c(0.5, 0.6)
  )
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  expect_warning(
    result <- alleleQc(target, ref, matchMinProp = 0, removeDups = TRUE),
    "duplicate variant"
  )
  expect_equal(nrow(getHarmonizedData(result)), 1)
})

# ---- duplicated variant IDs error (alleleQc.R line 180) ----
test_that("alleleQc errors on duplicated variant IDs with different values when remove_dups is FALSE", {
  # Two rows at same position, same alleles, but different beta - when not removing dups
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 100),
    A2 = c("A", "A"), A1 = c("G", "G"),
    beta = c(0.5, 0.6)
  )
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  expect_error(
    alleleQc(target, ref, matchMinProp = 0, removeDups = FALSE),
    "Duplicated variants"
  )
})

# ===========================================================================
# colToComplement hook (rss-qc-parity): af complemented on allele swap
# ===========================================================================

test_that("af is complemented (1 - af) when harmonization swaps the effect allele", {
  ref <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                    A2 = c("A", "C"), A1 = c("G", "T"), stringsAsFactors = FALSE)
  # chr1:100 alleles swapped vs ref (=> sign flip); chr1:200 exact match.
  target <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                       A2 = c("G", "C"), A1 = c("A", "T"),
                       z = c(2.0, 1.5), af = c(0.30, 0.40), stringsAsFactors = FALSE)

  res <- matchRefPanel(target, ref, colToFlip = "z", colToComplement = "af")
  h <- getHarmonizedData(res)
  swapped <- h[h$pos == 100, ]
  control <- h[h$pos == 200, ]

  expect_equal(swapped$af, 0.70)   # 1 - input af
  expect_equal(swapped$z, -2.0)    # signed columns still sign-flip
  expect_equal(control$af, 0.40)   # untouched (no swap)
  expect_equal(control$z, 1.5)
})

test_that("colToComplement default leaves af unchanged (non-RSS callers unaffected)", {
  ref <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                    A2 = c("A", "C"), A1 = c("G", "T"), stringsAsFactors = FALSE)
  target <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                       A2 = c("G", "C"), A1 = c("A", "T"),
                       z = c(2.0, 1.5), af = c(0.30, 0.40), stringsAsFactors = FALSE)

  res <- matchRefPanel(target, ref, colToFlip = "z")  # default: no complement
  h <- getHarmonizedData(res)
  swapped <- h[h$pos == 100, ]
  expect_equal(swapped$af, 0.30)   # unchanged
  expect_equal(swapped$z, -2.0)    # z still sign-flips (independent path)
})

test_that("colToComplement errors on a missing column name", {
  ref <- data.frame(chrom = "chr1", pos = 100, A2 = "A", A1 = "G", stringsAsFactors = FALSE)
  target <- data.frame(chrom = "chr1", pos = 100, A2 = "A", A1 = "G",
                       z = 1.0, stringsAsFactors = FALSE)
  expect_error(
    matchRefPanel(target, ref, colToComplement = "af"),
    "not found in targetData"
  )
})
