context("deprecated")

# All tests below exercise functions in R/deprecated.R. The functions
# emit .Deprecated() warnings on every call; we wrap each call in
# suppressWarnings() (or rely on the surrounding expect_warning when
# the test already asserts a warning) so the deprecation message is
# not surfaced as test output.
#
# Every test starts with `skip_on_covr()`: R/deprecated.R is listed in
# .covrignore (its lines are not measured for coverage), so running
# these tests under `covr::package_coverage()` / `covr::codecov()`
# contributes nothing to the coverage signal — they would just exercise
# the deprecation shims whose forwarding targets are tested elsewhere.
# Skipping under covr cuts the coverage-run wallclock without changing
# what's measured.

# ===========================================================================
# Helpers (moved from test_sumstatsQc.R / test_colocPipeline.R /
# test_qtlEnrichmentPipeline.R)
# ===========================================================================

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

.ep_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("v", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

.ep_makeFmEntry <- function(variant_ids = paste0("chr1:", 100*(1:5), ":A:G"),
                             n_eff = 2L) {
  pip <- seq(0.9, by = -0.15, length.out = length(variant_ids))
  tl <- data.frame(variant_id = variant_ids, pip = pip,
                   stringsAsFactors = FALSE)
  set.seed(1)
  fit <- list(
    alpha = matrix(1/length(variant_ids),
                   nrow = n_eff, ncol = length(variant_ids),
                   dimnames = list(NULL, variant_ids)),
    pip   = setNames(pip, variant_ids),
    V     = rep(0.05, n_eff),
    lbf_variable = matrix(rnorm(n_eff * length(variant_ids)),
                          nrow = n_eff, ncol = length(variant_ids),
                          dimnames = list(NULL, variant_ids)))
  FineMappingEntry(variantIds = variant_ids,
                   trimmedFit = fit,
                   topLoci    = tl)
}

.ep_mockColocBfBf <- function() {
  function(qLbf, gLbf, p1, p2, p12, ...) {
    list(summary = data.frame(
      idx1 = 1L, idx2 = 1L, nSnps = ncol(qLbf),
      PP.H0.abf = 0.1, PP.H1.abf = 0.2, PP.H2.abf = 0.2,
      PP.H3.abf = 0.2, PP.H4.abf = 0.3,
      p12_actual = p12,
      stringsAsFactors = FALSE))
  }
}

.ep_makeQtlFmr <- function(with_sketch = TRUE) {
  QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(.ep_makeFmEntry()),
    ldSketch = if (with_sketch) .ep_makeHandle() else NULL)
}

.ep_makeGwasFmr <- function(with_sketch = TRUE) {
  GwasFineMappingResult(
    study  = "G1", method = "susie",
    entry  = list(.ep_makeFmEntry()),
    ldSketch = if (with_sketch) .ep_makeHandle() else NULL)
}

.ep_makeGwasSumstats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = rnorm(5), N = rep(1000L, 5))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .ep_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

generate_mock_data <- function(seed=1, num_pips = 1000, num_susie_fits = 2) {
  # Simulate fake data for gwas_pip
  n_gwas_pip <- num_pips
  gwas_pip <- runif(n_gwas_pip)
  names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)
  gwas_fit <- list(pip=gwas_pip)

  # Simulate fake data for a single SuSiEFit object
  simulate_susiefit <- function(n, p) {
    pip <- runif(n)
    names(pip) <- paste0("snp", 1:n)
    alpha <- t(matrix(runif(n * p), nrow = n))
    alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
    list(
      pip = pip,
      alpha = alpha,
      prior_variance = runif(p)
    )
  }

  # Simulate multiple SuSiEFit objects
  n_susie_fits <- num_susie_fits
  susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 10), simplify = FALSE)
  # Add these fits to a list, providing names to each element
  names(susie_fits) <- paste0("fit", 1:length(susie_fits))
  return(list(gwas_fit=gwas_fit, susie_fits=susie_fits))
}

# ===========================================================================
# alleleQc (deprecated; replaced by summaryStatsQc())
# ===========================================================================

test_that("Check that we correctly remove stand ambiguous SNPs",{
  skip_on_covr()
  res <- create_allele_data(1, n=100, match_min_prop=0.8, ambiguous=TRUE)
  output <- suppressWarnings(alleleQc(
    res$target_data, res$ref_variants, colToFlip = "beta",
    matchMinProp = 0.2))
  expect_equal(nrow(output$harmonizedData), 80)
})

test_that("Check that we correctly remove non-ACTG coding SNPs",{
  skip_on_covr()
  res <- create_allele_data(1, n=100, match_min_prop=0.4, non_actg=TRUE)
  output <- suppressWarnings(alleleQc(
    res$target_data, res$ref_variants, colToFlip = "beta",
    matchMinProp = 0.2))
  expect_equal(nrow(output$harmonizedData), 40)
})

test_that("Check that execution stops if not enough variants are matched",{
  skip_on_covr()
  res <- create_allele_data(1, n=100, match_min_prop=0.1, ambiguous=TRUE)
  expect_error(suppressWarnings(alleleQc(
    res$target_data, res$ref_variants, colToFlip = "beta",
    matchMinProp = 0.2)), "Not enough variants have been matched.")
})

test_that("alleleQc matches exact alleles", {
  skip_on_covr()
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "C"), A1 = c("G", "T")
  )
  ref <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "C"), A1 = c("G", "T")
  )
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 2)
})

test_that("alleleQc detects sign flips", {
  skip_on_covr()
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G",
    z = 2.5
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "G", A1 = "A"
  )
  result <- suppressWarnings(alleleQc(target, ref, colToFlip = "z", matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 1)
  # z should be flipped
  expect_equal(result$harmonizedData$z, -2.5)
})

test_that("alleleQc handles string input format", {
  skip_on_covr()
  target <- c("1:100:A:G", "1:200:C:T")
  ref <- c("1:100:A:G", "1:200:C:T")
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 2)
})

test_that("alleleQc with chr prefix", {
  skip_on_covr()
  target <- c("chr1:100:A:G", "chr1:200:C:T")
  ref <- c("chr1:100:A:G", "chr1:200:C:T")
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 2)
})

test_that("alleleQc warns when too few matches", {
  skip_on_covr()
  target <- c("1:100:A:G")
  ref <- c("2:200:C:T", "2:300:A:G", "2:400:C:T", "2:500:A:G", "2:600:C:T")
  expect_warning(alleleQc(target, ref, matchMinProp = 0.5))
})

test_that("alleleQc with no matching positions returns empty", {
  skip_on_covr()
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
  expect_equal(nrow(result$harmonizedData), 0)
})

test_that("alleleQc preserves extra columns", {
  skip_on_covr()
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G",
    beta = 0.5, se = 0.1
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_true("beta" %in% colnames(result$harmonizedData))
  expect_true("se" %in% colnames(result$harmonizedData))
})

test_that("alleleQc with lowercase alleles", {
  skip_on_covr()
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "a", A1 = "g"
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 1)
})

# ---- sanitize_names edge cases (alleleQc.R lines 37, 42) ----
test_that("alleleQc handles data frame with NULL colnames after merge", {
  skip_on_covr()
  # Create a data frame where merge might produce empty names
  # by giving target_data a column with NA name
  target <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  colnames(target)[1] <- ""
  colnames(target) <- make.unique(colnames(target), sep = "_")
  # Restore chrom for the join
  colnames(target)[1] <- "chrom"
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 1)
})

# ---- target_data with redundant columns (alleleQc.R line 75) ----
test_that("alleleQc removes redundant columns from target_data before join", {
  skip_on_covr()
  target <- data.frame(
    chrom = 1, pos = 100, A2 = "A", A1 = "G",
    variant_id = "1:100:A:G", chromosome = "chr1", position = 100
  )
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  result <- suppressWarnings(alleleQc(target, ref, matchMinProp = 0))
  expect_equal(nrow(result$harmonizedData), 1)
  # The redundant columns should have been removed before the join
  expect_true("variant_id" %in% colnames(result$harmonizedData))
})

# ---- col_to_flip with nonexistent column (alleleQc.R line 130) ----
test_that("alleleQc errors when col_to_flip column does not exist", {
  skip_on_covr()
  target <- data.frame(chrom = 1, pos = 100, A2 = "G", A1 = "A")
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  expect_error(
    suppressWarnings(
      alleleQc(target, ref, colToFlip = "nonexistent_col", matchMinProp = 0)),
    "not found in targetData"
  )
})

# Duplicate-handling is no longer the responsibility of alleleQc /
# matchRefPanel: callers are expected to deduplicate (via MungeSumstats /
# summaryStatsQc) before harmonization. The new behavior is to error.
test_that("alleleQc errors on duplicate variants in target input", {
  skip_on_covr()
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 100),
    A2 = c("A", "A"), A1 = c("G", "G"),
    beta = c(0.5, 0.6)
  )
  ref <- data.frame(chrom = 1, pos = 100, A2 = "A", A1 = "G")
  expect_error(
    suppressWarnings(alleleQc(target, ref, matchMinProp = 0)),
    "Duplicated variant IDs"
  )
})

# ===========================================================================
# matchRefPanel (deprecated; replaced by summaryStatsQc())
# colToComplement hook (rss-qc-parity): af complemented on allele swap
# ===========================================================================

test_that("af is complemented (1 - af) when harmonization swaps the effect allele", {
  skip_on_covr()
  ref <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                    A2 = c("A", "C"), A1 = c("G", "T"), stringsAsFactors = FALSE)
  # chr1:100 alleles swapped vs ref (=> sign flip); chr1:200 exact match.
  target <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                       A2 = c("G", "C"), A1 = c("A", "T"),
                       z = c(2.0, 1.5), af = c(0.30, 0.40), stringsAsFactors = FALSE)

  res <- suppressWarnings(
    matchRefPanel(target, ref, colToFlip = "z", colToComplement = "af"))
  h <- res$harmonizedData
  swapped <- h[h$pos == 100, ]
  control <- h[h$pos == 200, ]

  expect_equal(swapped$af, 0.70)   # 1 - input af
  expect_equal(swapped$z, -2.0)    # signed columns still sign-flip
  expect_equal(control$af, 0.40)   # untouched (no swap)
  expect_equal(control$z, 1.5)
})

test_that("colToComplement default leaves af unchanged (non-RSS callers unaffected)", {
  skip_on_covr()
  ref <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                    A2 = c("A", "C"), A1 = c("G", "T"), stringsAsFactors = FALSE)
  target <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 200),
                       A2 = c("G", "C"), A1 = c("A", "T"),
                       z = c(2.0, 1.5), af = c(0.30, 0.40), stringsAsFactors = FALSE)

  res <- suppressWarnings(
    matchRefPanel(target, ref, colToFlip = "z"))  # default: no complement
  h <- res$harmonizedData
  swapped <- h[h$pos == 100, ]
  expect_equal(swapped$af, 0.30)   # unchanged
  expect_equal(swapped$z, -2.0)    # z still sign-flips (independent path)
})

test_that("colToComplement errors on a missing column name", {
  skip_on_covr()
  ref <- data.frame(chrom = "chr1", pos = 100, A2 = "A", A1 = "G", stringsAsFactors = FALSE)
  target <- data.frame(chrom = "chr1", pos = 100, A2 = "A", A1 = "G",
                       z = 1.0, stringsAsFactors = FALSE)
  expect_error(
    suppressWarnings(matchRefPanel(target, ref, colToComplement = "af")),
    "not found in targetData"
  )
})

# ===========================================================================
# enlocPipeline (deprecated; integrated into colocPipeline)
# Input-type validation
# ===========================================================================

test_that("enlocPipeline: rejects non-QtlFineMappingResult qtlFmr", {
  skip_on_covr()
  expect_error(
    suppressWarnings(enlocPipeline(qtlFineMappingResult = "no",
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment = data.frame(gwasStudy = "G1", qtlContext = "c1",
                                           enrichment = 2.0,
                                           stringsAsFactors = FALSE))),
    "must be a QtlFineMappingResult"
  )
})

test_that("enlocPipeline: rejects gwasInput that is neither GwasSumStats nor GwasFineMappingResult", {
  skip_on_covr()
  expect_error(
    suppressWarnings(enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = 42L,
                  enrichment = data.frame(gwasStudy = "G1", qtlContext = "c1",
                                           enrichment = 2.0,
                                           stringsAsFactors = FALSE))),
    "must be a GwasSumStats or a GwasFineMappingResult"
  )
})

test_that("enlocPipeline: enrichment must be a data.frame", {
  skip_on_covr()
  expect_error(
    suppressWarnings(enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = "not a df")),
    "must be a data.frame"
  )
})

test_that("enlocPipeline: enrichment missing required columns errors", {
  skip_on_covr()
  expect_error(
    suppressWarnings(enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = data.frame(gwasStudy = "G1"))),
    "missing column"
  )
})

test_that("enlocPipeline: un-QCd GwasSumStats input is rejected", {
  skip_on_covr()
  expect_error(
    suppressWarnings(enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasSumstats(qc = FALSE),
                  enrichment = data.frame(gwasStudy = "G1", qtlContext = "c1",
                                           enrichment = 2.0,
                                           stringsAsFactors = FALSE))),
    "has no QC record"
  )
})

# ===========================================================================
# enlocPipeline — Pair loop (runs end-to-end via the LBF + coloc.bf_bf path)
# ===========================================================================

test_that("enlocPipeline: pair loop produces one row per (QTL tuple, GWAS tuple) with adjusted p12", {
  skip_on_covr()
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = enr,
                  p12                  = 5e-6,
                  p12Max               = 1e-3))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_equal(out$enrichment, 2.0)
  # 5e-6 * (1 + 2.0) = 1.5e-5 < 1e-3, so capped value is the raw product.
  expect_equal(out$p12Used, 1.5e-5)
  expect_equal(out$p12_actual, 1.5e-5)
})

test_that("enlocPipeline: missing-enrichment pair falls back to baseline p12 with a warning", {
  skip_on_covr()
  # An enrichment frame that has no row for (G1, c1).
  enr <- data.frame(gwasStudy = "G_other", qtlContext = "c_other",
                    enrichment = 10.0, stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  expect_warning(
    out <- enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                          gwasInput            = .ep_makeGwasFmr(),
                          enrichment           = enr,
                          p12                  = 5e-6,
                          p12Max               = 1e-3),
    "no enrichment entry"
  )
  # Baseline p12 unchanged because the pair fell back (enRow = 0).
  expect_equal(out$p12Used, 5e-6)
})

test_that("enlocPipeline: p12Max caps the adjusted prior", {
  skip_on_covr()
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 1e6,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = enr,
                  p12                  = 5e-6,
                  p12Max               = 1e-4))
  expect_equal(out$p12Used, 1e-4)
})

test_that("enlocPipeline: coloc.bf_bf failures are caught and warned, pair skipped", {
  skip_on_covr()
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(
    coloc.bf_bf = function(q, g, ...) stop("boom"),
    .package = "coloc")
  expect_warning(
    out <- enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                          gwasInput            = .ep_makeGwasFmr(),
                          enrichment           = enr),
    "coloc.bf_bf failed"
  )
  expect_equal(nrow(out), 0L)
})

test_that("enlocPipeline: returnGwasFineMapping=TRUE attaches gwasFineMapping attr (non-empty result)", {
  skip_on_covr()
  enr <- data.frame(gwasStudy = "G1", qtlContext = "c1", enrichment = 2.0,
                    stringsAsFactors = FALSE)
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  # Use GwasFineMappingResult input: returnGwasFineMapping has no effect
  # for this branch (only GwasSumStats triggers attachment). To trigger
  # attachment we mock fineMappingPipeline so the GwasSumStats path
  # produces a usable FMR.
  fakeFmr <- .ep_makeGwasFmr()
  local_mocked_bindings(
    fineMappingPipeline = function(data, ...) fakeFmr,
    .package = "pecotmr")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult  = .ep_makeQtlFmr(),
                  gwasInput             = .ep_makeGwasSumstats(),
                  enrichment            = enr,
                  returnGwasFineMapping = TRUE))
  expect_true("gwasFineMapping" %in% names(attributes(out)))
  expect_s4_class(attr(out, "gwasFineMapping"), "GwasFineMappingResult")
})

test_that("enlocPipeline: qLbf NULL (QTL entry's LBF rows drop after priorTol) skips that QTL row", {
  skip_on_covr()
  # Build a QTL FMR whose lbf_variable is empty after the V > priorTol filter
  # (V = 0 < default priorTol 1e-9 -> drop all rows -> return NULL).
  emptyFit <- list(
    lbf_variable = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "v1")),
    V = 0.0)
  e <- FineMappingEntry(variantIds = "v1",
                        trimmedFit = emptyFit,
                        topLoci = data.frame(variant_id = "v1", pip = 0,
                                              stringsAsFactors = FALSE))
  qfmr <- QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(e),
    ldSketch = .ep_makeHandle())
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = .ep_makeGwasFmr(),
                  enrichment           = data.frame(
                    gwasStudy = "G1", qtlContext = "c1",
                    enrichment = 1.0, stringsAsFactors = FALSE)))
  expect_equal(nrow(out), 0L)
})

test_that("enlocPipeline: aligned NULL (disjoint variant sets) skips that pair", {
  skip_on_covr()
  # QTL fmr with variant ids that don't overlap the GWAS variant ids.
  qVids <- paste0("chr1:", 100*(1:5), ":A:G")
  gVids <- paste0("chr2:", 200*(1:5), ":A:G")
  qfmr <- QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(.ep_makeFmEntry(variant_ids = qVids)),
    ldSketch = .ep_makeHandle())
  gfmr <- GwasFineMappingResult(
    study   = "G1", method = "susie",
    entry   = list(.ep_makeFmEntry(variant_ids = gVids)),
    ldSketch = .ep_makeHandle())
  local_mocked_bindings(coloc.bf_bf = .ep_mockColocBfBf(), .package = "coloc")
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = qfmr,
                  gwasInput            = gfmr,
                  enrichment           = data.frame(
                    gwasStudy = "G1", qtlContext = "c1",
                    enrichment = 1.0, stringsAsFactors = FALSE)))
  expect_equal(nrow(out), 0L)
})

test_that("enlocPipeline: empty result schema includes enrichment + p12Used", {
  skip_on_covr()
  # Build a GWAS FMR whose entry has no usable LBF -> pre-extract returns empty.
  emptyFit <- list(alpha = matrix(0, 1, 1), pip = c(v1 = 0),
                   V = 0, lbf_variable = matrix(NA_real_, 1, 1))
  e <- FineMappingEntry(variantIds = "v1",
                        trimmedFit = emptyFit,
                        topLoci = data.frame(variant_id = "v1", pip = 0,
                                              stringsAsFactors = FALSE))
  gfmr <- GwasFineMappingResult(
    study = "G1", method = "susie",
    entry = list(e),
    ldSketch = .ep_makeHandle())
  out <- suppressWarnings(
    enlocPipeline(qtlFineMappingResult = .ep_makeQtlFmr(),
                  gwasInput            = gfmr,
                  enrichment           = data.frame(
                    gwasStudy = "G1", qtlContext = "c1",
                    enrichment = 1.0, stringsAsFactors = FALSE)))
  expect_equal(nrow(out), 0L)
  expect_true(all(c("enrichment", "p12Used") %in% colnames(out)))
})

# ===========================================================================
# computeQtlEnrichment (deprecated; renamed to qtlEnrichment())
# ===========================================================================

test_that("computeQtlEnrichment dummy data single-threaded works",{
  skip_on_covr()
  local_mocked_bindings(
      qtlEnrichmentRcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  expect_warning(
    computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, impN = 10, numThreads = 1),
    "numGwas is not provided. Estimating piGwas from the data. Note that this estimate may be biased if the input gwasPip does not contain genome-wide variants.")
  expect_warning(
    computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, impN = 10, numThreads = 1),
    "piQtl is not provided. Estimating piQtl from the data. Note that this estimate may be biased if either 1) the input susieQtlRegions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
  res <- suppressWarnings(computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, numGwas=5000, piQtl=0.49819, lambda = 1, impN = 10, numThreads = 1))
  expect_true(length(res) > 0)
})

test_that("computeQtlEnrichment dummy data single thread and multi-threaded are equivalent",{
  skip_on_covr()
  local_mocked_bindings(
      qtlEnrichmentRcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  res_single <- suppressWarnings(computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, numGwas=5000, piQtl=0.49819, lambda = 1, impN = 10, numThreads = 1))
  res_multi <- suppressWarnings(computeQtlEnrichment(input_data$gwas_fit$pip, input_data$susie_fits, numGwas=5000, piQtl=0.49819, lambda = 1, impN = 10, numThreads = 2))
  expect_equal(res_single, res_multi)
})

# ---- error paths (computeQtlEnrichment.R lines 86, 87, 91) ----
test_that("computeQtlEnrichment errors when pi_gwas is zero", {
  skip_on_covr()
  gwas_pip <- rep(0, 10)
  names(gwas_pip) <- paste0("snp", 1:10)
  susie_fits <- list(fit1 = list(pip = setNames(runif(10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(computeQtlEnrichment(gwas_pip, susie_fits, piQtl = 0.5)),
    "No association signal found in GWAS data"
  )
})

test_that("computeQtlEnrichment errors when pi_qtl is zero", {
  skip_on_covr()
  gwas_pip <- runif(10)
  names(gwas_pip) <- paste0("snp", 1:10)
  susie_fits <- list(fit1 = list(pip = setNames(rep(0, 10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(computeQtlEnrichment(gwas_pip, susie_fits, numGwas = 1000, piQtl = 0)),
    "No QTL associated"
  )
})

test_that("computeQtlEnrichment errors when gwas_pip has no names", {
  skip_on_covr()
  gwas_pip <- runif(10)  # no names
  susie_fits <- list(fit1 = list(pip = setNames(runif(10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(computeQtlEnrichment(gwas_pip, susie_fits, numGwas = 1000, piQtl = 0.5)),
    "Variant names are missing in gwasPip"
  )
})

# ---- real C++ qtlEnrichmentRcpp integration test ----
test_that("computeQtlEnrichment calls real C++ enrichment code and returns expected keys", {
  skip_on_covr()
  set.seed(42)
  n_snps <- 50
  variantNames <- paste0("1:", 1:n_snps, ":A:G")

  # GWAS PIPs: sparse signal
  gwas_pip <- rep(0.01, n_snps)
  gwas_pip[c(5, 20, 35)] <- c(0.8, 0.6, 0.9)
  names(gwas_pip) <- variantNames

  # SuSiE fit with 2 single effects over same variants
  L <- 2
  alpha <- matrix(1 / n_snps, nrow = L, ncol = n_snps)
  # Concentrate probability on causal variants
  alpha[1, ] <- 0.001; alpha[1, 5] <- 0.95; alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
  alpha[2, ] <- 0.001; alpha[2, 20] <- 0.95; alpha[2, ] <- alpha[2, ] / sum(alpha[2, ])
  pip <- colSums(alpha)
  names(pip) <- variantNames

  susie_fits <- list(
    fit1 = list(pip = pip, alpha = alpha, prior_variance = c(0.5, 0.3))
  )

  # Call without mocking - exercises the real C++ code
  res <- suppressWarnings(
    computeQtlEnrichment(gwas_pip, susie_fits,
                           numGwas = 5000, piQtl = 0.5,
                           lambda = 1, impN = 5, numThreads = 1)
  )
  expect_type(res, "list")
  # The enrichment results are in res[[1]] (the C++ output list)
  en <- res[[1]]
  expected_keys <- c("Intercept", "Enrichment (no shrinkage)", "Enrichment (w/ shrinkage)",
                     "sd (no shrinkage)", "sd (w/ shrinkage)",
                     "Alternative (coloc) p1", "Alternative (coloc) p2", "Alternative (coloc) p12")
  for (key in expected_keys) {
    expect_true(key %in% names(en), info = paste("Missing key:", key))
  }
  # All numeric and finite
  numeric_vals <- unlist(en[expected_keys])
  expect_true(all(is.finite(numeric_vals)))
})

# ---- unmatched variants tracking (computeQtlEnrichment.R line 102) ----
test_that("computeQtlEnrichment tracks unmatched QTL variants", {
  skip_on_covr()
  local_mocked_bindings(
    qtlEnrichmentRcpp = function(...) TRUE
  )
  gwas_pip <- runif(10)
  names(gwas_pip) <- paste0("1:", 1:10, ":A:G")
  # QTL has some variants not in GWAS
  qtl_pip <- runif(5)
  names(qtl_pip) <- c(paste0("1:", 1:3, ":A:G"), "1:999:A:G", "1:998:A:G")
  susie_fits <- list(fit1 = list(pip = qtl_pip,
                                  alpha = matrix(runif(5), 1, 5),
                                  prior_variance = 1))
  res <- suppressWarnings(
    computeQtlEnrichment(gwas_pip, susie_fits, numGwas = 1000, piQtl = 0.5)
  )
  expect_true("unused_xqtl_variants" %in% names(res))
})

# ===========================================================================
# coloc/xqtl pre-S4 entry points (deprecated; superseded by the S4 pipelines)
# These are pure no-op shims: each fires .Deprecated() and returns NULL.
# ===========================================================================

test_that("xqtlEnrichmentWrapper is a deprecated no-op", {
  skip_on_covr()
  expect_warning(
    res <- xqtlEnrichmentWrapper(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})

test_that("colocWrapper is a deprecated no-op", {
  skip_on_covr()
  expect_warning(
    res <- colocWrapper(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})

test_that("colocPostProcessor is a deprecated no-op", {
  skip_on_covr()
  expect_warning(
    res <- colocPostProcessor(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})
