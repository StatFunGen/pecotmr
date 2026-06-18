# Tests for R/sldscWrapper.R
#
# Fixture convention:
#   - 2 chromosomes (1, 2), 50 SNPs each → 100 total
#   - 2 target annotations: "annot_A" (binary), "annot_B" (continuous)
#   - 97 baseline annotations (baselineLD_0 .. baselineLD_96) in joint run
#   - 10 jackknife blocks
#   - Polyfun appends "_0" to target annotation names in .results

# =============================================================================
# Fixture generators
# =============================================================================

# Create a single .annot.gz file for one chromosome
# Real polyfun .annot.gz files have CHR, SNP, BP, CM + annotation columns only
# (no MAF/A1/A2 — those come from the .frq / PLINK files).
.make_annot_gz <- function(dir, chrom, nSnps = 50) {
  df <- data.frame(
    CHR = chrom,
    SNP = paste0("rs", (chrom - 1L) * 100L + seq_len(nSnps)),
    BP  = seq_len(nSnps) * 1000L,
    CM  = seq_len(nSnps) * 0.01,
    annot_A = sample(c(0L, 1L), nSnps, replace = TRUE),
    annot_B = rnorm(nSnps, 2, 0.5),
    stringsAsFactors = FALSE
  )
  path <- file.path(dir, sprintf("target.%d.annot.gz", chrom))
  gz <- gzfile(path, "wb")
  vroom::vroom_write(df, gz, delim = "\t")
  close(gz)
  invisible(df)
}

# Create a PLINK .frq file for one chromosome
.make_frq <- function(dir, chrom, plinkName ="ref_chr", nSnps = 50) {
  df <- data.frame(
    CHR = chrom,
    SNP = paste0("rs", (chrom - 1L) * 100L + seq_len(nSnps)),
    A1  = "A",
    A2  = "G",
    MAF = runif(nSnps, 0.01, 0.49),
    NCHROBS = 200L,
    stringsAsFactors = FALSE
  )
  path <- file.path(dir, sprintf("%s%d.frq", plinkName, chrom))
  vroom::vroom_write(df, path, delim = "\t")
  invisible(df)
}

# Create the three polyfun output files (.results, .log, .part_delete)
# for a single-target run. Real polyfun output includes baseline categories
# even in single-target mode, so we add 2 dummy baseline categories.
.make_polyfun_single <- function(dir, prefix, target_name, nBlocks = 10,
                                 h2g = 0.3, tau = 1e-7, enrichment = 2.5,
                                 n_baseline = 2) {
  target_cat <- paste0(target_name, "_0")
  baseline_cats <- paste0("baselineLD_", seq_len(n_baseline) - 1L)
  all_cats <- c(target_cat, baseline_cats)
  n_cats <- length(all_cats)

  taus_all <- c(tau, rep(1e-8, n_baseline))
  enrichments_all <- c(enrichment, rep(1.0, n_baseline))

  results <- data.frame(
    Category                = all_cats,
    Coefficient             = taus_all,
    Coefficient_std_error   = abs(taus_all) * 0.3,
    Enrichment              = enrichments_all,
    Enrichment_std_error    = enrichments_all * 0.2,
    Enrichment_p            = rep(0.01, n_cats),
    `Prop._h2`              = c(0.15, rep(0.425, n_baseline)),
    `Prop._SNPs`            = c(0.06, rep(0.47, n_baseline)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  vroom::vroom_write(results, paste0(prefix, ".results"), delim = "\t")

  writeLines(c(
    "Analysis started at 2024-01-01",
    sprintf("Total Observed scale h2: %g (0.05)", h2g),
    "Analysis finished"
  ), paste0(prefix, ".log"))

  blocks <- matrix(rnorm(nBlocks * n_cats,
                         mean = rep(taus_all, each = nBlocks),
                         sd = abs(rep(taus_all, each = nBlocks)) * 0.5),
                   nrow = nBlocks, ncol = n_cats)
  colnames(blocks) <- all_cats
  vroom::vroom_write(as.data.frame(blocks), paste0(prefix, ".part_delete"), delim = "\t")
  invisible(NULL)
}

# Create polyfun output files for a joint run (target + baseline annotations)
.make_polyfun_joint <- function(dir, prefix, target_names,
                                n_baseline = 3, nBlocks = 10, h2g = 0.3) {
  target_cats <- paste0(target_names, "_0")
  baseline_cats <- paste0("baselineLD_", seq_len(n_baseline) - 1L)
  all_cats <- c(target_cats, baseline_cats)
  n_cats <- length(all_cats)

  taus <- c(rep(1e-7, length(target_cats)), rep(1e-8, n_baseline))
  enrichments <- c(rep(2.0, length(target_cats)), rep(1.0, n_baseline))

  results <- data.frame(
    Category                = all_cats,
    Coefficient             = taus,
    Coefficient_std_error   = abs(taus) * 0.3,
    Enrichment              = enrichments,
    Enrichment_std_error    = enrichments * 0.2,
    Enrichment_p            = rep(0.05, n_cats),
    `Prop._h2`              = rep(1 / n_cats, n_cats),
    `Prop._SNPs`            = rep(1 / n_cats, n_cats),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  vroom::vroom_write(results, paste0(prefix, ".results"), delim = "\t")

  writeLines(c(
    "Analysis started at 2024-01-01",
    sprintf("Total Observed scale h2: %g (0.05)", h2g),
    "Analysis finished"
  ), paste0(prefix, ".log"))

  blocks <- matrix(rnorm(nBlocks * n_cats, mean = rep(taus, each = nBlocks),
                         sd = abs(rep(taus, each = nBlocks)) * 0.5),
                   nrow = nBlocks, ncol = n_cats)
  colnames(blocks) <- all_cats
  vroom::vroom_write(as.data.frame(blocks), paste0(prefix, ".part_delete"), delim = "\t")
  invisible(NULL)
}


# Build a complete fixture directory for the full pipeline
.make_sldsc_fixtures <- function(envir = parent.frame()) {
  base_dir <- withr::local_tempdir(.local_envir = envir)

  anno_dir <- file.path(base_dir, "annot")
  frq_dir  <- file.path(base_dir, "frq")
  out_dir  <- file.path(base_dir, "output")
  dir.create(anno_dir)
  dir.create(frq_dir)
  dir.create(out_dir)

  plink_name <- "ref_chr"

  # Annotation + freq files for 2 chromosomes
  for (chr in 1:2) {
    .make_annot_gz(anno_dir, chr)
    .make_frq(frq_dir, chr, plinkName =plink_name)
  }

  targets <- c("annot_A", "annot_B")

  # Single-target runs: 2 targets x 2 traits
  for (trait in c("traitX", "traitY")) {
    for (i in seq_along(targets)) {
      pref <- file.path(out_dir, sprintf("%s_single_%s", trait, targets[i]))
      .make_polyfun_single(out_dir, pref, targets[i], h2g = 0.3 + (i - 1) * 0.05)
    }
  }

  # Joint runs: 1 per trait

  for (trait in c("traitX", "traitY")) {
    pref <- file.path(out_dir, sprintf("%s_joint", trait))
    .make_polyfun_joint(out_dir, pref, targets, h2g = 0.3)
  }

  list(
    base_dir   = base_dir,
    anno_dir   = anno_dir,
    frq_dir    = frq_dir,
    out_dir    = out_dir,
    plinkName =plink_name,
    targets    = targets,
    trait_names = c("traitX", "traitY")
  )
}


# =============================================================================
# .sldscChromFromFilename
# =============================================================================

test_that(".sldscChromFromFilename parses chromosome number", {
  fn <- pecotmr:::.sldscChromFromFilename
  expect_equal(fn("target.1.annot.gz"), 1L)
  expect_equal(fn("target.22.annot.gz"), 22L)
  expect_true(is.na(fn("no_chrom_here.txt")))
  expect_true(is.na(fn("target.X.annot.gz")))
})


# =============================================================================
# .sldscDetectAnnotCols
# =============================================================================

test_that(".sldscDetectAnnotCols finds non-standard columns", {
  dir <- withr::local_tempdir()
  .make_annot_gz(dir, 1)
  f <- file.path(dir, "target.1.annot.gz")
  cols <- pecotmr:::.sldscDetectAnnotCols(f)
  expect_true("annot_A" %in% cols)
  expect_true("annot_B" %in% cols)
  expect_false("CHR" %in% cols)
  expect_false("SNP" %in% cols)
})


# =============================================================================
# readSldscTrait
# =============================================================================

test_that("readSldscTrait reads polyfun outputs correctly", {
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "test_trait")
  .make_polyfun_single(dir, prefix, "myannot", nBlocks = 5, h2g = 0.25)

  result <- readSldscTrait(prefix)
  expect_true(is.list(result))
  # 1 target + 2 baseline categories
  expect_true("myannot_0" %in% result$categories)
  expect_equal(result$h2g, 0.25)
  expect_equal(result$nBlocks, 5L)
  expect_equal(length(result$tau), 3L)
  expect_true("myannot_0" %in% names(result$tau))
  expect_true(is.matrix(result$tauBlocks))
  expect_equal(nrow(result$tauBlocks), 5L)
  expect_equal(ncol(result$tauBlocks), 3L)
})

test_that("readSldscTrait errors on missing files", {
  expect_error(readSldscTrait("/nonexistent/prefix"), "missing file")
})

test_that("readSldscTrait errors when h2 not in log", {
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "bad_log")
  .make_polyfun_single(dir, prefix, "a", nBlocks = 3)
  # Overwrite the log with no h2 line
  writeLines("No heritability here", paste0(prefix, ".log"))
  expect_error(readSldscTrait(prefix), "Total Observed scale h2")
})

test_that("readSldscTrait errors on column mismatch in part_delete", {
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "bad_delete")
  .make_polyfun_single(dir, prefix, "a", nBlocks = 3)
  # Overwrite part_delete with wrong number of columns (need != 3)
  vroom::vroom_write(data.frame(x = 1:3, y = 4:6, z = 7:9, w = 10:12),
                     paste0(prefix, ".part_delete"), delim = "\t")
  expect_error(readSldscTrait(prefix), "part_delete")
})


# =============================================================================
# computeSldscAnnotSd
# =============================================================================

test_that("computeSldscAnnotSd computes SDs with MAF filtering", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  frq_dir  <- file.path(dir, "frq")
  dir.create(anno_dir)
  dir.create(frq_dir)
  for (chr in 1:2) {
    .make_annot_gz(anno_dir, chr)
    .make_frq(frq_dir, chr, plinkName ="ref_chr")
  }

  sds <- computeSldscAnnotSd(anno_dir, frqfileDir =frq_dir,
                                 plinkName ="ref_chr", mafCutoff =0.05)
  expect_true(is.numeric(sds))
  expect_equal(length(sds), 2L)
  expect_named(sds, c("annot_A", "annot_B"))
  expect_true(all(sds > 0))
})

test_that("computeSldscAnnotSd works with mafCutoff =0", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  for (chr in 1:2) .make_annot_gz(anno_dir, chr)

  sds <- computeSldscAnnotSd(anno_dir, mafCutoff =0)
  expect_true(all(sds > 0))
})

test_that("computeSldscAnnotSd respects annot_cols argument (character)", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  for (chr in 1:2) .make_annot_gz(anno_dir, chr)

  sds <- computeSldscAnnotSd(anno_dir, mafCutoff =0,
                                 annotCols ="annot_A")
  expect_equal(length(sds), 1L)
  expect_named(sds, "annot_A")
})

test_that("computeSldscAnnotSd respects annot_cols argument (numeric)", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  for (chr in 1:2) .make_annot_gz(anno_dir, chr)

  sds <- computeSldscAnnotSd(anno_dir, mafCutoff =0, annotCols =2L)
  expect_equal(length(sds), 1L)
  expect_named(sds, "annot_B")
})

test_that("computeSldscAnnotSd errors on missing frqfile_dir when maf > 0", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  .make_annot_gz(anno_dir, 1)
  expect_error(computeSldscAnnotSd(anno_dir, frqfileDir =NULL, mafCutoff =0.05),
               "frqfileDir")
})

test_that("computeSldscAnnotSd errors on missing anno dir", {
  expect_error(computeSldscAnnotSd("/nonexistent/dir", mafCutoff =0),
               "does not exist")
})

test_that("computeSldscAnnotSd errors on empty anno dir", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "empty")
  dir.create(anno_dir)
  expect_error(computeSldscAnnotSd(anno_dir, mafCutoff =0), "no .annot.gz")
})

test_that("computeSldscAnnotSd errors on missing .frq file", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  frq_dir  <- file.path(dir, "frq")
  dir.create(anno_dir)
  dir.create(frq_dir)
  .make_annot_gz(anno_dir, 1)
  # No .frq file for chr1
  expect_error(computeSldscAnnotSd(anno_dir, frqfileDir =frq_dir,
                                       plinkName ="ref_chr", mafCutoff =0.05),
               "frq file not found")
})


# =============================================================================
# computeSldscMRef
# =============================================================================

test_that("computeSldscMRef counts SNPs from .frq files with MAF cutoff", {
  dir <- withr::local_tempdir()
  frq_dir <- file.path(dir, "frq")
  dir.create(frq_dir)
  for (chr in 1:2) .make_frq(frq_dir, chr, plinkName ="ref_chr")

  M <- computeSldscMRef(frqfileDir =frq_dir, plinkName ="ref_chr",
                            mafCutoff =0.05)
  expect_true(is.integer(M))
  expect_true(M > 0)
  expect_true(M <= 100L)  # 2 chroms x 50 SNPs max
})

test_that("computeSldscMRef counts all SNPs when mafCutoff =0", {
  dir <- withr::local_tempdir()
  frq_dir <- file.path(dir, "frq")
  dir.create(frq_dir)
  for (chr in 1:2) .make_frq(frq_dir, chr, plinkName ="ref_chr")

  M <- computeSldscMRef(frqfileDir =frq_dir, plinkName ="ref_chr",
                            mafCutoff =0)
  expect_equal(M, 100L)
})

test_that("computeSldscMRef falls back to .l2.ldscore with mafCutoff =0", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  # Create a fake .l2.ldscore.gz file with 40 rows
  df <- data.frame(CHR = 1, SNP = paste0("rs", 1:40), BP = 1:40, L2 = runif(40))
  gz <- gzfile(file.path(anno_dir, "scores.l2.ldscore.gz"), "wb")
  vroom::vroom_write(df, gz, delim = "\t")
  close(gz)

  M <- computeSldscMRef(targetAnnoDir =anno_dir, mafCutoff =0)
  expect_equal(M, 40L)
})

test_that("computeSldscMRef errors when maf > 0 and no frq dir", {
  expect_error(computeSldscMRef(mafCutoff =0.05), "frqfileDir")
})

test_that("computeSldscMRef errors with no dirs at all", {
  expect_error(computeSldscMRef(mafCutoff =0), "need frqfileDir")
})

test_that("computeSldscMRef warns on l2.ldscore fallback", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  df <- data.frame(CHR = 1, SNP = paste0("rs", 1:20), BP = 1:20, L2 = runif(20))
  gz <- gzfile(file.path(anno_dir, "scores.l2.ldscore.gz"), "wb")
  vroom::vroom_write(df, gz, delim = "\t")
  close(gz)

  expect_warning(
    computeSldscMRef(targetAnnoDir =anno_dir, mafCutoff =0),
    "UNDERCOUNTS"
  )
})

test_that("computeSldscMRef uses generic .frq glob when plink_name pattern fails", {
  dir <- withr::local_tempdir()
  frq_dir <- file.path(dir, "frq")
  dir.create(frq_dir)
  # Name doesn't match the plink_name pattern
  for (chr in 1:2) .make_frq(frq_dir, chr, plinkName ="other_chr")

  M <- computeSldscMRef(frqfileDir =frq_dir, plinkName ="nomatch_chr",
                            mafCutoff =0)
  expect_equal(M, 100L)
})


# =============================================================================
# isBinarySldscAnnot
# =============================================================================

test_that("isBinarySldscAnnot detects binary and continuous annotations", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  for (chr in 1:2) .make_annot_gz(anno_dir, chr)

  result <- isBinarySldscAnnot(anno_dir)
  expect_true(is.logical(result))
  expect_named(result, c("annot_A", "annot_B"))
  expect_true(result[["annot_A"]])   # binary (0/1)
  expect_false(result[["annot_B"]])  # continuous
})

test_that("isBinarySldscAnnot respects annot_cols (character)", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  for (chr in 1:2) .make_annot_gz(anno_dir, chr)

  result <- isBinarySldscAnnot(anno_dir, annotCols ="annot_A")
  expect_equal(length(result), 1L)
  expect_true(result[["annot_A"]])
})

test_that("isBinarySldscAnnot respects annot_cols (numeric)", {
  dir <- withr::local_tempdir()
  anno_dir <- file.path(dir, "annot")
  dir.create(anno_dir)
  for (chr in 1:2) .make_annot_gz(anno_dir, chr)

  result <- isBinarySldscAnnot(anno_dir, annotCols =2L)
  expect_equal(length(result), 1L)
  expect_named(result, "annot_B")
})

test_that("isBinarySldscAnnot errors on empty dir", {
  dir <- withr::local_tempdir()
  expect_error(isBinarySldscAnnot(dir), "no .annot.gz")
})


# =============================================================================
# standardizeSldscTrait
# =============================================================================

# Helper to build a trait_data list (as from readSldscTrait)
.make_trait_data <- function(cats = c("A_0", "B_0"), nBlocks = 10, h2g = 0.3) {
  n <- length(cats)
  taus <- rep(1e-7, n)
  blocks <- matrix(rnorm(nBlocks * n, mean = rep(taus, each = nBlocks),
                         sd = 1e-8),
                   nrow = nBlocks, ncol = n)
  colnames(blocks) <- cats

  list(
    categories     = cats,
    tau            = setNames(taus, cats),
    tauSe         = setNames(abs(taus) * 0.3, cats),
    enrichment     = setNames(rep(2.0, n), cats),
    enrichmentSe  = setNames(rep(0.4, n), cats),
    enrichmentP   = setNames(rep(0.01, n), cats),
    propH2        = setNames(rep(0.15, n), cats),
    propSnps      = setNames(rep(0.06, n), cats),
    h2g            = h2g,
    tauBlocks     = blocks,
    nBlocks       = nBlocks
  )
}

test_that("standardizeSldscTrait works in single mode", {
  td <- .make_trait_data()
  sd_annot <- c(A_0 = 0.5, B_0 = 1.2)
  M_ref <- 1000L

  result <- standardizeSldscTrait(td, sd_annot, M_ref, mode = "single")
  expect_true(is.list(result))
  expect_equal(result$mode, "single")
  expect_equal(result$h2g, 0.3)
  expect_equal(result$nBlocks, 10L)

  s <- result$summary
  expect_true(is.data.frame(s))
  expect_equal(nrow(s), 2L)
  expect_true(all(c("tauStar", "tauStarSe", "enrichment", "enrichmentSe",
                     "enrichmentP", "enrichstat", "enrichstatSe") %in% names(s)))

  # tauStar = tau * sd * M_ref / h2g
  expected_ts <- unname(td$tau) * c(0.5, 1.2) * 1000 / 0.3
  expect_equal(s$tauStar, expected_ts)

  # tau_star_blocks has correct dimensions
  expect_true(is.matrix(result$tau_star_blocks))
  expect_equal(dim(result$tau_star_blocks), c(10L, 2L))
})

test_that("standardizeSldscTrait works in joint mode", {
  td <- .make_trait_data()
  sd_annot <- c(A_0 = 0.5, B_0 = 1.2)
  M_ref <- 1000L

  result <- standardizeSldscTrait(td, sd_annot, M_ref, mode = "joint")
  expect_equal(result$mode, "joint")
  s <- result$summary
  # joint mode should NOT have enrichment columns
  expect_false("enrichment" %in% names(s))
  expect_true("tauStar" %in% names(s))
})

test_that("standardizeSldscTrait auto-detects target categories", {
  td <- .make_trait_data()
  sd_annot <- c(A_0 = 0.5)  # only one overlaps

  result <- standardizeSldscTrait(td, sd_annot, 1000L, mode = "joint")
  expect_equal(nrow(result$summary), 1L)
  expect_equal(result$summary$target, "A_0")
})

test_that("standardizeSldscTrait errors on empty categories", {
  td <- .make_trait_data()
  sd_annot <- c(X_0 = 0.5)  # no overlap
  expect_error(standardizeSldscTrait(td, sd_annot, 1000L), "no target categories")
})

test_that("standardizeSldscTrait errors on missing categories", {
  td <- .make_trait_data(cats = c("A_0"))
  sd_annot <- c(A_0 = 0.5, B_0 = 1.0)
  expect_error(
    standardizeSldscTrait(td, sd_annot, 1000L,
                            targetCategories =c("A_0", "B_0")),
    "missing categories"
  )
})

test_that("standardizeSldscTrait warns on zero sd", {
  td <- .make_trait_data(cats = "A_0")
  sd_annot <- c(A_0 = 0)
  expect_warning(
    standardizeSldscTrait(td, sd_annot, 1000L, mode = "joint"),
    "zero/NA sd"
  )
})

test_that("standardizeSldscTrait enrichstatSe handles p = 0", {
  td <- .make_trait_data(cats = "A_0")
  td$enrichmentP <- c(A_0 = 0)  # p = 0 → abs_z = Inf → se = 0
  sd_annot <- c(A_0 = 0.5)

  result <- standardizeSldscTrait(td, sd_annot, 1000L, mode = "single")
  # enrichstatSe should be NA when abs_z is infinite
  expect_true(is.na(result$summary$enrichstatSe))
})


# =============================================================================
# metaSldscRandom
# =============================================================================

# Helper to build per_trait_estimates for metaSldscRandom
.make_per_trait_meta <- function(nTraits = 3, category = "A_0",
                                 means = NULL, ses = NULL) {
  if (is.null(means)) means <- rnorm(nTraits, 1e-5, 1e-6)
  if (is.null(ses)) ses <- rep(1e-6, nTraits)
  per_trait <- list()
  for (i in seq_len(nTraits)) {
    per_trait[[paste0("trait", i)]] <- list(
      summary = data.frame(
        target      = category,
        tauStar    = means[i],
        tauStarSe = ses[i],
        enrichment    = means[i] * 100,
        enrichmentSe = ses[i] * 50,
        enrichstat    = means[i] * 10,
        enrichstatSe = ses[i] * 5,
        stringsAsFactors = FALSE
      )
    )
  }
  per_trait
}

test_that("metaSldscRandom works for tauStar", {
  pt <- .make_per_trait_meta(nTraits = 4)
  result <- metaSldscRandom(pt, "A_0", quantity = "tauStar")
  expect_true(is.list(result))
  expect_equal(result$nTraits, 4L)
  expect_true(is.numeric(result$mean))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$p))
  expect_true(result$se > 0)
  expect_equal(length(result$traitsUsed), 4L)
})

test_that("metaSldscRandom works for enrichment", {
  pt <- .make_per_trait_meta(nTraits = 3)
  result <- metaSldscRandom(pt, "A_0", quantity = "enrichment")
  expect_equal(result$nTraits, 3L)
  expect_true(is.finite(result$mean))
})

test_that("metaSldscRandom works for enrichstat", {
  pt <- .make_per_trait_meta(nTraits = 3)
  result <- metaSldscRandom(pt, "A_0", quantity = "enrichstat")
  expect_equal(result$nTraits, 3L)
  expect_true(is.finite(result$mean))
})

test_that("metaSldscRandom returns NA with < 2 traits", {
  pt <- .make_per_trait_meta(nTraits = 1)
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_true(is.na(result$mean))
  expect_true(is.na(result$se))
  expect_true(is.na(result$p))
  expect_equal(result$nTraits, 1L)
})

test_that("metaSldscRandom skips traits with missing category", {
  pt <- .make_per_trait_meta(nTraits = 3)
  # Change category in trait2
  pt$trait2$summary$target <- "other"
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_equal(result$nTraits, 2L)
  expect_equal(result$traitsUsed, c("trait1", "trait3"))
})

test_that("metaSldscRandom skips traits with NA or zero SE", {
  pt <- .make_per_trait_meta(nTraits = 3)
  pt$trait2$summary$tauStarSe <- NA
  pt$trait3$summary$tauStarSe <- 0
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_equal(result$nTraits, 1L)
  expect_true(is.na(result$mean))  # < 2 valid
})

test_that("metaSldscRandom skips NULL entries", {
  pt <- .make_per_trait_meta(nTraits = 3)
  pt$trait2 <- NULL
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_equal(result$nTraits, 2L)
})

test_that("metaSldscRandom generates names for unnamed list", {
  pt <- .make_per_trait_meta(nTraits = 2)
  names(pt) <- NULL
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_equal(result$traitsUsed, c("1", "2"))
})


# =============================================================================
# .sldscAssembleTraitSummary
# =============================================================================

test_that(".sldscAssembleTraitSummary combines single and joint", {
  fn <- pecotmr:::.sldscAssembleTraitSummary
  targets <- c("A_0", "B_0")
  is_bin <- c(A_0 = TRUE, B_0 = FALSE)

  single_df <- data.frame(
    target = targets,
    tau = c(1e-7, 2e-7), tauSe = c(3e-8, 4e-8),
    tauStar = c(0.01, 0.02), tauStarSe = c(0.003, 0.004),
    enrichment = c(2.0, 3.0), enrichmentSe = c(0.4, 0.6),
    enrichmentP = c(0.01, 0.05),
    enrichstat = c(0.001, 0.002), enrichstatSe = c(0.0003, 0.0004),
    stringsAsFactors = FALSE
  )

  joint_df <- data.frame(
    target = targets,
    tau = c(1.1e-7, 2.1e-7), tauSe = c(3.1e-8, 4.1e-8),
    tauStar = c(0.011, 0.021), tauStarSe = c(0.0031, 0.0041),
    stringsAsFactors = FALSE
  )

  result <- fn(single_df, joint_df, targets, is_bin)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2L)
  expect_true("isBinary" %in% names(result))
  expect_equal(result$isBinary, c(TRUE, FALSE))
  expect_true("tauStarSingle" %in% names(result))
  expect_true("tauStarJoint" %in% names(result))
  expect_equal(result$tauStarSingle, c(0.01, 0.02))
  expect_equal(result$tauStarJoint, c(0.011, 0.021))
})

test_that(".sldscAssembleTraitSummary handles NULL single", {
  fn <- pecotmr:::.sldscAssembleTraitSummary
  joint_df <- data.frame(target = "A_0", tauStar = 0.01, tauStarSe = 0.003,
                          stringsAsFactors = FALSE)
  is_bin <- c(A_0 = TRUE)
  result <- fn(NULL, joint_df, "A_0", is_bin)
  expect_equal(nrow(result), 1L)
  expect_true(all(is.na(result$tauStarSingle)))
  expect_equal(result$tauStarJoint, 0.01)
})

test_that(".sldscAssembleTraitSummary handles NULL joint", {
  fn <- pecotmr:::.sldscAssembleTraitSummary
  single_df <- data.frame(target = "A_0", tauStar = 0.01, tauStarSe = 0.003,
                           enrichment = 2.0, enrichmentSe = 0.4,
                           enrichmentP = 0.01, enrichstat = 0.001,
                           enrichstatSe = 0.0003,
                           stringsAsFactors = FALSE)
  is_bin <- c(A_0 = TRUE)
  result <- fn(single_df, NULL, "A_0", is_bin)
  expect_equal(result$tauStarSingle, 0.01)
  expect_true(all(is.na(result$tauStarJoint)))
})

test_that(".sldscAssembleTraitSummary handles both NULL", {
  fn <- pecotmr:::.sldscAssembleTraitSummary
  is_bin <- c(A_0 = TRUE)
  result <- fn(NULL, NULL, "A_0", is_bin)
  expect_equal(nrow(result), 1L)
  expect_equal(result$target, "A_0")
})


# =============================================================================
# .sldscViewForMeta
# =============================================================================

test_that(".sldscViewForMeta extracts single-mode columns", {
  fn <- pecotmr:::.sldscViewForMeta
  per_trait <- list(
    traitX = list(summary = data.frame(
      target = "A_0",
      tauStarSingle = 0.01, tauStarSeSingle = 0.003,
      enrichmentSingle = 2.0, enrichmentSeSingle = 0.4,
      stringsAsFactors = FALSE
    ))
  )
  view <- fn(per_trait, "single")
  expect_true(is.list(view))
  expect_equal(length(view), 1L)
  s <- view$traitX$summary
  expect_true("tauStar" %in% names(s))
  expect_true("tauStarSe" %in% names(s))
  expect_equal(s$tauStar, 0.01)
})

test_that(".sldscViewForMeta returns NULL for missing summary", {
  fn <- pecotmr:::.sldscViewForMeta
  per_trait <- list(traitX = list(summary = NULL))
  view <- fn(per_trait, "single")
  expect_null(view$traitX)
})

test_that(".sldscViewForMeta returns NULL when no matching columns", {
  fn <- pecotmr:::.sldscViewForMeta
  per_trait <- list(
    traitX = list(summary = data.frame(target = "A_0", other_col = 1,
                                        stringsAsFactors = FALSE))
  )
  view <- fn(per_trait, "single")
  expect_null(view$traitX)
})


# =============================================================================
# sldscPostprocessingPipeline (integration)
# =============================================================================

test_that("sldscPostprocessingPipeline runs end-to-end", {
  fix <- .make_sldsc_fixtures()

  trait_single_prefixes <- list(
    traitX = c(
      file.path(fix$out_dir, "traitX_single_annot_A"),
      file.path(fix$out_dir, "traitX_single_annot_B")
    ),
    traitY = c(
      file.path(fix$out_dir, "traitY_single_annot_A"),
      file.path(fix$out_dir, "traitY_single_annot_B")
    )
  )
  trait_joint_prefix <- c(
    traitX = file.path(fix$out_dir, "traitX_joint"),
    traitY = file.path(fix$out_dir, "traitY_joint")
  )

  result <- suppressMessages(sldscPostprocessingPipeline(
    traitSinglePrefixes =trait_single_prefixes,
    traitJointPrefix    = trait_joint_prefix,
    targetAnnoDir       =fix$anno_dir,
    frqfileDir           =fix$frq_dir,
    plinkName            =fix$plinkName,
    mafCutoff            =0.05
  ))

  expect_true(is.list(result))
  expect_named(result, c("per_trait", "meta", "params"))

  # per_trait
  expect_equal(length(result$per_trait), 2L)
  expect_named(result$per_trait, c("traitX", "traitY"))

  pt <- result$per_trait$traitX
  expect_true(is.data.frame(pt$summary))
  expect_true("isBinary" %in% names(pt$summary))
  expect_true(is.numeric(pt$h2g))

  # meta
  expect_named(result$meta, c("tauStar", "enrichment", "enrichstat"))
  for (nm in names(result$meta)) {
    m <- result$meta[[nm]]
    expect_true(is.data.frame(m))
    expect_true("target" %in% names(m))
    expect_true("isBinary" %in% names(m))
  }

  # params
  expect_equal(result$params$M_ref > 0, TRUE)
  expect_equal(result$params$maf_cutoff, 0.05)
  expect_equal(length(result$params$target_categories), 2L)
  expect_equal(result$params$trait_names, c("traitX", "traitY"))
})

test_that("sldscPostprocessingPipeline works without joint runs", {
  fix <- .make_sldsc_fixtures()

  trait_single_prefixes <- list(
    traitX = c(
      file.path(fix$out_dir, "traitX_single_annot_A"),
      file.path(fix$out_dir, "traitX_single_annot_B")
    ),
    traitY = c(
      file.path(fix$out_dir, "traitY_single_annot_A"),
      file.path(fix$out_dir, "traitY_single_annot_B")
    )
  )

  # No joint prefix
  result <- suppressMessages(sldscPostprocessingPipeline(
    traitSinglePrefixes =trait_single_prefixes,
    traitJointPrefix    = NULL,
    targetAnnoDir       =fix$anno_dir,
    frqfileDir           =fix$frq_dir,
    plinkName            =fix$plinkName,
    mafCutoff            =0.05
  ))

  expect_true(is.list(result))
  expect_true(all(is.na(result$meta$tauStar$jointMean)))
})

test_that("sldscPostprocessingPipeline applies target_labels", {
  fix <- .make_sldsc_fixtures()

  trait_single_prefixes <- list(
    traitX = c(
      file.path(fix$out_dir, "traitX_single_annot_A"),
      file.path(fix$out_dir, "traitX_single_annot_B")
    ),
    traitY = c(
      file.path(fix$out_dir, "traitY_single_annot_A"),
      file.path(fix$out_dir, "traitY_single_annot_B")
    )
  )
  trait_joint_prefix <- c(
    traitX = file.path(fix$out_dir, "traitX_joint"),
    traitY = file.path(fix$out_dir, "traitY_joint")
  )

  result <- suppressMessages(sldscPostprocessingPipeline(
    traitSinglePrefixes =trait_single_prefixes,
    traitJointPrefix    = trait_joint_prefix,
    targetAnnoDir       =fix$anno_dir,
    frqfileDir           =fix$frq_dir,
    plinkName            =fix$plinkName,
    mafCutoff            =0.05,
    targetLabels         =c("Pretty_A", "Pretty_B")
  ))

  # Check relabeling
  expect_equal(result$params$target_categories, c("Pretty_A", "Pretty_B"))
  expect_true(!is.null(result$params$target_categories_orig))
  expect_true(all(result$meta$tauStar$target %in% c("Pretty_A", "Pretty_B")))
  expect_true(all(result$per_trait$traitX$summary$target %in% c("Pretty_A", "Pretty_B")))
})

test_that("sldscPostprocessingPipeline errors on wrong target_labels length", {
  fix <- .make_sldsc_fixtures()

  trait_single_prefixes <- list(
    traitX = c(
      file.path(fix$out_dir, "traitX_single_annot_A"),
      file.path(fix$out_dir, "traitX_single_annot_B")
    )
  )
  trait_joint_prefix <- c(
    traitX = file.path(fix$out_dir, "traitX_joint")
  )

  expect_error(
    suppressMessages(sldscPostprocessingPipeline(
      traitSinglePrefixes =trait_single_prefixes,
      traitJointPrefix    = trait_joint_prefix,
      targetAnnoDir       =fix$anno_dir,
      frqfileDir           =fix$frq_dir,
      plinkName            =fix$plinkName,
      targetLabels         =c("only_one")
    )),
    "targetLabels"
  )
})

test_that("sldscPostprocessingPipeline errors on unnamed prefixes", {
  expect_error(
    sldscPostprocessingPipeline(
      traitSinglePrefixes =list(c("a", "b")),
      traitJointPrefix    = NULL,
      targetAnnoDir       ="."
    ),
    "named list"
  )
})
