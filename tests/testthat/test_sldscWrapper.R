# Tests for R/sldscWrapper.R
#
# File-based fixture builders (.make_annot_gz / .make_frq / .make_polyfun_single
# / .make_polyfun_joint / .make_sldsc_fixtures) and the in-memory SldscData
# builders (.sldscMkRun / .sldscMkData) live in helper-sldsc.R.
#
# The reader functions (readSldscTrait/readSldscAnnot/readSldscFrq) do file I/O
# and are tested against real fixture files. The compute functions
# (computeSldscMRef/computeSldscAnnotSd/isBinarySldscAnnot/standardizeSldscTrait)
# operate on an in-memory SldscData and are tested with in-memory fixtures.

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
  writeLines("No heritability here", paste0(prefix, ".log"))
  expect_error(readSldscTrait(prefix), "Total Observed scale h2")
})

test_that("readSldscTrait errors on column mismatch in part_delete", {
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "bad_delete")
  .make_polyfun_single(dir, prefix, "a", nBlocks = 3)
  vroom::vroom_write(data.frame(x = 1:3, y = 4:6, z = 7:9, w = 10:12),
                     paste0(prefix, ".part_delete"), delim = "\t")
  expect_error(readSldscTrait(prefix), "part_delete")
})

test_that("readSldscTrait errors when the h2g value is non-numeric", {
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "bad_h2")
  .make_polyfun_single(dir, prefix, "a", nBlocks = 3)
  writeLines(c("start", "Total Observed scale h2: abc (0.05)", "end"),
             paste0(prefix, ".log"))
  expect_error(readSldscTrait(prefix), "failed to parse h2g numeric")
})


# =============================================================================
# readSldscAnnot
# =============================================================================

test_that("readSldscAnnot stacks per-chromosome .annot.gz into one table", {
  dir <- withr::local_tempdir()
  for (chr in 1:2) .make_annot_gz(dir, chr)
  df <- readSldscAnnot(dir)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("CHR", "SNP", "annot_A", "annot_B") %in% names(df)))
  expect_equal(nrow(df), 100L)             # 2 chroms x 50 SNPs
  expect_setequal(unique(df$CHR), c(1, 2))
})

test_that("readSldscAnnot respects annotCols", {
  dir <- withr::local_tempdir()
  for (chr in 1:2) .make_annot_gz(dir, chr)
  df <- readSldscAnnot(dir, annotCols = "annot_A")
  expect_true("annot_A" %in% names(df))
  expect_false("annot_B" %in% names(df))
})

test_that("readSldscAnnot errors on missing dir / no files", {
  expect_error(readSldscAnnot("/nonexistent/dir"), "does not exist")
  empty <- withr::local_tempdir()
  expect_error(readSldscAnnot(empty), "no .annot.gz")
})

test_that("readSldscAnnot errors when annotCols resolves to nothing", {
  dir <- withr::local_tempdir()
  .make_annot_gz(dir, 1)
  expect_error(readSldscAnnot(dir, annotCols = character(0)),
               "no annotation columns")
})


# =============================================================================
# readSldscFrq
# =============================================================================

test_that("readSldscFrq stacks per-chromosome .frq into one table", {
  dir <- withr::local_tempdir()
  for (chr in 1:2) .make_frq(dir, chr, plinkName = "ref_chr")
  df <- readSldscFrq(dir, plinkName = "ref_chr")
  expect_s3_class(df, "data.frame")
  expect_true(all(c("CHR", "SNP", "MAF") %in% names(df)))
  expect_equal(nrow(df), 100L)
})

test_that("readSldscFrq falls back to a generic .frq glob when the prefix misses", {
  dir <- withr::local_tempdir()
  for (chr in 1:2) .make_frq(dir, chr, plinkName = "other_chr")
  df <- readSldscFrq(dir, plinkName = "nomatch_chr")
  expect_equal(nrow(df), 100L)
})

test_that("readSldscFrq errors on missing dir / no files", {
  expect_error(readSldscFrq("/nonexistent/dir"), "does not exist")
  empty <- withr::local_tempdir()
  expect_error(readSldscFrq(empty), "no .frq")
})


# =============================================================================
# computeSldscAnnotSd  (operates on SldscData)
# =============================================================================

test_that("computeSldscAnnotSd computes SDs with MAF filtering", {
  sds <- computeSldscAnnotSd(.sldscMkData(), mafCutoff = 0.05)
  expect_true(is.numeric(sds))
  expect_equal(length(sds), 2L)
  expect_named(sds, c("annot_A", "annot_B"))
  expect_true(all(sds > 0))
})

test_that("computeSldscAnnotSd works with mafCutoff = 0 (no frq needed)", {
  sds <- computeSldscAnnotSd(.sldscMkData(withFrq = FALSE), mafCutoff = 0)
  expect_true(all(sds > 0))
})

test_that("computeSldscAnnotSd respects annotCols (character)", {
  sds <- computeSldscAnnotSd(.sldscMkData(), mafCutoff = 0, annotCols = "annot_A")
  expect_equal(length(sds), 1L)
  expect_named(sds, "annot_A")
})

test_that("computeSldscAnnotSd respects annotCols (numeric)", {
  sds <- computeSldscAnnotSd(.sldscMkData(), mafCutoff = 0, annotCols = 2L)
  expect_equal(length(sds), 1L)
  expect_named(sds, "annot_B")
})

test_that("computeSldscAnnotSd errors when mafCutoff > 0 but no frq data", {
  expect_error(
    computeSldscAnnotSd(.sldscMkData(withFrq = FALSE), mafCutoff = 0.05),
    "requires frq data")
})

test_that("computeSldscAnnotSd errors when there are no annotation columns", {
  expect_error(
    computeSldscAnnotSd(.sldscMkData(), mafCutoff = 0, annotCols = character(0)),
    "no annotation columns")
})

test_that("computeSldscAnnotSd errors on non-SldscData input", {
  expect_error(computeSldscAnnotSd(list(a = 1)), "must be an SldscData")
})

test_that("computeSldscAnnotSd errors with zero degrees of freedom", {
  # One SNP per chromosome: after the per-chromosome split each block has
  # nrow <= 1 (the `next`), so no variance accumulates and den stays 0.
  annot <- data.frame(CHR = c(1, 2), SNP = c("rs1", "rs2"),
                      annot_A = c(1, 0), annot_B = c(2.1, 1.9),
                      stringsAsFactors = FALSE)
  frq <- data.frame(CHR = c(1, 2), SNP = c("rs1", "rs2"), MAF = c(0.2, 0.2),
                    stringsAsFactors = FALSE)
  sd <- SldscData(annot, frq, list())
  expect_error(computeSldscAnnotSd(sd, mafCutoff = 0.05),
               "zero degrees of freedom")
})


# =============================================================================
# computeSldscMRef  (operates on SldscData)
# =============================================================================

test_that("computeSldscMRef counts MAF > cutoff SNPs from the frq table", {
  M <- computeSldscMRef(.sldscMkData(), mafCutoff = 0.05)
  expect_true(is.integer(M))
  expect_equal(M, 6L)        # all 6 frq SNPs have MAF 0.2 > 0.05
})

test_that("computeSldscMRef counts all SNPs when mafCutoff = 0", {
  expect_equal(computeSldscMRef(.sldscMkData(), mafCutoff = 0), 6L)
})

test_that("computeSldscMRef applies the MAF cutoff", {
  annot <- data.frame(CHR = 1, SNP = paste0("rs", 1:4), a = c(1, 0, 1, 0),
                      stringsAsFactors = FALSE)
  frq <- data.frame(CHR = 1, SNP = paste0("rs", 1:4),
                    MAF = c(0.2, 0.01, 0.3, 0.02), stringsAsFactors = FALSE)
  sd <- SldscData(annot, frq, list())
  expect_equal(computeSldscMRef(sd, mafCutoff = 0.05), 2L)   # rs1, rs3
})

test_that("computeSldscMRef falls back to annot row count when frq absent and mafCutoff = 0", {
  expect_equal(computeSldscMRef(.sldscMkData(withFrq = FALSE), mafCutoff = 0), 6L)
})

test_that("computeSldscMRef errors when mafCutoff > 0 but no frq data", {
  expect_error(computeSldscMRef(.sldscMkData(withFrq = FALSE), mafCutoff = 0.05),
               "requires frq data")
})

test_that("computeSldscMRef errors on non-SldscData input", {
  expect_error(computeSldscMRef(list(a = 1)), "must be an SldscData")
})


# =============================================================================
# isBinarySldscAnnot  (operates on SldscData)
# =============================================================================

test_that("isBinarySldscAnnot detects binary and continuous annotations", {
  result <- isBinarySldscAnnot(.sldscMkData())
  expect_true(is.logical(result))
  expect_named(result, c("annot_A", "annot_B"))
  expect_true(result[["annot_A"]])    # binary (0/1)
  expect_false(result[["annot_B"]])   # continuous
})

test_that("isBinarySldscAnnot respects annotCols (character)", {
  result <- isBinarySldscAnnot(.sldscMkData(), annotCols = "annot_A")
  expect_equal(length(result), 1L)
  expect_true(result[["annot_A"]])
})

test_that("isBinarySldscAnnot respects annotCols (numeric)", {
  result <- isBinarySldscAnnot(.sldscMkData(), annotCols = 2L)
  expect_equal(length(result), 1L)
  expect_named(result, "annot_B")
})

test_that("isBinarySldscAnnot errors on non-SldscData input", {
  expect_error(isBinarySldscAnnot(list(a = 1)), "must be an SldscData")
})


# =============================================================================
# standardizeSldscTrait  (operates on SldscData via getTraitRun)
# =============================================================================

# Build a readSldscTrait-shaped run with specific values the tests assert on.
.make_trait_data <- function(cats = c("A_0", "B_0"), nBlocks = 10, h2g = 0.3) {
  n <- length(cats)
  taus <- rep(1e-7, n)
  blocks <- matrix(rnorm(nBlocks * n, mean = rep(taus, each = nBlocks), sd = 1e-8),
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
    nBlocks       = nBlocks)
}

# Wrap a run (and optional joint run) in a minimal SldscData so the new
# standardizeSldscTrait(sldscData, trait, mode, idx, ...) API can reach it.
.wrapRun <- function(run, joint = run) {
  SldscData(annot = data.frame(CHR = 1, SNP = "rs1", A = 1, stringsAsFactors = FALSE),
            traits = list(t = list(single = list(run), joint = joint)))
}

test_that("standardizeSldscTrait works in single mode", {
  td <- .make_trait_data()
  result <- standardizeSldscTrait(.wrapRun(td), "t", mode = "single", idx = 1L,
                                  sdAnnot = c(A_0 = 0.5, B_0 = 1.2), MRef = 1000L)
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
  expect_equal(s$tauStar, unname(td$tau) * c(0.5, 1.2) * 1000 / 0.3)
  expect_true(is.matrix(result$tau_star_blocks))
  expect_equal(dim(result$tau_star_blocks), c(10L, 2L))
})

test_that("standardizeSldscTrait works in joint mode", {
  td <- .make_trait_data()
  result <- standardizeSldscTrait(.wrapRun(td), "t", mode = "joint",
                                  sdAnnot = c(A_0 = 0.5, B_0 = 1.2), MRef = 1000L)
  expect_equal(result$mode, "joint")
  expect_false("enrichment" %in% names(result$summary))
  expect_true("tauStar" %in% names(result$summary))
})

test_that("standardizeSldscTrait auto-detects target categories", {
  td <- .make_trait_data()
  result <- standardizeSldscTrait(.wrapRun(td), "t", mode = "joint",
                                  sdAnnot = c(A_0 = 0.5), MRef = 1000L)
  expect_equal(nrow(result$summary), 1L)
  expect_equal(result$summary$target, "A_0")
})

test_that("standardizeSldscTrait errors on empty categories", {
  td <- .make_trait_data()
  expect_error(
    standardizeSldscTrait(.wrapRun(td), "t", mode = "single", idx = 1L,
                          sdAnnot = c(X_0 = 0.5), MRef = 1000L),
    "no target categories")
})

test_that("standardizeSldscTrait errors on missing categories", {
  td <- .make_trait_data(cats = "A_0")
  expect_error(
    standardizeSldscTrait(.wrapRun(td), "t", mode = "single", idx = 1L,
                          sdAnnot = c(A_0 = 0.5, B_0 = 1.0), MRef = 1000L,
                          targetCategories = c("A_0", "B_0")),
    "missing categories")
})

test_that("standardizeSldscTrait warns on zero sd", {
  td <- .make_trait_data(cats = "A_0")
  expect_warning(
    standardizeSldscTrait(.wrapRun(td), "t", mode = "joint",
                          sdAnnot = c(A_0 = 0), MRef = 1000L),
    "zero/NA sd")
})

test_that("standardizeSldscTrait enrichstatSe handles p = 0", {
  td <- .make_trait_data(cats = "A_0")
  td$enrichmentP <- c(A_0 = 0)   # p = 0 -> abs_z = Inf -> SE = NA
  result <- standardizeSldscTrait(.wrapRun(td), "t", mode = "single", idx = 1L,
                                  sdAnnot = c(A_0 = 0.5), MRef = 1000L)
  expect_true(is.na(result$summary$enrichstatSe))
})

test_that("standardizeSldscTrait errors when the requested run is absent", {
  td <- .make_trait_data()
  expect_error(
    standardizeSldscTrait(.wrapRun(td), "t", mode = "single", idx = 5L,
                          sdAnnot = c(A_0 = 0.5), MRef = 1000L),
    "no single run")
})

test_that("standardizeSldscTrait errors on non-SldscData input", {
  expect_error(
    standardizeSldscTrait(list(a = 1), "t", mode = "single", idx = 1L,
                          sdAnnot = c(A_0 = 0.5), MRef = 1000L),
    "must be an SldscData")
})


# =============================================================================
# metaSldscRandom  (unchanged: operates on standardized per-trait estimates)
# =============================================================================

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
  expect_true(is.na(result$mean))
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

test_that("metaSldscRandom skips NULL entries and NULL summaries", {
  valid <- .make_per_trait_meta(nTraits = 2, means = c(1e-5, 2e-5),
                                ses = c(1e-6, 1e-6))
  # list() keeps explicit NULLs (unlike `x$y <- NULL`), so trait3/trait4 model
  # a NULL entry and a NULL-summary entry the loop must skip.
  pt <- c(valid, list(trait3 = NULL, trait4 = list(summary = NULL)))
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_equal(result$nTraits, 2L)
})

test_that("metaSldscRandom skips a trait whose summary lacks the quantity columns", {
  valid <- .make_per_trait_meta(nTraits = 2, means = c(1e-5, 2e-5),
                                ses = c(1e-6, 1e-6))
  # trait3's row matches the category but has no tauStar/tauStarSe columns.
  pt <- c(valid, list(trait3 = list(
    summary = data.frame(target = "A_0", foo = 1, stringsAsFactors = FALSE))))
  result <- metaSldscRandom(pt, "A_0", "tauStar")
  expect_equal(result$nTraits, 2L)
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
  result <- fn(NULL, joint_df, "A_0", c(A_0 = TRUE))
  expect_equal(nrow(result), 1L)
  expect_true(all(is.na(result$tauStarSingle)))
  expect_equal(result$tauStarJoint, 0.01)
})

test_that(".sldscAssembleTraitSummary handles NULL joint", {
  fn <- pecotmr:::.sldscAssembleTraitSummary
  single_df <- data.frame(target = "A_0", tauStar = 0.01, tauStarSe = 0.003,
                           enrichment = 2.0, enrichmentSe = 0.4,
                           enrichmentP = 0.01, enrichstat = 0.001,
                           enrichstatSe = 0.0003, stringsAsFactors = FALSE)
  result <- fn(single_df, NULL, "A_0", c(A_0 = TRUE))
  expect_equal(result$tauStarSingle, 0.01)
  expect_true(all(is.na(result$tauStarJoint)))
})

test_that(".sldscAssembleTraitSummary handles both NULL", {
  fn <- pecotmr:::.sldscAssembleTraitSummary
  result <- fn(NULL, NULL, "A_0", c(A_0 = TRUE))
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
  view <- fn(list(traitX = list(summary = NULL)), "single")
  expect_null(view$traitX)
})

test_that(".sldscViewForMeta returns NULL when no matching columns", {
  fn <- pecotmr:::.sldscViewForMeta
  per_trait <- list(
    traitX = list(summary = data.frame(target = "A_0", other_col = 1,
                                        stringsAsFactors = FALSE)))
  view <- fn(per_trait, "single")
  expect_null(view$traitX)
})
