# Tests for R/sldscPostprocessingPipeline.R
# The pipeline is pure computation over an in-memory SldscData: no file I/O,
# no mocks. Fixtures come from helper-sldsc.R (.sldscMkData / .sldscMkRun).

test_that("pipeline runs end-to-end on a single + joint SldscData", {
  sd <- .sldscMkData()
  res <- suppressMessages(sldscPostprocessingPipeline(sd, mafCutoff = 0.05))

  expect_named(res, c("per_trait", "meta", "params"))
  expect_named(res$per_trait, c("traitX", "traitY"))
  expect_named(res$meta, c("tauStar", "enrichment", "enrichstat"))
  expect_true(is.data.frame(res$per_trait$traitX$summary))
  expect_true("isBinary" %in% names(res$per_trait$traitX$summary))

  # meta tauStar carries both single and joint channels
  expect_true(all(c("singleMean", "singleSe", "singleP",
                    "jointMean", "jointSe", "jointP") %in% names(res$meta$tauStar)))
  expect_equal(nrow(res$meta$tauStar), 2L)

  expect_gt(res$params$M_ref, 0)
  expect_equal(res$params$maf_cutoff, 0.05)
  expect_equal(res$params$target_categories, c("annot_A_0", "annot_B_0"))
  expect_equal(res$params$trait_names, c("traitX", "traitY"))
  # baseline annotation detected from the joint run
  expect_true("baselineLD_0" %in% res$params$baseline_categories)
})

test_that("pipeline without joint runs yields NA joint meta", {
  sd <- .sldscMkData(withJoint = FALSE)
  res <- suppressMessages(sldscPostprocessingPipeline(
    sd, targetCategories = c("annot_A_0", "annot_B_0")))
  expect_true(all(is.na(res$meta$tauStar$jointMean)))
  expect_equal(res$params$n_baseline, 0L)
})

test_that("pipeline applies targetLabels", {
  sd <- .sldscMkData()
  res <- suppressMessages(sldscPostprocessingPipeline(
    sd, targetLabels = c("Pretty_A", "Pretty_B")))
  expect_equal(res$params$target_categories, c("Pretty_A", "Pretty_B"))
  expect_false(is.null(res$params$target_categories_orig))
  expect_setequal(res$meta$tauStar$target, c("Pretty_A", "Pretty_B"))
  expect_setequal(res$per_trait$traitX$summary$target, c("Pretty_A", "Pretty_B"))
})

test_that("pipeline errors on wrong targetLabels length", {
  sd <- .sldscMkData()
  expect_error(
    suppressMessages(sldscPostprocessingPipeline(sd, targetLabels = c("only_one"))),
    "targetLabels")
})

test_that("pipeline errors on non-SldscData input", {
  expect_error(sldscPostprocessingPipeline(list(a = 1)),
               "must be an SldscData object")
})

test_that("pipeline errors when the SldscData has no traits", {
  sd <- SldscData(annot = data.frame(CHR = 1, SNP = "rs1", annot_A = 1))
  expect_error(suppressMessages(sldscPostprocessingPipeline(sd)),
               "no traits")
})

test_that("pipeline takes an explicit targetCategories (skips auto-detect)", {
  sd <- .sldscMkData()
  res <- suppressMessages(sldscPostprocessingPipeline(
    sd, targetCategories = c("annot_A_0")))
  expect_equal(res$params$target_categories, "annot_A_0")
  expect_equal(nrow(res$meta$tauStar), 1L)
})

test_that("pipeline falls back to positional rename when names don't match", {
  # annot columns -> annot_A_0/annot_B_0, but the runs' categories use the
  # polyfun --snp-list "L2" naming, so the intersect is empty and the pipeline
  # renames the first length(sdAnnot) .results rows positionally.
  annot <- data.frame(CHR = c(1, 1, 2, 2), SNP = paste0("rs", 1:4),
                      annot_A = c(1, 0, 1, 0), annot_B = c(2.1, 1.9, 2.4, 2.0),
                      stringsAsFactors = FALSE)
  frq <- data.frame(CHR = c(1, 1, 2, 2), SNP = paste0("rs", 1:4),
                    MAF = rep(0.2, 4), stringsAsFactors = FALSE)
  mkTrait <- function() list(
    single = list(.sldscMkRun(c("L2_1", "base1")), .sldscMkRun(c("L2_2", "base1"))),
    joint  = .sldscMkRun(c("L2_1", "L2_2", "base1")))
  sd <- SldscData(annot, frq, list(traitX = mkTrait(), traitY = mkTrait()))
  res <- suppressMessages(sldscPostprocessingPipeline(sd))
  # target categories were renamed positionally to the first 2 .results rows
  expect_equal(res$params$target_categories, c("L2_1", "L2_2"))
})

test_that("pipeline breaks the single loop when a trait has fewer runs than targets", {
  # traitX has only ONE single run but there are TWO target categories: the
  # i > length(singleRuns) break is hit on the 2nd target.
  annot <- data.frame(CHR = c(1, 1, 2, 2), SNP = paste0("rs", 1:4),
                      annot_A = c(1, 0, 1, 0), annot_B = c(2.1, 1.9, 2.4, 2.0),
                      stringsAsFactors = FALSE)
  frq <- data.frame(CHR = c(1, 1, 2, 2), SNP = paste0("rs", 1:4),
                    MAF = rep(0.2, 4), stringsAsFactors = FALSE)
  traits <- list(traitX = list(
    single = list(.sldscMkRun(c("annot_A_0", "baselineLD_0"))),  # only 1 run
    joint  = .sldscMkRun(c("annot_A_0", "annot_B_0", "baselineLD_0"))))
  sd <- SldscData(annot, frq, traits)
  res <- suppressMessages(sldscPostprocessingPipeline(
    sd, targetCategories = c("annot_A_0", "annot_B_0")))
  # The break stops before annot_B_0's single run; the single-keyed summary
  # therefore carries only annot_A_0.
  sm <- res$per_trait$traitX$summary
  expect_equal(sm$target, "annot_A_0")
  expect_false(is.na(sm$tauStarSingle))
})

test_that("pipeline warns and skips a single run that fails to standardize", {
  # traitX's 2nd single run lacks the target category annot_B_0, so
  # standardizeSldscTrait errors and the pipeline's tryCatch warns + skips it.
  annot <- data.frame(CHR = c(1, 1, 2, 2), SNP = paste0("rs", 1:4),
                      annot_A = c(1, 0, 1, 0), annot_B = c(2.1, 1.9, 2.4, 2.0),
                      stringsAsFactors = FALSE)
  frq <- data.frame(CHR = c(1, 1, 2, 2), SNP = paste0("rs", 1:4),
                    MAF = rep(0.2, 4), stringsAsFactors = FALSE)
  traits <- list(traitX = list(
    single = list(.sldscMkRun(c("annot_A_0", "baselineLD_0")),
                  .sldscMkRun(c("WRONG_0", "baselineLD_0"))),  # missing annot_B_0
    joint  = .sldscMkRun(c("annot_A_0", "annot_B_0", "baselineLD_0"))))
  sd <- SldscData(annot, frq, traits)
  expect_warning(
    res <- suppressMessages(sldscPostprocessingPipeline(
      sd, targetCategories = c("annot_A_0", "annot_B_0"))),
    "Failed to standardize single")
  # run1 (annot_A_0) succeeded; run2 (annot_B_0) failed to standardize and was
  # skipped, so the single-keyed summary carries only annot_A_0.
  sm <- res$per_trait$traitX$summary
  expect_equal(sm$target, "annot_A_0")
})

# annotation + frq with enough per-chromosome variance for the compute steps
# to succeed before the branch under test is reached.
.sldscBranchAnnotFrq <- function() {
  list(
    annot = data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
                       annot_A = c(1, 0, 1, 0, 1, 0),
                       annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0),
                       stringsAsFactors = FALSE),
    frq = data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
                     MAF = rep(0.2, 6), stringsAsFactors = FALSE))
}

test_that("pipeline errors when targetCategories cannot be auto-detected", {
  af <- .sldscBranchAnnotFrq()
  # The single trait has no joint run and an empty single list, so there is no
  # pivot run to auto-detect categories from.
  sd <- SldscData(af$annot, af$frq,
                  list(traitX = list(single = list(), joint = NULL)))
  expect_error(suppressMessages(sldscPostprocessingPipeline(sd)),
               "cannot auto-detect")
})

test_that("pipeline warns and skips a joint run that fails to standardize", {
  af <- .sldscBranchAnnotFrq()
  # The joint run lacks the target category annot_A_0, so its standardization
  # errors and the pipeline's joint-side tryCatch warns + skips it.
  traits <- list(traitX = list(
    single = list(.sldscMkRun(c("annot_A_0", "baselineLD_0"))),
    joint  = .sldscMkRun(c("WRONG_0", "baselineLD_0"))))
  sd <- SldscData(af$annot, af$frq, traits)
  expect_warning(
    res <- suppressMessages(sldscPostprocessingPipeline(
      sd, targetCategories = c("annot_A_0"))),
    "Failed to standardize joint")
  # the single side still produced an estimate
  expect_false(is.na(res$per_trait$traitX$summary$tauStarSingle[1]))
})
