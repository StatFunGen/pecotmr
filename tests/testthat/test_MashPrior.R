# Tests for R/MashPrior.R — the data-driven (mash) prior bundle: full-data prior
# (fullFit) + per-fold cross-validated priors (cvFits) consumed by
# twasWeightsPipeline (mr.mash) and, transitively, fineMappingPipeline (mvSuSiE).

test_that("MashPrior: construct + accessors (full + cv)", {
  U  <- list(U1 = diag(2))
  sp <- data.frame(Sample = paste0("s", 1:6), Fold = rep(1:3, each = 2),
                   stringsAsFactors = FALSE)
  perFold <- list(list(U = U, w = c(0.5, 0.5)),
                  list(U = U, w = c(0.5, 0.5)),
                  list(U = U, w = c(0.5, 0.5)))
  mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)),
                  cvFits = list(samplePartition = sp, perFoldFits = perFold))
  expect_s4_class(mp, "MashPrior")
  expect_identical(getFullFit(mp)$U, U)
  expect_length(getCvFits(mp)$perFoldFits, 3L)
  expect_identical(getCvFits(mp)$samplePartition, sp)
  expect_output(show(mp), "MashPrior")
})

test_that("MashPrior: full-only and cv-only bundles", {
  U <- list(U1 = diag(2))
  perFold <- list(list(U = U), list(U = U), list(U = U))
  mpFull <- MashPrior(fullFit = list(U = U))
  expect_null(getCvFits(mpFull))
  mpCv <- MashPrior(cvFits = list(perFoldFits = perFold))
  expect_null(getFullFit(mpCv))
  expect_length(getCvFits(mpCv)$perFoldFits, 3L)
})

test_that("MashPrior: validity rejects malformed bundles", {
  U <- list(U1 = diag(2))
  sp <- data.frame(Sample = paste0("s", 1:6), Fold = rep(1:3, each = 2),
                   stringsAsFactors = FALSE)
  # Empty bundle.
  expect_error(MashPrior(), "at least one")
  # perFoldFits count must match the partition's fold count.
  expect_error(
    MashPrior(cvFits = list(samplePartition = sp,
                            perFoldFits = list(list(U = U)))),
    "fold")
  # perFoldFits must be a non-empty list.
  expect_error(MashPrior(cvFits = list(perFoldFits = "nope")))
  # samplePartition must carry Sample + Fold columns.
  expect_error(
    MashPrior(cvFits = list(samplePartition = data.frame(x = 1),
                            perFoldFits = list(list(U = U)))),
    "Sample")
})

test_that("MashPrior: show reports 'cvFits: none' for a full-only bundle", {
  U  <- list(U1 = diag(2))
  mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)))  # cvFits NULL
  out <- capture.output(show(mp))
  expect_true(any(grepl("cvFits: none", out)))
  expect_true(any(grepl("fullFit: present", out)))
})
