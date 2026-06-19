context("twasWeights internal helpers (extra)")

# ===========================================================================
# .twasFineMappingFits
# ===========================================================================

.tw_makeFmEntry <- function(method_tag = "susie", n = 3) {
  FineMappingEntry(
    variantIds = paste0("v", seq_len(n)),
    trimmedFit = list(payload = method_tag),
    topLoci    = data.frame(variant_id = paste0("v", seq_len(n)),
                            pip = seq(0.9, by = -0.1, length.out = n),
                            stringsAsFactors = FALSE))
}

test_that(".twasFineMappingFits: NULL input returns an empty list", {
  out <- pecotmr:::.twasFineMappingFits(NULL,
                                        study = "s1", context = "c1",
                                        trait = "t1")
  expect_equal(out, list())
})

test_that(".twasFineMappingFits: non-FineMappingResult input errors", {
  expect_error(
    pecotmr:::.twasFineMappingFits("not_a_result",
                                    study = "s1", context = "c1", trait = "t1"),
    "must be a FineMappingResult or NULL"
  )
})

test_that(".twasFineMappingFits: pulls matching (study,context,trait) susie fit", {
  e_susie    <- .tw_makeFmEntry("susie")
  e_susieInf <- .tw_makeFmEntry("susieInf")
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susieInf"),
    entry   = list(e_susie, e_susieInf))
  out <- pecotmr:::.twasFineMappingFits(res,
                                         study = "s1", context = "c1",
                                         trait = "t1")
  expect_true("susie" %in% names(out))
  expect_true("susieInf" %in% names(out))
  expect_equal(out$susie$payload, "susie")
  expect_equal(out$susieInf$payload, "susieInf")
})

test_that(".twasFineMappingFits: returns empty list when no tuple matches", {
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(.tw_makeFmEntry("susie")))
  out <- pecotmr:::.twasFineMappingFits(res,
                                         study = "ghost", context = "c1",
                                         trait = "t1")
  expect_equal(out, list())
})

test_that(".twasFineMappingFits: ignores non-susie methods (e.g. lasso)", {
  e <- .tw_makeFmEntry("susie")
  e_lasso <- .tw_makeFmEntry("lasso")
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("susie", "lasso"),
    entry   = list(e, e_lasso))
  out <- pecotmr:::.twasFineMappingFits(res,
                                         study = "s1", context = "c1",
                                         trait = "t1")
  # Only susie should be picked up; lasso is not in the susie/susieInf/susieAsh
  # canonical menu.
  expect_equal(names(out), "susie")
})

# ===========================================================================
# .twasLdFromSketch (mocked extractBlockGenotypes)
# ===========================================================================

.tw_makeSketchHandle <- function(snp_n = 6L, n_samples = 30L) {
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

.tw_mockExtractor <- function(seed = 7, n_samples = 30L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds,
                                    handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges   = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx],
                                  width = 1L))
    S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
      SNP = handle@snpInfo$SNP[snpIdx],
      A1  = handle@snpInfo$A1[snpIdx],
      A2  = handle@snpInfo$A2[snpIdx])
    cd <- S4Vectors::DataFrame(sampleId = handle@sampleIds,
                               row.names = handle@sampleIds)
    dosage <- t(sub)
    rownames(dosage) <- handle@snpInfo$SNP[snpIdx]
    colnames(dosage) <- handle@sampleIds
    SummarizedExperiment::SummarizedExperiment(
      assays    = list(dosage = dosage),
      rowRanges = rr,
      colData   = cd)
  }
}

test_that(".twasLdFromSketch: rejects non-GenotypeHandle ldSketch", {
  expect_error(
    pecotmr:::.twasLdFromSketch("not_a_handle", c("v1", "v2")),
    "ldSketch must be a GenotypeHandle"
  )
})

test_that(".twasLdFromSketch: variants not in the panel error", {
  h <- .tw_makeSketchHandle()
  expect_error(
    pecotmr:::.twasLdFromSketch(h, c("v1", "ghost")),
    "variant id.*not present in the LD sketch"
  )
})

test_that(".twasLdFromSketch: returns a square LD matrix named by variantIds", {
  h <- .tw_makeSketchHandle()
  local_mocked_bindings(extractBlockGenotypes = .tw_mockExtractor(),
                        .package = "pecotmr")
  ids <- c("v2", "v4", "v5")
  R <- pecotmr:::.twasLdFromSketch(h, ids)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(rownames(R), ids)
  expect_equal(colnames(R), ids)
  # Symmetric, diagonal == 1 (sample correlation).
  expect_equal(unname(diag(R)), c(1, 1, 1), tolerance = 1e-12)
  expect_equal(R, t(R), tolerance = 1e-12)
})

# ===========================================================================
# .twasWeightsPipelineMatrix: susieFit pre-fit pass-through
# ===========================================================================

test_that(".twasWeightsPipelineMatrix: susieFit pre-fit is recorded in res", {
  set.seed(0)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:p)))
  y <- as.numeric(X %*% c(1.0, -0.5, 0, 0, 0) + rnorm(n, sd = 0.2))

  # Build a stub susie fit shape; the pipeline records its intermediate.
  fake_susie <- list(
    alpha = matrix(1/p, nrow = 2, ncol = p),
    mu    = matrix(0,    nrow = 2, ncol = p),
    X_column_scale_factors = rep(1, p),
    pip   = rep(0.1, p))

  # The intermediate-recording branch keys on snake_case `susie_weights`.
  res <- suppressMessages(
    pecotmr:::.twasWeightsPipelineMatrix(
      X = X, y = y,
      susieFit = fake_susie,
      cvFolds = 0,
      weightMethods = list(susie_weights = list()),
      estimatePi = FALSE,
      verbose = 0))
  expect_true("susieWeightsIntermediate" %in% names(res))
  expect_true("twasWeights" %in% names(res))
})

# ===========================================================================
# .twasWeightsPipelineMatrix: empirical pi path via mr.ash (mocked)
# ===========================================================================

test_that(".twasWeightsPipelineMatrix: empirical pi from mr.ash gets propagated", {
  set.seed(1)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:p)))
  y <- as.numeric(X %*% c(0.5, 0, 0, 0, 0) + rnorm(n, sd = 0.2))

  # Mock mrashWeights to return a fake matrix carrying a fit$pi attribute.
  local_mocked_bindings(
    mrashWeights = function(X, y, ...) {
      out <- matrix(rep(0.05, ncol(X)), ncol = 1)
      attr(out, "fit") <- list(pi = c(0.8, 0.1, 0.1))
      rownames(out) <- colnames(X)
      out
    },
    bayesCWeights = function(X, y, pi, ...) {
      # Capture the pi the pipeline injected.
      out <- matrix(pi, nrow = ncol(X), ncol = 1)
      rownames(out) <- colnames(X)
      out
    },
    .package = "pecotmr"
  )

  res <- suppressMessages(
    pecotmr:::.twasWeightsPipelineMatrix(
      X = X, y = y,
      cvFolds = 0,
      weightMethods = list(mrash_weights = list(),
                           bayes_c_weights = list()),
      estimatePi = TRUE,
      verbose = 0))
  expect_true("empiricalPi" %in% names(res))
  expect_equal(as.numeric(res$empiricalPi), 1 - 0.8, tolerance = 1e-12)
})
