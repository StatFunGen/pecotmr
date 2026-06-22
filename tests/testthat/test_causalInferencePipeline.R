context("causalInferencePipeline")

# ===========================================================================
# Strategy: mock extractBlockGenotypes so .cipLdFromSketch returns a real
# LD matrix on a small panel. Everything else (twasZ, MR, p-value combine)
# runs for real on the tiny fixture.
# ===========================================================================

.cip_makeHandle <- function(snp_n = 6L, n_samples = 30L,
                            sample_prefix = "s") {
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
    sampleIds = paste0(sample_prefix, seq_len(n_samples)),
    pgenPtr = NULL)
}

.cip_mockExtractor <- function(seed = 7, n_samples = 30L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds, handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges   = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx], width = 1L))
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

.cip_makeGwasSumstats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = c(2.0, -1.5, 0.5, 1.2, -0.8),
    N   = rep(1000L, 5),
    MAF = rep(0.3, 5))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .cip_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.cip_makeTwasWeights <- function(method = "susie",
                                  variant_ids = paste0("v", 1:5)) {
  entry <- TwasWeightsEntry(
    variantIds = variant_ids,
    weights    = c(0.1, 0.05, -0.2, 0.3, 0.0))
  TwasWeights(
    study    = "Q1", context = "c1", trait = "t1", method = method,
    entry    = list(entry),
    ldSketch = .cip_makeHandle())
}

.cip_makeQtlFmr <- function(variant_ids = paste0("v", 1:5)) {
  n <- length(variant_ids)
  tl <- data.frame(
    variant_id     = variant_ids,
    pip            = c(0.9, 0.05, 0.5, 0.8, 0.01),
    # posterior_mean / posterior_sd carry the "fine-mapped causal effect"
    # estimates that getTopLoci surfaces as beta / se in its projected
    # output (the column names downstream MR / TWAS code reads).
    posterior_mean = c(0.2, 0.05, -0.1, 0.3, 0.0),
    posterior_sd   = rep(0.05, n),
    stringsAsFactors = FALSE)
  e <- FineMappingEntry(variantIds = variant_ids,
                        susieFit   = list(),
                        topLoci    = tl)
  QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(e),
    ldSketch = .cip_makeHandle())
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("causalInferencePipeline: rejects non-GwasSumStats input", {
  expect_error(causalInferencePipeline(gwasSumStats = "no"),
               "must be a GwasSumStats")
})

test_that("causalInferencePipeline: rejects un-QCd GwasSumStats", {
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(qc = FALSE),
                            twasWeights  = .cip_makeTwasWeights()),
    "has no QC record"
  )
})

test_that("causalInferencePipeline: requires at least one of twasWeights/fineMappingResult", {
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats()),
    "at least one of"
  )
})

test_that("causalInferencePipeline: rejects non-TwasWeights twasWeights arg", {
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(),
                            twasWeights  = "not a TwasWeights"),
    "must be a TwasWeights"
  )
})

test_that("causalInferencePipeline: rejects GwasFineMappingResult for the QTL slot", {
  e <- FineMappingEntry(variantIds = "v1", susieFit = list(),
                        topLoci = data.frame(variant_id = "v1", pip = 0.1,
                                              stringsAsFactors = FALSE))
  gfmr <- GwasFineMappingResult(
    study = "G1", method = "susie", entry = list(e))
  expect_error(
    causalInferencePipeline(gwasSumStats      = .cip_makeGwasSumstats(),
                            fineMappingResult = gfmr),
    "does not accept GWAS-side fine"
  )
})

# ===========================================================================
# .cipRequireMatchingLdSketches branches
# ===========================================================================

test_that(".cipRequireMatchingLdSketches: NULL twas-side ldSketch is allowed", {
  # Build a TwasWeights without an ldSketch.
  twNoLd <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "lasso",
    entry = list(TwasWeightsEntry(variantIds = paste0("v", 1:5),
                                   weights = rep(0.1, 5))),
    ldSketch = NULL)
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats = .cip_makeGwasSumstats(),
    twasWeights  = twNoLd)
  expect_s4_class(out, "GRanges")
})

test_that(".cipRequireMatchingLdSketches: panel size mismatch errors", {
  bigSketch <- .cip_makeHandle(snp_n = 7L)
  twBig <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "lasso",
    entry = list(TwasWeightsEntry(variantIds = paste0("v", 1:5),
                                   weights = rep(0.1, 5))),
    ldSketch = bigSketch)
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(),
                            twasWeights  = twBig),
    "differ in size"
  )
})

# ===========================================================================
# Happy path: TwasWeights only
# ===========================================================================

test_that("causalInferencePipeline: returns GRanges with TWAS Z per tuple", {
  tw <- .cip_makeTwasWeights()
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(),
                                  twasWeights  = tw)
  expect_s4_class(out, "GRanges")
  expect_equal(length(out), 1L)
  mc <- S4Vectors::mcols(out)
  expect_equal(as.character(mc$qtlStudy[[1L]]), "Q1")
  expect_equal(as.character(mc$gwasStudy[[1L]]), "G1")
  expect_true(is.finite(mc$twasZ[[1L]]))
  expect_true(is.finite(mc$twasPval[[1L]]))
  # No FMR -> MR fields stay NA.
  expect_true(is.na(mc$waldRatio[[1L]]))
  expect_true(is.na(mc$mrPval[[1L]]))
})

# ===========================================================================
# Happy path: TwasWeights + matching FineMappingResult enables MR
# ===========================================================================

test_that("causalInferencePipeline: with both inputs, MR fields are populated", {
  tw <- .cip_makeTwasWeights()
  fmr <- .cip_makeQtlFmr()
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats       = .cip_makeGwasSumstats(),
    twasWeights        = tw,
    fineMappingResult  = fmr,
    mrPipCutoff        = 0.5)
  mc <- S4Vectors::mcols(out)
  # MR uses PIP > 0.5 variants from the FMR (v1, v3, v4).
  expect_true(is.finite(mc$waldRatio[[1L]]))
  expect_true(is.finite(mc$mrPval[[1L]]))
  expect_gt(mc$nIV[[1L]], 0L)
})

# ===========================================================================
# FineMappingResult-only path: weights come from topLoci$betahat
# ===========================================================================

test_that("causalInferencePipeline: FMR-only path extracts weights from topLoci$betahat", {
  fmr <- .cip_makeQtlFmr()
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats       = .cip_makeGwasSumstats(),
    fineMappingResult  = fmr)
  mc <- S4Vectors::mcols(out)
  expect_true(is.finite(mc$twasZ[[1L]]))
})

# ===========================================================================
# combineMethods integration
# ===========================================================================

test_that("causalInferencePipeline: combineMethods appends combined rows", {
  tw1 <- TwasWeights(
    study   = c("Q1", "Q1"), context = c("c1", "c1"),
    trait   = c("t1", "t1"), method = c("lasso", "enet"),
    entry = list(
      TwasWeightsEntry(variantIds = paste0("v", 1:5), weights = rep(0.1, 5)),
      TwasWeightsEntry(variantIds = paste0("v", 1:5), weights = rep(0.05, 5))),
    ldSketch = .cip_makeHandle())
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats   = .cip_makeGwasSumstats(),
    twasWeights    = tw1,
    combineMethods = "acat")
  # 2 per-tuple rows + 1 combined row = 3.
  expect_equal(length(out), 3L)
  mc <- S4Vectors::mcols(out)
  expect_true(any(grepl("^combined\\.", as.character(mc$method))))
})

# ===========================================================================
# .cipZToBeta / .cipZToSe fallbacks
# ===========================================================================

test_that(".cipZToBeta: falls back to z when maf/n are NA", {
  res <- pecotmr:::.cipZToBeta(z = c(1, 2), maf = NA, n = NA)
  expect_equal(res, c(1, 2))
})

test_that(".cipZToSe: falls back to vector of 1 when maf/n are NA", {
  res <- pecotmr:::.cipZToSe(z = c(1, 2), maf = NA, n = NA)
  expect_equal(res, c(1, 1))
})


context("twas: twasZ and harmonize deprecated wrappers")

# ===========================================================================
# Fixture builder: a small genotype + LD matrix.
# ===========================================================================

.tz_makeLd <- function(n = 100, p = 8, seed = 7) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, runif(p, 0.2, 0.8)), nrow = n, ncol = p)
  vid <- paste0("v", seq_len(p))
  colnames(X) <- vid
  af <- colMeans(X) / 2
  Xstd <- sweep(X, 2, 2 * af)
  Xstd <- sweep(Xstd, 2, sqrt(2 * af * (1 - af)), "/")
  R <- crossprod(Xstd) / (n - 1)
  rownames(R) <- colnames(R) <- vid
  list(X = X, Xstd = Xstd, R = R, vid = vid, n = n, p = p)
}

# ===========================================================================
# twasZ: input coercion and validity
# ===========================================================================

test_that("twasZ: numeric vector input is coerced to a single-method matrix", {
  d <- .tz_makeLd()
  w <- rnorm(d$p)
  names(w) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(w, z, R = d$R)
  expect_equal(dim(res$Z), c(1L, 2L))
  expect_equal(colnames(res$Z), c("Z", "pval"))
  expect_equal(rownames(res$Z), "method1")
})

test_that("twasZ: matrix without colnames gets method1/method2 names", {
  d <- .tz_makeLd()
  W <- cbind(rnorm(d$p), rnorm(d$p))
  rownames(W) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(W, z, R = d$R)
  expect_equal(rownames(res$Z), c("method1", "method2"))
})

test_that("twasZ: non-matrix non-numeric weights errors out", {
  expect_error(pecotmr:::twasZ(list(a = 1), z = 1, R = matrix(1)),
               "must be a numeric vector or a matrix")
})

test_that("twasZ: length mismatch between weights and z errors", {
  expect_error(
    pecotmr:::twasZ(c(0.1, 0.2, 0.3), z = c(1, 2)),
    "nrow\\(weights\\) must equal length\\(z\\)"
  )
})

test_that("twasZ: missing R/X/V triplet errors", {
  expect_error(
    pecotmr:::twasZ(c(0.1, 0.2), z = c(1, 2)),
    "provide R, X, or the .V, D, nSketch. SVD triplet"
  )
})

# ===========================================================================
# twasZ: R path
# ===========================================================================

test_that("twasZ: R path with named alignment matches the manual formula", {
  d <- .tz_makeLd()
  set.seed(1)
  z <- rnorm(d$p); names(z) <- d$vid
  w <- rnorm(d$p); names(w) <- d$vid
  res <- pecotmr:::twasZ(w, z, R = d$R)
  expected_stat <- sum(w * z)
  expected_denom <- sqrt(as.numeric(crossprod(w, d$R %*% w)))
  expect_equal(as.numeric(res$Z[, "Z"]),
               expected_stat / expected_denom,
               tolerance = 1e-10)
})

test_that("twasZ: R path realigns by rowname order", {
  d <- .tz_makeLd()
  z <- rnorm(d$p)
  w <- rnorm(d$p); names(w) <- d$vid
  # Permute R rows/columns; results must be identical to the un-permuted R.
  perm <- sample(seq_len(d$p))
  R_perm <- d$R[perm, perm]
  res_orig <- pecotmr:::twasZ(w, z, R = d$R)
  res_perm <- pecotmr:::twasZ(w, z, R = R_perm)
  expect_equal(res_perm$Z, res_orig$Z, tolerance = 1e-10)
})

test_that("twasZ: R path errors when R is missing rows named in weights", {
  d <- .tz_makeLd()
  w <- rnorm(d$p); names(w) <- d$vid
  R_short <- d$R[1:(d$p - 1), 1:(d$p - 1)]
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p), R = R_short),
    "R is missing rows for"
  )
})

test_that("twasZ: R path positional alignment errors on dim mismatch", {
  d <- .tz_makeLd()
  w <- rnorm(d$p)   # unnamed -> positional alignment
  R_short <- unname(d$R[1:(d$p - 1), 1:(d$p - 1)])
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p), R = R_short),
    "positional alignment requires nrow\\(R\\) == nrow\\(weights\\)"
  )
})

# ===========================================================================
# twasZ: X path (computes R via computeLd)
# ===========================================================================

test_that("twasZ: X path computes R via computeLd and matches that R explicitly", {
  d <- .tz_makeLd()
  set.seed(2)
  w <- rnorm(d$p); names(w) <- d$vid
  z <- rnorm(d$p)
  # The X path delegates to computeLd(X, method = "sample"), which uses
  # cor(X). Compare against the same correlation matrix passed explicitly.
  R_from_X <- cor(d$X)
  rownames(R_from_X) <- colnames(R_from_X) <- d$vid
  res_X   <- pecotmr:::twasZ(w, z, X = d$X)
  res_Rfx <- pecotmr:::twasZ(w, z, R = R_from_X)
  expect_equal(res_X$Z, res_Rfx$Z, tolerance = 1e-10)
})

# ===========================================================================
# twasZ: SVD path coverage (V-row missing error)
# ===========================================================================

test_that("twasZ: SVD path errors when V is missing rows named in weights", {
  d <- .tz_makeLd()
  s <- svd(d$Xstd)
  rownames(s$v) <- d$vid[seq_len(nrow(s$v))]
  w <- rnorm(d$p); names(w) <- c(d$vid[1:(d$p - 1)], "ghost")
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p),
                    V = s$v, D = s$d, nSketch = d$n),
    "V is missing rows for"
  )
})

test_that("twasZ: SVD path positional alignment errors on dim mismatch", {
  d <- .tz_makeLd()
  s <- svd(d$Xstd)
  # Pass an unnamed V with the wrong nrow.
  w <- rnorm(d$p)
  V_short <- s$v[1:(d$p - 1), , drop = FALSE]
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p),
                    V = V_short, D = s$d, nSketch = d$n),
    "positional alignment requires nrow\\(V\\) == nrow\\(weights\\)"
  )
})

# ===========================================================================
# twasZ: combineMethods integration
# ===========================================================================

test_that("twasZ: combineMethods K=1 returns the per-tuple p-value unchanged", {
  d <- .tz_makeLd()
  w <- rnorm(d$p); names(w) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(w, z, R = d$R, combineMethods = c("acat", "bonferroni"))
  expect_false(is.null(res$combined))
  expect_equal(res$combined$results$acat$pval, res$Z[1, "pval"])
  expect_equal(res$combined$results$bonferroni$pval, res$Z[1, "pval"])
  expect_equal(res$combined$input$nValid, 1L)
})

test_that("twasZ: combineMethods K>=2 forwards to combinePValues with correlation", {
  d <- .tz_makeLd()
  set.seed(11)
  W <- cbind(method_a = rnorm(d$p), method_b = rnorm(d$p))
  rownames(W) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(W, z, R = d$R, combineMethods = "acat")
  expect_false(is.null(res$combined))
  # ACAT combines two p-values; result must be in (0, 1).
  expect_true(res$combined$results$acat$pval >= 0 &&
              res$combined$results$acat$pval <= 1)
  expect_equal(res$combined$input$nValid, 2L)
})

# ===========================================================================
# Deprecated wrappers
# ===========================================================================

test_that("harmonizeTwas is a deprecated no-op", {
  expect_warning(res <- harmonizeTwas(), "has been removed", ignore.case = TRUE)
  expect_null(res)
})

test_that("harmonizeGwas is a deprecated no-op", {
  expect_warning(res <- harmonizeGwas(), "has been removed", ignore.case = TRUE)
  expect_null(res)
})

test_that("twasPipeline is a deprecated no-op", {
  expect_warning(res <- twasPipeline(), "has been removed", ignore.case = TRUE)
  expect_null(res)
})

# ===========================================================================
# Tests migrated from test_twas.R (twasZ behaviour)
# ===========================================================================

test_that("twasZ errors when weights and z lengths differ", {
  expect_error(twasZ(c(1, 2), c(1, 2, 3)), "must equal")
})


test_that("twasZ: single weight and single z-score returns valid result", {
  weights <- 0.7
  z <- 2.5
  R <- matrix(1, nrow = 1, ncol = 1)
  result <- twasZ(weights, z, R = R)
  expect_true(is.list(result))
  expect_equal(colnames(result$Z), c("Z", "pval"))
  # With single variant: stat = 0.7 * 2.5 = 1.75, denom = 0.7 * 1 * 0.7 = 0.49
  # zscore = 1.75 / sqrt(0.49) = 1.75 / 0.7 = 2.5
  expect_equal(as.numeric(result$Z[, "Z"]), 2.5, tolerance = 1e-10)
  expect_true(result$Z[, "pval"] > 0 && result$Z[, "pval"] < 1)
})


test_that("twasZ: all-zero weights produce NaN z-score", {
  weights <- c(0, 0, 0)
  z <- c(1.5, -0.5, 2.0)
  R <- diag(3)
  result <- twasZ(weights, z, R = R)
  # stat = 0, denom = 0, so zscore = 0/0 = NaN
  expect_true(is.nan(as.numeric(result$Z[, "Z"])))
})


test_that("twasZ: very large z-scores still produce finite results", {
  set.seed(42)
  p <- 5
  weights <- rnorm(p)
  z <- rep(1e6, p)
  R <- diag(p)
  result <- twasZ(weights, z, R = R)
  expect_true(is.finite(as.numeric(result$Z[, "Z"])))
  # p-value should be extremely small for large z
  expect_true(result$Z[, "pval"] < 1e-10 || result$Z[, "pval"] == 0)
})


test_that("twasZ: identical z-scores with equal weights gives proportional result", {
  p <- 5
  weights <- rep(1, p)
  z <- rep(3.0, p)
  R <- diag(p)
  result <- twasZ(weights, z, R = R)
  # stat = sum(weights * z) = 5 * 3 = 15
  # denom = t(w) %*% I %*% w = 5
  # zscore = 15 / sqrt(5) = 6.7082...
  expect_equal(as.numeric(result$Z[, "Z"]), 15 / sqrt(5), tolerance = 1e-10)
})


test_that("twasZ: negative weights flip the sign of the z-score", {
  weights_pos <- c(0.5, 0.3)
  weights_neg <- c(-0.5, -0.3)
  z <- c(2.0, 1.0)
  R <- diag(2)
  result_pos <- twasZ(weights_pos, z, R = R)
  result_neg <- twasZ(weights_neg, z, R = R)
  # z-score should have opposite sign but same p-value
  expect_equal(as.numeric(result_pos$Z[, "Z"]),
               -as.numeric(result_neg$Z[, "Z"]),
               tolerance = 1e-10)
  expect_equal(as.numeric(result_pos$Z[, "pval"]),
               as.numeric(result_neg$Z[, "pval"]),
               tolerance = 1e-10)
})


test_that("twasZ: off-diagonal correlation in R changes the result", {
  weights <- c(0.5, 0.5)
  z <- c(2.0, 2.0)
  R_identity <- diag(2)
  R_correlated <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  result_identity <- twasZ(weights, z, R = R_identity)
  result_correlated <- twasZ(weights, z, R = R_correlated)
  # Same stat but different denominators, so different z-scores
  expect_false(isTRUE(all.equal(
    as.numeric(result_identity$Z[, "Z"]),
    as.numeric(result_correlated$Z[, "Z"]))))
})


test_that("twasZ: computing R from X matches providing R directly", {
  set.seed(123)
  n <- 20
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("SNP", 1:p)
  R <- cor(X)
  weights <- rnorm(p)
  z <- rnorm(p)
  result_with_R <- twasZ(weights, z, R = R)
  result_with_X <- twasZ(weights, z, X = X)
  expect_equal(as.numeric(result_with_R$Z[, "Z"]),
               as.numeric(result_with_X$Z[, "Z"]),
               tolerance = 1e-6)
  expect_equal(as.numeric(result_with_R$Z[, "pval"]),
               as.numeric(result_with_X$Z[, "pval"]),
               tolerance = 1e-6)
})


test_that("twasZ: p-value is always between 0 and 1 for random inputs", {
  set.seed(999)
  for (i in 1:5) {
    p <- sample(2:10, 1)
    weights <- rnorm(p)
    z <- rnorm(p)
    R <- diag(p) # use identity to avoid singularity
    result <- twasZ(weights, z, R = R)
    pval <- as.numeric(result$Z[, "pval"])
    expect_true(pval >= 0 && pval <= 1,
      info = paste("Iteration", i, "p-value out of range:", pval))
  }
})


test_that("twasZ: single very large weight with tiny z gives moderate result", {
  weights <- c(1e6, 0, 0)
  z <- c(1e-6, 5.0, 5.0)
  R <- diag(3)
  result <- twasZ(weights, z, R = R)
  # stat = 1e6 * 1e-6 + 0 + 0 = 1.0
  # denom = (1e6)^2 * 1 = 1e12
  # zscore = 1 / sqrt(1e12) = 1e-6
  expect_equal(as.numeric(result$Z[, "Z"]), 1e-6, tolerance = 1e-10)
})

# ===========================================================================
# twasZ: more mathematical edge cases
# ===========================================================================


test_that("twasZ: perfectly correlated R matrix (all ones off-diagonal)", {
  p <- 4
  R <- matrix(1, nrow = p, ncol = p) # perfectly correlated
  weights <- c(0.25, 0.25, 0.25, 0.25)
  z <- c(2, 2, 2, 2)
  result <- twasZ(weights, z, R = R)
  # stat = sum(0.25 * 2) = 2
  # denom = t(w) %*% ones_matrix %*% w = (sum(w))^2 = 1
  # zscore = 2 / sqrt(1) = 2
  expect_equal(as.numeric(result$Z[, "Z"]), 2.0, tolerance = 1e-10)
})


test_that("twasZ: sparse weights (only one non-zero) extracts single SNP signal", {
  p <- 5
  weights <- c(0, 0, 1, 0, 0)
  z <- c(1, 2, 3, 4, 5)
  R <- diag(p)
  result <- twasZ(weights, z, R = R)
  # With identity R and sparse weight: zscore = w3 * z3 / sqrt(w3^2) = 3
  expect_equal(as.numeric(result$Z[, "Z"]), 3.0, tolerance = 1e-10)
})


test_that("twasZ: z-scores of zero give zero TWAS z-score", {
  weights <- c(0.5, 0.3, 0.2)
  z <- c(0, 0, 0)
  R <- diag(3)
  result <- twasZ(weights, z, R = R)
  expect_equal(as.numeric(result$Z[, "Z"]), 0.0, tolerance = 1e-10)
  # p-value for z=0 should be 1
  expect_equal(as.numeric(result$Z[, "pval"]), 1.0, tolerance = 1e-10)
})

# ===========================================================================
# twasZ: edge case with near-singular R matrix
# ===========================================================================


test_that("twasZ: near-singular R matrix still produces a result", {
  set.seed(777)
  p <- 3
  # Create a nearly singular R by making two rows almost identical
  R <- matrix(c(1, 0.999, 0.999, 0.999, 1, 0.999, 0.999, 0.999, 1), nrow = 3)
  weights <- c(0.3, 0.4, 0.3)
  z <- c(2.0, 2.5, 1.8)
  result <- twasZ(weights, z, R = R)
  expect_true(is.list(result))
  expect_true(is.finite(as.numeric(result$Z[, "Z"])))
})


test_that("twasZ: length-one input vectors produce correct scalar output", {
  result <- twasZ(1.0, 3.0, R = matrix(1, 1, 1))
  expect_equal(as.numeric(result$Z[, "Z"]), 3.0, tolerance = 1e-10)
  expect_true(result$Z[, "pval"] > 0 && result$Z[, "pval"] < 1)
})


test_that("twasZ: R=NULL and X=NULL still errors consistently", {
  # When neither R nor X is provided, twasZ should error
  expect_error(twasZ(c(1, 2), c(3, 4)))
})

# ===========================================================================
# twasZ: matrix weights path -- multiple methods/conditions
# ===========================================================================


test_that("twasZ: matrix weights produce one Z row per column", {
  set.seed(10)
  p <- 5
  k <- 3
  weights <- matrix(rnorm(p * k), nrow = p, ncol = k)
  rownames(weights) <- paste0("SNP", 1:p)
  colnames(weights) <- paste0("Cond", 1:k)
  z <- rnorm(p)
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("SNP", 1:p)
  result <- twasZ(weights, z, R = R)
  expect_equal(nrow(result$Z), k)
  expect_equal(rownames(result$Z), paste0("Cond", 1:k))
  expect_equal(colnames(result$Z), c("Z", "pval"))
  # combineMethods omitted -> combined is NULL
  expect_null(result$combined)
})


test_that("twasZ: combineMethods returns combined p-value summary", {
  skip_if_not_installed("ACAT")
  set.seed(11)
  p <- 4
  k <- 2
  weights <- matrix(rnorm(p * k), nrow = p, ncol = k,
                    dimnames = list(paste0("SNP", 1:p), paste0("Cond", 1:k)))
  z <- rnorm(p)
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("SNP", 1:p)
  result <- twasZ(weights, z, R = R, combineMethods = "ACAT")
  expect_false(is.null(result$combined))
})

# ===========================================================================
# Tests from test_twas_predict.R
# ===========================================================================

# ---- twasPredict ----


# ===========================================================================
# Tests migrated from test_twasSketch.R (twasZ SVD path)
# ===========================================================================

test_that("twasZ: SVD path matches R path for full-rank genotype matrix", {
  set.seed(42)
  n <- 100 # samples (sketch size)
  p <- 20 # variants

  # Generate genotype-like matrix (dosages 0/1/2)
  X <- matrix(rbinom(n * p, 2, runif(p, 0.1, 0.9)), nrow = n, ncol = p)

  # HWE-based standardization
  af <- colMeans(X) / 2
  X_std <- sweep(X, 2, 2 * af)
  X_std <- sweep(X_std, 2, sqrt(2 * af * (1 - af)), "/")

  # Compute R from standardized X
  R <- crossprod(X_std) / (n - 1)

  # SVD of standardized X
  svd_result <- svd(X_std)

  # Random weights and z-scores
  weights <- rnorm(p)
  z <- rnorm(p)

  # R path
  result_R <- pecotmr:::twasZ(weights, z, R = R)

  # SVD path
  result_SVD <- pecotmr:::twasZ(weights, z, V = svd_result$v, D = svd_result$d, nSketch = n)

  expect_equal(as.numeric(result_SVD$z), as.numeric(result_R$z), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD$pval), as.numeric(result_R$pval), tolerance = 1e-10)
})


test_that("twasZ: SVD path matches R path for rank-deficient matrix (n < p)", {
  set.seed(123)
  n <- 15 # fewer samples than variants

  p <- 30

  X <- matrix(rbinom(n * p, 2, runif(p, 0.15, 0.85)), nrow = n, ncol = p)

  af <- colMeans(X) / 2
  # Remove monomorphic columns
  keep <- af > 0 & af < 1
  X <- X[, keep, drop = FALSE]
  af <- af[keep]
  p <- ncol(X)

  X_std <- sweep(X, 2, 2 * af)
  X_std <- sweep(X_std, 2, sqrt(2 * af * (1 - af)), "/")

  R <- crossprod(X_std) / (n - 1)
  svd_result <- svd(X_std)

  weights <- rnorm(p)
  z <- rnorm(p)

  result_R <- pecotmr:::twasZ(weights, z, R = R)
  result_SVD <- pecotmr:::twasZ(weights, z, V = svd_result$v, D = svd_result$d, nSketch = n)

  expect_equal(as.numeric(result_SVD$z), as.numeric(result_R$z), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD$pval), as.numeric(result_R$pval), tolerance = 1e-10)
})


test_that("twasZ: error when weights and z have different lengths", {
  expect_error(
    pecotmr:::twasZ(rnorm(5), rnorm(3), V = matrix(1, 5, 2), D = c(1, 1), nSketch = 10),
    "nrow\\(weights\\) must equal length\\(z\\)"
  )
})

# Phase 2: loadLdSketch() and standardize_genotype_hwe()


