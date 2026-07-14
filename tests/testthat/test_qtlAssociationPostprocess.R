context("qtlAssociationPostprocess")

# A deterministic per-gene QtlSumStats: one row per gene carrying the regional
# summary as row columns (p_beta / beta shapes / n_variants), each entry a
# per-variant GRanges (P / af / tss_distance / tes_distance). Deterministic so
# qvalue::qvalue sees a stable p-value distribution.
#   drop       row columns to omit (e.g. "n_variants", "p_beta", beta shapes),
#              for the required-column / missing-input error paths.
#   withQvalue attach a per-variant `qvalue` mcol (needed by the qvalue SNP
#              significance query), with v1 tracking the gene signal.
#   emptyIdx   gene indices whose entry is an empty GRanges (no p-values), for
#              the "skip genes with no variants" branch.
#   naShapeIdx gene indices given an NA beta_shape1 (so qbeta -> NA nominal
#              threshold), for the "skip genes with an NA threshold" branch.
# Column counts and the p_beta distribution are unchanged from the base
# fixture, so qvalue::qvalue still estimates pi0 on its default path.
.qapFixture <- function(nVar = 50L, nVarFilt = 30L, drop = character(),
                        withQvalue = FALSE, emptyIdx = integer(0),
                        naShapeIdx = integer(0)) {
  # 6 strong signal p_beta + a deterministic uniform null (ppoints) so
  # qvalue::qvalue estimates pi0 on its default path (no fallback needed).
  pBeta <- c(10^(-c(8, 7, 6, 5, 4, 3)), stats::ppoints(54))
  G <- length(pBeta)
  entries <- lapply(seq_len(G), function(i) {
    if (i %in% emptyIdx) return(GenomicRanges::GRanges())  # no variants
    k  <- 4L
    pv <- c(pBeta[i] / 5, 0.2, 0.5, 0.8)             # min tracks the gene signal
    gr <- GenomicRanges::GRanges("chr1",
            IRanges::IRanges(start = seq(1000L, by = 50L, length.out = k), width = 1L))
    mc <- S4Vectors::DataFrame(
      SNP = paste0("g", i, "_v", seq_len(k)), A1 = "A", A2 = "G",
      P = pv, af = c(0.30, 0.20, 0.005, 0.40),        # v3 fails MAF 0.01
      tss_distance = c(0L, 500L, 900000L, 2000000L),  # v4 outside 1Mb cis
      tes_distance = c(0L, 500L, 900000L, 2000000L))
    if (withQvalue) mc$qvalue <- c(pBeta[i] / 5, 0.30, 0.60, 0.90)
    S4Vectors::mcols(gr) <- mc
    gr
  })
  shape1 <- rep(1, G); shape1[naShapeIdx] <- NA_real_
  args <- list(
    study = rep("s", G), context = rep("brain", G),
    trait = paste0("g", seq_len(G)), entry = entries, genome = "hg19",
    n_variants = rep(nVar, G), n_variants_filtered = rep(nVarFilt, G),
    p_beta = pBeta, beta_shape1 = shape1, beta_shape2 = rep(200, G))
  do.call(QtlSumStats, args[setdiff(names(args), drop)])
}

test_that("qtlAssociationPostprocess enriches with package-computed columns", {
  x <- .qapFixture()
  r <- qtlAssociationPostprocess(x, fdrThreshold = 0.1,
                                 mafCutoff = 0.01, cisWindow = 1e6)

  expect_s4_class(r, "QtlSumStats")
  expect_false(is.null(r@qcInfo$associationPostprocess))          # recipe stashed

  # Bonferroni original == min over variants of p.adjust(P, "bonferroni", n).
  expP <- vapply(seq_len(nrow(x)), function(i)
    min(stats::p.adjust(S4Vectors::mcols(x$entry[[i]])$P, "bonferroni", n = 50)),
    numeric(1))
  expect_equal(as.numeric(r$p_bonferroni_min_original), expP)
  expect_equal(as.numeric(r$fdr_bonferroni_min_original),
               stats::p.adjust(expP, "fdr"))
  expect_true(all(r$q_bonferroni_min_original >= 0 &
                  r$q_bonferroni_min_original <= 1))

  # Filtered flavour present and uses the smaller test count (n_variants_filtered).
  expect_true(all(c("p_bonferroni_min_filtered", "fdr_bonferroni_min_filtered",
                    "q_bonferroni_min_filtered") %in% names(r)))
})

test_that("q-values come from qvalue::qvalue, not a hand-rolled fallback", {
  skip_if_not_installed("qvalue")
  x <- .qapFixture()
  r <- qtlAssociationPostprocess(x, methods = "permutation")
  expect_equal(as.numeric(r$q_beta),
               qvalue::qvalue(as.numeric(x$p_beta))$qvalues)
  expect_equal(as.numeric(r$fdr_beta),
               stats::p.adjust(as.numeric(x$p_beta), "fdr"))
})

test_that("permutation nominal threshold == stats::qbeta of the empirical cutoff", {
  x  <- .qapFixture()
  r  <- qtlAssociationPostprocess(x, fdrThreshold = 0.1, methods = "permutation")
  qb <- as.numeric(r$q_beta)
  pb <- as.numeric(x$p_beta)
  lb <- sort(pb[qb <= 0.1]); ub <- sort(pb[qb > 0.1])
  expect_gt(length(lb), 0)
  thr <- if (length(ub) > 0) (lb[length(lb)] + ub[1]) / 2 else lb[length(lb)]
  expect_equal(as.numeric(r$p_nominal_threshold),
               stats::qbeta(thr, rep(1, nrow(x)), rep(200, nrow(x))))
})

test_that("getSignificantQtls (bonferroni) matches the derived threshold rule", {
  x <- .qapFixture()
  r <- qtlAssociationPostprocess(x, mafCutoff = 0.01, cisWindow = 1e6)
  sig <- getSignificantQtls(r, "bonferroni_original", threshold = 0.5)
  expect_s4_class(sig, "GRanges")
  expect_true("trait" %in% names(S4Vectors::mcols(sig)))

  # Reproduce the rule: variants with P*n_variants <= max(p_bonf_min over FDR-sig genes).
  fdr <- as.numeric(r$fdr_bonferroni_min_original)
  sigGenes <- which(fdr < 0.5)
  expect_gt(length(sigGenes), 0)
  varThr <- max(as.numeric(r$p_bonferroni_min_original)[sigGenes])
  expN <- sum(vapply(seq_len(nrow(r)), function(i)
    sum(pmin(1, S4Vectors::mcols(r$entry[[i]])$P * 50) <= varThr), integer(1)))
  expect_equal(length(sig), expN)
})

test_that("getSumStats(annotateSignificance=) adds a derived logical mcol", {
  x   <- .qapFixture()
  r   <- qtlAssociationPostprocess(x, methods = "permutation")
  gr  <- getSumStats(r, study = "s", context = "brain", trait = "g1",
                     annotateSignificance = "permutation")
  expect_true("significant" %in% names(S4Vectors::mcols(gr)))
  expect_type(S4Vectors::mcols(gr)$significant, "logical")
  # matches the permutation rule P < p_nominal_threshold[gene 1]
  thr1 <- as.numeric(r$p_nominal_threshold)[1]
  exp1 <- if (is.na(thr1)) rep(FALSE, length(gr)) else
    S4Vectors::mcols(gr)$P < thr1
  expect_equal(S4Vectors::mcols(gr)$significant, exp1)
})

test_that("significance accessors require a postprocessed object", {
  x <- .qapFixture()                                   # not postprocessed
  expect_error(getSignificantQtls(x, "bonferroni_original"),
               "qtlAssociationPostprocess")
})

# ---- qvalue-native edge-case retries (the ONLY substitute is qvalue itself) --
# The fixture carries p_beta but no q_beta, so permutation postprocessing routes
# q_beta through .qapSafeQvalue -> qvalue::qvalue. Mock qvalue to raise each
# qvalue-native error on the primary (no-lambda) call and succeed on the retry
# (lambda supplied), then assert the retry value flows through unchanged.

test_that(".qapSafeQvalue retries with lambda=0 on qvalue 'missing or infinite'", {
  skip_if_not_installed("qvalue")
  testthat::local_mocked_bindings(
    qvalue = function(p, lambda, ...) {
      if (missing(lambda)) stop("missing or infinite values in the smoother")
      list(qvalues = rep(0.111, length(p)))
    }, .package = "qvalue")
  r <- qtlAssociationPostprocess(.qapFixture(), methods = "permutation")
  expect_true(all(as.numeric(r$q_beta) == 0.111))
})

test_that(".qapSafeQvalue retries with bootstrap pi0 on qvalue 'pi0 <= 0'", {
  skip_if_not_installed("qvalue")
  testthat::local_mocked_bindings(
    qvalue = function(p, lambda, ...) {
      if (missing(lambda)) stop("ERROR: The estimated pi0 <= 0.")
      list(qvalues = rep(0.222, length(p)))
    }, .package = "qvalue")
  r <- qtlAssociationPostprocess(.qapFixture(), methods = "permutation")
  expect_true(all(as.numeric(r$q_beta) == 0.222))
})

test_that(".qapSafeQvalue re-raises a non-native qvalue error (no hand-rolled q)", {
  skip_if_not_installed("qvalue")
  testthat::local_mocked_bindings(
    qvalue = function(p, ...) stop("some unrelated numerical failure"),
    .package = "qvalue")
  expect_error(
    qtlAssociationPostprocess(.qapFixture(), methods = "permutation"),
    "qvalue::qvalue failed")
})

test_that("permutation nominal threshold is NA for every gene when none pass FDR", {
  # An impossibly small FDR threshold leaves no q_beta-significant gene, so the
  # empirical p_beta bracketing is empty and every nominal threshold is NA.
  r <- qtlAssociationPostprocess(.qapFixture(), fdrThreshold = 1e-300,
                                 methods = "permutation")
  expect_true(all(is.na(as.numeric(r$p_nominal_threshold))))
})

# ---- getSignificantQtls: one branch of the significance mask per method -------

test_that("getSignificantQtls (permutation) uses each gene's nominal threshold", {
  r   <- qtlAssociationPostprocess(.qapFixture(), fdrThreshold = 0.1,
                                   methods = "permutation")
  sig <- getSignificantQtls(r, "permutation")          # threshold defaults to recipe
  expect_s4_class(sig, "GRanges")
  expect_true("trait" %in% names(S4Vectors::mcols(sig)))
  # Reproduce: per gene, variants with P < p_nominal_threshold[gene].
  thr  <- as.numeric(r$p_nominal_threshold)
  expN <- sum(vapply(seq_len(nrow(r)), function(i) {
    if (is.na(thr[i])) return(0L)
    sum(S4Vectors::mcols(r$entry[[i]])$P < thr[i])
  }, integer(1)))
  expect_gt(expN, 0)
  expect_equal(length(sig), expN)
})

test_that("getSignificantQtls (permutation) skips genes with an NA threshold", {
  # Gene 1's NA beta_shape1 makes its qbeta nominal threshold NA, so the mask
  # loop must `next` past it while still processing the non-NA genes.
  r   <- qtlAssociationPostprocess(.qapFixture(naShapeIdx = 1L),
                                   fdrThreshold = 0.1, methods = "permutation")
  thr <- as.numeric(r$p_nominal_threshold)
  expect_true(is.na(thr[1]))              # NA-threshold gene present ...
  expect_true(any(!is.na(thr)))           # ... alongside non-NA genes
  sig <- getSignificantQtls(r, "permutation")
  expect_s4_class(sig, "GRanges")
  # Gene 1 contributes nothing; total is the sum over the non-NA-threshold genes.
  expN <- sum(vapply(seq_len(nrow(r)), function(i) {
    if (is.na(thr[i])) return(0L)
    sum(S4Vectors::mcols(r$entry[[i]])$P < thr[i])
  }, integer(1)))
  expect_equal(length(sig), expN)
})

test_that("getSignificantQtls (permutation) errors without a nominal threshold", {
  # No beta shapes -> postprocess produces no p_nominal_threshold column.
  r <- qtlAssociationPostprocess(
    .qapFixture(drop = c("beta_shape1", "beta_shape2")), methods = "permutation")
  expect_null(r$p_nominal_threshold)
  expect_error(getSignificantQtls(r, "permutation"), "p_nominal_threshold")
})

test_that("getSignificantQtls (bonferroni_filtered) applies the MAF/cis keep filter", {
  r   <- qtlAssociationPostprocess(.qapFixture(), mafCutoff = 0.01, cisWindow = 1e6)
  sig <- getSignificantQtls(r, "bonferroni_filtered", threshold = 0.5)
  expect_s4_class(sig, "GRanges")
  expect_gt(length(sig), 0)
  # Every returned variant must satisfy the filtered keep rule (MAF & cis window).
  af <- S4Vectors::mcols(sig)$af
  expect_true(all(pmin(af, 1 - af) > 0.01))
  expect_true(all(S4Vectors::mcols(sig)$tss_distance >= -1e6 &
                  S4Vectors::mcols(sig)$tes_distance <=  1e6))
})

test_that("getSignificantQtls (bonferroni) errors when its columns are absent", {
  r <- qtlAssociationPostprocess(.qapFixture(), methods = "permutation")  # no bonferroni
  expect_error(getSignificantQtls(r, "bonferroni_original"), "columns absent")
})

test_that("getSignificantQtls (qvalue) selects variants by their qvalue mcol", {
  skip_if_not_installed("qvalue")
  r   <- qtlAssociationPostprocess(.qapFixture(withQvalue = TRUE),
                                   fdrThreshold = 0.1, methods = "permutation")
  sig <- getSignificantQtls(r, "qvalue", threshold = 0.1)
  expect_s4_class(sig, "GRanges")
  # Reproduce: within genes whose event q_beta < 0.1, variants with qvalue < 0.1.
  qb       <- as.numeric(r$q_beta)
  sigGenes <- which(qb < 0.1)
  expect_gt(length(sigGenes), 0)
  expN <- sum(vapply(sigGenes, function(i)
    sum(S4Vectors::mcols(r$entry[[i]])$qvalue < 0.1), integer(1)))
  expect_gt(expN, 0)
  expect_equal(length(sig), expN)
})

test_that("getSignificantQtls (qvalue) errors when no event q column exists", {
  # Drop p_beta so permutation adds no q_beta, and skip bonferroni so there is
  # no q_bonferroni_min_original either -- yet the object is still postprocessed.
  r <- qtlAssociationPostprocess(.qapFixture(drop = "p_beta"), methods = "permutation")
  expect_null(r$q_beta)
  expect_null(r$q_bonferroni_min_original)
  expect_error(getSignificantQtls(r, "qvalue"), "event q-value column")
})

test_that("getSignificantQtls returns an empty GRanges when nothing is significant", {
  r   <- qtlAssociationPostprocess(.qapFixture(), mafCutoff = 0.01, cisWindow = 1e6)
  sig <- getSignificantQtls(r, "bonferroni_original", threshold = 1e-300)
  expect_s4_class(sig, "GRanges")
  expect_length(sig, 0)
})

# ---- qtlAssociationPostprocess: required-column + empty-entry guards ----------

test_that("qtlAssociationPostprocess requires n_variants for Bonferroni", {
  expect_error(
    qtlAssociationPostprocess(.qapFixture(drop = "n_variants"), methods = "bonferroni"),
    "required for the Bonferroni")
})

test_that("qtlAssociationPostprocess requires n_variants_filtered when filtering", {
  expect_error(
    qtlAssociationPostprocess(.qapFixture(drop = "n_variants_filtered"),
                              methods = "bonferroni", mafCutoff = 0.01),
    "n_variants_filtered")
})

test_that("qtlAssociationPostprocess skips genes whose entry has no variants", {
  r     <- qtlAssociationPostprocess(.qapFixture(emptyIdx = 1L), methods = "bonferroni")
  pmin1 <- as.numeric(r$p_bonferroni_min_original)
  expect_true(is.na(pmin1[1]))          # empty gene -> NA Bonferroni min
  expect_true(all(!is.na(pmin1[-1])))   # every other gene is finite
})

test_that("FILTERED Bonferroni drops MAF/cis-failing variants + uses n_variants_filtered", {
  # Fixture entries: af = (.30, .20, .005, .40), tss/tes_distance = (0, 500, 9e5, 2e6).
  # With mafCutoff 0.01 + cisWindow 1e6 the filtered set keeps only v1, v2
  # (v3 fails MAF, v4 is outside the cis window); the filtered count is 30.
  x <- .qapFixture()
  r <- qtlAssociationPostprocess(x, mafCutoff = 0.01, cisWindow = 1e6,
                                 methods = "bonferroni")
  P <- lapply(seq_len(nrow(x)), function(i) S4Vectors::mcols(x$entry[[i]])$P)
  expFilt <- vapply(P, function(p) min(pmin(1, p[1:2] * 30)), numeric(1))   # v1,v2 @ n=30
  expOrig <- vapply(P, function(p) min(pmin(1, p * 50)), numeric(1))         # all @ n=50
  expect_equal(as.numeric(r$p_bonferroni_min_filtered), expFilt)
  expect_equal(as.numeric(r$p_bonferroni_min_original), expOrig)
  # A smaller test count makes the filtered flavour no less significant.
  expect_true(all(r$p_bonferroni_min_filtered <= r$p_bonferroni_min_original + 1e-12))
})

test_that("getSignificantQtls(bonferroni_filtered) applies the derived rule on the filtered set", {
  x <- .qapFixture()
  r <- qtlAssociationPostprocess(x, mafCutoff = 0.01, cisWindow = 1e6,
                                 methods = "bonferroni")
  sig <- getSignificantQtls(r, "bonferroni_filtered", threshold = 0.5)
  expect_s4_class(sig, "GRanges")
  fdr <- as.numeric(r$fdr_bonferroni_min_filtered)
  sigGenes <- which(!is.na(fdr) & fdr < 0.5)
  expect_gt(length(sigGenes), 0)                          # signal genes pass
  varThr <- max(as.numeric(r$p_bonferroni_min_filtered)[sigGenes])
  # only MAF/cis-passing variants (v1,v2) with P*n_filtered <= threshold qualify
  expN <- sum(vapply(seq_len(nrow(r)), function(i) {
    p <- S4Vectors::mcols(r$entry[[i]])$P[1:2]
    sum(pmin(1, p * 30) <= varThr)
  }, integer(1)))
  expect_equal(length(sig), expN)
})
