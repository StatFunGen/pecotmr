context("fsusieAccessors")

.fsa_makeFit <- function() {
  set.seed(1); n <- 150L; p <- 24L; J <- 16L
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("chr1:", (seq_len(p)) * 100, ":A:G")))
  b1 <- sin(seq(0, 2 * pi, length.out = J))
  b2 <- cos(seq(0, pi, length.out = J))
  Y <- X[, 3] %o% b1 + X[, 10] %o% b2 + matrix(rnorm(n * J, sd = 0.3), n, J)
  colnames(Y) <- paste0("f", seq_len(J))
  suppressWarnings(fsusieR::susiF(X = X, Y = Y, pos = seq_len(J), L = 5,
                                  post_processing = "none", verbose = FALSE))
}
.fsa_entry <- function(fit, purityVal = 0.85) {
  vids <- names(fit$csd_X)
  pipv <- if (!is.null(fit$pip)) as.numeric(fit$pip) else rep(0.1, length(vids))
  # Canonical CS labels + purity from the fit's CS membership, mirroring what
  # fsusieWrapper stamps into a real FMR's topLoci (cs_95 / cs_95_purity).
  cs95 <- rep("fsusie_0", length(vids))
  for (l in seq_along(fit$cs)) {
    idx <- fit$cs[[l]]; if (is.numeric(idx)) idx <- as.integer(idx)
    cs95[idx] <- paste0("fsusie_", l)
  }
  tl <- data.frame(variant_id = vids, pip = pipv, cs_95 = cs95,
                   cs_95_purity = ifelse(cs95 != "fsusie_0", purityVal, 0),
                   stringsAsFactors = FALSE)
  FineMappingEntry(variantIds = vids, susieFit = fit, topLoci = tl)
}

test_that("fsusieCredibleBand returns a long effect + band table (lower <= effect <= upper)", {
  skip_if_not_installed("fsusieR"); skip_if_not_installed("wavethresh")
  cb <- fsusieCredibleBand(.fsa_entry(.fsa_makeFit()))
  expect_named(cb, c("cs", "chrom", "pos", "effect", "lower", "upper"))
  expect_gt(nrow(cb), 0)
  expect_true(all(cb$lower <= cb$effect + 1e-9 & cb$effect <= cb$upper + 1e-9))
  expect_true(all(grepl("^fsusie_", cb$cs)))
})

test_that("fsusieAffectedRegions returns a GRanges with cs / purity / direction", {
  skip_if_not_installed("fsusieR"); skip_if_not_installed("wavethresh")
  gr <- fsusieAffectedRegions(.fsa_entry(.fsa_makeFit(), purityVal = 0.85))
  expect_s4_class(gr, "GRanges")
  expect_gt(length(gr), 0)
  expect_true(all(c("cs", "purity", "direction") %in% names(S4Vectors::mcols(gr))))
  expect_true(all(S4Vectors::mcols(gr)$direction %in% c("pos", "neg", NA)))
  # purity sourced from the entry's cs_95_purity (matched by CS variant membership),
  # NOT the (nonexistent) fit$purity slot
  expect_true(all(S4Vectors::mcols(gr)$purity == 0.85))
  expect_true(all(grepl("^fsusie_", S4Vectors::mcols(gr)$cs)))
})

test_that("fsusie accessors degrade to empty for a non-fSuSiE / trimmed fit (no wavelet slots)", {
  e <- FineMappingEntry(variantIds = "chr1:100:A:G", susieFit = list(pip = 0.5),
                        topLoci = data.frame(variant_id = "chr1:100:A:G", pip = 0.5,
                                             stringsAsFactors = FALSE))
  expect_equal(nrow(fsusieCredibleBand(e)), 0L)
  expect_equal(length(fsusieAffectedRegions(e)), 0L)
})

test_that("fsusieCredibleBand + fsusieAffectedRegions aggregate across a collection", {
  skip_if_not_installed("fsusieR"); skip_if_not_installed("wavethresh")
  fit <- .fsa_makeFit()                       # deterministic; reuse for both rows
  res <- QtlFineMappingResult(
    study = c("s", "s"), context = c("brain", "blood"),
    trait = c("g", "g"), method = c("fsusie", "fsusie"),
    entry = list(.fsa_entry(fit), .fsa_entry(fit)))

  cb <- fsusieCredibleBand(res)
  expect_true(all(c("study", "context", "trait", "method",
                    "cs", "effect", "lower", "upper") %in% names(cb)))
  expect_gt(nrow(cb), 0)
  expect_setequal(unique(cb$context), c("brain", "blood"))

  gr <- fsusieAffectedRegions(res)
  expect_s4_class(gr, "GRanges")
  expect_gt(length(gr), 0)
  expect_true(all(c("study", "context", "trait", "cs", "purity", "direction")
                  %in% names(S4Vectors::mcols(gr))))
  expect_setequal(unique(S4Vectors::mcols(gr)$context), c("brain", "blood"))
})

test_that("fsusieAffectedRegions on a collection of non-fSuSiE entries is an empty GRanges", {
  e <- FineMappingEntry(variantIds = "chr1:100:A:G", susieFit = list(pip = 0.5),
                        topLoci = data.frame(variant_id = "chr1:100:A:G", pip = 0.5,
                                             stringsAsFactors = FALSE))
  res <- QtlFineMappingResult(study = "s", context = "c", trait = "t",
                              method = "susie", entry = list(e))
  expect_equal(length(fsusieAffectedRegions(res)), 0L)
})

test_that(".fsusieChrom returns NA when the fit has no named variants", {
  expect_true(is.na(pecotmr:::.fsusieChrom(list(csd_X = c(1, 2, 3)))))  # unnamed
  expect_true(is.na(pecotmr:::.fsusieChrom(list())))                    # no csd_X
})

test_that(".fsusieCsMapFromTopLoci falls back to index labels when topLoci lacks CS columns", {
  fit <- list(cs = list(1L, 2L),
              csd_X = stats::setNames(1:3, paste0("chr1:", 1:3 * 100, ":A:G")))
  m <- pecotmr:::.fsusieCsMapFromTopLoci(fit, topLoci = NULL)
  expect_equal(unname(m$label), c("fsusie_1", "fsusie_2"))
  expect_true(all(is.na(m$purity)))
  # a topLoci missing the required cs_95 / cs_95_purity columns takes the same path
  m2 <- pecotmr:::.fsusieCsMapFromTopLoci(fit, topLoci = data.frame(variant_id = "x"))
  expect_equal(unname(m2$label), c("fsusie_1", "fsusie_2"))
})

# A minimal object that passes .isFsusieFit (the wavelet-slot gate) so the band /
# affected-region degradation paths can be exercised with the (untrimmed) band
# computation replaced by a mock returning a controlled, degenerate fit.
.fsa_fakeFit <- function()
  list(fitted_wc2 = 1, fitted_func = list(1), outing_grid = 1:2, alpha = 1)
.fsa_bareEntry <- function()
  FineMappingEntry(variantIds = "chr1:100:A:G", susieFit = .fsa_fakeFit(),
                   topLoci = data.frame(variant_id = "chr1:100:A:G", pip = 0.5,
                                        stringsAsFactors = FALSE))

test_that("fsusieCredibleBand skips NULL effect bands and degrades to empty", {
  skip_if_not_installed("fsusieR")
  testthat::local_mocked_bindings(
    .fsusiePopulateCredibleBand = function(fit) list(
      outing_grid = c(100, 200),
      csd_X = stats::setNames(1:2, c("chr1:100:A:G", "chr1:200:A:G")),
      cred_band = list(NULL, NULL), fitted_func = list(NULL, NULL)),
    .package = "pecotmr")
  expect_equal(nrow(fsusieCredibleBand(.fsa_bareEntry())), 0L)
})

test_that("fsusieAffectedRegions degrades to empty when affected_reg finds nothing", {
  skip_if_not_installed("fsusieR")
  testthat::local_mocked_bindings(
    .fsusiePopulateCredibleBand = function(fit) list(
      outing_grid = c(100, 200),
      csd_X = stats::setNames(1:2, c("chr1:100:A:G", "chr1:200:A:G")),
      fitted_func = list(c(1, 1))),
    .package = "pecotmr")
  testthat::local_mocked_bindings(affected_reg = function(...) NULL, .package = "fsusieR")
  expect_equal(length(fsusieAffectedRegions(.fsa_bareEntry())), 0L)
})

test_that("fsusieAffectedRegions yields NA direction for a region outside the grid", {
  skip_if_not_installed("fsusieR")
  testthat::local_mocked_bindings(
    .fsusiePopulateCredibleBand = function(fit) list(
      outing_grid = c(100, 200),
      csd_X = stats::setNames(1:2, c("chr1:100:A:G", "chr1:200:A:G")),
      fitted_func = list(c(1, 1))),
    .package = "pecotmr")
  # Region [500, 600] lies outside the grid [100, 200]: no grid points fall in it,
  # so the per-region effect-direction summary is NA.
  testthat::local_mocked_bindings(
    affected_reg = function(...) data.frame(CS = 1L, Start = 500, End = 600),
    .package = "fsusieR")
  gr <- fsusieAffectedRegions(.fsa_bareEntry())
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 1L)
  expect_true(is.na(S4Vectors::mcols(gr)$direction))
})
