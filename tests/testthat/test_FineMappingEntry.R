# Tests for FineMappingEntry (S4 class)

# Helper for adjustPips tests, migrated from test_dataStructures.R

.makeAdjustEntry <- function(vids, L = 2L) {
  p <- length(vids)
  set.seed(11L)
  lbf <- matrix(rnorm(L * p), nrow = L, ncol = p)
  colnames(lbf) <- vids
  alpha <- lbfToAlpha(lbf)
  pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
  FineMappingEntry(
    variantIds = vids,
    susieFit = list(
      pip          = pip,
      alpha        = alpha,
      lbf_variable = lbf,
      mu           = matrix(0, L, p),
      X_column_scale_factors = rep(1, p)
    ),
    topLoci = data.frame(
      variant_id = vids,
      pip        = pip,
      betahat    = rep(0, p),
      sebetahat  = rep(1, p),
      stringsAsFactors = FALSE
    )
  )
}


# ===========================================================================
# Tests migrated from test_dataStructures.R (getTopLoci, adjustPips)
# ===========================================================================

test_that("getTopLoci(type='GRanges') converts topLoci data.frame to GRanges", {
  tl <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.9, 0.1),
    betahat = c(0.5, -0.2),
    sebetahat = c(0.1, 0.2),
    cs = c(1L, 0L),
    method = "susie",
    stringsAsFactors = FALSE
  )
  ent <- FineMappingEntry(variantIds = tl$variant_id,
                          susieFit = list(), topLoci = tl)
  gr <- getTopLoci(ent, type = "GRanges")
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(S4Vectors::mcols(gr)$pip, c(0.9, 0.1))
})


test_that("getTopLoci(type='GRanges') handles empty input", {
  ent <- FineMappingEntry(variantIds = character(0),
                          susieFit = list(),
                          topLoci = data.frame())
  gr <- getTopLoci(ent, type = "GRanges")
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 0)
})


test_that("getTopLoci defaults to data.frame", {
  tl <- data.frame(
    variant_id = "1:100:A:G", pip = 0.9,
    betahat = 0.5, sebetahat = 0.1, cs = 1L,
    stringsAsFactors = FALSE
  )
  ent <- FineMappingEntry(variantIds = tl$variant_id,
                          susieFit = list(), topLoci = tl)
  expect_s3_class(getTopLoci(ent), "data.frame")
})

# =============================================================================
# extractBlockGenotypes returns RSE
# =============================================================================


test_that("adjustPips renormalizes PIPs on a kept FineMappingEntry subset", {
  vids <- paste0("chr1:", 1:6, ":A:G")
  entry <- .makeAdjustEntry(vids)
  keep <- vids[2:5]
  adj <- adjustPips(entry, keep)
  expect_s4_class(adj, "FineMappingEntry")
  expect_equal(adj@variantIds, keep)
  expect_equal(ncol(adj@susieFit$lbf_variable), 4)
  # Renormalized: each effect's alpha row sums to 1 (when row has any signal)
  expect_true(all(abs(rowSums(adj@susieFit$alpha) - 1) < 1e-10))
  # PIPs match topLoci
  expect_equal(adj@topLoci$pip, adj@susieFit$pip)
  # PIPs change under renormalization
  origPips <- getPip(entry)
  expect_false(identical(unname(origPips[keep]), adj@susieFit$pip))
})


test_that("adjustPips errors when the intersection is empty", {
  vids <- paste0("chr1:", 1:4, ":A:G")
  entry <- .makeAdjustEntry(vids)
  expect_error(
    adjustPips(entry, paste0("chr2:", 1:4, ":A:G")),
    "intersection.*empty"
  )
})


test_that("adjustPips tolerates a chr-prefix difference between entry and keepVariants", {
  vids <- paste0("chr1:", 1:6, ":A:G")
  entry <- .makeAdjustEntry(vids)
  keep <- paste0("1:", 2:5, ":A:G")            # same variants, no "chr" prefix
  adj <- adjustPips(entry, keep)
  expect_s4_class(adj, "FineMappingEntry")
  expect_equal(adj@variantIds, vids[2:5])      # entry keeps its own labels
  expect_equal(ncol(adj@susieFit$lbf_variable), 4)
})


test_that("adjustPips on a FineMappingResultBase collection renormalizes each entry", {
  vidsA <- paste0("chr1:", 1:6, ":A:G")
  vidsB <- paste0("chr1:", 3:8, ":A:G")
  entryA <- .makeAdjustEntry(vidsA)
  entryB <- .makeAdjustEntry(vidsB)
  fmr <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("g1", "g1"),
    method  = c("susie", "susie"),
    entry   = list(entryA, entryB))
  # Keep only variants shared by both entries' raw sets.
  keep <- intersect(vidsA, vidsB)
  adj <- adjustPips(fmr, keep)
  expect_s4_class(adj, "QtlFineMappingResult")
  expect_equal(nrow(adj), 2L)
  expect_equal(adj@listData$entry[[1L]]@variantIds, keep)
  expect_equal(adj@listData$entry[[2L]]@variantIds, keep)
})


# === Tests migrated from test_s4Constructors.R (FineMappingEntry) ===

test_that("FineMappingEntry: constructor stores slots and accessors return them", {
  tl <- .sc_makeTopLoci(3)
  tl$variant_id <- c("a", "b", "c")
  entry <- FineMappingEntry(
    variantIds = c("a", "b", "c"),
    susieFit   = list(payload = 1L),
    topLoci    = tl)
  expect_s4_class(entry, "FineMappingEntry")
  expect_equal(getVariantIds(entry), c("a", "b", "c"))
  expect_equal(getSusieFit(entry), list(payload = 1L))
  # getTopLoci returns the projected posterior view, not the raw slot
  out <- getTopLoci(entry, signalCutoff = 0)
  expect_equal(out$variant_id, c("a", "b", "c"))
})


test_that("FineMappingEntry: getPip returns named pip vector keyed by variant_id", {
  entry <- .sc_makeFineMappingEntry(3)
  pip <- getPip(entry)
  expect_equal(length(pip), 3L)
  expect_equal(names(pip),
               paste0("chr1:", 100 * 1:3, ":A:G"))
})

test_that("FineMappingEntry: resolveWeights returns topLoci posterior effect aligned to variant_id", {
  entry <- .sc_makeFineMappingEntry(3)   # topLoci posterior_mean = 0.05 for all
  wr <- resolveWeights(entry)
  expect_equal(wr$variantIds, paste0("chr1:", 100 * 1:3, ":A:G"))
  expect_equal(wr$weights, rep(0.05, 3))
  expect_equal(length(wr$variantIds), length(wr$weights))
  # empty topLoci -> empty pair
  empty <- FineMappingEntry(variantIds = character(0), susieFit = list(),
    topLoci = data.frame(variant_id = character(0), pip = numeric(0),
                         stringsAsFactors = FALSE))
  expect_length(resolveWeights(empty)$variantIds, 0L)
})


test_that("FineMappingEntry: getPip returns numeric(0) when topLoci is empty", {
  entry <- FineMappingEntry(
    variantIds = character(0),
    susieFit = list(),
    topLoci    = data.frame(variant_id = character(0), pip = numeric(0),
                            stringsAsFactors = FALSE))
  expect_equal(getPip(entry), numeric(0))
})


test_that("FineMappingEntry: getCs filters to rows in any credible set", {
  entry <- .sc_makeFineMappingEntry(3)  # last row has cs_95 = "susie_0"
  res <- getCs(entry)
  expect_equal(nrow(res), 2L)
})


test_that("FineMappingEntry: getCs/getTopLoci surface the directional af column", {
  # Regression: the posterior view must carry the topLoci `af` (effect-allele
  # frequency) through to getCs() / getTopLoci(), not drop it to NA. The
  # value is directional (0.87 > 0.5), so a folded MAF would be a bug.
  entry <- .sc_makeFineMappingEntry(3)
  cs <- getCs(entry)
  expect_true("af" %in% names(cs))
  expect_equal(unname(cs$af), rep(0.87, nrow(cs)))
  expect_false(anyNA(cs$af))

  tl <- getTopLoci(entry, signalCutoff = 0)
  expect_true("af" %in% names(tl))
  expect_equal(unname(tl$af), rep(0.87, nrow(tl)))
})


test_that("FineMappingEntry: validity errors when topLoci is missing required cols", {
  expect_error(
    FineMappingEntry(
      variantIds = "v1",
      susieFit   = list(),
      topLoci    = data.frame(other = 1, stringsAsFactors = FALSE)),
    "topLoci missing required columns"
  )
})


# === cvResult slot (cross-validation payload) ===

test_that("FineMappingEntry stores and returns a cvResult payload", {
  tl <- data.frame(variant_id = c("v1", "v2"), pip = c(0.8, 0.2),
                   stringsAsFactors = FALSE)
  cv <- list(samplePartition = data.frame(Sample = c("s1", "s2"), Fold = c(1L, 2L)),
             prediction = list(susie_predicted = matrix(0, 2, 1)),
             performance = list(susie_performance = matrix(0, 1, 6)))
  e <- FineMappingEntry(variantIds = tl$variant_id, susieFit = list(),
                        topLoci = tl, cvResult = cv)
  expect_identical(getCvResult(e), cv)
})

test_that("FineMappingEntry cvResult defaults to NULL and rejects non-list", {
  tl <- data.frame(variant_id = "v1", pip = 0.5, stringsAsFactors = FALSE)
  e <- FineMappingEntry(variantIds = "v1", susieFit = list(), topLoci = tl)
  expect_null(getCvResult(e))
  expect_error(
    FineMappingEntry(variantIds = "v1", susieFit = list(), topLoci = tl,
                     cvResult = 1:3),
    "cvResult must be NULL or a list")
})

# ===========================================================================
# TwasWeightsEntry
# ===========================================================================


test_that("QtlFineMappingResult: validity rejects non-FineMappingEntry rows", {
  expect_error(
    QtlFineMappingResult(
      study = "s1", context = "c1", trait = "t1", method = "susie",
      entry = list("not_an_entry")),
    "every element of the `entry` column must be a FineMappingEntry"
  )
})



# === Tests migrated from test_showMethods.R (FineMappingEntry) ===

test_that("show.FineMappingEntry reports variant count and CS count", {
  e_with_cs <- .sh_makeFmEntry(n = 3, with_cs = TRUE)  # 2 distinct cs > 0
  out <- capture.output(show(e_with_cs))
  expect_true(any(grepl("FineMappingEntry: 3 variants.*1 credible sets", out)))

  # No cs column -> 0 credible sets reported.
  tl <- data.frame(variant_id = c("a", "b"), pip = c(0.1, 0.2),
                   stringsAsFactors = FALSE)
  e_no_cs <- FineMappingEntry(variantIds = c("a", "b"),
                              susieFit = list(), topLoci = tl)
  out_no <- capture.output(show(e_no_cs))
  expect_true(any(grepl("0 credible sets", out_no)))
})


# === getMarginalEffects maxPval filter ===

test_that("FineMappingEntry: getMarginalEffects applies the maxPval filter", {
  tl <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    pip        = c(0.9, 0.5, 0.1),
    marginal_p = c(0.001, 0.5, NA_real_),
    stringsAsFactors = FALSE)
  entry <- FineMappingEntry(variantIds = tl$variant_id,
                            susieFit = list(), topLoci = tl)
  out <- getMarginalEffects(entry, maxPval = 0.01)
  # Drops the p = 0.5 row and the NA-p row; keeps only v1.
  expect_equal(nrow(out), 1L)
  expect_equal(out$variant_id, "v1")
})


# === getCs empty / cs-less topLoci projections ===

test_that("FineMappingEntry: getCs returns empty posterior view when topLoci is empty", {
  entry <- FineMappingEntry(
    variantIds = character(0),
    susieFit = list(),
    topLoci = data.frame(variant_id = character(0), pip = numeric(0),
                         stringsAsFactors = FALSE))
  res <- getCs(entry)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
  expect_true(all(c("variant_id", "pip") %in% names(res)))
})


test_that("FineMappingEntry: getCs returns empty posterior view when the cs column is absent", {
  tl <- data.frame(variant_id = c("a", "b"), pip = c(0.1, 0.2),
                   stringsAsFactors = FALSE)
  entry <- FineMappingEntry(variantIds = c("a", "b"),
                            susieFit = list(), topLoci = tl)
  res <- getCs(entry)  # no cs_95 column -> empty view
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})


# === adjustPips: missing lbf_variable + mu2-driven posterior recompute ===

test_that("FineMappingEntry: adjustPips errors when susieFit lacks lbf_variable", {
  tl <- data.frame(variant_id = c("v1", "v2"), pip = c(0.6, 0.4),
                   stringsAsFactors = FALSE)
  entry <- FineMappingEntry(variantIds = c("v1", "v2"),
                            susieFit = list(), topLoci = tl)
  expect_error(adjustPips(entry, c("v1", "v2")),
               "no `lbf_variable` matrix")
})


test_that("FineMappingEntry: adjustPips subsets mu2 and recomputes posterior_sd", {
  vids <- paste0("chr1:", 1:5, ":A:G")
  p <- length(vids); L <- 2L
  set.seed(101L)
  lbf <- matrix(rnorm(L * p), nrow = L, ncol = p)
  colnames(lbf) <- vids
  alpha <- lbfToAlpha(lbf)
  pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
  mu  <- matrix(rnorm(L * p), L, p)
  mu2 <- mu^2 + 1                      # plausible second moment (>= mean^2)
  entry <- FineMappingEntry(
    variantIds = vids,
    susieFit = list(pip = pip, alpha = alpha, lbf_variable = lbf,
                    mu = mu, mu2 = mu2,
                    X_column_scale_factors = rep(1, p)),
    topLoci = data.frame(variant_id = vids, pip = pip,
                         stringsAsFactors = FALSE))
  keep <- vids[2:4]
  adj <- adjustPips(entry, keep)
  expect_s4_class(adj, "FineMappingEntry")
  # mu2 carried through the variant subsetting alongside lbf/mu.
  expect_equal(ncol(adj@susieFit$mu2), 3L)
  # posterior_mean / posterior_sd recomputed from the subset alpha/mu/mu2.
  expect_equal(nrow(adj@topLoci), 3L)
  expect_true("posterior_sd" %in% names(adj@topLoci))
  expect_equal(length(adj@topLoci$posterior_sd), 3L)
  expect_true(all(adj@topLoci$posterior_sd >= 0))
})

test_that("getCs / getTopLoci minPurity filters CS by purity, independent of coverage + pip", {
  # Two 0.95 credible sets: CS1 (v1,v2) pure (0.9), CS2 (v3,v4) impure (0.3); v5 non-CS.
  vn <- paste0("chr1:", (1:5) * 100, ":A:G")
  L <- 2L; P <- 5L
  alpha <- matrix(0.1, L, P); alpha[1, 1:2] <- c(0.5, 0.4); alpha[2, 3:4] <- c(0.5, 0.4)
  fit <- list(alpha = alpha, mu = matrix(0.3, L, P), mu2 = matrix(1.2, L, P),
              pip = c(0.6, 0.5, 0.4, 0.4, 0.02))
  class(fit) <- "susie"
  cst <- list(list(sets = list(cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
                               purity = data.frame(min.abs.corr = c(0.9, 0.3)))),
              list(sets = list(cs = list())), list(sets = list(cs = list())))
  attr(cst, "coverage") <- c(0.95, 0.70, 0.50)
  tl <- buildTopLoci(fit, cst, variantNames = vn, method = "susie")
  e  <- FineMappingEntry(variantIds = vn, susieFit = fit, topLoci = tl)

  # getCs: minPurity is orthogonal to coverage -> keeps only the pure CS members
  expect_equal(nrow(getCs(e, coverage = 0.95)), 4L)
  expect_equal(getCs(e, coverage = 0.95, minPurity = 0.8)$variant_id, vn[1:2])
  # getTopLoci: minPurity is orthogonal to the pip cutoff -> drops impure-CS
  # variants (v3,v4), keeps pure-CS (v1,v2) and non-CS (v5)
  expect_setequal(getTopLoci(e, signalCutoff = 0, minPurity = 0.8)$variant_id, vn[c(1, 2, 5)])
  expect_equal(nrow(getTopLoci(e, signalCutoff = 0)), 5L)  # default (NULL) unchanged
})


