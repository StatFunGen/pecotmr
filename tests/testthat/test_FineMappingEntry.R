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


