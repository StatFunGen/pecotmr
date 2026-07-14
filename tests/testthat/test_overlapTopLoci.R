context("overlapTopLoci")

# Minimal top-loci frame (the FineMappingEntry slot schema).
.ot_tl <- function(vids, A1, A2, pm, af, cs = "susie_1") {
  n <- length(vids)
  p <- pecotmr:::parseVariantId(vids)
  data.frame(
    variant_id = vids, chrom = p$chrom, pos = as.integer(p$pos),
    A1 = A1, A2 = A2, N = rep(1000, n), af = af,
    marginal_beta = NA_real_, marginal_se = NA_real_,
    marginal_z = NA_real_, marginal_p = NA_real_,
    pip = rep(0.9, n), posterior_mean = pm, posterior_sd = rep(0.02, n),
    cs_95 = rep(cs, n), cs_95_purity = rep(0.9, n),
    stringsAsFactors = FALSE)
}
.ot_entry <- function(tl) FineMappingEntry(variantIds = tl$variant_id,
                                           susieFit = list(fake = TRUE), topLoci = tl)

test_that("overlapTopLoci joins QTL x GWAS with allele-aware matching + sign-flip", {
  # QTL: two contexts, variants v1 (chr1:100:A:G) and v2 (chr1:200:A:G).
  q1 <- .ot_entry(.ot_tl(c("chr1:100:A:G", "chr1:200:A:G"),
                         A1 = c("G", "G"), A2 = c("A", "A"), pm = c(0.30, 0.40),
                         af = c(0.3, 0.4)))
  q2 <- .ot_entry(.ot_tl(c("chr1:100:A:G", "chr1:200:A:G"),
                         A1 = c("G", "G"), A2 = c("A", "A"), pm = c(0.20, 0.25),
                         af = c(0.3, 0.4)))
  qtl <- QtlFineMappingResult(study = c("s", "s"), context = c("ctx1", "ctx2"),
                              trait = c("t", "t"), method = c("susie", "susie"),
                              entry = list(q1, q2))
  # GWAS: v1 same orientation; v2 REF/ALT-SWAPPED (chr1:200:G:A).
  g1 <- .ot_entry(.ot_tl(c("chr1:100:A:G", "chr1:200:G:A"),
                         A1 = c("G", "A"), A2 = c("A", "G"), pm = c(0.60, 0.50),
                         af = c(0.2, 0.7)))
  gwas <- GwasFineMappingResult(study = "gwas", method = "susie", entry = list(g1))

  ov <- overlapTopLoci(qtl, gwas, signalCutoff = 0)

  # variant key kept once; every other column prefixed qtl_ / gwas_
  expect_true(all(c("variant_id", "chrom", "pos", "A1", "A2") %in% names(ov)))
  expect_true(any(grepl("^qtl_", names(ov))) && any(grepl("^gwas_", names(ov))))
  expect_false(any(grepl("^qtl_(chrom|pos|A1|A2)$", names(ov))))
  # 2 shared variants x 2 QTL contexts x 1 GWAS study = 4 rows (wide cross-product)
  expect_equal(nrow(ov), 4L)
  expect_setequal(unique(ov$variant_id), c("chr1:100:A:G", "chr1:200:A:G"))
  # per-context QTL rows preserved
  expect_setequal(unique(ov$qtl_context), c("ctx1", "ctx2"))

  # v1 unflipped: gwas beta unchanged (+0.6); v2 swapped: sign-flipped (-0.5)
  v1 <- ov[ov$variant_id == "chr1:100:A:G", ]
  v2 <- ov[ov$variant_id == "chr1:200:A:G", ]
  expect_true(all(abs(v1$gwas_beta - 0.60) < 1e-9))
  expect_true(all(abs(v2$gwas_beta + 0.50) < 1e-9))
  # af complemented on the swapped variant (0.7 -> 0.3), untouched otherwise
  expect_true(all(abs(v2$gwas_af - 0.30) < 1e-9))
  expect_true(all(abs(v1$gwas_af - 0.20) < 1e-9))
})

test_that("overlapTopLoci returns zero rows when no variants overlap", {
  q <- .ot_entry(.ot_tl("chr1:100:A:G", A1 = "G", A2 = "A", pm = 0.3, af = 0.3))
  g <- .ot_entry(.ot_tl("chr2:999:C:T", A1 = "T", A2 = "C", pm = 0.6, af = 0.2))
  qtl  <- QtlFineMappingResult(study = "s", context = "c", trait = "t",
                               method = "susie", entry = list(q))
  gwas <- GwasFineMappingResult(study = "g", method = "susie", entry = list(g))
  ov <- overlapTopLoci(qtl, gwas, signalCutoff = 0)
  expect_equal(nrow(ov), 0L)
  expect_true(any(grepl("^qtl_", names(ov))) && any(grepl("^gwas_", names(ov))))
})

test_that("overlapTopLoci returns an empty merge when a side has no signal (df + GRanges)", {
  q <- .ot_entry(.ot_tl("chr1:100:A:G", A1 = "G", A2 = "A", pm = 0.3, af = 0.3))
  g <- .ot_entry(.ot_tl("chr1:100:A:G", A1 = "G", A2 = "A", pm = 0.6, af = 0.2))
  qtl  <- QtlFineMappingResult(study = "s", context = "c", trait = "t",
                               method = "susie", entry = list(q))
  gwas <- GwasFineMappingResult(study = "g", method = "susie", entry = list(g))
  # signalCutoff above every pip (0.9) -> both top-loci tables empty -> emptyMerge.
  df <- overlapTopLoci(qtl, gwas, signalCutoff = 0.99)
  expect_equal(nrow(df), 0L)
  expect_true(any(grepl("^qtl_", names(df))) && any(grepl("^gwas_", names(df))))
  # same empty result routed through the GRanges branch (.overlapToGRanges empty guard).
  gr <- overlapTopLoci(qtl, gwas, signalCutoff = 0.99, type = "GRanges")
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 0L)
})

test_that("overlapTopLoci type='GRanges' returns a GRanges of the shared variants", {
  q <- .ot_entry(.ot_tl(c("chr1:100:A:G", "chr1:200:A:G"),
                        A1 = c("G", "G"), A2 = c("A", "A"), pm = c(0.3, 0.4),
                        af = c(0.3, 0.4)))
  g <- .ot_entry(.ot_tl("chr1:100:A:G", A1 = "G", A2 = "A", pm = 0.6, af = 0.2))
  qtl  <- QtlFineMappingResult(study = "s", context = "c", trait = "t",
                               method = "susie", entry = list(q))
  gwas <- GwasFineMappingResult(study = "g", method = "susie", entry = list(g))
  gr <- overlapTopLoci(qtl, gwas, signalCutoff = 0, type = "GRanges")
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 1L)
  expect_true("gwas_beta" %in% names(S4Vectors::mcols(gr)))
})
