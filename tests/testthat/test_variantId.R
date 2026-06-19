context("variantId")

# ===========================================================================
# parseVariantId
# ===========================================================================

test_that("parseVariantId: canonical colon format with chr prefix", {
  res <- parseVariantId(c("chr1:100:A:G", "chr2:200:T:C"))
  expect_s3_class(res, "data.frame")
  expect_equal(res$chrom, c(1L, 2L))
  expect_equal(res$pos, c(100L, 200L))
  expect_equal(res$A2, c("A", "T"))
  expect_equal(res$A1, c("G", "C"))
})

test_that("parseVariantId: canonical colon format without chr prefix", {
  res <- parseVariantId(c("1:100:A:G", "2:200:T:C"))
  expect_equal(res$chrom, c(1L, 2L))
  expect_equal(res$pos, c(100L, 200L))
  expect_equal(attr(res, "convention")$hasChr, FALSE)
})

test_that("parseVariantId: underscore separator format", {
  res <- parseVariantId(c("chr1_100_A_G", "1_200_T_C"))
  expect_equal(res$chrom, c(1L, 1L))
  expect_equal(res$pos, c(100L, 200L))
  expect_equal(res$A2, c("A", "T"))
  expect_equal(res$A1, c("G", "C"))
})

test_that("parseVariantId: mixed colon-underscore format", {
  res <- parseVariantId("chr1:100_A_G")
  expect_equal(res$chrom, 1L)
  expect_equal(res$pos, 100L)
  expect_equal(res$A2, "A")
  expect_equal(res$A1, "G")
  expect_equal(attr(res, "convention")$alleleSep, "_")
})

test_that("parseVariantId: strips build suffix", {
  res <- parseVariantId(c("chr1:100:A:G:b38", "chr1:200:T:C_b37"))
  expect_equal(res$pos, c(100L, 200L))
  expect_equal(res$A1, c("G", "C"))
  expect_equal(attr(res, "convention")$hasBuild, TRUE)
})

test_that("parseVariantId: convention attribute records hasChr", {
  with_chr <- parseVariantId("chr1:100:A:G")
  no_chr   <- parseVariantId("1:100:A:G")
  expect_equal(attr(with_chr, "convention")$hasChr, TRUE)
  expect_equal(attr(no_chr, "convention")$hasChr, FALSE)
})

test_that("parseVariantId: data.frame input with chrom/pos/A2/A1 returns as-is", {
  df <- data.frame(chrom = "chr1", pos = 100L, A2 = "A", A1 = "G",
                   stringsAsFactors = FALSE)
  res <- parseVariantId(df)
  expect_equal(res$chrom, 1L)
  expect_equal(res$pos, 100L)
  expect_equal(attr(res, "convention")$hasChr, TRUE)
})

test_that("parseVariantId: data.frame input with >=4 columns assigns positional names", {
  df <- data.frame(c1 = "1", c2 = "100", c3 = "A", c4 = "G", extra = "x",
                   stringsAsFactors = FALSE)
  res <- parseVariantId(df)
  expect_equal(res$chrom, 1L)
  expect_equal(res$pos, 100L)
  expect_equal(res$A2, "A")
  expect_equal(res$A1, "G")
})

# ===========================================================================
# normalizeVariantId
# ===========================================================================

test_that("normalizeVariantId: canonical output adds chr prefix by default", {
  res <- normalizeVariantId(c("1:100:A:G", "2:200:T:C"))
  expect_equal(res, c("chr1:100:A:G", "chr2:200:T:C"))
})

test_that("normalizeVariantId: chrPrefix = FALSE strips the chr prefix", {
  res <- normalizeVariantId(c("chr1:100:A:G", "chr2:200:T:C"),
                            chrPrefix = FALSE)
  expect_equal(res, c("1:100:A:G", "2:200:T:C"))
})

test_that("normalizeVariantId: convention argument preserves the input format", {
  ids <- c("chr1:100_A_G", "chr2:200_T_C")
  conv <- attr(parseVariantId(ids), "convention")
  res <- normalizeVariantId(c("1:100:A:G", "2:200:T:C"), convention = conv)
  expect_equal(res, c("chr1:100_A_G", "chr2:200_T_C"))
})

test_that("normalizeVariantId: round-trips underscore-only input", {
  res <- normalizeVariantId("1_100_A_G")
  expect_equal(res, "chr1:100:A:G")
})

test_that("normalizeVariantId: strips build suffix and re-emits canonical", {
  res <- normalizeVariantId("chr1:100:A:G:b38")
  expect_equal(res, "chr1:100:A:G")
})

# ===========================================================================
# parseRegion
# ===========================================================================

test_that("parseRegion: canonical chr:start-end string", {
  res <- parseRegion("chr1:100-200")
  expect_s3_class(res, "data.frame")
  expect_equal(res$chrom, "1")
  expect_equal(res$start, 100L)
  expect_equal(res$end, 200L)
})

test_that("parseRegion: handles X chromosome", {
  res <- parseRegion("chrX:1000-2000")
  expect_equal(res$chrom, "X")
  expect_equal(res$start, 1000L)
  expect_equal(res$end, 2000L)
})

test_that("parseRegion: rejects malformed input", {
  expect_error(parseRegion("notARegion"),
               "format must be 'chr:start-end'")
  expect_error(parseRegion("1:100:200"),
               "format must be 'chr:start-end'")
})

test_that("parseRegion: returns non-character input unchanged", {
  df_in <- data.frame(chrom = "1", start = 100L, end = 200L,
                      stringsAsFactors = FALSE)
  expect_identical(parseRegion(df_in), df_in)
  expect_identical(parseRegion(42L), 42L)
})

# ===========================================================================
# regionToDf
# ===========================================================================

test_that("regionToDf: parses chrom_start_end LD region IDs", {
  res <- regionToDf(c("1_100_200", "2_300_400"))
  expect_equal(res$chrom, c(1L, 2L))
  expect_equal(res$start, c(100L, 300L))
  expect_equal(res$end,   c(200L, 400L))
})

test_that("regionToDf: strips chr prefix from chromosome", {
  res <- regionToDf("chr5_1000_2000")
  expect_equal(res$chrom, 5L)
})

test_that("regionToDf: accepts colon/dash separators", {
  res <- regionToDf("chr1:100-200")
  expect_equal(res$chrom, 1L)
  expect_equal(res$start, 100L)
  expect_equal(res$end, 200L)
})

test_that("regionToDf: honours custom column names", {
  res <- regionToDf("1_10_20", colnames = c("seq", "from", "to"))
  expect_equal(colnames(res), c("seq", "from", "to"))
})

# ===========================================================================
# regionsOverlap
# ===========================================================================

test_that("regionsOverlap: TRUE when regions share a base pair", {
  expect_true(regionsOverlap("chr1:100-200", "chr1:150-250"))
  expect_true(regionsOverlap("chr1:100-200", "chr1:200-300"))  # touching
})

test_that("regionsOverlap: FALSE when regions are disjoint", {
  expect_false(regionsOverlap("chr1:100-200", "chr1:300-400"))
})

test_that("regionsOverlap: FALSE across different chromosomes", {
  # Bioconductor warns when comparing GRanges with disjoint seqlevels —
  # the semantics of "no overlap" is exactly what we want here.
  expect_false(suppressWarnings(regionsOverlap("chr1:100-200", "chr2:100-200")))
})

test_that("regionsOverlap: accepts data.frame input", {
  a <- data.frame(chrom = "1", start = 100L, end = 200L,
                  stringsAsFactors = FALSE)
  b <- data.frame(chrom = "1", start = 150L, end = 300L,
                  stringsAsFactors = FALSE)
  expect_true(regionsOverlap(a, b))
})

# ===========================================================================
# findOverlappingRegions
# ===========================================================================

test_that("findOverlappingRegions: returns the indices of overlapping targets", {
  targets <- c("chr1:100-200", "chr1:300-400", "chr1:150-250", "chr2:100-200")
  res <- findOverlappingRegions("chr1:175-225", targets)
  expect_equal(sort(res), c(1L, 3L))
})

test_that("findOverlappingRegions: empty result when no overlap", {
  res <- findOverlappingRegions("chr1:500-600",
                                c("chr1:100-200", "chr1:300-400"))
  expect_equal(res, integer(0))
})

test_that("findOverlappingRegions: deduplicates target hits", {
  res <- findOverlappingRegions("chr1:150-160",
                                c("chr1:100-200", "chr1:300-400"))
  expect_equal(res, 1L)
})

# ===========================================================================
# classifyVariantType
# ===========================================================================

test_that("classifyVariantType: identifies SNPs", {
  expect_equal(
    classifyVariantType(c("chr1:100:A:G", "chr1:200:C:T")),
    c("SNP", "SNP"))
})

test_that("classifyVariantType: identifies insertions and deletions", {
  expect_equal(classifyVariantType("chr1:100:A:ATG"), "insertion")
  expect_equal(classifyVariantType("chr1:100:ATG:A"), "deletion")
})

test_that("classifyVariantType: identifies MNPs (equal-length multi-base)", {
  expect_equal(classifyVariantType("chr1:100:AT:GC"), "MNP")
})

test_that("classifyVariantType: accepts data.frame input with A2/A1 columns", {
  df <- data.frame(A2 = c("A", "AT"), A1 = c("G", "GC"),
                   stringsAsFactors = FALSE)
  expect_equal(classifyVariantType(df), c("SNP", "MNP"))
})

test_that("classifyVariantType: errors when given an unsupported input", {
  expect_error(classifyVariantType(list(a = 1)),
               "character vector of variant IDs or a data.frame")
  expect_error(classifyVariantType(data.frame(foo = 1)),
               "A2 and A1 columns")
})
