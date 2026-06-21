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

# ===========================================================================
# Tests migrated from test_misc.R (variant-id helpers)
# ===========================================================================

test_that("Test formatVariantId",{
    expect_equal(formatVariantId(c(1, 1), c(123, 132), c("G", "A"), c("C", "T")), c("chr1:123:G:C", "chr1:132:A:T"))
})


test_that("formatVariantId uses convention parameter automatically", {
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G"), "chr1:100:A:G")

  conv_mixed <- list(hasChr = TRUE, alleleSep = "_")
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", convention = conv_mixed), "chr1:100_A_G")

  conv_nochr <- list(hasChr = FALSE, alleleSep = "_")
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", convention = conv_nochr), "1:100_A_G")

  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", chrPrefix = FALSE, convention = conv_mixed), "chr1:100_A_G")
})


test_that("formatVariantId constructs canonical IDs", {
  expect_equal(pecotmr:::formatVariantId(c(1, 2), c(100, 200), c("A", "C"), c("G", "T")),
               c("chr1:100:A:G", "chr2:200:C:T"))
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", chrPrefix = FALSE), "1:100:A:G")
  expect_equal(pecotmr:::formatVariantId("chr1", 100, "A", "G"), "chr1:100:A:G")
})

# =============================================================================
# =============================================================================
# pvalHmp
# =============================================================================


test_that("parseRegion parses valid region string", {
  result <- parseRegion("chr1:100-200")
  expect_s3_class(result, "data.frame")
  expect_equal(result$chrom, "1")
  expect_equal(result$start, 100L)
  expect_equal(result$end, 200L)
})


test_that("parseRegion handles X chromosome", {
  result <- parseRegion("chrX:500-1000")
  expect_equal(result$chrom, "X")
  expect_equal(result$start, 500L)
  expect_equal(result$end, 1000L)
})


test_that("parseRegion errors on invalid format", {
  expect_error(parseRegion("1:100-200"), "format must be")
  expect_error(parseRegion("chr1-100-200"), "format must be")
  expect_error(parseRegion("chr1:abc-200"), "format must be")
})


test_that("parseRegion returns non-string input unchanged", {
  df <- data.frame(chrom = 1, start = 100, end = 200)
  result <- parseRegion(df)
  expect_identical(result, df)
})


test_that("parseRegion returns non-single-string input unchanged", {
  input <- c("chr1:100-200", "chr2:300-400")
  result <- parseRegion(input)
  expect_identical(result, input)
})

# =============================================================================
# parseVariantId
# =============================================================================


test_that("parseVariantId parses single variant with chr prefix", {
  result <- parseVariantId("chr1:12345:A:G")
  expect_equal(result$chrom, 1L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_true(conv$hasChr)
  expect_equal(conv$alleleSep, ":")
})


test_that("parseVariantId parses single variant without chr prefix", {
  result <- parseVariantId("5:12345:A:G")
  expect_equal(result$chrom, 5L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_false(conv$hasChr)
})


test_that("parseVariantId parses multiple variants", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T", "chr3:300:G:A")
  result <- parseVariantId(ids)
  expect_equal(nrow(result), 3)
  expect_equal(result$chrom, c(1L, 2L, 3L))
  expect_equal(result$pos, c(100L, 200L, 300L))
  expect_equal(result$A2, c("A", "C", "G"))
  expect_equal(result$A1, c("G", "T", "A"))
})

# =============================================================================
# detectVariantConvention
# =============================================================================


test_that("detectVariantConvention detects chr prefix and allele separators", {
  conv <- pecotmr:::detectVariantConvention(c("chr1:100:A:G", "chr2:200:C:T"))
  expect_true(conv$hasChr)
  expect_equal(conv$alleleSep, ":")
  expect_false(conv$hasBuild)

  conv2 <- pecotmr:::detectVariantConvention(c("1_100_A_G", "2_200_C_T"))
  expect_false(conv2$hasChr)
  expect_equal(conv2$alleleSep, "_")

  conv3 <- pecotmr:::detectVariantConvention(c("chr1:100:A:G:b38"))
  expect_true(conv3$hasBuild)

  conv4 <- pecotmr:::detectVariantConvention(c("chr1:100_A_G"))
  expect_true(conv4$hasChr)
  expect_equal(conv4$alleleSep, "_")

  conv5 <- pecotmr:::detectVariantConvention(c("1:100_A_G"))
  expect_false(conv5$hasChr)
  expect_equal(conv5$alleleSep, "_")
})

# =============================================================================
# normalizeVariantId
# =============================================================================


test_that("normalizeVariantId normalizes various formats", {
  expect_equal(normalizeVariantId("1_100_A_G"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100:A:G"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("1:100:A:G"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100:A:G", chrPrefix = FALSE), "1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100:A:G:b38"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100_A_G"), "chr1:100:A:G")
  conv <- pecotmr:::detectVariantConvention(c("chr1:100_A_G"))
  expect_equal(normalizeVariantId("1:200:C:T", convention = conv), "chr1:200_C_T")
})

# =============================================================================
# variantIdToDf
# =============================================================================


test_that("variantIdToDf handles colon-separated format", {
  ids <- c("1:100:A:G", "2:200:C:T")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
  expect_equal(result$A2, c("A", "C"))
  expect_equal(result$A1, c("G", "T"))
})


test_that("variantIdToDf handles underscore-separated format", {
  ids <- c("1:100_A_G", "2:200_C_T")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$A2, c("A", "C"))
})


test_that("variantIdToDf strips chr prefix", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(result$chrom, c(1L, 2L))
})


test_that("variantIdToDf handles data.frame input with named columns", {
  df <- data.frame(chrom = c("chr1", "2"), pos = c(100, 200),
                   A2 = c("A", "C"), A1 = c("G", "T"))
  suppressWarnings(result <- pecotmr:::variantIdToDf(df))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
})


test_that("variantIdToDf handles 5-part IDs with build suffix", {
  ids <- c("chr1:100:A:G:b38", "chr2:200:T:C")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("chrom", "pos", "A2", "A1"))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$A2, c("A", "T"))
  expect_equal(result$A1, c("G", "C"))
})


test_that("variantIdToDf handles mixed 4/5-part IDs", {
  ids <- c("1:100:A:G", "chr2:200:T:C:b38", "3:300:G:A:b37")
  suppressWarnings(result <- pecotmr:::variantIdToDf(ids))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
  expect_equal(result$A1, c("G", "C", "A"))
})

# =============================================================================
# =============================================================================
# =============================================================================
# computeLd — uncovered branches
# =============================================================================


test_that("detectVariantConvention returns defaults for all-NA input", {
  result <- detectVariantConvention(c(NA, NA, NA))
  expect_false(result$hasChr)
  expect_equal(result$alleleSep, ":")
  expect_false(result$hasBuild)
  expect_true(is.na(result$example))
})

# =============================================================================
# parseVariantId — uncovered line 622
# =============================================================================


test_that("parseVariantId handles data.frame with generic column names", {
  df <- data.frame(
    col1 = c("chr1", "chr2"),
    col2 = c(100, 200),
    col3 = c("A", "T"),
    col4 = c("G", "C"),
    stringsAsFactors = FALSE
  )
  result <- parseVariantId(df)
  expect_equal(names(result)[1:4], c("chrom", "pos", "A2", "A1"))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
  expect_equal(result$A2, c("A", "T"))
  expect_equal(result$A1, c("G", "C"))
})

# =============================================================================

# =============================================================================
# regionsOverlap
# =============================================================================


test_that("regionsOverlap detects overlapping regions on same chromosome", {
  expect_true(regionsOverlap("chr1:100-300", "chr1:200-400"))
})


test_that("regionsOverlap returns FALSE for non-overlapping same-chr regions", {
  expect_false(regionsOverlap("chr1:100-200", "chr1:300-400"))
})


test_that("regionsOverlap returns FALSE for different chromosomes", {
  expect_false(regionsOverlap("chr1:100-300", "chr2:100-300"))
})


test_that("regionsOverlap detects touching boundaries", {
  expect_true(regionsOverlap("chr1:100-200", "chr1:200-300"))
})


test_that("regionsOverlap works with underscore-separated IDs", {
  expect_true(regionsOverlap("1_100_300", "1_200_400"))
  expect_false(regionsOverlap("1_100_200", "2_100_200"))
})


test_that("regionsOverlap works with data.frame input", {
  df_a <- data.frame(chrom = 1, start = 100, end = 300)
  df_b <- data.frame(chrom = 1, start = 200, end = 400)
  expect_true(regionsOverlap(df_a, df_b))
})

# =============================================================================
# findOverlappingRegions
# =============================================================================


test_that("findOverlappingRegions returns correct indices", {
  query <- "chr1:100-300"
  targets <- c("chr1:200-400", "chr2:100-200", "chr1:50-150")
  result <- findOverlappingRegions(query, targets)
  expect_true(1 %in% result)
  expect_true(3 %in% result)
  expect_false(2 %in% result)
})


test_that("findOverlappingRegions returns empty vector for no matches", {
  query <- "chr1:100-200"
  targets <- c("chr2:100-200", "chr3:100-200")
  result <- findOverlappingRegions(query, targets)
  expect_length(result, 0)
})


test_that("findOverlappingRegions works with data.frame targets", {
  query <- "chr1:100-300"
  targets <- data.frame(chrom = c(1, 2, 1), start = c(200, 100, 50), end = c(400, 200, 150))
  result <- findOverlappingRegions(query, targets)
  expect_true(1 %in% result)
  expect_true(3 %in% result)
  expect_false(2 %in% result)
})

# =============================================================================
# classifyVariantType
# =============================================================================


test_that("classifyVariantType identifies SNPs", {
  expect_equal(classifyVariantType("chr1:100:A:G"), "SNP")
})


test_that("classifyVariantType identifies insertions", {
  expect_equal(classifyVariantType("chr1:100:A:ATG"), "insertion")
})


test_that("classifyVariantType identifies deletions", {
  expect_equal(classifyVariantType("chr1:100:ATG:A"), "deletion")
})


test_that("classifyVariantType identifies MNPs", {
  expect_equal(classifyVariantType("chr1:100:AT:GC"), "MNP")
})


test_that("classifyVariantType handles vector input", {
  ids <- c("chr1:100:A:G", "chr1:200:ATG:A", "chr1:300:A:ATG", "chr1:400:AT:GC")
  result <- classifyVariantType(ids)
  expect_equal(result, c("SNP", "deletion", "insertion", "MNP"))
})


test_that("classifyVariantType accepts data.frame input", {
  df <- data.frame(A2 = c("A", "ATG"), A1 = c("G", "A"))
  result <- classifyVariantType(df)
  expect_equal(result, c("SNP", "deletion"))
})

# =============================================================================
# ensureChrMatch
# =============================================================================


test_that("ensureChrMatch returns unchanged when both have chr prefix", {
  idsA <- c("chr1:100:A:G", "chr1:200:C:T")
  idsB <- c("chr1:150:A:G", "chr1:250:C:T")
  result <- pecotmr:::ensureChrMatch(idsA, idsB)
  expect_equal(result$idsA, idsA)
  expect_equal(result$idsB, idsB)
})


test_that("ensureChrMatch normalizes when prefixes mismatch", {
  idsA <- c("chr1:100:A:G", "chr1:200:C:T")
  idsB <- c("1:150:A:G", "1:250:C:T")
  result <- pecotmr:::ensureChrMatch(idsA, idsB)
  expect_true(all(grepl("^chr", result$idsA)))
  expect_true(all(grepl("^chr", result$idsB)))
})


test_that("ensureChrMatch returns unchanged when both lack chr prefix", {
  idsA <- c("1:100:A:G", "1:200:C:T")
  idsB <- c("1:150:A:G", "1:250:C:T")
  result <- pecotmr:::ensureChrMatch(idsA, idsB)
  # Both already match (no prefix), so returned unchanged
  expect_equal(result$idsA, idsA)
  expect_equal(result$idsB, idsB)
})

