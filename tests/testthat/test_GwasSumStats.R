# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (GwasSumStats) ===

test_that("GwasSumStats(df) errors when required mcols are missing", {
  df <- data.frame(SNP = "rs1", CHR = "1", BP = 100,
                   A1 = "A",
                   stringsAsFactors = FALSE)
  expect_error(
    makeGwasSumStatsFromDf(df),
    "Missing required columns"
  )
})


test_that("GwasSumStats valid object passes with all required mcols", {
  set.seed(1)
  df <- data.frame(
    SNP = paste0("rs", 1:5),
    CHR = rep("1", 5),
    BP = 1:5,
    A1 = rep("A", 5),
    A2 = rep("G", 5),
    Z = rnorm(5),
    N = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  obj <- makeGwasSumStatsFromDf(df)
  expect_true(methods::validObject(obj))
})


test_that("makeGwasSumStatsFromDf() constructor creates object from data.frame", {
  df <- make_test_sumstats_df(20)
  obj <- makeGwasSumStatsFromDf(df, traitName = "height", genome = "hg38")

  expect_s4_class(obj, "GwasSumStats")
  expect_equal(as.character(obj$study)[[1L]], "height")
  expect_equal(getGenome(obj), "hg38")
  expect_equal(length(getSumStats(obj)), 20)
})


test_that("makeGwasSumStatsFromDf() normalizes chr prefix", {
  df <- make_test_sumstats_df(5)
  # Input has CHR = "1" (no prefix)
  obj <- makeGwasSumStatsFromDf(df)
  chrs <- as.character(GenomicRanges::seqnames(getSumStats(obj)))
  expect_true(all(startsWith(chrs, "chr")))

  # Input already has "chr" prefix
  df2 <- df
  df2$CHR <- "chr1"
  obj2 <- makeGwasSumStatsFromDf(df2)
  chrs2 <- as.character(GenomicRanges::seqnames(getSumStats(obj2)))
  # Should not double-prefix
  expect_true(all(chrs2 == "chr1"))
  expect_false(any(grepl("^chrchr", chrs2)))
})


test_that("makeGwasSumStatsFromDf() errors on missing columns", {
  df <- data.frame(SNP = "rs1", CHR = "1", BP = 100)
  expect_error(makeGwasSumStatsFromDf(df), "Missing required columns")
})


test_that("makeGwasSumStatsFromDf() removes rows with NA in required columns", {
  df <- make_test_sumstats_df(10)
  df$Z[1] <- NA
  df$N[3] <- NA
  expect_message(obj <- makeGwasSumStatsFromDf(df), "Removed.*SNPs with missing")
  expect_equal(length(getSumStats(obj)), 8)
})


test_that("getz() returns correct Z vector", {
  set.seed(99)
  df <- make_test_sumstats_df(5)
  obj <- makeGwasSumStatsFromDf(df)
  z <- getZ(obj)
  expect_type(z, "double")
  expect_equal(length(z), 5)
})


test_that("getn() returns correct N vector", {
  df <- make_test_sumstats_df(5)
  obj <- makeGwasSumStatsFromDf(df)
  n <- getN(obj)
  expect_equal(length(n), 5)
  expect_true(all(n == 10000))
})


test_that("getmaf() returns MAF when present, NULL when absent", {
  df <- make_test_sumstats_df(5)
  obj_no_maf <- makeGwasSumStatsFromDf(df)
  expect_null(getMaf(obj_no_maf))

  df$MAF <- runif(5, 0.01, 0.5)
  obj_with_maf <- makeGwasSumStatsFromDf(df)
  maf <- getMaf(obj_with_maf)
  expect_type(maf, "double")
  expect_equal(length(maf), 5)
})


test_that("nSnps() returns correct count", {
  df <- make_test_sumstats_df(30)
  obj <- makeGwasSumStatsFromDf(df)
  expect_equal(nSnps(obj), 30)
})


test_that("subsetchr() filters correctly", {
  df <- make_test_sumstats_df(10)
  df$CHR <- c(rep("1", 6), rep("2", 4))
  obj <- makeGwasSumStatsFromDf(df)

  chr1 <- subsetChr(obj, "1")
  expect_equal(nSnps(chr1), 6)

  # Also works with "chr" prefix
  chr2 <- subsetChr(obj, "chr2")
  expect_equal(nSnps(chr2), 4)
})


test_that("getvary() returns var_y and NULL cases", {
  df <- make_test_sumstats_df(5)

  obj_null <- makeGwasSumStatsFromDf(df, varY = NULL)
  expect_null(getVarY(obj_null))

  obj_vy <- makeGwasSumStatsFromDf(df, varY = 4.5)
  expect_equal(getVarY(obj_vy), 4.5)
})


test_that("as.data.frame.makeGwasSumStatsFromDf() round-trips", {
  df_in <- make_test_sumstats_df(15)
  obj <- makeGwasSumStatsFromDf(df_in)
  df_out <- as.data.frame(obj)

  expect_true(is.data.frame(df_out))
  expect_true(all(c("SNP", "CHR", "BP", "A1", "A2", "Z", "N") %in%
                    names(df_out)))
  expect_equal(nrow(df_out), 15)
  expect_equal(df_out$SNP, df_in$SNP)
  # BP should round-trip
  expect_equal(df_out$BP, as.integer(df_in$BP))
})


# =============================================================================
# AnnotationMatrix (h2Annotations.R)
# =============================================================================



# === Tests migrated from test_showMethods.R (GwasSumStats) ===

test_that("show.GwasSumStats prints nrow and genome build", {
  ss <- GwasSumStats(
    study = c("g1", "g2"),
    entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
    genome = "hg19",
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(ss))
  expect_true(any(grepl("GwasSumStats: 2 studies, genome build hg19", out)))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})

# === Tests migrated from test_h2ClassesSumstats.R (showMethods) ===


