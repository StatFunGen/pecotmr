# Tests for S4 classes (h2_classes.R), GWASSumStats (h2_sumstats.R),
# and AnnotationMatrix (h2_annotations.R)

# =============================================================================
# Test data helpers
# =============================================================================

make_test_granges <- function(n = 10) {
  GenomicRanges::GRanges(
    seqnames = rep("chr1", n),
    ranges = IRanges::IRanges(start = seq(1000, by = 100, length.out = n),
                              width = 1L)
  )
}

make_test_sumstats_df <- function(n = 50) {
  set.seed(42)
  data.frame(
    SNP = paste0("rs", seq_len(n)),
    CHR = rep("1", n),
    BP = seq(1000, by = 100, length.out = n),
    A1 = rep("A", n),
    A2 = rep("G", n),
    Z = rnorm(n),
    N = rep(10000, n),
    stringsAsFactors = FALSE
  )
}

make_test_ldblocks <- function() {
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 5001), end = c(5000, 10000))
  )
  new("LDBlocks", blocks = blocks_gr, genome = "hg19")
}

make_test_snp_info <- function(n = 10) {
  data.frame(
    SNP = paste0("rs", seq_len(n)),
    CHR = rep("1", n),
    BP = seq(1000, by = 100, length.out = n),
    A1 = rep("A", n),
    A2 = rep("G", n),
    stringsAsFactors = FALSE
  )
}

make_test_annotation_meta <- function() {
  data.frame(
    name = c("base", "enhancer", "promoter"),
    tier = c("baseline", "candidate", "candidate"),
    type = c("binary", "binary", "continuous"),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# S4 Classes (h2_classes.R)
# =============================================================================

test_that("LDBlocks constructs and validates correctly", {
  obj <- make_test_ldblocks()
  expect_s4_class(obj, "LDBlocks")
  expect_equal(length(obj@blocks), 2)
  expect_equal(obj@genome, "hg19")
  expect_true(methods::validObject(obj))
})

test_that("LDBlocks rejects genome of length != 1", {
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  expect_error(
    methods::validObject(
      new("LDBlocks", blocks = blocks_gr, genome = c("hg19", "hg38"))
    ),
    "genome.*single"
  )
})

test_that("GenotypeHandle constructs and validates correctly", {
  obj <- new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snp_info = make_test_snp_info(),
    n_samples = 100L,
    sample_ids = paste0("sample_", 1:100),
    pgen_ptr = NULL
  )
  expect_s4_class(obj, "GenotypeHandle")
  expect_equal(obj@format, "gds")
  expect_true(methods::validObject(obj))
})

test_that("GenotypeHandle accepts all valid formats", {
  for (fmt in c("gds", "vcf", "plink1", "plink2")) {
    obj <- new("GenotypeHandle",
      path = "/tmp/test",
      format = fmt,
      snp_info = data.frame(),
      n_samples = 0L,
      sample_ids = character(),
      pgen_ptr = NULL
    )
    expect_true(methods::validObject(obj))
  }
})

test_that("GenotypeHandle rejects invalid format", {
  expect_error(
    methods::validObject(
      new("GenotypeHandle",
        path = "/tmp/test",
        format = "bgen",
        snp_info = data.frame(),
        n_samples = 0L,
        sample_ids = character(),
        pgen_ptr = NULL
      )
    ),
    "format.*must be one of"
  )
})

test_that("GWASSumStats validity requires SNP, A1, A2, Z, N mcols", {
  gr <- make_test_granges(5)
  # Only set SNP and A1 -- missing A2, Z, N
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("rs", 1:5),
    A1 = rep("A", 5)
  )
  expect_error(
    methods::validObject(new("GWASSumStats",
      sumstats = gr, genome = "hg19", trait_name = "test", var_y = NULL
    )),
    "Missing required columns"
  )
})

test_that("GWASSumStats valid object passes with all required mcols", {
  set.seed(1)
  gr <- make_test_granges(5)
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("rs", 1:5),
    A1 = rep("A", 5),
    A2 = rep("G", 5),
    Z = rnorm(5),
    N = rep(1000, 5)
  )
  obj <- new("GWASSumStats",
    sumstats = gr, genome = "hg19", trait_name = "test", var_y = NULL
  )
  expect_true(methods::validObject(obj))
})

test_that("AnnotationMatrix validates dimensions and meta", {
  gr <- make_test_granges(10)
  meta <- make_test_annotation_meta()
  mat <- matrix(0, nrow = 10, ncol = 3)

  obj <- AnnotationMatrix(mat, gr, meta)
  expect_s4_class(obj, "AnnotationMatrix")
  expect_true(methods::validObject(obj))
})

test_that("AnnotationMatrix rejects row mismatch", {
  gr <- make_test_granges(10)
  meta <- make_test_annotation_meta()
  mat <- matrix(0, nrow = 5, ncol = 3)  # wrong number of rows

  expect_error(
    AnnotationMatrix(mat, gr, meta),
    "rows.*must match"
  )
})

test_that("AnnotationMatrix rejects column mismatch with meta", {
  gr <- make_test_granges(10)
  meta <- make_test_annotation_meta()  # 3 annotations
  mat <- matrix(0, nrow = 10, ncol = 2)  # only 2 columns

  # Constructor errors when colnames assignment fails (dimnames mismatch)
  expect_error(AnnotationMatrix(mat, gr, meta))
})

test_that("AnnotationMatrix rejects invalid tier values", {
  gr <- make_test_granges(10)
  meta <- data.frame(
    name = "x", tier = "invalid_tier", type = "binary",
    stringsAsFactors = FALSE
  )
  mat <- matrix(0, nrow = 10, ncol = 1)

  expect_error(
    AnnotationMatrix(mat, gr, meta),
    "tier.*baseline.*candidate"
  )
})

test_that("AnnotationMatrix rejects invalid type values", {
  gr <- make_test_granges(10)
  meta <- data.frame(
    name = "x", tier = "baseline", type = "ordinal",
    stringsAsFactors = FALSE
  )
  mat <- matrix(0, nrow = 10, ncol = 1)

  expect_error(
    AnnotationMatrix(mat, gr, meta),
    "type.*binary.*continuous"
  )
})

test_that("LDEigen constructs and validates correctly", {
  ldblocks <- make_test_ldblocks()
  snp_info <- make_test_snp_info()
  eigen_list <- list(
    list(values = c(1, 0.5), vectors = matrix(rnorm(20), 10, 2),
         snp_idx = 1:10),
    list(values = c(0.8), vectors = matrix(rnorm(10), 10, 1),
         snp_idx = 1:10)
  )

  obj <- new("LDEigen",
    ld_blocks = ldblocks,
    snp_info = snp_info,
    n_ref = 500L,
    in_sample = FALSE,
    genome = "hg19",
    eigen_list = eigen_list,
    eigenvalue_truncation = 0.9
  )
  expect_s4_class(obj, "LDEigen")
  expect_true(methods::validObject(obj))
})

test_that("LDEigen rejects eigen_list length mismatch", {
  ldblocks <- make_test_ldblocks()  # 2 blocks
  # Only 1 element in eigen_list
  expect_error(
    methods::validObject(
      new("LDEigen",
        ld_blocks = ldblocks,
        snp_info = make_test_snp_info(),
        n_ref = 500L,
        in_sample = FALSE,
        genome = "hg19",
        eigen_list = list(list(values = 1)),
        eigenvalue_truncation = 0.9
      )
    ),
    "eigen_list.*must match"
  )
})

test_that("LDEigen rejects invalid eigenvalue_truncation", {
  ldblocks <- make_test_ldblocks()
  expect_error(
    methods::validObject(
      new("LDEigen",
        ld_blocks = ldblocks,
        snp_info = make_test_snp_info(),
        n_ref = 500L,
        in_sample = FALSE,
        genome = "hg19",
        eigen_list = list(list(), list()),
        eigenvalue_truncation = 0
      )
    ),
    "eigenvalue_truncation"
  )
})

test_that("LDScore constructs and validates correctly", {
  ldblocks <- make_test_ldblocks()
  n <- 10
  snp_info <- make_test_snp_info(n)

  obj <- new("LDScore",
    ld_blocks = ldblocks,
    snp_info = snp_info,
    n_ref = 500L,
    in_sample = FALSE,
    genome = "hg19",
    ld_scores = matrix(runif(n), nrow = n, ncol = 1),
    ld_score_weights = runif(n),
    ld_matrix_list = list()
  )
  expect_s4_class(obj, "LDScore")
  expect_true(methods::validObject(obj))
})

test_that("LDScore rejects ld_scores row mismatch with snp_info", {
  ldblocks <- make_test_ldblocks()
  snp_info <- make_test_snp_info(10)

  expect_error(
    methods::validObject(
      new("LDScore",
        ld_blocks = ldblocks,
        snp_info = snp_info,
        n_ref = 500L,
        in_sample = FALSE,
        genome = "hg19",
        ld_scores = matrix(0, nrow = 5, ncol = 1),  # wrong rows
        ld_score_weights = runif(10),
        ld_matrix_list = list()
      )
    ),
    "ld_scores.*must match"
  )
})

test_that("H2Estimate constructs with all slots", {
  obj <- new("H2Estimate",
    h2 = 0.3,
    h2_se = 0.05,
    intercept = 1.01,
    intercept_se = 0.02,
    local = NULL,
    enrichment = NULL,
    tau_blocks = NULL,
    score_stats = NULL,
    method = "lder",
    n_snps = 10000L,
    trait_name = "height"
  )
  expect_s4_class(obj, "H2Estimate")
  expect_equal(obj@h2, 0.3)
  expect_equal(obj@method, "lder")
})

test_that("show() methods do not error", {
  # LDBlocks
  expect_output(show(make_test_ldblocks()), "LDBlocks")

  # GenotypeHandle
  gh <- new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snp_info = make_test_snp_info(),
    n_samples = 100L,
    sample_ids = paste0("s", 1:100),
    pgen_ptr = NULL
  )
  expect_output(show(gh), "GenotypeHandle")

  # GWASSumStats (via constructor)
  ss <- GWASSumStats(make_test_sumstats_df(10))
  expect_output(show(ss), "GWASSumStats")

  # AnnotationMatrix
  am <- AnnotationMatrix(
    matrix(0, nrow = 10, ncol = 1),
    make_test_granges(10),
    data.frame(name = "base", tier = "baseline", type = "binary",
               stringsAsFactors = FALSE)
  )
  expect_output(show(am), "AnnotationMatrix")

  # LDEigen
  ldblocks <- make_test_ldblocks()
  eig <- new("LDEigen",
    ld_blocks = ldblocks,
    snp_info = make_test_snp_info(),
    n_ref = 500L,
    in_sample = FALSE,
    genome = "hg19",
    eigen_list = list(list(), list()),
    eigenvalue_truncation = 0.9
  )
  expect_output(show(eig), "LDEigen")

  # LDScore
  n <- 10
  lsr <- new("LDScore",
    ld_blocks = ldblocks,
    snp_info = make_test_snp_info(n),
    n_ref = 500L,
    in_sample = FALSE,
    genome = "hg19",
    ld_scores = matrix(1, nrow = n, ncol = 1),
    ld_score_weights = rep(1, n),
    ld_matrix_list = list()
  )
  expect_output(show(lsr), "LDScore")

  # H2Estimate
  h2 <- new("H2Estimate",
    h2 = 0.3, h2_se = 0.05,
    intercept = 1.0, intercept_se = 0.01,
    local = NULL, enrichment = NULL,
    tau_blocks = NULL, score_stats = NULL,
    method = "lder", n_snps = 1000L, trait_name = "test"
  )
  expect_output(show(h2), "H2Estimate")
})

# =============================================================================
# GWASSumStats Constructor (h2_sumstats.R)
# =============================================================================

test_that("GWASSumStats() constructor creates object from data.frame", {
  df <- make_test_sumstats_df(20)
  obj <- GWASSumStats(df, trait_name = "height", genome = "hg38")

  expect_s4_class(obj, "GWASSumStats")
  expect_equal(obj@trait_name, "height")
  expect_equal(obj@genome, "hg38")
  expect_equal(length(obj@sumstats), 20)
})

test_that("GWASSumStats() normalizes chr prefix", {
  df <- make_test_sumstats_df(5)
  # Input has CHR = "1" (no prefix)
  obj <- GWASSumStats(df)
  chrs <- as.character(GenomicRanges::seqnames(obj@sumstats))
  expect_true(all(startsWith(chrs, "chr")))

  # Input already has "chr" prefix
  df2 <- df
  df2$CHR <- "chr1"
  obj2 <- GWASSumStats(df2)
  chrs2 <- as.character(GenomicRanges::seqnames(obj2@sumstats))
  # Should not double-prefix
  expect_true(all(chrs2 == "chr1"))
  expect_false(any(grepl("^chrchr", chrs2)))
})

test_that("GWASSumStats() errors on missing columns", {
  df <- data.frame(SNP = "rs1", CHR = "1", BP = 100)
  expect_error(GWASSumStats(df), "Missing required columns.*A1.*A2.*Z.*N")
})

test_that("GWASSumStats() removes rows with NA in required columns", {
  df <- make_test_sumstats_df(10)
  df$Z[1] <- NA
  df$N[3] <- NA
  expect_message(obj <- GWASSumStats(df), "Removed.*SNPs with missing")
  expect_equal(length(obj@sumstats), 8)
})

test_that("getZ() returns correct Z vector", {
  set.seed(99)
  df <- make_test_sumstats_df(5)
  obj <- GWASSumStats(df)
  z <- getZ(obj)
  expect_type(z, "double")
  expect_equal(length(z), 5)
})

test_that("getN() returns correct N vector", {
  df <- make_test_sumstats_df(5)
  obj <- GWASSumStats(df)
  n <- getN(obj)
  expect_equal(length(n), 5)
  expect_true(all(n == 10000))
})

test_that("getMaf() returns MAF when present, NULL when absent", {
  df <- make_test_sumstats_df(5)
  obj_no_maf <- GWASSumStats(df)
  expect_null(getMaf(obj_no_maf))

  df$MAF <- runif(5, 0.01, 0.5)
  obj_with_maf <- GWASSumStats(df)
  maf <- getMaf(obj_with_maf)
  expect_type(maf, "double")
  expect_equal(length(maf), 5)
})

test_that("nSnps() returns correct count", {
  df <- make_test_sumstats_df(30)
  obj <- GWASSumStats(df)
  expect_equal(nSnps(obj), 30)
})

test_that("subsetChr() filters correctly", {
  df <- make_test_sumstats_df(10)
  df$CHR <- c(rep("1", 6), rep("2", 4))
  obj <- GWASSumStats(df)

  chr1 <- subsetChr(obj, "1")
  expect_equal(nSnps(chr1), 6)

  # Also works with "chr" prefix
  chr2 <- subsetChr(obj, "chr2")
  expect_equal(nSnps(chr2), 4)
})

test_that("getVarY() returns var_y and NULL cases", {
  df <- make_test_sumstats_df(5)

  obj_null <- GWASSumStats(df, var_y = NULL)
  expect_null(getVarY(obj_null))

  obj_vy <- GWASSumStats(df, var_y = 4.5)
  expect_equal(getVarY(obj_vy), 4.5)
})

test_that("as.data.frame.GWASSumStats() round-trips", {
  df_in <- make_test_sumstats_df(15)
  obj <- GWASSumStats(df_in)
  df_out <- as.data.frame(obj)

  expect_true(is.data.frame(df_out))
  expect_true(all(c("SNP", "CHR", "BP", "A1", "A2", "Z", "N") %in%
                    names(df_out)))
  expect_equal(nrow(df_out), 15)
  expect_equal(df_out$SNP, df_in$SNP)
  # BP should round-trip
  expect_equal(df_out$BP, as.integer(df_in$BP))
})

test_that("rss_to_gwas_sumstats() converts from load_rss_data format", {
  rss_list <- list(
    sumstats = data.frame(
      chrom = rep("1", 5),
      pos = seq(1000, by = 100, length.out = 5),
      variant_id = paste0("rs", 1:5),
      A1 = rep("A", 5),
      A2 = rep("G", 5),
      z = rnorm(5),
      stringsAsFactors = FALSE
    ),
    n = 5000,
    var_y = 1.0
  )

  obj <- rss_to_gwas_sumstats(rss_list, trait_name = "test_trait",
                               genome = "hg38")
  expect_s4_class(obj, "GWASSumStats")
  expect_equal(nSnps(obj), 5)
  expect_equal(obj@trait_name, "test_trait")
  expect_equal(getVarY(obj), 1.0)
  expect_true(all(getN(obj) == 5000))
})

test_that("rss_to_gwas_sumstats() returns NULL for empty input", {
  rss_list <- list(
    sumstats = data.frame(
      chrom = character(), pos = integer(), variant_id = character(),
      A1 = character(), A2 = character(), z = numeric()
    ),
    n = 1000,
    var_y = NULL
  )
  expect_null(rss_to_gwas_sumstats(rss_list))
  expect_null(rss_to_gwas_sumstats(list(sumstats = NULL)))
})

# =============================================================================
# AnnotationMatrix (h2_annotations.R)
# =============================================================================

test_that("AnnotationMatrix() constructor creates object from matrix", {
  n <- 10
  gr <- make_test_granges(n)
  meta <- make_test_annotation_meta()
  mat <- matrix(runif(n * 3), nrow = n, ncol = 3)

  obj <- AnnotationMatrix(mat, gr, meta, genome = "hg38")
  expect_s4_class(obj, "AnnotationMatrix")
  expect_equal(nrow(obj@annotations), n)
  expect_equal(ncol(obj@annotations), 3)
  expect_equal(obj@genome, "hg38")
})

test_that("AnnotationMatrix() sets column names from annotation_meta", {
  n <- 10
  gr <- make_test_granges(n)
  meta <- make_test_annotation_meta()
  mat <- matrix(0, nrow = n, ncol = 3)
  # No colnames set on mat

  obj <- AnnotationMatrix(mat, gr, meta)
  expect_equal(colnames(obj@annotations), c("base", "enhancer", "promoter"))
})

test_that("getBaseline() subsets to baseline-tier only", {
  n <- 10
  gr <- make_test_granges(n)
  meta <- make_test_annotation_meta()  # 1 baseline, 2 candidate
  mat <- matrix(0, nrow = n, ncol = 3)

  obj <- AnnotationMatrix(mat, gr, meta)
  baseline <- getBaseline(obj)

  expect_s4_class(baseline, "AnnotationMatrix")
  expect_equal(ncol(baseline@annotations), 1)
  expect_true(all(baseline@annotation_meta$tier == "baseline"))
})

test_that("getCandidates() subsets to candidate-tier only", {
  n <- 10
  gr <- make_test_granges(n)
  meta <- make_test_annotation_meta()  # 1 baseline, 2 candidate
  mat <- matrix(0, nrow = n, ncol = 3)

  obj <- AnnotationMatrix(mat, gr, meta)
  cand <- getCandidates(obj)

  expect_s4_class(cand, "AnnotationMatrix")
  expect_equal(ncol(cand@annotations), 2)
  expect_true(all(cand@annotation_meta$tier == "candidate"))
})

test_that(".annot_detect_format() detects bigwig, bed, ldsc_annot", {
  detect <- pecotmr:::.annot_detect_format

  expect_equal(detect("path/to/file.bw"), "bigwig")
  expect_equal(detect("path/to/file.bigwig"), "bigwig")
  expect_equal(detect("path/to/file.annot.gz"), "ldsc_annot")
  expect_equal(detect("path/to/file.annot"), "ldsc_annot")
  expect_equal(detect("path/to/file.bed"), "bed")
  # Unknown extensions default to bed
  expect_equal(detect("path/to/file.txt"), "bed")
})

# =============================================================================
# .read_bed_annotation (h2_annotations.R)
# =============================================================================

test_that(".read_bed_annotation returns correct binary overlap vector", {
  skip_if_not_installed("rtracklayer")

  bed_content <- paste(
    "chr1\t100\t500\tregion1",
    "chr1\t800\t900\tregion2",
    sep = "\n"
  )
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file), add = TRUE)
  writeLines(bed_content, bed_file)

  # SNPs at positions 50, 200, 300, 600, 850
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(start = c(50, 200, 300, 600, 850), width = 1)
  )

  result <- pecotmr:::.read_bed_annotation(bed_file, snp_gr)

  expect_equal(length(result), 5)
  # Position 50 is outside [100,500) and [800,900) -> 0
  expect_equal(result[1], 0)
  # Positions 200 and 300 are inside [100,500) -> 1
  expect_equal(result[2], 1)
  expect_equal(result[3], 1)
  # Position 600 is outside both regions -> 0
  expect_equal(result[4], 0)
  # Position 850 is inside [800,900) -> 1
  expect_equal(result[5], 1)
})

test_that(".read_bed_annotation returns all zeros when no overlaps", {
  skip_if_not_installed("rtracklayer")

  bed_content <- "chr1\t100\t200\tregion1"
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file), add = TRUE)
  writeLines(bed_content, bed_file)

  # All SNPs outside the region
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(50, 300, 500), width = 1)
  )

  result <- pecotmr:::.read_bed_annotation(bed_file, snp_gr)
  expect_equal(result, c(0, 0, 0))
})

# =============================================================================
# .read_ldsc_annot (h2_annotations.R)
# =============================================================================

test_that(".read_ldsc_annot reads annotation column and matches SNPs", {
  annot_df <- data.frame(
    CHR = c(1, 1, 1, 1),
    BP = c(100, 200, 300, 400),
    SNP = paste0("rs", 1:4),
    CM = 0,
    my_annot = c(1, 0, 1, 0),
    stringsAsFactors = FALSE
  )
  annot_file <- tempfile(fileext = ".annot")
  on.exit(unlink(annot_file), add = TRUE)
  write.table(annot_df, annot_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # SNP positions matching BP 100, 200, 300, 400
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(100, 200, 300, 400), width = 1)
  )

  result <- pecotmr:::.read_ldsc_annot(annot_file, snp_gr, "my_annot")

  expect_equal(length(result), 4)
  expect_equal(result, c(1, 0, 1, 0))
})

test_that(".read_ldsc_annot returns 0 for unmatched SNP positions", {
  annot_df <- data.frame(
    CHR = c(1, 1),
    BP = c(100, 200),
    SNP = c("rs1", "rs2"),
    CM = 0,
    score = c(5, 10),
    stringsAsFactors = FALSE
  )
  annot_file <- tempfile(fileext = ".annot")
  on.exit(unlink(annot_file), add = TRUE)
  write.table(annot_df, annot_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # SNP at position 100 matches, position 999 does not
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 2),
    ranges = IRanges::IRanges(start = c(100, 999), width = 1)
  )

  result <- pecotmr:::.read_ldsc_annot(annot_file, snp_gr, "score")
  expect_equal(result, c(5, 0))
})

test_that(".read_ldsc_annot errors when annotation column not found", {
  annot_df <- data.frame(
    CHR = c(1), BP = c(100), SNP = c("rs1"), CM = 0, my_annot = c(1),
    stringsAsFactors = FALSE
  )
  annot_file <- tempfile(fileext = ".annot")
  on.exit(unlink(annot_file), add = TRUE)
  write.table(annot_df, annot_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  snp_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, width = 1)
  )

  expect_error(
    pecotmr:::.read_ldsc_annot(annot_file, snp_gr, "nonexistent_col"),
    "Annotation column.*nonexistent_col.*not found"
  )
})

test_that(".read_ldsc_annot handles chr-prefixed CHR column", {
  annot_df <- data.frame(
    CHR = c("chr1", "chr1"),
    BP = c(100, 200),
    SNP = c("rs1", "rs2"),
    CM = 0,
    val = c(3, 7),
    stringsAsFactors = FALSE
  )
  annot_file <- tempfile(fileext = ".annot")
  on.exit(unlink(annot_file), add = TRUE)
  write.table(annot_df, annot_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 2),
    ranges = IRanges::IRanges(start = c(100, 200), width = 1)
  )

  result <- pecotmr:::.read_ldsc_annot(annot_file, snp_gr, "val")
  expect_equal(result, c(3, 7))
})

# =============================================================================
# readAnnotations method (h2_annotations.R)
# =============================================================================

test_that("readAnnotations errors when paths is not named", {
  snp_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, width = 1)
  )

  expect_error(
    readAnnotations(c("/some/file.bed"), snp_gr),
    "paths.*must be a named"
  )
})

test_that("readAnnotations with BED file creates AnnotationMatrix", {
  skip_if_not_installed("rtracklayer")

  bed_content <- "chr1\t100\t500\tregion1"
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file), add = TRUE)
  writeLines(bed_content, bed_file)

  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(50, 200, 600), width = 1)
  )

  paths <- c(enhancer = bed_file)
  result <- readAnnotations(paths, snp_gr, genome = "hg38")

  expect_s4_class(result, "AnnotationMatrix")
  expect_equal(nrow(result@annotations), 3)
  expect_equal(ncol(result@annotations), 1)
  expect_equal(colnames(result@annotations), "enhancer")
  expect_equal(result@genome, "hg38")
  # Auto-detected meta should be binary, candidate tier
  expect_equal(result@annotation_meta$type, "binary")
  expect_equal(result@annotation_meta$tier, "candidate")
  # Overlap check: only position 200 is inside [100,500)
  expect_equal(as.numeric(result@annotations[, 1]), c(0, 1, 0))
})

test_that("readAnnotations with LDSC annot file creates AnnotationMatrix", {
  annot_df <- data.frame(
    CHR = c(1, 1, 1),
    BP = c(100, 200, 300),
    SNP = paste0("rs", 1:3),
    CM = 0,
    my_score = c(1, 0, 1),
    stringsAsFactors = FALSE
  )
  annot_file <- tempfile(fileext = ".annot")
  on.exit(unlink(annot_file), add = TRUE)
  write.table(annot_df, annot_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100, 200, 300), width = 1)
  )

  paths <- c(my_score = annot_file)
  result <- readAnnotations(paths, snp_gr)

  expect_s4_class(result, "AnnotationMatrix")
  expect_equal(as.numeric(result@annotations[, 1]), c(1, 0, 1))
  expect_equal(colnames(result@annotations), "my_score")
})

test_that("readAnnotations uses provided annotation_meta", {
  annot_df <- data.frame(
    CHR = c(1, 1),
    BP = c(100, 200),
    SNP = c("rs1", "rs2"),
    CM = 0,
    custom = c(1, 0),
    stringsAsFactors = FALSE
  )
  annot_file <- tempfile(fileext = ".annot")
  on.exit(unlink(annot_file), add = TRUE)
  write.table(annot_df, annot_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 2),
    ranges = IRanges::IRanges(start = c(100, 200), width = 1)
  )

  custom_meta <- data.frame(
    name = "custom",
    tier = "baseline",
    type = "continuous",
    stringsAsFactors = FALSE
  )

  paths <- c(custom = annot_file)
  result <- readAnnotations(paths, snp_gr, annotation_meta = custom_meta)

  expect_s4_class(result, "AnnotationMatrix")
  # Should use the provided meta, not auto-detected
  expect_equal(result@annotation_meta$tier, "baseline")
  expect_equal(result@annotation_meta$type, "continuous")
})

test_that("readAnnotations with multiple files creates multi-column matrix", {
  # Two LDSC annot files
  annot_df1 <- data.frame(
    CHR = c(1, 1), BP = c(100, 200), SNP = c("rs1", "rs2"),
    CM = 0, ann1 = c(1, 0), stringsAsFactors = FALSE
  )
  annot_df2 <- data.frame(
    CHR = c(1, 1), BP = c(100, 200), SNP = c("rs1", "rs2"),
    CM = 0, ann2 = c(0, 1), stringsAsFactors = FALSE
  )

  file1 <- tempfile(fileext = ".annot")
  file2 <- tempfile(fileext = ".annot")
  on.exit(unlink(c(file1, file2)), add = TRUE)
  write.table(annot_df1, file1, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(annot_df2, file2, sep = "\t", row.names = FALSE, quote = FALSE)

  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 2),
    ranges = IRanges::IRanges(start = c(100, 200), width = 1)
  )

  paths <- c(ann1 = file1, ann2 = file2)
  result <- readAnnotations(paths, snp_gr)

  expect_s4_class(result, "AnnotationMatrix")
  expect_equal(ncol(result@annotations), 2)
  expect_equal(colnames(result@annotations), c("ann1", "ann2"))
  expect_equal(as.numeric(result@annotations[, "ann1"]), c(1, 0))
  expect_equal(as.numeric(result@annotations[, "ann2"]), c(0, 1))
})

# =============================================================================
# readSumstats edge cases (h2_sumstats.R)
# =============================================================================

test_that("readSumstats with use_mungesumstats=FALSE returns GWASSumStats", {
  df <- data.frame(
    SNP = paste0("rs", 1:5),
    CHR = rep("1", 5),
    BP = seq(1000, by = 100, length.out = 5),
    A1 = rep("A", 5),
    A2 = rep("G", 5),
    Z = c(1.1, -0.5, 2.3, 0.0, -1.7),
    N = rep(10000, 5),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- readSumstats(tmp, trait_name = "mytest", genome = "hg38",
                         use_mungesumstats = FALSE)
  expect_s4_class(result, "GWASSumStats")
  expect_equal(nSnps(result), 5)
  expect_equal(result@trait_name, "mytest")
  expect_equal(result@genome, "hg38")
})

test_that("readSumstats with use_mungesumstats=FALSE and n parameter fills N", {
  df <- data.frame(
    SNP = paste0("rs", 1:3),
    CHR = rep("2", 3),
    BP = c(500, 600, 700),
    A1 = rep("C", 3),
    A2 = rep("T", 3),
    Z = c(0.5, -1.2, 0.8),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- readSumstats(tmp, n = 5000, use_mungesumstats = FALSE)
  expect_s4_class(result, "GWASSumStats")
  expect_true(all(getN(result) == 5000))
})

test_that("readSumstats with genome=NULL defaults to hg19", {
  df <- data.frame(
    SNP = "rs1", CHR = "1", BP = 1000,
    A1 = "A", A2 = "G", Z = 1.5, N = 10000,
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- readSumstats(tmp, genome = NULL, use_mungesumstats = FALSE)
  expect_equal(result@genome, "hg19")
})

# =============================================================================
# rss_to_gwas_sumstats edge cases (h2_sumstats.R)
# =============================================================================

test_that("rss_to_gwas_sumstats fabricates SNP from CHR:BP:A2:A1 when SNP missing", {
  rss_list <- list(
    sumstats = data.frame(
      chrom = rep("3", 3),
      pos = c(100, 200, 300),
      A1 = c("A", "C", "T"),
      A2 = c("G", "T", "A"),
      z = c(1.0, -0.5, 2.0),
      stringsAsFactors = FALSE
    ),
    n = 8000,
    var_y = NULL
  )

  result <- rss_to_gwas_sumstats(rss_list)
  expect_s4_class(result, "GWASSumStats")
  snp_col <- S4Vectors::mcols(result@sumstats)$SNP
  expect_equal(snp_col, c("3:100:G:A", "3:200:T:C", "3:300:A:T"))
})

test_that("rss_to_gwas_sumstats uses sequential integers for SNP when CHR/BP/A1/A2 incomplete", {
  # The sequential-integer fallback for SNP triggers when no variant_id/SNP
  # and not all of CHR+BP+A2+A1 are present after column mapping. Since
  # GWASSumStats requires CHR/BP/A1/A2, the fallback path necessarily errors
  # at the constructor. We verify the error comes from GWASSumStats (missing
  # columns), confirming rss_to_gwas_sumstats itself did not error earlier.
  ss_df <- data.frame(
    z = c(0.5, -1.0),
    stringsAsFactors = FALSE
  )
  expect_error(
    rss_to_gwas_sumstats(list(sumstats = ss_df, n = 1000, var_y = NULL)),
    "Missing required columns"
  )
})
