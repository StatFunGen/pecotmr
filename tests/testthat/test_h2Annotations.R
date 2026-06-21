# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (h2Annotations) ===

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


test_that("getbaseline() subsets to baseline-tier only", {
  n <- 10
  gr <- make_test_granges(n)
  meta <- make_test_annotation_meta()  # 1 baseline, 2 candidate
  mat <- matrix(0, nrow = n, ncol = 3)

  obj <- AnnotationMatrix(mat, gr, meta)
  baseline <- getBaseline(obj)

  expect_s4_class(baseline, "AnnotationMatrix")
  expect_equal(ncol(baseline@annotations), 1)
  expect_true(all(baseline@annotationMeta$tier == "baseline"))
})


test_that("getcandidates() subsets to candidate-tier only", {
  n <- 10
  gr <- make_test_granges(n)
  meta <- make_test_annotation_meta()  # 1 baseline, 2 candidate
  mat <- matrix(0, nrow = n, ncol = 3)

  obj <- AnnotationMatrix(mat, gr, meta)
  cand <- getCandidates(obj)

  expect_s4_class(cand, "AnnotationMatrix")
  expect_equal(ncol(cand@annotations), 2)
  expect_true(all(cand@annotationMeta$tier == "candidate"))
})


test_that(".annot_detect_format() detects bigwig, bed, ldsc_annot", {
  detect <- pecotmr:::.annotDetectFormat

  expect_equal(detect("path/to/file.bw"), "bigwig")
  expect_equal(detect("path/to/file.bigwig"), "bigwig")
  expect_equal(detect("path/to/file.annot.gz"), "ldsc_annot")
  expect_equal(detect("path/to/file.annot"), "ldsc_annot")
  expect_equal(detect("path/to/file.bed"), "bed")
  # Unknown extensions default to bed
  expect_equal(detect("path/to/file.txt"), "bed")
})

# =============================================================================
# .read_bed_annotation (h2Annotations.R)
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

  result <- pecotmr:::.readBedAnnotation(bed_file, snp_gr)

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

  result <- pecotmr:::.readBedAnnotation(bed_file, snp_gr)
  expect_equal(result, c(0, 0, 0))
})

# =============================================================================
# .read_ldsc_annot (h2Annotations.R)
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

  result <- pecotmr:::.readLdscAnnot(annot_file, snp_gr, "my_annot")

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

  result <- pecotmr:::.readLdscAnnot(annot_file, snp_gr, "score")
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
    pecotmr:::.readLdscAnnot(annot_file, snp_gr, "nonexistent_col"),
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

  result <- pecotmr:::.readLdscAnnot(annot_file, snp_gr, "val")
  expect_equal(result, c(3, 7))
})

# =============================================================================
# readannotations method (h2Annotations.R)
# =============================================================================


test_that("readannotations errors when paths is not named", {
  snp_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, width = 1)
  )

  expect_error(
    readAnnotations(c("/some/file.bed"), snp_gr),
    "paths.*must be a named"
  )
})


test_that("readannotations with BED file creates AnnotationMatrix", {
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
  expect_equal(result@annotationMeta$type, "binary")
  expect_equal(result@annotationMeta$tier, "candidate")
  # Overlap check: only position 200 is inside [100,500)
  expect_equal(as.numeric(result@annotations[, 1]), c(0, 1, 0))
})


test_that("readannotations with LDSC annot file creates AnnotationMatrix", {
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


test_that("readannotations uses provided annotation_meta", {
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
  result <- readAnnotations(paths, snp_gr, annotationMeta = custom_meta)

  expect_s4_class(result, "AnnotationMatrix")
  # Should use the provided meta, not auto-detected
  expect_equal(result@annotationMeta$tier, "baseline")
  expect_equal(result@annotationMeta$type, "continuous")
})


test_that("readannotations with multiple files creates multi-column matrix", {
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
# readsumstats edge cases (h2_sumstats.R)
# =============================================================================


# =============================================================================

