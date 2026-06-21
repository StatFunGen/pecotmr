# Tests for S4 classes (h2_classes.R), GwasSumStats (h2_sumstats.R),
# and AnnotationMatrix (h2Annotations.R)

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
  new("LdBlocks", blocks = blocks_gr, genome = "hg19")
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

# Bridge helper: turn the legacy per-study data.frame shape into a
# single-row GwasSumStats collection using the new API. Keeps the bulk of
# the per-accessor tests below readable.
.testGenotypeHandle <- function() {
  new("GenotypeHandle",
    path = "/tmp/test.gds", format = "gds",
    snpInfo = data.frame(), nSamples = 0L,
    sampleIds = character(), pgenPtr = NULL)
}

.dfToSumstatsGr <- function(df) {
  chrs <- as.character(df$CHR)
  if (!all(grepl("^chr", chrs))) chrs <- paste0("chr", sub("^chr", "", chrs))
  gr <- GenomicRanges::GRanges(
    seqnames = chrs,
    ranges = IRanges::IRanges(start = as.integer(df$BP), width = 1L)
  )
  mc <- df[, setdiff(colnames(df), c("CHR", "BP")), drop = FALSE]
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mc)
  gr
}

makeGwasSumStatsFromDf <- function(df, traitName = "test",
                                   genome = "hg19", varY = NA_real_) {
  required <- c("SNP", "A1", "A2", "Z", "N")
  missingCols <- setdiff(required, colnames(df))
  if (length(missingCols) > 0L)
    stop("Missing required columns: ",
         paste(missingCols, collapse = ", "))
  keep <- stats::complete.cases(df[, required, drop = FALSE])
  if (!all(keep))
    message(sprintf("Removed %d SNPs with missing required-column values.",
                    sum(!keep)))
  df <- df[keep, , drop = FALSE]
  if (is.null(varY)) varY <- NA_real_
  GwasSumStats(
    study    = traitName,
    entry    = list(.dfToSumstatsGr(df)),
    genome   = genome,
    ldSketch = .testGenotypeHandle(),
    varY     = varY)
}


