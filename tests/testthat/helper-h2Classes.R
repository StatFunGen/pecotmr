# Tests for S4 classes (h2_classes.R), GwasSumStats (h2_sumstats.R),
# and AnnotationMatrix (h2Annotations.R)

# =============================================================================
# Test data helpers
# =============================================================================

make_test_granges <- function(n = 10) {
    GenomicRanges::GRanges(
        seqnames = rep("chr1", n),
        ranges = IRanges::IRanges(
            start = seq(1000, by = 100, length.out = n),
            width = 1L
        )
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
    blocks_gr
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
    new(
        "GenotypeHandle",
        path = "/tmp/test.gds",
        format = "gds",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
    )
}

.dfToSumstatsGr <- function(df) {
    chrs <- as.character(df$CHR)
    if (!all(grepl("^chr", chrs))) {
        chrs <- paste0("chr", sub("^chr", "", chrs))
    }
    gr <- GenomicRanges::GRanges(
        seqnames = chrs,
        ranges = IRanges::IRanges(start = as.integer(df$BP), width = 1L)
    )
    mc <- df[, setdiff(colnames(df), c("CHR", "BP")), drop = FALSE]
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mc)
    gr
}

makeGwasSumStatsFromDf <- function(
    df,
    traitName = "test",
    genome = "hg19",
    varY = NA_real_
) {
    required <- c("SNP", "A1", "A2", "Z", "N")
    missingCols <- setdiff(required, colnames(df))
    if (length(missingCols) > 0L) {
        stop("Missing required columns: ", paste(missingCols, collapse = ", "))
    }
    keep <- stats::complete.cases(df[, required, drop = FALSE])
    if (!all(keep)) {
        message(sprintf(
            "Removed %d SNPs with missing required-column values.",
            sum(!keep)
        ))
    }
    df <- df[keep, , drop = FALSE]
    if (is.null(varY)) {
        varY <- NA_real_
    }
    GwasSumStats(
        study = traitName,
        entry = list(.dfToSumstatsGr(df)),
        genome = genome,
        ldSketch = .testGenotypeHandle(),
        varY = varY
    )
}

# A handle over a real (fake-path) variant set, so tests can build the
# per-block LD panels a block-parallel pipeline produces: same genotype source,
# each narrowed to its own variants.
.blockGenotypeHandle <- function(bp, path = "/tmp/test.pgen") {
    new(
        "GenotypeHandle",
        path = path,
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("chr1:", bp, ":A:G"),
            CHR = rep("1", length(bp)),
            BP = as.integer(bp),
            A1 = rep("G", length(bp)),
            A2 = rep("A", length(bp)),
            fileIdx = seq_along(bp),
            stringsAsFactors = FALSE
        ),
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10)),
        pgenPtr = NULL
    )
}

# One LD block's worth of GWAS summary statistics: `n` variants starting at
# `start`, carrying a QC record and an LD panel narrowed to just those
# variants (what a per-block pipeline step writes out).
makeGwasBlock <- function(blockId, start, n = 5L, genome = "hg19",
                          qcOptions = list(mafCutoff = 0.01)) {
    bp <- seq(start, by = 100L, length.out = n)
    df <- data.frame(
        SNP = paste0("chr1:", bp, ":A:G"),
        CHR = rep("1", n),
        BP = bp,
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq_len(n) / 10,
        N = rep(1000L, n),
        stringsAsFactors = FALSE
    )
    obj <- GwasSumStats(
        study = "trait1",
        entry = list(.dfToSumstatsGr(df)),
        genome = genome,
        ldSketch = .blockGenotypeHandle(bp),
        blockId = blockId,
        qcInfo = list(
            timestamp = NA_character_,
            options = qcOptions,
            entryAudit = list(list(variantsIn = n, block = blockId))
        )
    )
    obj
}
