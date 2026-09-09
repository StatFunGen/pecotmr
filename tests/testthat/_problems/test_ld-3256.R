# Extracted from test_ld.R:3256

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "pecotmr", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("LD")
library(tidyverse)
make_test_ld_data <- function(variant_ids, R = NULL, blockMetadata = NULL) {
    if (is.null(R)) {
        p <- length(variant_ids)
        R <- diag(p)
        rownames(R) <- colnames(R) <- variant_ids
    }
    ref_panel <- pecotmr:::parseVariantId(variant_ids)
    ref_panel$variant_id <- variant_ids
    variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
    if (is.null(blockMetadata)) {
        blockMetadata <- data.frame(
            blockId = 1L,
            chrom = as.character(ref_panel$chrom[1]),
            blockStart = min(ref_panel$pos),
            blockEnd = max(ref_panel$pos),
            size = length(variant_ids),
            startIdx = 1L,
            endIdx = length(variant_ids),
            stringsAsFactors = FALSE
        )
    }
    LdData(
        correlation = R,
        variants = variants_gr,
        blockMetadata = blockMetadata
    )
}
generate_dummy_data <- function() {
    region <- data.frame(
        chrom = "chr1",
        start = c(1000),
        end = c(1190)
    )
    meta_df <- data.frame(
        chrom = "chr1",
        start = c(1000, 1200, 1400, 1600, 1800),
        end = c(1200, 1400, 1600, 1800, 2000),
        path = c(
            "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
            "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
            "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
            "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
            "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
        )
    )
    return(list(region = region, meta = meta_df))
}
generate_multi_block_data <- function() {
    region <- data.frame(
        chrom = "chr1",
        start = c(1000),
        end = c(1500)
    )
    meta_df <- data.frame(
        chrom = "chr1",
        start = c(1000, 1200, 1400, 1600, 1800),
        end = c(1200, 1400, 1600, 1800, 2000),
        path = c(
            "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
            "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
            "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
            "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
            "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
        )
    )
    return(list(region = region, meta = meta_df))
}
geno_test_data_dir <- test_path("test_data")
geno_region_all <- "chr21:17513228-17592874"
test_data_dir <- test_path("test_data")
library(testthat)
.lds_makeHandle <- function(snpN = 6L, nSamples = 30L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(snpN)),
            CHR = rep("1", snpN),
            BP = seq(100L, by = 100L, length.out = snpN),
            A1 = rep("A", snpN),
            A2 = rep("G", snpN),
            stringsAsFactors = FALSE
        ),
        nSamples = nSamples,
        sampleIds = sprintf("s%d", seq_len(nSamples)),
        pgenPtr = NULL
    )
}
.lds_mockExtractor <- function(seed = 7, nSamples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        nSnp <- nrow(getSnpInfo(handle))
        panel <- matrix(
            rbinom(nSamples * nSnp, 2, 0.3),
            nrow = nSamples,
            ncol = nSnp,
            dimnames = list(getSampleIds(handle), getSnpInfo(handle)$SNP)
        )
        rr <- GenomicRanges::GRanges(
            seqnames = str_c("chr", getSnpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = getSnpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = getSnpInfo(handle)$SNP[snpIdx],
            A1 = getSnpInfo(handle)$A1[snpIdx],
            A2 = getSnpInfo(handle)$A2[snpIdx]
        )
        dosage <- t(panel[, snpIdx, drop = FALSE])
        dimnames(dosage) <- list(
            getSnpInfo(handle)$SNP[snpIdx],
            getSampleIds(handle)
        )
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = S4Vectors::DataFrame(
                sampleId = getSampleIds(handle),
                row.names = getSampleIds(handle)
            )
        )
    }
}

# test -------------------------------------------------------------------------
h <- .lds_makeHandle()
local_mocked_bindings(
    extractBlockGenotypes = .lds_mockExtractor(),
    .package = "pecotmr"
)
same <- c("chr1:200:A:G", "chr1:400:A:G", "chr1:500:A:G")
flipped <- c("chr1:200:G:A", "chr1:400:A:G", "chr1:500:A:G")
rSame <- pecotmr:::.ldFromSketch(h, same)
rFlip <- pecotmr:::.ldFromSketch(h, flipped)
sgn <- c(-1, 1, 1)
expect_equal(
    unname(rFlip),
    unname(rSame) * outer(sgn, sgn),
    tolerance = 1e-12
)
