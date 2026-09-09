# Extracted from test_ctwasPipeline.R:468

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "pecotmr", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("ctwasPipeline")
.ctp_makeHandle <- function(snp_n = 12L, n_samples = 30L) {
    # Use a per-process tempfile so .ctwasLdPanelKey's file.exists check
    # succeeds against the fixture handle (real LD-sketch payloads exist
    # by construction; mock fixtures need an equivalent on-disk anchor).
    gdsPath <- file.path(tempdir(), "ctp_sketch.gds")
    if (!file.exists(gdsPath)) {
        file.create(gdsPath)
    }
    positions <- seq(100L, by = 100L, length.out = snp_n)
    # SNP IDs follow the canonical chr:pos:A2:A1 layout so allele
    # harmonization inside .ctwasBuildWeights / .ctwasHarmonizeWeights can
    # parse them via parseVariantId().
    snpIds <- sprintf("chr1:%d:G:A", positions)
    new(
        "GenotypeHandle",
        path = gdsPath,
        format = "gds",
        snpInfo = data.frame(
            SNP = snpIds,
            CHR = rep("1", snp_n),
            BP = positions,
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = n_samples,
        sampleIds = paste0("s", seq_len(n_samples)),
        pgenPtr = NULL
    )
}
.ctp_snpId <- function(i) sprintf("chr1:%d:G:A", 100L * i)
.ctp_mockExtractor <- function(seed = 5, n_samples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(getSnpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(getSnpInfo(handle)),
            dimnames = list(getSampleIds(handle), getSnpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", getSnpInfo(handle)$CHR[snpIdx]),
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
        cd <- S4Vectors::DataFrame(
            sampleId = getSampleIds(handle),
            row.names = getSampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- getSnpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- getSampleIds(handle)
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = cd
        )
    }
}
.ctp_blockVariants <- function(b) {
    idx <- seq_len(6L) + 6L * (b - 1L)
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100L * idx, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = vapply(idx, .ctp_snpId, character(1)),
        A1 = rep("A", 6),
        A2 = rep("G", 6),
        Z = rnorm(6),
        N = rep(1000L, 6)
    )
    gr
}
.ctp_makeGwasSumstats <- function(qc = TRUE, blockIds = c("block1", "block2")) {
    n <- length(blockIds)
    GwasSumStats(
        study = rep("G1", n),
        entry = map(seq_len(n), .ctp_blockVariants),
        genome = "hg19",
        ldSketch = .ctp_makeHandle(snp_n = 6L * n),
        blockId = blockIds,
        qcInfo = if (qc) list(step1 = "ok") else list()
    )
}
.ctp_makeTwasWeights <- function(variantIdx = 1:5) {
    e <- twasWeightsRow(
        variantIds = vapply(variantIdx, .ctp_snpId, character(1)),
        weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
    )
    TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e),
        ldSketch = .ctp_makeHandle()
    )
}
.ctp_makeFmrWeightSource <- function() {
    vids <- vapply(1:5, .ctp_snpId, character(1))
    w <- c(0.1, 0.05, -0.2, 0.3, 0.0)
    tl <- data.frame(
        variant_id = vids,
        chrom = rep("1", 5),
        pos = 100L * (1:5),
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        N = rep(100, 5),
        af = rep(0.3, 5),
        pip = rep(0.5, 5),
        posterior_mean = w,
        posterior_sd = rep(0.1, 5),
        cs_95 = rep("susie_1", 5),
        stringsAsFactors = FALSE
    )
    fe <- fineMappingRow(variantIds = vids, susieFit = list(), topLoci = tl)
    QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fe),
        ldSketch = .ctp_makeHandle()
    )
}
.ctp_makeMultiBlockInputs <- function(qc = TRUE) {
    ss <- .ctp_makeGwasSumstats(qc = qc)
    list(
        gwasSumStats = ss,
        # ctwas fine-maps one region at a time, so each block's gene has to
        # carry weights on variants that block actually has LD for: block 1
        # holds variants 1..6, block 2 holds 7..12.
        twasWeights = list(
            block1 = .ctp_makeTwasWeights(),
            block2 = .ctp_makeTwasWeights(7:11)
        )
    )
}

# test -------------------------------------------------------------------------
ss <- .ctp_makeGwasSumstats(blockIds = "block1")
gwasIds <- vapply(1:6, .ctp_snpId, character(1))
panelIds <- gwasIds
panelIds[2] <- "chr1:200:A:G"
plain <- pecotmr:::.ctwasBuildZSnp(ss, gwasIds)
flipped <- pecotmr:::.ctwasBuildZSnp(ss, panelIds)
i <- match("chr1:200:G:A", plain$id)
expect_equal(flipped$id[i], "chr1:200:A:G")
